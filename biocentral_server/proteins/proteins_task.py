import os
import shutil
import tempfile
from typing import Callable, Dict, List

from pymmseqs.commands import easy_cluster, easy_linclust
from ..server_management import TaskInterface, TaskDTO, TaskStatus

class ClusterSequencesTask(TaskInterface):
    """
    Task to cluster sequences via pymmseqs by converting an input dictionary
    to a temporary FASTA file, running the alignment, and clean-up.
    """

    def __init__(
        self,
        sequence_data: Dict[str, str],
        sequence_identity_threshold: float,
        use_linear_clustering: bool = False,
    ):
        self.sequence_data = sequence_data
        self.sequence_identity_threshold = sequence_identity_threshold
        self.use_linear_clustering = use_linear_clustering

    def run_task(self, update_dto_callback: Callable) -> TaskDTO:
        # Set task status to RUNNING
        update_dto_callback(TaskDTO(status=TaskStatus.RUNNING))

        temp_input_path = None
        temp_output_prefix = None

        try:
            # 1. Open Context Manager for temporary FASTA file
            with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta", delete=False) as temp_input:
                for seq_id, seq in self.sequence_data.items():
                    temp_input.write(f">{seq_id}\n{seq}\n")
                temp_input_path = temp_input.name

                # 2. Create a unique prefix for MMseqs2 output files within the context
                temp_output_prefix = os.path.join(tempfile.gettempdir(), f"mmseqs_out_{os.getpid()}")

                print("Running pymmseqs command from temporary FASTA file...")

                # 3. Execute MMseqs2 clustering inside the context scope
                if self.use_linear_clustering:
                    easy_linclust(
                        temp_input_path,
                        temp_output_prefix,
                        tempfile.gettempdir(),
                        min_seq_id=self.sequence_identity_threshold,
                    )
                else:
                    easy_cluster(
                        temp_input_path,
                        temp_output_prefix,
                        tempfile.gettempdir(),
                        min_seq_id=self.sequence_identity_threshold,
                    )

            # 4. Parse the generated TSV file to map representatives to their cluster members
            tsv_file = f"{temp_output_prefix}_cluster.tsv"
            clustered_results: Dict[str, List[str]] = {}

            if not os.path.exists(tsv_file):
                raise FileNotFoundError("MMseqs2 did not generate the expected TSV cluster file.")

            with open(tsv_file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        rep_id, member_id = parts[0], parts[1]
                        if rep_id not in clustered_results:
                            clustered_results[rep_id] = []
                        clustered_results[rep_id].append(member_id)

            # Return the finished DTO with the mapped cluster IDs
            return TaskDTO(status=TaskStatus.FINISHED, clustered_data=clustered_results)

        except Exception as e:
            print(f"Error during clustering execution: {str(e)}")
            return TaskDTO.errored(f"Clustering failed: {str(e)}")

        finally:
            # 5. CLEANUP: Delete temporary input file
            if temp_input_path and os.path.exists(temp_input_path):
                try:
                    os.unlink(temp_input_path)
                except Exception:
                    pass

            # Clean up generated MMseqs output files
            if temp_output_prefix:
                for suffix in ["_rep_seq.fasta", "_all_seqs.fasta", "_cluster.tsv", ""]:
                    file_to_delete = f"{temp_output_prefix}{suffix}"
                    if os.path.exists(file_to_delete):
                        try:
                            if os.path.isdir(file_to_delete):
                                shutil.rmtree(file_to_delete)
                            else:
                                os.unlink(file_to_delete)
                        except Exception:
                            pass