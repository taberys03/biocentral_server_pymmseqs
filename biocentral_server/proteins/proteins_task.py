import tempfile
from pathlib import Path
from typing import List, Callable
from pymmseqs.commands import easy_cluster, easy_linclust

# 1. Natively import his read_FASTA parser and record structure
from biotrainer.input_files.fasta import read_FASTA
from biotrainer.input_files import BiotrainerSequenceRecord

from ..server_management import TaskInterface, TaskDTO, TaskStatus, EmbeddingProgress

class ClusterSequencesTask(TaskInterface):
    """Runs MMseqs2 sequence clustering asynchronously in the background queue"""

    def __init__(
        self,
        file_path: str,
        sequence_identity_threshold: float = 0.5,
        use_linear_clustering: bool = False
    ):
        self.file_path = file_path
        self.sequence_identity_threshold = sequence_identity_threshold
        self.use_linear_clustering = use_linear_clustering

    def run_task(self, update_dto_callback: Callable) -> TaskDTO:
        try:
            # 2. Use Sebastian's read_FASTA function to parse the file
            sequence_input = read_FASTA(self.file_path)
            
            # Update status to RUNNING immediately
            update_dto_callback(
                TaskDTO(
                    status=TaskStatus.RUNNING,
                    embedding_progress=EmbeddingProgress(current=0, total=len(sequence_input))
                )
            )

            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                input_fasta = tmp_path / "clustering_input.fasta"
                output_prefix = str(tmp_path / "cluster_out")
                mmseqs_tmp = str(tmp_path / "mmseqs_tmp")

                # 3. Write input sequences using seq_record.seq_id to match fasta.py
                id_map = {}
                with open(input_fasta, "w") as f:
                    for seq_record in sequence_input:
                        record_id = seq_record.seq_id  # <-- Aligned to his exact attribute
                        id_map[record_id] = seq_record.seq
                        f.write(f">{record_id}\n{str(seq_record.seq)}\n")

                # 4. Execute the corresponding mmseqs algorithm
                if self.use_linear_clustering:
                    cluster_results = easy_linclust(
                        str(input_fasta), output_prefix, mmseqs_tmp,
                        min_seq_id=self.sequence_identity_threshold
                    )
                else:
                    cluster_results = easy_cluster(
                        str(input_fasta), output_prefix, mmseqs_tmp,
                        min_seq_id=self.sequence_identity_threshold
                    )

                # 5. Parse out the unique cluster representatives
                clustered_results_dict = {}
                for cluster in cluster_results.to_gen():
                    rep_id = cluster.get("rep")
                    if rep_id and rep_id in id_map:
                        clustered_results_dict[rep_id] = id_map[rep_id]

            # 6. Push final metrics update
            update_dto_callback(
                TaskDTO(
                    status=TaskStatus.RUNNING,
                    embedding_progress=EmbeddingProgress(
                        current=len(clustered_results_dict), 
                        total=len(sequence_input)
                    )
                )
            )

            # 7. Return finished DTO
            return TaskDTO(status=TaskStatus.FINISHED, clustered_data=clustered_results_dict)

        except Exception as e:
            return TaskDTO.errored(f"MMseqs2 clustering failed: {str(e)}")