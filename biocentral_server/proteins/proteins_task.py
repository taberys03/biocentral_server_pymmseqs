import os
import tempfile
from typing import Callable, Dict


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
        # Status auf RUNNING setzen
        update_dto_callback(TaskDTO(status=TaskStatus.RUNNING))

        temp_input_path = None
        temp_output_prefix = None

        try:
            # 1. Dictionary in temporäre FASTA-Datei schreiben
            with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta", delete=False) as temp_input:
                for seq_id, seq in self.sequence_data.items():
                    temp_input.write(f">{seq_id}\n{seq}\n")
                temp_input_path = temp_input.name

            # 2. Ein flüchtiges Präfix für die MMseqs2-Ausgabedateien im Standard-Temp-Pfad erzeugen
            # Wir erstellen keinen Ordner, sondern nutzen einen eindeutigen Datei-Präfix
            temp_output_prefix = os.path.join(tempfile.gettempdir(), f"mmseqs_out_{os.getpid()}")

            print("Running pymmseqs command from temporary FASTA file...")

            # 3. MMseqs2-Clustering ausführen
            if self.use_linear_clustering:
                cluster_results = easy_linclust(
                    temp_input_path,
                    temp_output_prefix,
                    tempfile.gettempdir(),
                    min_seq_id=self.sequence_identity_threshold,
                )
            else:
                cluster_results = easy_cluster(
                    temp_input_path,
                    temp_output_prefix,
                    tempfile.gettempdir(),
                    min_seq_id=self.sequence_identity_threshold,
                )

            # 4. Resultate direkt über pymmseqs parsen (Kein manuelles Suchen nach Dateien nötig!)
            clustered_results = {}
            for cluster in cluster_results.to_gen():
                rep_id = cluster.get("rep")
                # Wenn der Representative-Key in unseren Ursprungsdaten liegt, übernehmen wir ihn
                if rep_id and rep_id in self.sequence_data:
                    clustered_results[rep_id] = self.sequence_data[rep_id]

            # Rückgabe der finalen repräsentativen Sequenzen
            return TaskDTO(status=TaskStatus.FINISHED, clustered_data=clustered_results)

        except Exception as e:
            print(f"Error during clustering execution: {str(e)}")
            return TaskDTO.errored(f"Clustering failed: {str(e)}")

        finally:
            # 5. AUFRÄUMEN: Temporäre Eingabedatei löschen
            if temp_input_path and os.path.exists(temp_input_path):
                try:
                    os.unlink(temp_input_path)
                except Exception:
                    pass

            # Generierte MMseqs-Dateien (z.B. _rep_seq.fasta, _all_seqs.fasta, etc.) aufräumen
            if temp_output_prefix:
                for suffix in ["_rep_seq.fasta", "_all_seqs.fasta", "_cluster.tsv", ""]:
                    file_to_delete = f"{temp_output_prefix}{suffix}"
                    if os.path.exists(file_to_delete):
                        try:
                            if os.path.isdir(file_to_delete):
                                import shutil
                                shutil.rmtree(file_to_delete)
                            else:
                                os.unlink(file_to_delete)
                        except Exception:
                            pass