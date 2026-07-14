import sys
from pathlib import Path
from typing import List, Dict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .base import PipelineStep

class MergedReporter(PipelineStep):
    """
    Compiles and writes merged Bakta-style TSV report files.
    """
    def run(self, all_tsv_hits: List[Dict], output_tsv: Path) -> bool:
        # Create output directories if needed
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
            
        # Write merged TSV
        try:
            with open(output_tsv, "w") as tsv_file:
                tsv_file.write("Genome_ID\tContig_ID\tStart\tStop\tStrand\tBit_Score\tSequence_Length\tPredicted_protein\tLocation\tStpA_cov\tStpA_identity\tHNS_cov\tHNS_identity\n")
                for hit in all_tsv_hits:
                    tsv_file.write(
                        f"{hit['Genome_ID']}\t{hit['Contig_ID']}\t{hit['Start']}\t{hit['Stop']}\t"
                        f"{hit['Strand']}\t{hit['Bit_Score']}\t{hit['Sequence_Length']}\t"
                        f"{hit.get('Predicted_protein', 'Unknown')}\t{hit.get('Location', 'Unknown')}\t"
                        f"{hit.get('StpA_cov', '0')}\t{hit.get('StpA_identity', '0')}\t"
                        f"{hit.get('HNS_cov', '0')}\t{hit.get('HNS_identity', '0')}\n"
                    )
            print(f"Saved merged Bakta-style TSV report to: {output_tsv}")
            return True
        except Exception as e:
            print(f"Error writing output TSV: {e}", file=sys.stderr)
            raise
