import subprocess
import sys
from pathlib import Path
from .base import PipelineStep
from constants import SHARED_APPTAINER_DIR, CD_HIT_SIF_NAME, FALLBACK_CD_HIT_URI

class CDHitDeduplicator(PipelineStep):
    """
    Runs CD-HIT via Apptainer to cluster proteins at 100% identity.
    """
    def run(self, input_fasta: Path, output_unique_fasta: Path, threads: int = 4, mem_mb: int = 4000) -> bool:
        # Ensure the output directory exists
        output_unique_fasta.parent.mkdir(parents=True, exist_ok=True)
        
        print("Running CD-HIT dereplication via Apptainer...")
        
        shared_sif = SHARED_APPTAINER_DIR / CD_HIT_SIF_NAME
        project_dir = Path(__file__).resolve().parent.parent.parent
        local_sif = project_dir / CD_HIT_SIF_NAME
        
        if shared_sif.exists():
            image_uri = str(shared_sif)
        elif local_sif.exists():
            image_uri = str(local_sif)
        else:
            image_uri = FALLBACK_CD_HIT_URI
            
        cmd = [
            "apptainer", "exec", 
            image_uri, 
            "cd-hit",
            "-i", str(input_fasta),
            "-o", str(output_unique_fasta),
            "-c", "1.0",
            "-d", "0",
            "-T", str(threads),
            "-M", str(mem_mb)
        ]
        
        print(f"Executing command: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, check=True, text=True, capture_output=True)
            print(result.stdout)
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error running CD-HIT Apptainer container: {e}", file=sys.stderr)
            print(f"Stdout:\n{e.stdout}", file=sys.stderr)
            print(f"Stderr:\n{e.stderr}", file=sys.stderr)
            raise
