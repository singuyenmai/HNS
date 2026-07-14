import os
import sys
import subprocess
import pandas as pd
from pathlib import Path
from typing import List, Dict, Tuple, Optional

from .base import PipelineStep
from constants import PLATON_DB_PATH, SHARED_APPTAINER_DIR, PLATON_SIF_NAME, FALLBACK_PLATON_URI

class PlatonClassifier(PipelineStep):
    """
    Executes Platon on contig FASTA files to predict origin (Plasmid vs Chromosome).
    """
    def run(self, tsv_hits: List[Dict], unique_genomes: List[Tuple[str, str]], output_dir: str, cpus: int) -> List[Dict]:
        print(f"Starting Platon Contig Origin Classification on {len(unique_genomes)} genomes...")
        
        samples_with_hits = set()
        for hit in tsv_hits:
            samples_with_hits.add(hit["Genome_ID"])
            
        out_path = Path(output_dir).resolve()
        db_path = PLATON_DB_PATH.resolve()
        
        # Platon SIF image location
        shared_sif = SHARED_APPTAINER_DIR / PLATON_SIF_NAME
        project_dir = Path(__file__).resolve().parent.parent.parent
        local_sif = project_dir / PLATON_SIF_NAME
        
        if shared_sif.exists():
            image_uri = str(shared_sif)
        elif local_sif.exists():
            image_uri = str(local_sif)
        else:
            image_uri = FALLBACK_PLATON_URI
            
        import concurrent.futures
        import hashlib
        
        def run_platon_for_sample(sample_name: str) -> Tuple[Optional[str], bool]:
            contigs_fasta = out_path / f"{sample_name}.contigs.fasta"
            platon_out_dir = out_path / f"{sample_name}_platon"
            platon_tsv = platon_out_dir / f"{sample_name}.contigs.tsv"
            
            # Check if Platon output already exists and is non-empty
            if platon_tsv.exists() and platon_tsv.stat().st_size > 0:
                return None, True
                
            if not contigs_fasta.exists():
                return None, False
                
            os.makedirs(platon_out_dir, exist_ok=True)
            
            # Execute Platon via Apptainer (using 1 thread per job for parallel efficiency)
            cmd = [
                "apptainer", "exec",
                "--bind", f"{db_path}:/db",
                "--bind", f"{out_path}:{out_path}",
                image_uri,
                "platon",
                "--db", "/db",
                "-t", "1",
                "-o", str(platon_out_dir),
                str(contigs_fasta)
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                return None, False
            except subprocess.CalledProcessError as e:
                return f"Error running Platon for {sample_name}: {e.stderr}", False

        # Cap Platon workers based on available Slurm memory to prevent OOM
        # Assuming each Platon instance requires ~2048 MB (2 GB)
        slurm_mem = os.getenv("SLURM_MEM_PER_NODE")
        if slurm_mem:
            try:
                available_mem_mb = int(slurm_mem)
                max_mem_workers = max(1, available_mem_mb // 2048)
                if max_mem_workers < cpus:
                    print(f"Capping Platon parallel workers from {cpus} to {max_mem_workers} due to Slurm memory limit ({slurm_mem} MB).")
                    cpus = max_mem_workers
            except ValueError:
                pass

        # Run Platon in parallel using ThreadPoolExecutor
        print(f"Running Platon concurrently for genomes with hits using up to {cpus} workers...")
        from tqdm import tqdm
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=cpus) as executor:
            futures = {
                executor.submit(run_platon_for_sample, sample_name): sample_name
                for sample_name, _ in unique_genomes
                if sample_name in samples_with_hits
            }
            
            with tqdm(total=len(futures), desc="Platon Classifier", unit="genome") as pbar:
                for future in concurrent.futures.as_completed(futures):
                    sample_name = futures[future]
                    pbar.update(1)
                    try:
                        err_msg, is_cached = future.result()
                        
                        task_hash = hashlib.md5(f"{sample_name}_platon".encode()).hexdigest()
                        prefix = f"[{task_hash[:2]}/{task_hash[2:8]}]"
                        
                        if err_msg:
                            tqdm.write(err_msg)
                        elif is_cached:
                            tqdm.write(f"{prefix} Cached process > PlatonClassifier ({sample_name})")
                        else:
                            tqdm.write(f"{prefix} Completed process > PlatonClassifier ({sample_name})")
                    except Exception as e:
                        tqdm.write(f"Exception running Platon for {sample_name}: {e}")
                
        print("Parsing Platon outputs and merging into report...")
        platon_dict = {}
        
        for sample_name in samples_with_hits:
            platon_dir = out_path / f"{sample_name}_platon"
            
            # Read chromosome FASTA contig IDs
            chrom_fasta = platon_dir / f"{sample_name}.contigs.chromosome.fasta"
            if chrom_fasta.exists():
                try:
                    with open(chrom_fasta, "r") as f:
                        for line in f:
                            if line.startswith(">"):
                                contig_id = line[1:].strip().split()[0]
                                platon_dict[contig_id] = "Chromosome"
                except Exception as e:
                    print(f"Error reading Platon chromosome FASTA for {sample_name}: {e}", file=sys.stderr)
                    
            # Read plasmid FASTA contig IDs
            plasmid_fasta = platon_dir / f"{sample_name}.contigs.plasmid.fasta"
            if plasmid_fasta.exists():
                try:
                    with open(plasmid_fasta, "r") as f:
                        for line in f:
                            if line.startswith(">"):
                                contig_id = line[1:].strip().split()[0]
                                platon_dict[contig_id] = "Plasmid"
                except Exception as e:
                    print(f"Error reading Platon plasmid FASTA for {sample_name}: {e}", file=sys.stderr)
                    
        # Update tsv_hits with Location
        for hit in tsv_hits:
            contig_id = hit["Contig_ID"]
            hit["Location"] = platon_dict.get(contig_id, "Unknown")
            
        return tsv_hits
