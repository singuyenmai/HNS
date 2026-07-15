import os
import sys
import subprocess
import pandas as pd
from pathlib import Path
from typing import List, Dict, Tuple, Optional

from .base import PipelineStep
from constants import (
    PLATON_DB_PATH, 
    SHARED_APPTAINER_DIR, 
    PLATON_SIF_NAME, 
    FALLBACK_PLATON_URI,
    MLPLASMIDS_SIF_NAME,
    FALLBACK_MLPLASMIDS_URI,
    GENOMAD_SIF_NAME,
    FALLBACK_GENOMAD_URI,
    GENOMAD_DB_PATH
)

class LocationClassifier(PipelineStep):
    """
    Executes genomic location origin classification (Plasmid vs Chromosome) 
    using mlplasmids (default) or Platon.
    """
    def run(
        self, 
        tsv_hits: List[Dict], 
        unique_genomes: List[Tuple[str, str]], 
        output_dir: str, 
        cpus: int,
        classifier_engine: str = "mlplasmids",
        species: str = "Klebsiella pneumoniae"
    ) -> List[Dict]:
        print(f"Starting Contig Origin Classification ({classifier_engine}) on {len(unique_genomes)} genomes...")
        
        samples_with_hits = set()
        for hit in tsv_hits:
            samples_with_hits.add(hit["Genome_ID"])
            
        out_path = Path(output_dir).resolve()
        
        # Resolve Apptainer image based on classifier engine
        if classifier_engine == "mlplasmids":
            shared_sif = SHARED_APPTAINER_DIR / MLPLASMIDS_SIF_NAME
            project_dir = Path(__file__).resolve().parent.parent.parent
            local_sif = project_dir / MLPLASMIDS_SIF_NAME
            if shared_sif.exists():
                image_uri = str(shared_sif)
            elif local_sif.exists():
                image_uri = str(local_sif)
            else:
                image_uri = FALLBACK_MLPLASMIDS_URI
        elif classifier_engine == "genomad":
            db_path = GENOMAD_DB_PATH.resolve()
            shared_sif = SHARED_APPTAINER_DIR / GENOMAD_SIF_NAME
            project_dir = Path(__file__).resolve().parent.parent.parent
            local_sif = project_dir / GENOMAD_SIF_NAME
            if shared_sif.exists():
                image_uri = str(shared_sif)
            elif local_sif.exists():
                image_uri = str(local_sif)
            else:
                image_uri = FALLBACK_GENOMAD_URI
        else:
            db_path = PLATON_DB_PATH.resolve()
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
        
        # Platon classification function
        def run_platon_for_sample(sample_name: str) -> Tuple[Optional[str], bool]:
            contigs_fasta = out_path / f"{sample_name}.contigs.fasta"
            platon_out_dir = out_path / f"{sample_name}_platon"
            platon_tsv = platon_out_dir / f"{sample_name}.contigs.tsv"
            
            if platon_tsv.exists() and platon_tsv.stat().st_size > 0:
                return None, True
                
            if not contigs_fasta.exists():
                return None, False
                
            os.makedirs(platon_out_dir, exist_ok=True)
            
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

        # geNomad classification function
        def run_genomad_for_sample(sample_name: str) -> Tuple[Optional[str], bool]:
            contigs_fasta = out_path / f"{sample_name}.contigs.fasta"
            genomad_out_dir = out_path / f"{sample_name}_genomad"
            
            plasmid_summary = genomad_out_dir / f"{sample_name}.contigs_summary" / f"{sample_name}.contigs_plasmid_summary.tsv"
            virus_summary = genomad_out_dir / f"{sample_name}.contigs_summary" / f"{sample_name}.contigs_virus_summary.tsv"
            
            if (plasmid_summary.exists() and plasmid_summary.stat().st_size > 0) or (virus_summary.exists() and virus_summary.stat().st_size > 0):
                return None, True
                
            if not contigs_fasta.exists():
                return None, False
                
            os.makedirs(genomad_out_dir, exist_ok=True)
            
            cmd = [
                "apptainer", "exec",
                "--bind", f"{db_path}:/db",
                "--bind", f"{out_path}:{out_path}",
                image_uri,
                "genomad", "end-to-end",
                "--cleanup",
                "--splits", "1",
                str(contigs_fasta),
                str(genomad_out_dir),
                "/db"
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                return None, False
            except subprocess.CalledProcessError as e:
                return f"Error running geNomad for {sample_name}: {e.stderr}", False

        # mlplasmids classification function
        def run_mlplasmids_for_sample(sample_name: str) -> Tuple[Optional[str], bool]:
            contigs_fasta = out_path / f"{sample_name}.contigs.fasta"
            ml_out_dir = out_path / f"{sample_name}_mlplasmids"
            ml_tsv = ml_out_dir / f"{sample_name}_mlplasmids.tsv"
            
            if ml_tsv.exists() and ml_tsv.stat().st_size > 0:
                return None, True
                
            if not contigs_fasta.exists():
                return None, False
                
            os.makedirs(ml_out_dir, exist_ok=True)
            
            r_cmd = (
                f"library(mlplasmids); "
                f"res <- plasmid_classification(path_input_file = '{contigs_fasta}', species = '{species}'); "
                f"write.table(res, file = '{ml_tsv}', sep = '\\t', quote = FALSE, row.names = FALSE)"
            )
            
            cmd = [
                "apptainer", "exec",
                "--bind", f"{out_path}:{out_path}",
                image_uri,
                "Rscript", "-e", r_cmd
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                return None, False
            except subprocess.CalledProcessError as e:
                return f"Error running mlplasmids for {sample_name}: {e.stderr}", False

        # Cap parallel workers based on available Slurm memory to prevent OOM
        slurm_mem = os.getenv("SLURM_MEM_PER_NODE")
        if slurm_mem:
            try:
                available_mem_mb = int(slurm_mem)
                max_mem_workers = max(1, available_mem_mb // 2048)
                if max_mem_workers < cpus:
                    print(f"Capping {classifier_engine} parallel workers from {cpus} to {max_mem_workers} due to Slurm memory limit ({slurm_mem} MB).")
                    cpus = max_mem_workers
            except ValueError:
                pass

        print(f"Running {classifier_engine} concurrently for genomes with hits using up to {cpus} workers...")
        from tqdm import tqdm
        
        if classifier_engine == "mlplasmids":
            target_fn = run_mlplasmids_for_sample
            classifier_name = "MlplasmidsClassifier"
        elif classifier_engine == "genomad":
            target_fn = run_genomad_for_sample
            classifier_name = "GeNomadClassifier"
        else:
            target_fn = run_platon_for_sample
            classifier_name = "PlatonClassifier"
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=cpus) as executor:
            futures = {
                executor.submit(target_fn, sample_name): sample_name
                for sample_name, _ in unique_genomes
                if sample_name in samples_with_hits
            }
            
            with tqdm(total=len(futures), desc=f"{classifier_engine.capitalize()} Classifier", unit="genome") as pbar:
                for future in concurrent.futures.as_completed(futures):
                    sample_name = futures[future]
                    pbar.update(1)
                    try:
                        err_msg, is_cached = future.result()
                        
                        task_hash = hashlib.md5(f"{sample_name}_{classifier_engine}".encode()).hexdigest()
                        prefix = f"[{task_hash[:2]}/{task_hash[2:8]}]"
                        
                        if err_msg:
                            tqdm.write(err_msg)
                        elif is_cached:
                            tqdm.write(f"{prefix} Cached process > {classifier_name} ({sample_name})")
                        else:
                            tqdm.write(f"{prefix} Completed process > {classifier_name} ({sample_name})")
                    except Exception as e:
                        tqdm.write(f"Exception running {classifier_engine} for {sample_name}: {e}")
                
        print(f"Parsing {classifier_engine} outputs and merging into report...")
        location_dict = {}
        virus_intervals = {} # Key: (sample_name, contig_id), Value: List of tuples (start, end) or "Entire"
        
        for sample_name in samples_with_hits:
            if classifier_engine == "mlplasmids":
                ml_tsv = out_path / f"{sample_name}_mlplasmids" / f"{sample_name}_mlplasmids.tsv"
                if ml_tsv.exists():
                    try:
                        df = pd.read_csv(ml_tsv, sep='\t')
                        for _, row in df.iterrows():
                            contig_id = str(row['Contig_name'])
                            prediction = str(row['Prediction']).capitalize() # e.g. "Plasmid" or "Chromosome"
                            location_dict[(sample_name, contig_id)] = prediction
                    except Exception as e:
                        print(f"Error reading mlplasmids output for {sample_name}: {e}", file=sys.stderr)
            elif classifier_engine == "genomad":
                genomad_dir = out_path / f"{sample_name}_genomad" / f"{sample_name}.contigs_summary"
                
                # Read plasmid summary
                plasmid_tsv = genomad_dir / f"{sample_name}.contigs_plasmid_summary.tsv"
                if plasmid_tsv.exists():
                    try:
                        df = pd.read_csv(plasmid_tsv, sep='\t')
                        if 'seq_name' in df.columns:
                            for _, row in df.iterrows():
                                contig_id = str(row['seq_name']).split('|')[0].split()[0]
                                location_dict[(sample_name, contig_id)] = "Plasmid"
                    except Exception as e:
                        print(f"Error reading geNomad plasmid summary for {sample_name}: {e}", file=sys.stderr)
                        
                # Read virus summary
                virus_tsv = genomad_dir / f"{sample_name}.contigs_virus_summary.tsv"
                if virus_tsv.exists():
                    try:
                        df = pd.read_csv(virus_tsv, sep='\t')
                        if 'seq_name' in df.columns:
                            for _, row in df.iterrows():
                                contig_id = str(row['seq_name']).split('|')[0].split()[0]
                                topology = str(row.get('topology', ''))
                                coordinates = str(row.get('coordinates', ''))
                                
                                key = (sample_name, contig_id)
                                if topology == "Provirus" and coordinates and coordinates != 'nan' and '-' in coordinates:
                                    try:
                                        start_end = coordinates.split('-')
                                        v_start = int(start_end[0])
                                        v_end = int(start_end[1])
                                        if key not in virus_intervals or virus_intervals[key] == "Entire":
                                            virus_intervals[key] = []
                                        virus_intervals[key].append((v_start, v_end))
                                    except Exception:
                                        virus_intervals[key] = "Entire"
                                else:
                                    virus_intervals[key] = "Entire"
                    except Exception as e:
                        print(f"Error reading geNomad virus summary for {sample_name}: {e}", file=sys.stderr)
            else:
                platon_dir = out_path / f"{sample_name}_platon"
                
                # Read chromosome FASTA contig IDs
                chrom_fasta = platon_dir / f"{sample_name}.contigs.chromosome.fasta"
                if chrom_fasta.exists():
                    try:
                        with open(chrom_fasta, "r") as f:
                            for line in f:
                                if line.startswith(">"):
                                    contig_id = line[1:].strip().split()[0]
                                    location_dict[(sample_name, contig_id)] = "Chromosome"
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
                                    location_dict[(sample_name, contig_id)] = "Plasmid"
                    except Exception as e:
                        print(f"Error reading Platon plasmid FASTA for {sample_name}: {e}", file=sys.stderr)
                        
        # Update tsv_hits with Location
        for hit in tsv_hits:
            sample_name = hit["Genome_ID"]
            contig_id = hit["Contig_ID"]
            key = (sample_name, contig_id)
            
            default_loc = "Chromosome" if classifier_engine == "genomad" else "Unknown"
            
            if classifier_engine == "genomad":
                if key in location_dict:
                    hit["Location"] = location_dict[key]
                elif key in virus_intervals:
                    intervals = virus_intervals[key]
                    if intervals == "Entire":
                        hit["Location"] = "Virus"
                    else:
                        try:
                            hit_start = int(hit["Start"])
                            hit_stop = int(hit["Stop"])
                            is_in_provirus = False
                            for v_start, v_end in intervals:
                                if max(hit_start, v_start) <= min(hit_stop, v_end):
                                    is_in_provirus = True
                                    break
                            hit["Location"] = "Virus" if is_in_provirus else "Chromosome"
                        except (ValueError, TypeError):
                            hit["Location"] = "Virus"
                else:
                    hit["Location"] = default_loc
            else:
                hit["Location"] = location_dict.get(key, default_loc)
            
        return tsv_hits
