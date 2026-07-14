import os
import sys
import concurrent.futures
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyrodigal
import pyhmmer
from Bio import Align

from .base import PipelineStep

def process_genome_worker(sample_name: str, fasta_path_str: str, hmm_path_str: str, bit_score_threshold: float, out_dir: str, hns_ref: str, stpa_ref: str) -> Tuple[Optional[str], List[Dict], List[Dict], bool]:
    """
    Worker function to process a single genome assembly.
    Runs inside a child process.
    """
    import json
    import hashlib
    
    done_file = Path(out_dir) / f"{sample_name}.done"
    hits_json = Path(out_dir) / f"{sample_name}_hits.json"
    
    if done_file.exists():
        fasta_hits = []
        tsv_hits = []
        if hits_json.exists():
            try:
                with open(hits_json, "r") as f:
                    data = json.load(f)
                    fasta_hits = data.get("fasta_hits", [])
                    tsv_hits = data.get("tsv_hits", [])
            except Exception:
                pass
        return None, fasta_hits, tsv_hits, True
        
    try:
        fasta_path = Path(fasta_path_str)
        hmm_path = Path(hmm_path_str)
        fasta_hits = []
        tsv_hits = []
        
        # Load HMM
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            hmm = hmm_file.read()
            if hmm is None:
                return f"No HMM profiles found in {hmm_path_str}", [], []
            hmms = [hmm]
            
        # Initialize ORF finder
        orf_finder = pyrodigal.GeneFinder(meta=True)
        
        # Read the assembly FASTA file in memory and keep contig records
        proteins = []
        seq_dict = {}
        contig_dict = {}
        
        with open(fasta_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                contig_dict[record.id] = record
                genes = orf_finder.find_genes(str(record.seq))
                for i, gene in enumerate(genes):
                    strand_char = "+" if gene.strand == 1 else "-"
                    seq_id = f"{sample_name}||{record.id}||{gene.begin}||{gene.end}||{strand_char}||{i+1}"
                    aa_seq = gene.translate()
                    
                    # Store in dictionary to look up original sequence later
                    seq_dict[seq_id] = aa_seq
                    
                    # Digitized sequence for pyhmmer
                    pyhmmer_seq = pyhmmer.easel.TextSequence(name=seq_id.encode(), sequence=aa_seq).digitize(hmm.alphabet)
                    proteins.append(pyhmmer_seq)
        
        if not proteins:
            return None, [], [] # No ORFs predicted
            
        # Scan predicted proteins against H-NS HMM profile using 1 CPU core
        seq_block = pyhmmer.easel.DigitalSequenceBlock(hmm.alphabet, proteins)
        
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        
        for hits in pyhmmer.hmmsearch(hmms, seq_block, cpus=1):
            for hit in hits:
                if hit.score >= bit_score_threshold:
                    seq_id = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                    parts = seq_id.split("||")
                    if len(parts) >= 5:
                        g_id, c_id, start, stop, strand = parts[0], parts[1], parts[2], parts[3], parts[4]
                    else:
                        g_id, c_id, start, stop, strand = sample_name, seq_id, "unknown", "unknown", "unknown"
                        
                    aa_seq = seq_dict.get(seq_id, "")
                    
                    # Paralog Classification
                    hns_score, stpa_score = 0, 0
                    if hns_ref and aa_seq:
                        hns_score = aligner.align(aa_seq.replace('*', ''), hns_ref.replace('*', ''))[0].score
                    if stpa_ref and aa_seq:
                        stpa_score = aligner.align(aa_seq.replace('*', ''), stpa_ref.replace('*', ''))[0].score
                        
                    predicted = "H-NS" if hns_score >= stpa_score else "StpA"
                    
                    hns_cov = (hns_score / len(hns_ref.replace('*', ''))) * 100 if hns_ref else 0
                    stpa_cov = (stpa_score / len(stpa_ref.replace('*', ''))) * 100 if stpa_ref else 0
                    
                    fasta_hits.append({
                        "id": seq_id,
                        "seq": aa_seq,
                        "desc": f"Score={hit.score} Type={predicted}"
                    })
                    
                    tsv_hits.append({
                        "Genome_ID": g_id,
                        "Contig_ID": c_id,
                        "Start": start,
                        "Stop": stop,
                        "Strand": strand,
                        "Bit_Score": f"{hit.score:.1f}",
                        "Sequence_Length": str(len(aa_seq)),
                        "Predicted_protein": predicted,
                        "StpA_cov": f"{stpa_cov:.1f}",
                        "StpA_identity": f"{stpa_cov:.1f}",
                        "HNS_cov": f"{hns_cov:.1f}",
                        "HNS_identity": f"{hns_cov:.1f}"
                    })
        
        # If H-NS hits are detected for this genome, write separate output files
        if fasta_hits:
            os.makedirs(out_dir, exist_ok=True)
            
            # A. Write individual amino acid sequences (.faa)
            faa_path = Path(out_dir) / f"{sample_name}.faa"
            records = [
                SeqRecord(Seq(hit["seq"]), id=hit["id"], description=hit["desc"])
                for hit in fasta_hits
            ]
            SeqIO.write(records, faa_path, "fasta")
            
            # B. Write nucleotide contigs carrying detected H-NS (.contigs.fasta)
            contigs_path = Path(out_dir) / f"{sample_name}.contigs.fasta"
            contigs_with_hits = {hit["id"].split("||")[1] for hit in fasta_hits}
            matching_contig_records = [
                contig_dict[cid] for cid in sorted(contigs_with_hits) if cid in contig_dict
            ]
            SeqIO.write(matching_contig_records, contigs_path, "fasta")
            
        # Save hits to cache
        os.makedirs(out_dir, exist_ok=True)
        hits_data = {
            "fasta_hits": fasta_hits,
            "tsv_hits": tsv_hits
        }
        with open(hits_json, "w") as f:
            json.dump(hits_data, f)
            
        # Write done file
        done_file.touch()
        
        return None, fasta_hits, tsv_hits, False
        
    except Exception as e:
        import traceback
        err_msg = f"Error: {e}\n{traceback.format_exc()}"
        return err_msg, [], [], False

class HNSScanner(PipelineStep):
    """
    Orchestrates process-level parallel scanning of multiple genome assemblies
    using pyrodigal and pyhmmer.
    """
    def run(self, unique_genomes: List[Tuple[str, str]], hmm_path: Path, bit_score_threshold: float, cpus: int, output_dir: str, hns_ref: str = "", stpa_ref: str = "") -> Tuple[List[Dict], List[Dict], Dict[str, str]]:
        all_fasta_hits = []
        all_tsv_hits = []
        errors = {}
        
        from tqdm import tqdm
        import hashlib
        
        print(f"Starting parallel scan on {len(unique_genomes)} genomes using {cpus} workers...")
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
            future_to_sample = {
                executor.submit(process_genome_worker, sample_name, fasta_path, str(hmm_path), bit_score_threshold, output_dir, hns_ref, stpa_ref): sample_name
                for sample_name, fasta_path in unique_genomes
            }
            
            # Print initial submitted logs for non-cached files immediately
            # Actually, we can log when they are completed/cached
            with tqdm(total=len(unique_genomes), desc="Scanning Genomes", unit="genome") as pbar:
                for future in concurrent.futures.as_completed(future_to_sample):
                    sample_name = future_to_sample[future]
                    pbar.update(1)
                    try:
                        err_msg, fasta_hits, tsv_hits, is_cached = future.result()
                        
                        task_hash = hashlib.md5(sample_name.encode()).hexdigest()
                        prefix = f"[{task_hash[:2]}/{task_hash[2:8]}]"
                        
                        if err_msg:
                            tqdm.write(f"{prefix} Error processing sample {sample_name}: {err_msg}")
                            errors[sample_name] = err_msg
                        else:
                            all_fasta_hits.extend(fasta_hits)
                            all_tsv_hits.extend(tsv_hits)
                            if is_cached:
                                tqdm.write(f"{prefix} Cached process > ScanGenome ({sample_name})")
                            else:
                                tqdm.write(f"{prefix} Completed process > ScanGenome ({sample_name})")
                    except Exception as e:
                        tqdm.write(f"Exception processing sample {sample_name}: {e}")
                        errors[sample_name] = str(e)
                        
        return all_fasta_hits, all_tsv_hits, errors
