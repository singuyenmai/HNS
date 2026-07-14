import os
import sys
import argparse
from pathlib import Path

# Add the current directory to sys.path to ensure modules package can be imported 
# when run directly as a script.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from modules.downloader import HMMDownloader
from modules.extractor import HNSScanner
from modules.reporter import MergedReporter
from modules.deduplicator import CDHitDeduplicator
from modules.classifier import PlatonClassifier
from Bio import SeqIO
from constants import DEFAULT_HMM_FILE, DEFAULT_BIT_SCORE_THRESHOLD

def main():
    parser = argparse.ArgumentParser(description="Modular H-NS Extraction Pipeline using pyrodigal and pyhmmer.")
    parser.add_argument("--hmm", "-m", type=str, default=DEFAULT_HMM_FILE, help=f"Path to H-NS HMM profile (default: {DEFAULT_HMM_FILE}).")
    parser.add_argument("--download-hmm", action="store_true", help="Download the H-NS HMM profile directly from InterPro API and exit.")
    parser.add_argument("--input-dir", "-i", type=str, help="Directory containing genome assembly FASTA files.")
    parser.add_argument("--input-list", "-l", type=str, help="Text file (simple list of paths, or TSV with sample_name and fasta_path).")

    parser.add_argument("--output-tsv", "-t", type=str, help="Path for the output Bakta-style TSV report.")
    parser.add_argument("--bit-score", "-s", type=float, default=DEFAULT_BIT_SCORE_THRESHOLD, help=f"HMMER bit-score threshold (default: {DEFAULT_BIT_SCORE_THRESHOLD}).")
    parser.add_argument("--cpus", "-c", type=int, default=os.cpu_count(), help="Number of parallel processes (default: available CPU cores).")
    parser.add_argument("--output-dir", "-d", type=str, default="extracted_hns", help="Output directory for individual genome .faa and .contigs.fasta files.")
    parser.add_argument("--output-unique", "-u", action="store_true", help="Enable CD-HIT unique variants clustering per sample.")
    parser.add_argument("--mem-limit", type=int, default=4000, help="CD-HIT memory limit in Megabytes (default: 4000).")
    parser.add_argument("fasta_files", nargs="*", type=str, help="List of explicit FASTA files to process.")
    
    args = parser.parse_args()
    
    # Step 1: Handle HMM downloader execution
    if args.download_hmm:
        downloader = HMMDownloader()
        downloader.run(args.hmm)
        sys.exit(0)
        
    if not args.output_tsv:
        parser.error("the following arguments are required: --output-tsv/-t")
        
    # Gather genomes to process as (sample_name, fasta_path)
    genomes = []
    
    # 1. Positional arguments
    for f in args.fasta_files:
        path = Path(f)
        genomes.append((path.stem, str(path)))
        
    # 2. Input directory
    if args.input_dir:
        input_path = Path(args.input_dir)
        extensions = {".fasta", ".fa", ".fna"}
        if input_path.is_dir():
            for f in sorted(input_path.iterdir()):
                if f.suffix in extensions:
                    genomes.append((f.stem, str(f)))
        else:
            print(f"Error: input directory {args.input_dir} does not exist.", file=sys.stderr)
            sys.exit(1)
            
    # 3. Input list file (can be simple list, or TSV format with sample_name\tfasta_path)
    if args.input_list:
        list_path = Path(args.input_list)
        if list_path.is_file():
            with open(list_path, "r") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        if "\t" in line:
                            parts = line.split("\t")
                            sample_name = parts[0].strip()
                            fasta_path = parts[1].strip()
                            genomes.append((sample_name, fasta_path))
                        else:
                            path = Path(line)
                            genomes.append((path.stem, line))
        else:
            print(f"Error: input list file {args.input_list} does not exist.", file=sys.stderr)
            sys.exit(1)
            
    # Remove duplicates while preserving order
    seen = set()
    unique_genomes = []
    for sample_name, fasta_path in genomes:
        if fasta_path not in seen:
            seen.add(fasta_path)
            unique_genomes.append((sample_name, fasta_path))
            
    if not unique_genomes:
        print("Error: No FASTA files or sample lists specified to process.", file=sys.stderr)
        sys.exit(1)
        
    # Step 2: Ensure HMM profile exists, download if missing
    downloader = HMMDownloader()
    hmm_path = downloader.run(args.hmm)
    
    # Load References
    def read_ref(filepath: str) -> str:
        try:
            with open(filepath, "r") as f:
                for record in SeqIO.parse(f, "fasta"):
                    return str(record.seq)
        except Exception:
            return ""
        return ""
        
    db_dir = Path(__file__).resolve().parent / "databases"
    hns_ref = read_ref(str(db_dir / "HNS.faa"))
    stpa_ref = read_ref(str(db_dir / "StpA.faa"))
    
    # Step 3: Run H-NS extraction step
    scanner = HNSScanner()
    fasta_hits, tsv_hits, errors = scanner.run(
        unique_genomes=unique_genomes,
        hmm_path=hmm_path,
        bit_score_threshold=args.bit_score,
        cpus=args.cpus,
        output_dir=args.output_dir,
        hns_ref=hns_ref,
        stpa_ref=stpa_ref
    )
    
    # Step 3.5: Run Platon Contig Classification
    classifier = PlatonClassifier()
    tsv_hits = classifier.run(
        tsv_hits=tsv_hits,
        unique_genomes=unique_genomes,
        output_dir=args.output_dir,
        cpus=args.cpus
    )
    
    # Step 4: Write merged/aggregated outputs
    reporter = MergedReporter()
    reporter.run(
        all_tsv_hits=tsv_hits,
        output_tsv=Path(args.output_tsv)
    )
    
    if args.output_unique:
        samples_with_hits = {hit["Genome_ID"] for hit in tsv_hits}
        if samples_with_hits:
            deduplicator = CDHitDeduplicator()
            for sample_name, _ in unique_genomes:
                if sample_name not in samples_with_hits:
                    continue
                sample_faa = Path(args.output_dir) / f"{sample_name}.faa"
                sample_unique_faa = Path(args.output_dir) / f"{sample_name}_unique.faa"
                if sample_faa.exists():
                    deduplicator.run(
                        input_fasta=sample_faa,
                        output_unique_fasta=sample_unique_faa,
                        threads=args.cpus,
                        mem_mb=args.mem_limit
                    )
        else:
            print("No hits found; skipping CD-HIT clustering.")
            
    if errors:
        print(f"\nPipeline completed with errors in {len(errors)}/{len(unique_genomes)} genomes.", file=sys.stderr)
        sys.exit(2)
    else:
        print("\nPipeline completed successfully!")
        sys.exit(0)

if __name__ == "__main__":
    main()
