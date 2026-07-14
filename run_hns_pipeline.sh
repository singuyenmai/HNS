#!/bin/bash
#SBATCH --partition=big
#SBATCH --nodelist=s0063
#SBATCH --cpus-per-task=14
#SBATCH --mem=54G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=hns_pipeline
#SBATCH --output=hns_pipeline_%j.log

# Exit immediately if a command exits with a non-zero status
set -e

###################
SUMMARY_CSV="/home/nguyenmts/slurm_history.csv"
NOW=$(date +%Y%m%d_%H%M)
RUN_LOG="/home/nguyenmts/slurm_logs/${NOW}_${SLURM_JOB_ID}.log"
# Everything from here on goes to log
exec > "$RUN_LOG" 2>&1

echo "=== H-NS Extraction Pipeline Slurm Job Started ==="
echo "Node: ${SLURMD_NODENAME}"
echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
echo "JobID: ${SLURM_JOB_ID}"
###################


# 1. Dependency Bootstrapping
echo "Loading python-base/3.11.3 module..."
module load python-base/3.11.3

# Create local virtual environment if not already present
if [ ! -d "venv" ]; then
    echo "Creating virtual environment in venv/..."
    python3 -m venv venv
fi

echo "Activating virtual environment..."
source venv/bin/activate

# Install uv inside the virtual environment if not already installed
if ! command -v uv &> /dev/null; then
    echo "Installing uv package manager in virtualenv..."
    pip install uv
fi

echo "Ensuring required packages are installed using uv..."
uv pip install biopython pyrodigal pyhmmer pandas tqdm

# 2. Dynamic Resource Mapping
# Determine CPU threads
THREADS="${SLURM_CPUS_PER_TASK:-4}"
echo "Using THREADS: $THREADS (derived from SLURM_CPUS_PER_TASK)"

# Determine memory for cd-hit in Megabytes (default: 4000 MB)
MEM_MB=4000
if [ -n "$SLURM_MEM_PER_NODE" ]; then
    # Parse memory from SLURM_MEM_PER_NODE (e.g., 12G, 12288M, or just a number in MB)
    MEM_VAL=$(echo "$SLURM_MEM_PER_NODE" | sed 's/[Mm]//g')
    if [[ "$MEM_VAL" =~ [Gg]$ ]]; then
        MEM_VAL=$(echo "$MEM_VAL" | sed 's/[Gg]//g')
        MEM_MB=$((MEM_VAL * 1024))
    else
        MEM_MB=$MEM_VAL
    fi
    echo "Derived MEM_MB: $MEM_MB (from SLURM_MEM_PER_NODE)"
elif [ -n "$SLURM_MEM_PER_CPU" ] && [ -n "$SLURM_CPUS_PER_TASK" ]; then
    MEM_CPU_VAL=$(echo "$SLURM_MEM_PER_CPU" | sed 's/[Mm]//g')
    if [[ "$MEM_CPU_VAL" =~ [Gg]$ ]]; then
        MEM_CPU_VAL=$(echo "$MEM_CPU_VAL" | sed 's/[Gg]//g')
        MEM_CPU_VAL=$((MEM_CPU_VAL * 1024))
    fi
    MEM_MB=$((MEM_CPU_VAL * SLURM_CPUS_PER_TASK))
    echo "Derived MEM_MB: $MEM_MB (from SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK)"
else
    echo "Using default memory limit MEM_MB: $MEM_MB"
fi

# Inputs and outputs can be parameterized, defaulting to test values if not specified
INPUT_PATH="${1:-test_input/test_samples.tsv}"
OUTPUT_TSV="${2:-all_hns_variants_report.tsv}"
OUTPUT_DIR="${3:-extracted_hns}"

echo "Input Path: $INPUT_PATH"
echo "Output TSV: $OUTPUT_TSV"
echo "Output Directory: $OUTPUT_DIR"

# 3. Execution Sequence

HMM_PROFILE="PF00816.hmm"

# Ensure Apptainer temp, cache, and output directories exist
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-/project/GROUP-MOLEPI/apptainer/nguyenmts_tmp}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-/project/GROUP-MOLEPI/apptainer/nguyenmts_cache}"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

# Execute hns_extractor.py (which manages scanning, reporting, and CD-HIT clustering)
if [ -f "$INPUT_PATH" ]; then
    echo "Running hns_extractor.py with file list: $INPUT_PATH"
    python3 src/hns_extractor.py \
        --hmm "$HMM_PROFILE" \
        --input-list "$INPUT_PATH" \
        --output-tsv "$OUTPUT_TSV" \
        --output-dir "$OUTPUT_DIR" \
        --output-unique \
        --mem-limit "$MEM_MB" \
        --cpus "$THREADS"
else
    echo "Running hns_extractor.py with input directory: $INPUT_PATH"
    python3 src/hns_extractor.py \
        --hmm "$HMM_PROFILE" \
        --input-dir "$INPUT_PATH" \
        --output-tsv "$OUTPUT_TSV" \
        --output-dir "$OUTPUT_DIR" \
        --output-unique \
        --mem-limit "$MEM_MB" \
        --cpus "$THREADS"
fi

echo "=== H-NS Extraction Pipeline Completed Successfully ==="
