import os
from pathlib import Path

# HMM Settings
DEFAULT_HMM_FILE = "PF00816.hmm"
DEFAULT_BIT_SCORE_THRESHOLD = 25.0

# Shared Apptainer SIF files
SHARED_APPTAINER_DIR = Path("/project/GROUP-MOLEPI/apptainer")
PLATON_SIF_NAME = "platon.sif"
CD_HIT_SIF_NAME = "cd-hit.sif"

# Fallback Docker URIs
FALLBACK_PLATON_URI = "docker://quay.io/biocontainers/platon:1.6--pyhdfd78af_1"
FALLBACK_CD_HIT_URI = "docker://quay.io/biocontainers/cd-hit:4.8.1--h5ca1c30_13"

# Platon Settings
PLATON_DB_PATH = Path("/project/GROUP-MOLEPI/databases/platon_db/db")
