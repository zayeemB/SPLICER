"""Project-wide settings: thresholds and file locations.

All paths are resolved from the project root, so nothing depends on the working directory.
"""
from pathlib import Path

# Statistical and effect size thresholds
RMATS_FDR_THRESHOLD = 0.05       # rMATS q-value (FDR) cutoff
MAJIQ_PROB_THRESHOLD = 0.90      # MAJIQ Probability of Changing cutoff (>= 90% confidence)
MIN_DELTA_PSI = 0.10             # Minimum absolute magnitude change (10%)
MAX_BOUNDARY_DRIFT = 2           # Max bp difference at each junction end (tools differ by ~1 bp by convention)

# Inputs
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
LEAFCUTTER_FILE = DATA_DIR / "leafcutter" / "leafcutter_ds_effect_sizes.csv"
MAJIQ_FILE = DATA_DIR / "majiq" / "tsv_f.csv"
RMATS_DIR = DATA_DIR / "rmats"

# Outputs: tables and reports go in datasets/, figures in plots/
OUTPUT_DIR = BASE_DIR / "output"
DATASETS_DIR = OUTPUT_DIR / "datasets"
PLOTS_DIR = OUTPUT_DIR / "plots"
