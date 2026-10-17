from pathlib import Path

import pandas as pd

from src.parsers.leafcutter_parser import parse_leafcutter
from src.parsers.majiq_parser import parse_majiq_voila
from src.parsers.rmats_parser import (
    parse_rmats_se,
    parse_rmats_a3ss_a5ss,
    parse_rmats_mxe,
    parse_rmats_ri,
)
from src.utils import merge_overlap_genomic

# Statistical and effect size thresholds
RMATS_FDR_THRESHOLD = 0.05       # rMATS q-value (FDR) cutoff
MAJIQ_PROB_THRESHOLD = 0.90      # MAJIQ Probability of Changing cutoff (>= 90% confidence)
MIN_DELTA_PSI = 0.10             # Minimum absolute magnitude change (10%)
MAX_BOUNDARY_DRIFT = 2           # Max bp difference at each junction end (tools differ by ~1 bp by convention)

# All paths are resolved relative to this file, so the script works from any working directory.
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

LEAFCUTTER_FILE = DATA_DIR / "leafcutter" / "leafcutter_ds_effect_sizes.csv"
MAJIQ_FILE = DATA_DIR / "majiq" / "tsv_f.csv"
RMATS_DIR = DATA_DIR / "rmats"
OUTPUT_FILE = OUTPUT_DIR / "multitool_streamlined_filtered.csv"
SUMMARY_FILE = OUTPUT_DIR / "multitool_summary.csv"

# All tools report DeltaPSI as Hypoxia - Normoxia after parsing.
# rMATS (group 1 = hypoxia) is natively Hyp - Norm; LeafCutter and MAJIQ are natively Norm - Hyp.

# 1. Load LeafCutter data
try:
    lc_data = parse_leafcutter(LEAFCUTTER_FILE, flip_sign=True)
    print(f"Successfully parsed {len(lc_data)} introns from LeafCutter.")
except FileNotFoundError:
    lc_data = pd.DataFrame()
    print(f"Skipping LeafCutter: File {LEAFCUTTER_FILE} not found.")

# 2. Load MAJIQ data
try:
    majiq_data = parse_majiq_voila(MAJIQ_FILE, flip_sign=True)
    print(f"Successfully parsed and mapped {len(majiq_data)} flat introns from MAJIQ.")
except FileNotFoundError:
    majiq_data = pd.DataFrame()
    print(f"Skipping MAJIQ: File {MAJIQ_FILE} not found.")

# 3. Load and parse all rMATS event files
event_files = {
    'SE': (RMATS_DIR / "SE.csv", parse_rmats_se),
    'A3SS': (RMATS_DIR / "A3SS.csv", lambda p: parse_rmats_a3ss_a5ss(p, 'A3SS')),
    'A5SS': (RMATS_DIR / "A5SS.csv", lambda p: parse_rmats_a3ss_a5ss(p, 'A5SS')),
    'RI': (RMATS_DIR / "RI.csv", parse_rmats_ri),
    'MXE': (RMATS_DIR / "MXE.csv", parse_rmats_mxe),
}

rmats_dfs = []
for event_type, (filename, parser_func) in event_files.items():
    try:
        df_parsed = parser_func(filename)
        rmats_dfs.append(df_parsed)
        print(f"Successfully parsed {len(df_parsed)} introns from rMATS {event_type}")
    except FileNotFoundError:
        print(f"Skipping rMATS {event_type}: File {filename} not found.")

if rmats_dfs:
    all_rmats_introns = pd.concat(rmats_dfs, ignore_index=True)
    print(f"Total rMATS Introns Across All Event Types: {len(all_rmats_introns)}")
    
    # 5. Multi-Tool Spatial Reconciliation (LeafCutter <-> rMATS <-> MAJIQ)
    overlap_matches = pd.DataFrame()
    
    if not lc_data.empty and not all_rmats_introns.empty:
        # Step A: Find overlap between LeafCutter and rMATS (using '_rMATS' suffix)
        lc_rmats_overlap = merge_overlap_genomic(lc_data, all_rmats_introns, max_boundary_drift=MAX_BOUNDARY_DRIFT, suffix='_rMATS')
        
        # Step B: Intersect the LeafCutter-rMATS consensus with MAJIQ (using '_MAJIQ' suffix)
        if not majiq_data.empty and not lc_rmats_overlap.empty:
            overlap_matches = merge_overlap_genomic(lc_rmats_overlap, majiq_data, max_boundary_drift=MAX_BOUNDARY_DRIFT, suffix='_MAJIQ')
        else:
            overlap_matches = lc_rmats_overlap

if not overlap_matches.empty:
    # Streamlined filter: rMATS FDR, MAJIQ Probability, and Delta PSI magnitudes across tools
    filtered_consensus = overlap_matches[
        # 1. rMATS Statistical Significance (FDR / q-value)
        (overlap_matches.get('FDR', overlap_matches.get('FDR_rMATS', 1.0)) < RMATS_FDR_THRESHOLD) & 
        
        # 2. MAJIQ Statistical Significance (Probability Changing)
        (overlap_matches.get('ProbabilityChanging', overlap_matches.get('ProbabilityChanging_MAJIQ', 0.0)) >= MAJIQ_PROB_THRESHOLD) &
        
        # 3. Effect Size Magnitudes (Checking Delta PSI / IncLevelDifference equivalents for all tools)
        (overlap_matches.get('IncLevelDifference', overlap_matches.get('deltapsi', 0.0)).abs() >= MIN_DELTA_PSI) &
        (overlap_matches.get('IncLevelDifference_rMATS', 0.0).abs() >= MIN_DELTA_PSI) &
        (overlap_matches.get('IncLevelDifference_MAJIQ', 0.0).abs() >= MIN_DELTA_PSI)
    ].copy()

    filtered_consensus.to_csv(OUTPUT_FILE, index=False)
    
    print(f"Successfully filtered down to {len(filtered_consensus)} high-confidence consensus events.")
    print(f"Saved results to '{OUTPUT_FILE}'")

    # Condensed view: coordinates, chromosome, gene, DeltaPSIs (Hyp - Norm) and significance per tool
    summary = pd.DataFrame({
        'Chromosome': filtered_consensus['Chromosome'],
        'Strand': filtered_consensus['Strand'],
        'Gene': filtered_consensus['geneSymbol'].fillna(filtered_consensus['geneSymbol_MAJIQ']),
        'EventType': filtered_consensus['EventType'],
        'LeafCutter_Start': filtered_consensus['Start'],
        'LeafCutter_End': filtered_consensus['End'],
        'rMATS_Start': filtered_consensus['Start_rMATS'],
        'rMATS_End': filtered_consensus['End_rMATS'],
        'MAJIQ_Start': filtered_consensus['Start_MAJIQ'],
        'MAJIQ_End': filtered_consensus['End_MAJIQ'],
        'LeafCutter_DeltaPSI': filtered_consensus['IncLevelDifference'],
        'rMATS_DeltaPSI': filtered_consensus['IncLevelDifference_rMATS'],
        'MAJIQ_DeltaPSI': filtered_consensus['IncLevelDifference_MAJIQ'],
        'rMATS_FDR': filtered_consensus['FDR'],
        'MAJIQ_ProbabilityChanging': filtered_consensus['ProbabilityChanging'],
    })
    summary.to_csv(SUMMARY_FILE, index=False)
    print(f"Saved condensed summary ({len(summary)} rows) to '{SUMMARY_FILE}'")
else:
    print("No overlapping consensus events found across the tools.")