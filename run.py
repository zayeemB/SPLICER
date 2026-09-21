import pandas as pd
from src.leafcutter_parser import parse_leafcutter
from src.rmats_parser import (
    parse_rmats_se, 
    parse_rmats_a3ss_a5ss, 
    parse_rmats_mxe, 
    parse_rmats_ri
)
from src.majiq_parser import parse_majiq_voila
from src.utils import merge_overlap_genomic, load_gene_chromosome_map

# Define streamlined statistical and effect size thresholds
RMATS_FDR_THRESHOLD = 0.05       # rMATS q-value (FDR) cutoff
MAJIQ_PROB_THRESHOLD = 0.90      # MAJIQ Probability of Changing cutoff (>= 90% confidence)
MIN_DELTA_PSI = 0.10             # Minimum absolute magnitude change (10%)

DATA_DIR = "./data"
OUTPUT_DIR = "./output"

gtf_file = "gencode.v49.annotation.gtf" 
majiq_file = f"{DATA_DIR}/majiq/tsv_f.csv"
leafcutter_file = f"{DATA_DIR}/leafcutter/leafcutter_ds_effect_sizes.csv"

# 1. Build the GTF lookup map as a fallback (used if MAJIQ lacks a native chromosome column)
try:
    gene_chrom_map = load_gene_chromosome_map(gtf_file)
except FileNotFoundError:
    gene_chrom_map = None
    print(f"Warning: GTF file {gtf_file} not found. MAJIQ files without native chromosome columns will be skipped.")

# 2. Load LeafCutter data
try:
    lc_data = parse_leafcutter(leafcutter_file)
    print(f"Successfully parsed {len(lc_data)} introns from LeafCutter.")
except FileNotFoundError:
    lc_data = pd.DataFrame()
    print(f"Skipping LeafCutter: File {leafcutter_file} not found.")

# 3. Load and parse MAJIQ data (auto-checks for native chromosomes or uses GTF map)
try:
    majiq_data = parse_majiq_voila(majiq_file, gene_to_chrom_map=gene_chrom_map)
    print(f"Successfully parsed and mapped {len(majiq_data)} flat introns from MAJIQ.")
except FileNotFoundError:
    majiq_data = pd.DataFrame()
    print(f"Skipping MAJIQ: File {majiq_file} not found.")

# 4. Load and parse all rMATS event files
rmats_dfs = []
event_files = {
    'SE': (f'{DATA_DIR}/rmats/SE.csv', parse_rmats_se),
    'A3SS': (f'{DATA_DIR}/rmats/A3SS.csv', lambda p: parse_rmats_a3ss_a5ss(p, 'A3SS')),
    'A5SS': (f'{DATA_DIR}/rmats/A5SS.csv', lambda p: parse_rmats_a3ss_a5ss(p, 'A5SS')),
    'RI': (f'{DATA_DIR}/rmats/RI.csv', parse_rmats_ri),
    'MXE': (f'{DATA_DIR}/rmats/MXE.csv', parse_rmats_mxe)
}

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
        lc_rmats_overlap = merge_overlap_genomic(lc_data, all_rmats_introns, max_boundary_drift=30, suffix='_rMATS')
        
        # Step B: Intersect the LeafCutter-rMATS consensus with MAJIQ (using '_MAJIQ' suffix)
        if not majiq_data.empty and not lc_rmats_overlap.empty:
            overlap_matches = merge_overlap_genomic(lc_rmats_overlap, majiq_data, max_boundary_drift=30, suffix='_MAJIQ')
        else:
            overlap_matches = lc_rmats_overlap

if not overlap_matches.empty:
        print(f"Total rows before filtering: {len(overlap_matches)}")
        print("Columns available in overlap_matches:", overlap_matches.columns.tolist())

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

    output_file = f"{OUTPUT_DIR}/multitool_streamlined_filtered.csv"
    filtered_consensus.to_csv(output_file, index=False)
    
    print(f"Successfully filtered down to {len(filtered_consensus)} high-confidence consensus events.")
    print(f"Saved results to '{output_file}'")
else:
    print("No overlapping consensus events found across the tools.")