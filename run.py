from src.leafcutter_parser import parse_leafcutter
import src.rmats_parser
import src.utils  


DATA_DIR = "./data"
OUTPUT_DIR = "./output"

leafcutter_file = f"{DATA_DIR}/leafcutter/leafcutter_ds_effect_sizes.csv"

# Load LeafCutter data
lc_data = parse_leafcutter(leafcutter_file)

# Load and parse all rMATS event files
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
        print(f"Successfully parsed {len(df_parsed)} introns from {event_type}")
    except FileNotFoundError:
        print(f"Skipping {event_type}: File {filename} not found.")
        
if rmats_dfs:
    all_rmats_introns = pd.concat(rmats_dfs, ignore_index=True)
    print(f"Total rMATS Introns Across All Event Types: {len(all_rmats_introns)}")
    
    # Run drift-tolerant spatial join against LeafCutter
    overlap_matches = merge_overlap_genomic(lc_data, all_rmats_introns, max_boundary_drift=30)

    # Define your significance and effect size thresholds
    PVALUE_THRESHOLD = 0.01      # Adjust as needed (e.g., 0.01 or 0.05)
    FDR_THRESHOLD = 0.05         # Standard multiple testing threshold
    MIN_DELTA_PSI = 0.10         # Minimum absolute change in inclusion level (10%)

    # 1. P-Value Filtered Dataset
    # Filters for events meeting the p-value cutoff AND showing a real magnitude change
    pvalue_filtered = overlap_matches[
        (overlap_matches['PValue'] < PVALUE_THRESHOLD) & 
        (overlap_matches['IncLevelDifference'].abs() >= MIN_DELTA_PSI)
    ].copy()

    pvalue_filtered.to_csv(f"{OUTPUT_DIR}/rmats_leafcutter_pvalue_filtered.csv", index=False)
    print(f"Saved {len(pvalue_filtered)} rows to 'rmats_leafcutter_pvalue_filtered.csv'")


    # 2. FDR Filtered Dataset
    # Filters for events meeting the False Discovery Rate cutoff AND showing a real magnitude change
    fdr_filtered = overlap_matches[
        (overlap_matches['FDR'] < FDR_THRESHOLD) & 
        (overlap_matches['IncLevelDifference'].abs() >= MIN_DELTA_PSI)
    ].copy()

    fdr_filtered.to_csv(f"{OUTPUT_DIR}/rmats_leafcutter_fdr_filtered.csv", index=False)
    print(f"Saved {len(fdr_filtered)} rows to 'rmats_leafcutter_fdr_filtered.csv'")

    print(f"Total validated matches across all event types: {len(overlap_matches)}")