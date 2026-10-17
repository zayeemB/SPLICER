import pandas as pd

from src.config import DATASETS_DIR, MAX_BOUNDARY_DRIFT
from src.pipeline import load_filtered_tools
from src.analysis import rmats_rate_comparison, unmerged_junctions
from src.junction_table import build_junction_table
from src.report import write_report
from src.utils import merge_overlap_genomic

DATASETS_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = DATASETS_DIR / "multitool_streamlined_filtered.csv"
SUMMARY_FILE = DATASETS_DIR / "multitool_summary.csv"
JUNCTION_FILE = DATASETS_DIR / "multitool_junctions.csv"
REPORT_FILE = DATASETS_DIR / "multitool_report.txt"
DISAGREEMENTS_FILE = DATASETS_DIR / "sign_disagreements.csv"

# 1. Load each tool and apply its filters (shared with the visualisation scripts)
lc_data, all_rmats_introns, majiq_data, filter_counts = load_filtered_tools()

# 2. Multi-tool spatial reconciliation (LeafCutter <-> rMATS <-> MAJIQ)
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
    # Filters were already applied per tool before merging, so every merged row is a consensus event.
    filtered_consensus = overlap_matches.copy()

    filtered_consensus.to_csv(OUTPUT_FILE, index=False)
    
    print(f"Found {len(filtered_consensus)} high-confidence consensus events.")
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

    # One row per junction (event-level detail summarised, not discarded)
    junction_table = build_junction_table(filtered_consensus)
    junction_table.to_csv(JUNCTION_FILE, index=False)
    print(f"Saved junction-level table ({len(junction_table)} junctions) to '{JUNCTION_FILE}'")

    # 6. Sign-agreement report
    # Filtered junctions that never made it into the merged data
    unmerged = unmerged_junctions(
        {'LeafCutter': lc_data, 'rMATS': all_rmats_introns, 'MAJIQ': majiq_data},
        filtered_consensus, MAX_BOUNDARY_DRIFT)
    print(unmerged.to_string(index=False))

    # Were the unmerged rMATS junctions less strongly expressed?
    tool_data = {'LeafCutter': lc_data, 'rMATS': all_rmats_introns, 'MAJIQ': majiq_data}
    rate_summary, rate_pvalues = rmats_rate_comparison(
        all_rmats_introns, tool_data, filtered_consensus, MAX_BOUNDARY_DRIFT)
    print(rate_summary.to_string(index=False))

    write_report(summary, filter_counts, REPORT_FILE, DISAGREEMENTS_FILE, MAX_BOUNDARY_DRIFT, unmerged,
                 rate_summary, rate_pvalues)
    print(f"Saved sign-agreement report to '{REPORT_FILE}'")
else:
    print("No overlapping consensus events found across the tools.")
