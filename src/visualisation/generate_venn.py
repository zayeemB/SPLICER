import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import pandas as pd
import pyranges as pr

# Import your custom modular parser functions
from SPLICER.src.parsers.leafcutter_parser import parse_leafcutter
from SPLICER.src.parsers.majiq_parser import parse_majiq_voila
from SPLICER.src.parsers.rmats_parser import (
    parse_rmats_a3ss_a5ss,
    parse_rmats_mxe,
    parse_rmats_ri,
    parse_rmats_se,
)

OUTPUT_DIR = "./output"
print(
    "Processing tool outputs using modular parsers and significance filters..."
)

# 1. LOAD & FILTER rMATS OUTPUTS (All Types)
RMATS_DIR = "/Users/zaiem/Desktop/BTP Tool/SPLICER/data/rmats"
rmats_dfs = []

rmats_tasks = [
    ("SE.csv", parse_rmats_se, None),
    ("RI.csv", parse_rmats_ri, None),
    ("MXE.csv", parse_rmats_mxe, None),
    ("A3SS.csv", parse_rmats_a3ss_a5ss, "A3SS"),
    ("A5SS.csv", parse_rmats_a3ss_a5ss, "A5SS"),
]

for filename, parser_func, event_type in rmats_tasks:
  file_path = os.path.join(RMATS_DIR, filename)
  if os.path.exists(file_path):
    try:
      if event_type:
        df = parser_func(file_path, event_type)
      else:
        df = parser_func(file_path)
      if df is not None and not df.empty:
        rmats_dfs.append(df)
    except Exception as e:
      print(f"Skipping rMATS file {filename} due to error: {e}")

if rmats_dfs:
  rmats_combined = pd.concat(rmats_dfs, ignore_index=True)
  # Apply rMATS significance filter
  rmats_filtered = rmats_combined[rmats_combined["FDR"] < 0.05].copy()
else:
  rmats_filtered = pd.DataFrame(columns=["Chromosome", "Start", "End"])

pr_rmats = (
    pr.PyRanges(rmats_filtered[["Chromosome", "Start", "End"]])
    if not rmats_filtered.empty
    else None
)
print(f"Significant rMATS events (FDR < 0.05): {len(rmats_filtered)}")

# 2. LOAD & FILTER MAJIQ OUTPUT
MAJIQ_FILE = "/Users/zaiem/Desktop/BTP Tool/SPLICER/data/majiq/tsv_f.csv"
majiq_df = pd.DataFrame(columns=["Chromosome", "Start", "End"])

if os.path.exists(MAJIQ_FILE):
  try:
    raw_majiq = parse_majiq_voila(MAJIQ_FILE)
    # Apply MAJIQ significance filter
    if "ProbabilityChanging" in raw_majiq.columns:
      majiq_df = raw_majiq[raw_majiq["ProbabilityChanging"] >= 0.8].copy()
    else:
      majiq_df = raw_majiq.copy()
  except Exception as e:
    print(f"Could not load MAJIQ file: {e}")

pr_majiq = (
    pr.PyRanges(majiq_df[["Chromosome", "Start", "End"]])
    if not majiq_df.empty
    else None
)
print(
    "Significant MAJIQ events (ProbabilityChanging >= 0.8):"
    f" {len(majiq_df)}"
)

# 3. LOAD LEAFCUTTER OUTPUT
LC_FILE = (
    "/Users/zaiem/Desktop/BTP Tool/SPLICER/data/leafcutter/leafcutter_ds_effect"
    "_sizes.csv"
)
lc_df = pd.DataFrame(columns=["Chromosome", "Start", "End"])

if os.path.exists(LC_FILE):
  try:
    lc_df = parse_leafcutter(LC_FILE)
    # If LeafCutter needs a p-value filter, add it here (e.g., p.adjust < 0.05)
  except Exception as e:
    print(f"Could not load LeafCutter file: {e}")

pr_lc = (
    pr.PyRanges(lc_df[["Chromosome", "Start", "End"]])
    if not lc_df.empty
    else None
)
print(f"Total LeafCutter events: {len(lc_df)}")

# ==========================================
# 4. CALCULATE VENN SUBSETS SAFELY (Clean Joins)
# ==========================================
n_lc = len(lc_df) if not lc_df.empty else 0
n_rmats = len(rmats_filtered) if not rmats_filtered.empty else 0
n_majiq = len(majiq_df) if not majiq_df.empty else 0

if (
    pr_lc is not None
    and pr_rmats is not None
    and pr_majiq is not None
    and n_lc > 0
    and n_rmats > 0
    and n_majiq > 0
):
  # Pairwise intersections with column cleanup to prevent suffix collision errors
  res_lc_rmats = pr_lc.join(pr_rmats)
  lc_rmats_clean = pr.PyRanges(
      res_lc_rmats.df[["Chromosome", "Start", "End"]].drop_duplicates()
  )
  lc_rmats_any = len(lc_rmats_clean)

  res_lc_majiq = pr_lc.join(pr_majiq)
  lc_majiq_clean = pr.PyRanges(
      res_lc_majiq.df[["Chromosome", "Start", "End"]].drop_duplicates()
  )
  lc_majiq_any = len(lc_majiq_clean)

  res_rmats_majiq = pr_rmats.join(pr_majiq)
  rmats_majiq_clean = pr.PyRanges(
      res_rmats_majiq.df[["Chromosome", "Start", "End"]].drop_duplicates()
  )
  rmats_majiq_any = len(rmats_majiq_clean)

  # Three-way intersection using the cleaned intermediate object
  res_all = lc_rmats_clean.join(pr_majiq)
  all_three = len(
      res_all.df[["Chromosome", "Start", "End"]].drop_duplicates()
  )
else:
  all_three = 0
  lc_rmats_any = 0
  lc_majiq_any = 0
  rmats_majiq_any = 0

lc_rmats_only = max(0, lc_rmats_any - all_three)
lc_majiq_only = max(0, lc_majiq_any - all_three)
rmats_majiq_only = max(0, rmats_majiq_any - all_three)

only_lc = max(0, n_lc - lc_rmats_only - lc_majiq_only - all_three)
only_rmats = max(0, n_rmats - lc_rmats_only - rmats_majiq_only - all_three)
only_majiq = max(0, n_majiq - lc_majiq_only - rmats_majiq_only - all_three)

# Order for venn3: (LeafCutter only, rMATS only, LC&rMATS, MAJIQ only, LC&MAJIQ, rMATS&MAJIQ, All Three)
subsets = (
    only_lc,
    only_rmats,
    lc_rmats_only,
    only_majiq,
    lc_majiq_only,
    rmats_majiq_only,
    all_three,
)
print(f"Calculated Genomic Interval Venn Subsets: {subsets}")

# ==========================================
# 5. PLOT VENN DIAGRAM
# ==========================================
plt.figure(figsize=(8, 8))
venn = venn3(subsets=subsets, set_labels=("LeafCutter", "rMATS", "MAJIQ"))

for subset_id in ["100", "010", "001", "110", "101", "011", "111"]:
  patch = venn.get_patch_by_id(subset_id)
  if patch:
    patch.set_alpha(0.5)

plt.title(
    "Multi-Tool Genomic Interval Overlap (Venn Diagram)",
    fontsize=14,
    fontweight="bold",
    pad=20,
)
os.makedirs(OUTPUT_DIR, exist_ok=True)
output_path = f"{OUTPUT_DIR}/multitool_venn_diagram.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()

print(
    f"Successfully generated genomic interval Venn diagram and saved to"
    f" '{output_path}'"
)