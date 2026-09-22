import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import pandas as pd
import pyranges as pr

# Import your custom modular parser functions
from leafcutter_parser import parse_leafcutter
from majiq_parser import parse_majiq_voila
from rmats_parser import (
    parse_rmats_a3ss_a5ss,
    parse_rmats_mxe,
    parse_rmats_ri,
    parse_rmats_se,
)

OUTPUT_DIR = "./output"
print("Processing tool outputs with robust coordinate & strand-agnostic joining...")


def clean_chromosomes(df):
  """Standardizes chromosome names to ensure consistent merging (adds 'chr' if missing)."""
  if df is not None and not df.empty and "Chromosome" in df.columns:
    df["Chromosome"] = df["Chromosome"].astype(str)
    # If chromosomes don't start with 'chr', prepend it for uniformity
    mask = ~df["Chromosome"].str.startswith("chr")
    df.loc[mask, "Chromosome"] = "chr" + df.loc[mask, "Chromosome"]
  return df


# ==========================================
# 1. LOAD & PARSE rMATS OUTPUTS (All Types)
# ==========================================
RMATS_DIR = "/Users/zaiem/Desktop/BTP Tool/SPLICER/data/rmats"
rmats_dfs = []

rmats_tasks = [
    ("SE.MATS.JC.csv", parse_rmats_se),
    ("RI.MATS.JC.csv", parse_rmats_ri),
    ("MXE.MATS.JC.csv", parse_rmats_mxe),
    ("A3SS.MATS.JC.csv", parse_rmats_a3ss_a5ss),
    ("A5SS.MATS.JC.csv", parse_rmats_a3ss_a5ss),
]

for filename, parser_func in rmats_tasks:
  file_path = os.path.join(RMATS_DIR, filename)
  if os.path.exists(file_path):
    try:
      df = parser_func(file_path)
      if df is not None and not df.empty:
        rmats_dfs.append(clean_chromosomes(df))
    except Exception as e:
      print(f"Skipping rMATS file {filename} due to error: {e}")

if rmats_dfs:
  rmats_combined = pd.concat(rmats_dfs, ignore_index=True).drop_duplicates(
      subset=["Chromosome", "Start", "End"]
  )
else:
  rmats_combined = pd.DataFrame(columns=["Chromosome", "Start", "End"])

pr_rmats = pr.PyRanges(rmats_combined) if not rmats_combined.empty else None
print(f"Total unique significant rMATS events: {len(rmats_combined)}")

# ==========================================
# 2. LOAD & PARSE MAJIQ OUTPUT
# ==========================================
MAJIQ_FILE = "/Users/zaiem/Desktop/BTP Tool/SPLICER/data/majiq/tsv_f.csv"
majiq_df = pd.DataFrame(columns=["Chromosome", "Start", "End"])

if os.path.exists(MAJIQ_FILE):
  try:
    majiq_df = parse_majiq_voila(MAJIQ_FILE)
    majiq_df = clean_chromosomes(majiq_df)
  except Exception as e:
    print(f"Could not load MAJIQ file: {e}")

pr_majiq = pr.PyRanges(majiq_df) if not majiq_df.empty else None
print(f"Total unique significant MAJIQ events: {len(majiq_df)}")

# ==========================================
# 3. LOAD & PARSE LEAFCUTTER OUTPUT
# ==========================================
LC_FILE = (
    "/Users/zaiem/Desktop/BTP Tool/SPLICER/data/leafcutter/leafcutter_ds_effect"
    "_sizes.csv"
)
lc_df = pd.DataFrame(columns=["Chromosome", "Start", "End"])

if os.path.exists(LC_FILE):
  try:
    lc_df = parse_leafcutter(LC_FILE)
    lc_df = clean_chromosomes(lc_df)
  except Exception as e:
    print(f"Could not load LeafCutter file: {e}")

pr_lc = pr.PyRanges(lc_df) if not lc_df.empty else None
print(f"Total unique significant LeafCutter events: {len(lc_df)}")

# ==========================================
# 4. CALCULATE VENN SUBSETS SAFELY (Strand-Agnostic)
# ==========================================
n_lc = len(lc_df) if not lc_df.empty else 0
n_rmats = len(rmats_combined) if not rmats_combined.empty else 0
n_majiq = len(majiq_df) if not majiq_df.empty else 0

if (
    pr_lc is not None
    and pr_rmats is not None
    and pr_majiq is not None
    and n_lc > 0
    and n_rmats > 0
    and n_majiq > 0
):
  # Use stranded=False to prevent strand-mismatch filtering bugs
  all_three = len(
      pr_lc.join(pr_rmats, stranded=False).join(pr_majiq, stranded=False)
  )
  lc_rmats_any = len(pr_lc.join(pr_rmats, stranded=False))
  lc_majiq_any = len(pr_lc.join(pr_majiq, stranded=False))
  rmats_majiq_any = len(pr_rmats.join(pr_majiq, stranded=False))
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