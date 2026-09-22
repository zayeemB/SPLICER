import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Load your integrated dataset
INTEGRATED_FILE = "/Users/zaiem/Desktop/BTP Tool/SPLICER/output/multitool_streamlined_filtered.csv"
df = pd.read_csv(INTEGRATED_FILE)

# Define columns for each tool's DeltaPSI
col_lc = "IncLevelDifference"  # LeafCutter
col_rmats = "IncLevelDifference_rMATS"  # rMATS
col_majiq = "IncLevelDifference_MAJIQ"  # MAJIQ

# 1. Clean data: Drop rows missing any tool's DeltaPSI or gene symbol
df_clean = df.dropna(subset=[col_lc, col_rmats, col_majiq, "geneSymbol"]).copy()

# 2. Aggregate per gene (mean DeltaPSI if a gene has multiple events)
gene_aggregated = (
    df_clean.groupby("geneSymbol")[[col_lc, col_rmats, col_majiq]]
    .mean()
    .reset_index()
)

# 3. Select top 20 genes with the highest overall absolute effect size magnitude across tools
gene_aggregated["mean_abs_dpsi"] = gene_aggregated[
    [col_lc, col_rmats, col_majiq]
].abs().mean(axis=1)
top_genes = gene_aggregated.sort_values(
    by="mean_abs_dpsi", ascending=False
).head(20)

# 4. Melt the dataframe into long format for Seaborn grouping
melted_df = top_genes.melt(
    id_vars=["geneSymbol"],
    value_vars=[col_lc, col_rmats, col_majiq],
    var_name="Tool",
    value_name="DeltaPSI",
)

# Map clean tool names for the legend
tool_name_map = {col_lc: "LeafCutter", col_rmats: "rMATS", col_majiq: "MAJIQ"}
melted_df["Tool"] = melted_df["Tool"].map(tool_name_map)

# ==========================================
# 5. GENERATE GROUPED BAR PLOT
# ==========================================
plt.figure(figsize=(14, 7))
sns.barplot(
    data=melted_df,
    x="geneSymbol",
    y="DeltaPSI",
    hue="Tool",
    palette="Set2",
    alpha=0.9,
)

plt.axhline(0, color="black", linestyle="--", linewidth=1)
plt.title(
    "Top 20 Genes: Multi-Tool DeltaPSI Comparison",
    fontsize=14,
    fontweight="bold",
    pad=15,
)
plt.xlabel("Gene Symbol", fontsize=12, fontweight="bold")
plt.ylabel("Mean DeltaPSI", fontsize=12, fontweight="bold")
plt.xticks(rotation=45, ha="right", fontsize=10)
plt.legend(title="Tool", fontsize=10, title_fontsize=11)
plt.tight_layout()

output_path = "./output/top_genes_deltapsi_comparison.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()

print(
    f"Successfully generated gene-level comparison plot and saved to"
    f" '{output_path}'"
)