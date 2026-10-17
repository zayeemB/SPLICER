"""Grouped bar chart of the genes with the largest splicing change in the merged data.

Run from the project root:  python -m src.visualisation.gene_splice_variations

For each gene the junction with the largest mean |DeltaPSI| across the three tools is used (a gene
can have several junctions with opposite signs, so averaging them would cancel real changes out).
The 30 genes with the largest values are plotted, with one bar per tool. DeltaPSI is Hyp - Norm.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.pipeline import OUTPUT_DIR

MERGED_FILE = OUTPUT_DIR / "multitool_summary.csv"
OUTPUT_PATH = OUTPUT_DIR / "top_genes_deltapsi_comparison.png"
TOP_N = 30
TOOLS = ['LeafCutter', 'rMATS', 'MAJIQ']
COLORS = {'LeafCutter': '#e8605d', 'rMATS': '#5aa469', 'MAJIQ': '#5b6ee8'}


def top_genes(df, n=TOP_N):
    """Returns the n genes with the largest mean |DeltaPSI| across tools (strongest junction each)."""
    delta_cols = [f'{t}_DeltaPSI' for t in TOOLS]
    df = df.dropna(subset=['Gene'] + delta_cols).copy()
    df['mean_abs_dpsi'] = df[delta_cols].abs().mean(axis=1)
    best = df.sort_values('mean_abs_dpsi', ascending=False).drop_duplicates('Gene')
    return best.head(n)


def main():
    df = pd.read_csv(MERGED_FILE)
    top = top_genes(df)
    print(f"{len(df)} merged rows -> {df['Gene'].nunique()} genes; plotting top {len(top)}")

    labels = [f"{g}\n({e})" for g, e in zip(top['Gene'], top['EventType'])]
    x = np.arange(len(top))
    width = 0.27

    fig, ax = plt.subplots(figsize=(18, 7))
    for i, tool in enumerate(TOOLS):
        ax.bar(x + (i - 1) * width, top[f'{tool}_DeltaPSI'], width, label=tool, color=COLORS[tool])

    ax.axhline(0, color='black', linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=60, ha='right', fontsize=9)
    ax.set_ylabel("DeltaPSI (Hyp - Norm)", fontsize=12, fontweight='bold')
    ax.set_title(f"Top {len(top)} genes by splicing change (strongest junction per gene)",
                 fontsize=14, fontweight='bold', pad=15)
    ax.legend(title="Tool")
    ax.margins(x=0.01)
    fig.tight_layout()

    OUTPUT_DIR.mkdir(exist_ok=True)
    fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved bar chart to '{OUTPUT_PATH}'")


if __name__ == "__main__":
    main()
