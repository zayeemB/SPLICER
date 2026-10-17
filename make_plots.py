"""Generates every plot. Run after run.py:  python make_plots.py

Plots are written to output/plots/. To add a plot, give its module a main() and list it below.
"""
import sys

from src.config import DATASETS_DIR, PLOTS_DIR
from src.visualisation import event_type_comparison, gene_splice_variations, generate_venn

PLOTS = {
    "Venn diagram": generate_venn.main,
    "Top genes bar chart": gene_splice_variations.main,
    "Event type comparison": event_type_comparison.main,
}

MERGED_FILE = DATASETS_DIR / "multitool_summary.csv"

if not MERGED_FILE.exists():
    sys.exit(f"{MERGED_FILE} not found. Run run.py first.")

PLOTS_DIR.mkdir(parents=True, exist_ok=True)
failed = []
for name, make_plot in PLOTS.items():
    print(f"\n=== {name} ===")
    try:
        make_plot()
    except Exception as exc:  # keep going so one broken plot doesn't block the rest
        print(f"FAILED: {exc!r}")
        failed.append(name)

print(f"\n{len(PLOTS) - len(failed)}/{len(PLOTS)} plots generated in '{PLOTS_DIR}'")
if failed:
    sys.exit(f"Failed: {', '.join(failed)}")
