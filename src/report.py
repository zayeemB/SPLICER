import numpy as np

TOOLS = ['LeafCutter', 'rMATS', 'MAJIQ']


def _pct(n, total):
    return f"{n} ({100 * n / total:.1f}%)" if total else "0"


def sign_agreement(summary):
    """Returns a boolean Series: True where all three tools report the same DeltaPSI sign."""
    signs = np.sign(summary[[f'{t}_DeltaPSI' for t in TOOLS]])
    return signs.nunique(axis=1) == 1


def write_report(summary, filter_counts, report_path, disagreements_path, max_boundary_drift):
    """Writes a short plain-text sign-agreement report and a CSV of the disagreeing rows."""
    delta = {t: summary[f'{t}_DeltaPSI'] for t in TOOLS}
    signs = {t: np.sign(s) for t, s in delta.items()}
    agree = sign_agreement(summary)
    total = len(summary)

    lines = ["SPLICER multi-tool report", "========================", ""]

    lines.append("Rows passing per-tool filters (before merging)")
    for tool, (n_in, n_out, rule) in filter_counts.items():
        lines.append(f"  {tool:<11} {n_out:>6} / {n_in:<7} {rule}")
    lines += [f"  Matching: strand-aware, within {max_boundary_drift} bp, closest match only",
              f"  Merged consensus rows: {total}", ""]

    lines.append("DeltaPSI sign agreement (Hyp - Norm)")
    lines.append(f"  All three agree:  {_pct(int(agree.sum()), total)}")
    lines.append(f"  Any disagreement: {_pct(int((~agree).sum()), total)}")
    for i, a in enumerate(TOOLS):
        for b in TOOLS[i + 1:]:
            same = int((signs[a] == signs[b]).sum())
            lines.append(f"  {a} vs {b}: {_pct(same, total)}")
    lines.append("")

    lines.append("Odd one out (disagreements where only one tool differs)")
    for t in TOOLS:
        a, b = [u for u in TOOLS if u != t]
        odd = (signs[a] == signs[b]) & (signs[t] != signs[a])
        lines.append(f"  {t:<11} {int(odd.sum())}")
    lines.append("")

    lines.append("By event type (rows / sign disagreements)")
    for ev, idx in summary.groupby('EventType').groups.items():
        lines.append(f"  {ev:<5} {len(idx):>5} / {int((~agree.loc[idx]).sum())}")
    lines.append("")

    lines.append("Spearman correlation of DeltaPSI")
    for i, a in enumerate(TOOLS):
        for b in TOOLS[i + 1:]:
            lines.append(f"  {a} vs {b}: {delta[a].corr(delta[b], method='spearman'):.2f}")
    lines += ["", f"Rows with disagreeing signs are listed in {disagreements_path.name}"]

    report_path.write_text("\n".join(lines) + "\n")
    summary[~agree].to_csv(disagreements_path, index=False)
