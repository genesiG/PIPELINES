#!/usr/bin/env python3
"""qc_library_complexity_plot.py -- library complexity figure, drawn only from library_complexity.tsv.

Called by QC_LIBCOMPLEXITY_PLOT (modules/local/qc_libcomplexity_plot.nf):
    qc_library_complexity_plot.py --table library_complexity.tsv --theme my_theme_py.py

Arguments:
  --table   library_complexity.tsv (from qc_library_complexity_table.py)
  --theme   my_theme_py.py (palette and rcParams); a minimal built-in style is used when absent
  --out     output SVG (default library_complexity.svg)
Output: library_complexity.svg, four panels, one bar per sample, coloured by antibody (colours cycle through the
  theme palette in order of appearance in the data): NRF (guides 0.8 and 0.9), PBC1 (0.8 and 0.9), PBC2 (3 and 10),
  Picard EstimateLibraryComplexity % duplication.
"""
import argparse
import importlib.util
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

FALLBACK = ["#86187f", "#6b93ad", "#425b74", "#a7b4c2", "#1d3557"]


def load_theme(path):
    if path and os.path.exists(path):
        spec = importlib.util.spec_from_file_location("my_theme", path)
        th = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(th)
        th.apply()
        return [th.PURPLE, th.STEEL, th.SLATE, th.GREY, th.NAVY], th.INK
    matplotlib.rcParams.update({"svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False})
    return FALLBACK, "#222222"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True)
    ap.add_argument("--theme", default=None)
    ap.add_argument("--out", default="library_complexity.svg")
    a = ap.parse_args()
    palette, ink = load_theme(a.theme)

    d = pd.read_csv(a.table, sep="\t", na_values=["NA", ""])
    if d.empty:
        sys.stderr.write("ERROR: empty table %s\n" % a.table)
        sys.exit(3)
    ctrl = d["antibody"].astype(str).str.fullmatch(r"(?i)igg|input")
    d = d.assign(_ctrl=ctrl).sort_values(["_ctrl", "antibody", "group", "sample_id"], kind="stable").reset_index(drop=True)
    abs_ = list(dict.fromkeys(d["antibody"]))
    col = {ab: palette[i % len(palette)] for i, ab in enumerate(abs_)}
    colors = [col[x] for x in d["antibody"]]

    panels = [
        ("NRF", "NRF (distinct / total fragments)", [(0.9, "0.9"), (0.8, "0.8")], (0, 1.05)),
        ("PBC1", "PBC1 (single-count / distinct)", [(0.9, "0.9"), (0.8, "0.8")], (0, 1.05)),
        ("PBC2", "PBC2 (single-count / double-count)", [(10, "10"), (3, "3")], None),
        ("elc_pct_duplication", "Picard ELC duplication (%)", [], (0, 100)),
    ]
    n = len(d)
    w = max(7, 0.32 * n + 2)
    fig, axes = plt.subplots(len(panels), 1, figsize=(w, 2.6 * len(panels) + 1), sharex=True)
    for ax, (colname, ylabel, guides, ylim) in zip(axes, panels):
        vals = d[colname] if colname in d else pd.Series([float("nan")] * n)
        if vals.isna().all():
            ax.text(0.5, 0.5, "not available", ha="center", va="center", transform=ax.transAxes, color="#5f6670")
            ax.set_yticks([])
        else:
            ax.bar(range(n), vals.fillna(0), color=colors, width=0.75)
            for y, lab in guides:
                ax.axhline(y, ls="--", lw=0.8, color=ink)
                ax.text(n - 0.4, y, lab, va="bottom", ha="right", fontsize=8, color=ink)
            if ylim:
                ax.set_ylim(*ylim)
        ax.set_ylabel(ylabel, fontsize=9)
    axes[-1].set_xticks(range(n))
    axes[-1].set_xticklabels(d["sample_id"], rotation=45, ha="right", fontsize=7)
    handles = [plt.Rectangle((0, 0), 1, 1, color=col[ab]) for ab in abs_]
    axes[0].legend(handles, abs_, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(a.out, format="svg")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
