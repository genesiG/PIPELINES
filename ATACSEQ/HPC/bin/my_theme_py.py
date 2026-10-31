#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
my_theme_py.py - Shared figure style for the paper (retention-axis figures and beyond)

Canonical copy: /path/to/utils/my_theme.py
R equivalent:   /path/to/utils/my_theme.R (same palette and sizes)

Use from any project (do not copy; import the canonical file so figures stay consistent):

    import sys
    sys.path.insert(0, "/path/to/utils")
    import my_theme as th
    th.apply()                                  # rcParams: fonts, line widths, SVG text
    fig, ax = plt.subplots()
    th.fold_axis(ax, values)                    # log2 y-axis (fold over the expected 5%), 0.25x ... 4x
    th.panel_border(ax)                         # full frame (ggplot2 panel.border)
    th.panel_letter(ax, "A")
    th.save(fig, out_dir, "my_figure")          # .svg + .png

Colour roles (validated for colour-vision deficiency with the dataviz palette validator):
  - Diverging retention / share (fold over expected 5%) scale: NAVY (least retained) -> near-white -> PURPLE (most retained)
  - Treatments: TREATMENT_COLORS (DMSO is the neutral reference; always label it)
"""

import os

import matplotlib
from matplotlib.colors import LinearSegmentedColormap, Normalize

### Palette
PURPLE   = "#86187f"   # H3K27M line, EZH2i, most retained end of the diverging scale
STEEL    = "#6b93ad"   # untreated / DMSO reference, IDH-mutant lines
SLATE    = "#425b74"   # K27M-KO DMSO genome track
GREY     = "#a7b4c2"   # K27M-KO EZH2i genome track
NAVY     = "#1d3557"   # least retained end of the diverging scale
INK      = "#222222"   # text, panel borders
MUTED    = "#5f6670"   # secondary text
RULE     = "#9aa1a9"   # guides, baselines
UNSCORED = "#dfe3e7"   # bins excluded from scoring (schematics)

# Treatment comparisons (e.g. GG_EED_CUTANDRUN). Validated all-pairs: worst CVD dE 14.3, contrast >= 3:1.
# Inhibitor colours match the inhibitor viability figure (user, 2026-10-09): EZH2i = EPZ6438 black,
# EEDi 1 = EEDi-5285 dark purple, EEDi 2 = MAK683 blue-lilac; DMSO neutral grey (dashed, always labelled)
TREATMENT_COLORS = {"DMSO": "#808080", "EZH2i": "#000000", "EEDi 1": "#76069a", "EEDi 2": "#8d8dff"}
TREATMENT_LINESTYLE = {"DMSO": "--", "EZH2i": "-", "EEDi 1": "-", "EEDi 2": "-"}

CMAP  = LinearSegmentedColormap.from_list("share_vs_even", [NAVY, "#f7f7f7", PURPLE])
CNORM = Normalize(vmin=-2.3, vmax=2.3)            # log2(share fold over expected 5%); ~0.2x .. 4.9x
###

### Labels and sizes
N_GROUPS        = 20
XLABEL          = "Retention Group"
# A group is 5% of the scored bins, so evenly spread signal puts 5% of it in each group;
# a value of 2x means the group holds 10% of the library's signal.
FOLD_DESC       = "fold over expected 5%"
FOLD_LABEL      = f"Share of signal\n({FOLD_DESC})"
EXPECTED_LABEL  = "Expected share (5%)"
END_LABELS      = ("Least retained", "Most retained")
FOLD_TICKS      = [0.25, 0.5, 1, 2, 4]
FOLD_TICKLABELS = ["0.25x", "0.5x", "1x", "2x", "4x"]
SIZE = dict(base=11, axis_title=13, legend=9, annotation=8, heatmap_value=8, panel_letter=16,
            colorbar_label=9, tick=11)
LINEWIDTH = dict(axis=0.6, border=0.8, data_main=2.5, data_other=1.8, reference=0.8)
MARKER_SIZE = 4

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": SIZE["base"],
    "svg.fonttype": "none",              # keep text editable in Illustrator
    "pdf.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": LINEWIDTH["axis"],
    "xtick.major.width": LINEWIDTH["axis"], "ytick.major.width": LINEWIDTH["axis"],
    "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
}
###


def apply() -> None:
    """Set the shared rcParams (call once, after importing matplotlib)."""
    matplotlib.rcParams.update(RC)


def panel_border(ax, lw: float = LINEWIDTH["border"]) -> None:
    """Full rectangular frame around the plotting area (ggplot2 theme(panel.border = ...))."""
    for side, sp in ax.spines.items():
        sp.set_visible(True)
        sp.set_color(INK)
        sp.set_linewidth(lw)


def fold_axis(ax, values=None, lo: float = 0.4, hi: float = 4.5) -> None:
    """log2 y-axis for shares as fold over the expected 5% per group, with a dashed 1x reference line."""
    ax.set_yscale("log", base=2)
    ax.set_yticks(FOLD_TICKS)
    ax.set_yticklabels(FOLD_TICKLABELS)
    ax.minorticks_off()
    if values is not None:
        lo = min(lo, min(values) * 0.85)
        hi = max(hi, max(values) * 1.1)
    ax.set_ylim(lo, hi)
    ax.axhline(1, color="#333", ls="--", lw=LINEWIDTH["reference"])
    ax.text(1, 1.05, EXPECTED_LABEL, fontsize=SIZE["annotation"])


def group_axis(ax, label: bool = True, ends: bool = True) -> None:
    """x-axis over the 20 retention groups, with optional title and least/most retained end labels."""
    ax.set_xticks([1, 5, 10, 15, 20])
    ax.set_xlim(0.5, N_GROUPS + 0.5)
    if label:
        ax.set_xlabel(XLABEL, fontsize=SIZE["axis_title"])
    if ends:
        for xa, lab, ha in ((0, END_LABELS[0], "left"), (1, END_LABELS[1], "right")):
            ax.annotate(lab, xy=(xa, 0), xycoords="axes fraction", xytext=(0, -30),
                        textcoords="offset points", fontsize=SIZE["annotation"], va="top", ha=ha)


def colorbar_fold(fig, cax, label: str = None):
    """Colour bar for the diverging share heatmaps (fold over the expected 5% per group)."""
    label = FOLD_LABEL if label is None else label
    import matplotlib.pyplot as plt
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=CNORM, cmap=CMAP), cax=cax)
    cb.set_ticks([-2, -1, 0, 1, 2])
    cb.set_ticklabels(FOLD_TICKLABELS)
    cb.set_label(label, fontsize=SIZE["colorbar_label"])
    cb.outline.set_visible(False)
    cax.tick_params(labelsize=SIZE["annotation"])
    return cb


def panel_letter(ax, letter: str, x: float = -0.12, y: float = 1.02) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=SIZE["panel_letter"], fontweight="bold")


def save(fig, folder: str, name: str, exts=("svg", "png"), dpi: int = 300) -> None:
    os.makedirs(folder, exist_ok=True)
    for ext in exts:
        fig.savefig(os.path.join(folder, f"{name}.{ext}"), dpi=dpi)
    print(f"Saved: {os.path.join(folder, name)}.{'/.'.join(exts)}", flush=True)
