#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
Minor allele frequency (MAF) category distribution visualization
===============================================================================

Description
-----------
This script generates a publication-quality bar plot summarizing the
distribution of genetic variants across predefined minor allele frequency
(MAF) categories.

Variants are classified into:
    - Rare variants
    - Common variants

according to standard allele frequency thresholds commonly used in
population genetics and rare variant association studies.

The figure includes:
    - Percentage of variants within each MAF interval,
    - Category labels,
    - Rare/common variant group annotations,
    - Percentage labels above each bar,
    - Publication-ready formatting.

This type of visualization is commonly used in:
    - Whole-exome sequencing studies,
    - Rare variant analyses,
    - Population genetics,
    - Variant quality control reports,
    - Genomic cohort characterization.

Input
-----
A PLINK-style frequency file:

    Frecuencias_allpopulation.frq

Required column:
    - MAF

Output
------
The script generates:
    - MAF_categories_publication.png
    - MAF_categories_publication.pdf
    - MAF_categories_publication.svg

Usage
-----
Run directly from the command line:

    python plot_maf_categories.py

Dependencies
------------
    pandas
    numpy
    matplotlib

Author
------
Prepared for genomic allele frequency visualization and sequencing studies.

===============================================================================
"""

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FuncFormatter
import matplotlib.gridspec as gridspec

# ============================================================
# Configuration
# ============================================================

# Global matplotlib style configuration
matplotlib.rcParams.update({

    # Font settings
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial"],

    # Font sizes
    "font.size":          16,
    "axes.labelsize":     22,
    "axes.titlesize":     22,
    "xtick.labelsize":    18,
    "ytick.labelsize":    18,

    # Axis line styling
    "axes.linewidth":     1.0,

    # Tick styling
    "xtick.major.width":  1.0,
    "ytick.major.width":  1.0,
    "xtick.minor.width":  0.6,
    "ytick.minor.width":  0.6,

    # Tick direction
    "xtick.direction":    "out",
    "ytick.direction":    "out",

    # Figure resolution
    "figure.dpi":         150,
    "savefig.dpi":        300,

    # Font embedding for vector graphics
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
})

# ============================================================
# 1. Load and clean data
# ============================================================

# Load allele frequency data
df = pd.read_csv(
    "Frecuencias_allpopulation.frq",
    delim_whitespace=True
)

# Extract MAF column
maf = df["MAF"]

# Remove monomorphic variants and invalid frequencies
maf = maf[(maf > 0) & (maf < 0.5)]

# ============================================================
# 2. Define non-overlapping frequency bins
# ============================================================

# MAF interval limits
bins = [0, 0.001, 0.01, 0.05, 0.5]

# Category labels
labels = [
    r"$\!<\!0,\!001$",
    r"$0,\!001\!\leq\!\mathrm{MAF}\!<\!0,\!01$",
    r"$0,\!01\!\leq\!\mathrm{MAF}\!<\!0,\!05$",
    r"$\!\geq\!0,\!05$",
]

# Assign variants to MAF categories
maf_cat = pd.cut(
    maf,
    bins=bins,
    labels=labels,
    right=False,
    include_lowest=True,
)

# Count variants per category
counts = maf_cat.value_counts().sort_index()

# Calculate proportions
proportions = counts / counts.sum()

# Total number of variants
n_total = counts.sum()

# ============================================================
# Formatting utilities
# ============================================================

def format_thousands(x):
    """
    Format integer values using thin-space thousands separator.
    """
    return f"{int(x):,}".replace(",", "\u2009")

# ============================================================
# Print summary table
# ============================================================

print("=" * 48)
print(f"{'Categoría':<18} {'n':>8}  {'%':>7}")
print("-" * 48)

for lbl, cnt, prop in zip(labels, counts, proportions):

    print(
        f"{lbl:<18} "
        f"{format_thousands(cnt):>8}  "
        f"{str(f'{prop*100:>6.2f}').replace('.', ',')}%"
    )

print("-" * 48)

print(
    f"{'TOTAL':<18} "
    f"{format_thousands(n_total):>8}  "
    f"{'100,00%':>7}"
)

print("=" * 48)

# ============================================================
# 3. Colour palette
# ============================================================

# Rare variant colors
RARE_DARK  = "#008B8B"
RARE_LIGHT = "#008B8B"

# Common variant colors
COMMON_DARK  = "lightgrey"
COMMON_LIGHT = "lightgrey"

# Bar colors
bar_colors = [
    RARE_DARK,
    RARE_LIGHT,
    COMMON_DARK,
    COMMON_LIGHT
]

# ============================================================
# 4. Figure layout
# ============================================================

# Create figure
fig = plt.figure(figsize=(15, 5.8))

# Grid specification
gs = gridspec.GridSpec(
    1, 1,
    top=0.88,
    bottom=0.18,
    left=0.10,
    right=0.98
)

# Create subplot
ax = fig.add_subplot(gs[0, 0])

# X-axis positions
x = np.arange(len(labels))

# Bar width
bw = 0.55

# Create bar plot
bars = ax.bar(
    x,
    proportions.values * 100,
    width=bw,
    color=bar_colors,
    edgecolor="black",
    linewidth=0.7,
    zorder=3,
)

# ============================================================
# Axis formatting
# ============================================================

def integer_formatter(x, pos):
    """
    Format integer tick labels using thin-space separator.
    """
    return format_thousands(x)

# Apply formatter to y-axis
ax.yaxis.set_major_formatter(
    FuncFormatter(integer_formatter)
)

# Define y-axis tick intervals
ax.yaxis.set_major_locator(
    MultipleLocator(20)
)

# Add horizontal grid
ax.grid(
    axis="y",
    which="major",
    linestyle="--",
    linewidth=0.6,
    color="grey",
    alpha=0.4,
    zorder=0
)

# ============================================================
# Axis labels and styling
# ============================================================

# X-axis labels
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=19)

# Y-axis label
ax.set_ylabel("Variantes (%)", fontsize=22, labelpad=8)

# Empty x-axis label placeholder
ax.set_xlabel("", labelpad=6)

# Axis limits
ax.set_xlim(-0.5, len(labels) - 0.5)

ax.set_ylim(
    0,
    max(proportions.values * 100) * 1.24
)

# Remove upper and right borders
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Tick styling
ax.tick_params(
    axis="both",
    which="both",
    top=False,
    right=False,
    labelsize=18
)

# ============================================================
# Percentage labels above bars
# ============================================================

for bar, prop, cnt in zip(bars, proportions.values, counts):

    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + max(proportions.values * 100) * 0.014,
        f"{str(f'{prop * 100:.2f}').replace('.', ',')}%",
        ha="center",
        va="bottom",
        fontsize=19,
        color="#222222",
    )

# ============================================================
# 5. Group brackets
# ============================================================

def draw_bracket(
    ax,
    fig,
    x_left_bar,
    x_right_bar,
    label,
    color,
    y_frac=0.85
):
    """
    Draw annotation brackets grouping related bars.
    """

    trans = ax.transData + ax.transAxes.inverted()

    ax_left, _ = trans.transform(
        (x_left_bar - 0.5 * 0.55 * 0.88, 0)
    )

    ax_right, _ = trans.transform(
        (x_right_bar + 0.5 * 0.55 * 0.88, 0)
    )

    ax_mid = (ax_left + ax_right) / 2

    fig_trans = ax.transAxes + fig.transFigure.inverted()

    fl, _ = fig_trans.transform((ax_left, 0))
    fr, _ = fig_trans.transform((ax_right, 0))

    fm = (fl + fr) / 2

    y_brac = y_frac
    y_tick = y_frac - 0.012
    y_text = y_frac + 0.012

    line_kw = dict(
        transform=fig.transFigure,
        color=color,
        linewidth=1.0,
        clip_on=False
    )

    fig.add_artist(
        matplotlib.lines.Line2D(
            [fl, fr],
            [y_brac, y_brac],
            **line_kw
        )
    )

    fig.add_artist(
        matplotlib.lines.Line2D(
            [fl, fl],
            [y_tick, y_brac],
            **line_kw
        )
    )

    fig.add_artist(
        matplotlib.lines.Line2D(
            [fr, fr],
            [y_tick, y_brac],
            **line_kw
        )
    )

    fig.text(
        fm,
        y_text,
        label,
        ha="center",
        va="bottom",
        fontsize=18,
        color=color,
        transform=fig.transFigure
    )

# Draw bracket for rare variants
draw_bracket(
    ax,
    fig,
    x_left_bar=0,
    x_right_bar=1,
    label="Variantes raras",
    color="black"
)

# Draw bracket for common variants
draw_bracket(
    ax,
    fig,
    x_left_bar=2,
    x_right_bar=3,
    label="Variantes comunes",
    color="black"
)

# ============================================================
# Total sample size annotation
# ============================================================

ax.text(
    0.98,
    0.85,
    f"Nº Total = {format_thousands(n_total)}",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=19,
    color="black",
)

# ============================================================
# Figure title
# ============================================================

fig.suptitle(
    "",
    x=0.545,
    y=0.985,
    ha="center",
    va="top",
    fontsize=18,
    fontweight="bold",
)

# ============================================================
# 6. Save figures
# ============================================================

for fmt in ("png", "pdf", "svg"):

    fig.savefig(
        f"MAF_categories_publication.{fmt}",
        dpi=300,
        bbox_inches="tight"
    )

    print(f"Guardado → MAF_categories_publication.{fmt}")

# Display figure
plt.show()
