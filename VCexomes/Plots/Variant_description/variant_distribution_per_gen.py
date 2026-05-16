#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
Distribution of variants per gene visualization
===============================================================================

Description
-----------
This script generates a publication-quality histogram summarizing the
distribution of the number of variants detected per gene in a genomic dataset.

Genes are grouped into predefined bins according to their number of variants,
allowing visualization of the overall mutational burden distribution across
genes.

Intergenic regions are excluded from the analysis.

The resulting figure is commonly used in:
    - Whole-exome sequencing studies,
    - Variant burden analyses,
    - Gene-level variant summaries,
    - Genomic quality control reports,
    - Rare variant studies.

The histogram displays:
    - Number of genes within each variant-count category,
    - Absolute counts above each bar,
    - Thousand-separated formatting for improved readability.

Input
-----
A tab-delimited file named:

    variants_per_gene.tsv

Required columns:
    - gene
    - n_variants

Output
------
The script generates:
    - figures/hist_variants_per_gene_final.png
    - figures/hist_variants_per_gene_final.pdf

The PDF version is vectorized and suitable for manuscripts and doctoral theses.

Usage
-----
Run directly from the command line:

    python plot_variants_per_gene.py

Dependencies
------------
    pandas
    matplotlib
    numpy

Author
------
Prepared for genomic variant distribution visualization.

===============================================================================
"""

#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# =========================
# LOAD AND FILTER DATA
# =========================

# Load gene-level variant counts
df = pd.read_csv("variants_per_gene.tsv", sep="\t")

# Remove intergenic regions
df = df[df["gene"] != "INTERGENIC"].copy()

# Extract number of variants per gene
n_variants = df["n_variants"]

# =========================
# GENERAL STYLE SETTINGS
# =========================

plt.rcParams.update({

    # Font configuration
    "font.family":       "DejaVu Sans",
    "font.size":         13,

    # Axis styling
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    0.8,
    "axes.edgecolor":    "black",

    # Tick styling
    "xtick.direction":   "out",
    "ytick.direction":   "out",
    "xtick.major.size":  4,
    "ytick.major.size":  4,
    "xtick.color":       "black",
    "ytick.color":       "black",

    # Figure resolution
    "figure.dpi":        300,
})

# Main bar color
COLOR_MAIN = "#008B8B"

# =========================
# THOUSANDS FORMATTER
# =========================

def format_thousands(x):
    """
    Format integer values using thin-space thousands separator.
    """
    return f"{int(x):,}".replace(",", "\u2009")

# =========================
# DEFINE HISTOGRAM BINS
# =========================

# Variant count intervals
bins = [0, 5, 10, 20, 50, 100, 200, 500, np.inf]

# Bin labels
labels = [
    "1 – 5",
    "6 – 10",
    "11 – 20",
    "21 – 50",
    "51 – 100",
    "101 – 200",
    "201 – 500",
    "> 500"
]

# Count number of genes within each interval
counts = (
    pd.cut(
        n_variants,
        bins=bins,
        labels=labels,
        right=True
    )
    .value_counts()
    .reindex(labels)
)

# =========================
# CREATE FIGURE
# =========================

# Create figure
fig, ax = plt.subplots(figsize=(15, 5))

# Generate bar plot
bars = ax.bar(
    range(len(labels)),
    counts.values,
    color=COLOR_MAIN,
    edgecolor="black",
    linewidth=0.4,
    width=0.65,
    zorder=2
)

# =========================
# VALUE LABELS ABOVE BARS
# =========================

# Add counts above each bar
for bar, val in zip(bars, counts.values):

    if val > 0:

        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + counts.max() * 0.012,
            format_thousands(val),
            ha="center",
            va="bottom",
            fontsize=20,
            color="black"
        )

# =========================
# AXIS CONFIGURATION
# =========================

# X-axis ticks and labels
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels)

# Tick label size
ax.tick_params(axis="both", labelsize=20)

# Axis labels
ax.set_xlabel(
    "Número de variantes",
    fontsize=23,
    labelpad=11
)

ax.set_ylabel(
    "Número de genes",
    fontsize=23,
    labelpad=11
)

# Empty title placeholder
ax.set_title(
    "",
    fontsize=14,
    fontweight="bold",
    pad=12
)

# Format y-axis labels using thousands separator
ax.yaxis.set_major_formatter(
    ticker.FuncFormatter(
        lambda x, _: format_thousands(x)
    )
)

# Add horizontal grid
ax.grid(
    axis="y",
    linewidth=0.4,
    color="#DDDDDD",
    zorder=1
)

# Place grid below bars
ax.set_axisbelow(True)

# =========================
# SAVE FIGURES
# =========================

# Optimize layout
plt.tight_layout()

# Save PNG version
plt.savefig(
    "figures/hist_variants_per_gene_final.png",
    dpi=300,
    bbox_inches="tight",
    facecolor="white"
)

# Save vector PDF version
plt.savefig(
    "figures/hist_variants_per_gene_final.pdf",
    bbox_inches="tight",
    facecolor="white"
)

# Close figure
plt.close()

# =========================
# FINAL STATUS MESSAGE
# =========================

print("Guardado: PNG y PDF en figures/")
