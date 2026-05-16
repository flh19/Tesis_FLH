#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
Variant-level quality control (QC) metrics visualization
===============================================================================

Description
-----------
This script generates publication-quality histograms comparing variant-level
quality control (QC) metrics before and after filtering procedures in genomic
variant datasets.

The figure contains three panels:

    (a) Mean sequencing depth per variant (DP)
    (b) Mean genotype quality per variant (GQ)
    (c) Mapping quality per variant (MQ)

For each metric, the script overlays:
    - Pre-filtering distributions
    - Post-filtering distributions
    - Threshold lines used for quality filtering

The resulting visualization allows assessment of how filtering improves the
overall quality profile of detected variants.

This type of figure is commonly used in:
    - Whole-exome sequencing studies,
    - Variant QC pipelines,
    - Rare variant association studies,
    - GWAS preprocessing workflows.

Input files
-----------
The script requires six text files containing numeric values:

    DP_mean_per_variant_before.txt
    DP_mean_per_variant.txt

    GQ_mean_per_variant_before.txt
    GQ_mean_per_variant.txt

    MQ_variant_clean_before.txt
    MQ_variant_clean.txt

Each file should contain one numeric value per line.

Output
------
The script generates:
    - QC_variants_DP_GQ_MQ_hist_before_after.png
    - QC_variants_DP_GQ_MQ_hist_before_after.pdf

Output files are stored inside:
    figuras_filtros/

The script also prints summary QC statistics to the console.

Usage
-----
Run directly from the command line:

    python plot_variant_qc_metrics.py

Dependencies
------------
    numpy
    matplotlib

Author
------
Prepared for genomic variant quality control visualization.

===============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.ticker import FuncFormatter
from matplotlib.gridspec import GridSpec

# ===================== CONFIGURATION ===================== #

# Output directory
OUTPUT_DIR = Path("figuras_filtros")
OUTPUT_DIR.mkdir(exist_ok=True)

# Color configuration
COLOR_AFTER  = "#008080"
COLOR_BEFORE = "#B0BEC5"
COLOR_EDGE   = "black"
COLOR_THRESH = "#C62828"

# Transparency settings
ALPHA_BEFORE = 0.5
ALPHA_AFTER  = 0.65

# Output resolution
DPI = 600

# Global matplotlib style settings
plt.rcParams.update({
    "font.family": "sans-serif",
    "axes.linewidth": 1.0,
    "axes.spines.top": False,
    "axes.spines.right": False
})

# Formatter for y-axis labels using thin-space thousands separator
def comma_formatter(x, pos):
    return f"{int(x):,}".replace(",", "\u2009")

# ===================== LOAD DATA ===================== #

# Mean depth per variant
dp_before = np.loadtxt("DP_mean_per_variant_before.txt")
dp_after  = np.loadtxt("DP_mean_per_variant.txt")

# Mean genotype quality per variant
gq_before = np.loadtxt("GQ_mean_per_variant_before.txt")
gq_after  = np.loadtxt("GQ_mean_per_variant.txt")

# Mapping quality per variant
mq_before = np.loadtxt("MQ_variant_clean_before.txt")
mq_after  = np.loadtxt("MQ_variant_clean.txt")

# ===================== QC SUMMARY METRICS ===================== #

def resumen_metricas(valores, umbral):
    """
    Compute summary statistics for a QC metric.

    Parameters
    ----------
    valores : array-like
        Input metric values.
    umbral : float
        QC threshold value.

    Returns
    -------
    tuple
        Mean, median, number of variants passing threshold,
        and total number of variants.
    """

    media   = np.mean(valores)
    mediana = np.median(valores)
    pasadas = np.sum(valores >= umbral)
    total   = len(valores)

    return media, mediana, pasadas, total

# Compute QC summaries for each metric
dp_b = resumen_metricas(dp_before, 20)
dp_a = resumen_metricas(dp_after, 20)

gq_b = resumen_metricas(gq_before, 30)
gq_a = resumen_metricas(gq_after, 30)

mq_b = resumen_metricas(mq_before, 40)
mq_a = resumen_metricas(mq_after, 40)

# ===================== FIGURE LAYOUT ===================== #

# Create figure
fig = plt.figure(figsize=(18, 7))

# Define panel layout
gs = GridSpec(
    1, 5,
    figure=fig,
    width_ratios=[1.4, 0.04, 1.4, 0.18, 0.9],
    wspace=0.30
)

# Create subplots
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[2])
ax_c = fig.add_subplot(gs[4])

# ======================
# (a) Mean depth per variant (DP)
# ======================

ax = ax_a

# Pre-filter distribution
ax.hist(
    dp_before,
    bins=100,
    histtype="stepfilled",
    color=COLOR_BEFORE,
    alpha=ALPHA_BEFORE,
    edgecolor=COLOR_EDGE,
    linewidth=0.8
)

# Post-filter distribution
ax.hist(
    dp_after,
    bins=100,
    histtype="stepfilled",
    color=COLOR_AFTER,
    alpha=ALPHA_AFTER,
    edgecolor=COLOR_EDGE,
    linewidth=0.8
)

# Threshold line
ax.axvline(20, color=COLOR_THRESH, linestyle="--", linewidth=3.0)

# Axis configuration
ax.set_xlim(0, 80)
ax.set_yscale("log")

ax.set_xlabel("Profundidad media de cobertura", fontsize=23)
ax.set_ylabel("Nº variantes (log)", fontsize=22)

ax.tick_params(axis="both", labelsize=22)

# Panel label
ax.text(
    -0.21, 1.10, "(a)",
    transform=ax.transAxes,
    fontsize=24,
    fontweight="bold",
    va="top"
)

# ======================
# (b) Mean genotype quality per variant (GQ)
# ======================

ax = ax_b

# Pre-filter distribution
ax.hist(
    gq_before,
    bins=100,
    histtype="stepfilled",
    color=COLOR_BEFORE,
    alpha=ALPHA_BEFORE,
    edgecolor=COLOR_EDGE,
    linewidth=0.8
)

# Post-filter distribution
ax.hist(
    gq_after,
    bins=100,
    histtype="stepfilled",
    color=COLOR_AFTER,
    alpha=ALPHA_AFTER,
    edgecolor=COLOR_EDGE,
    linewidth=0.8
)

# Threshold line
ax.axvline(30, color=COLOR_THRESH, linestyle="--", linewidth=3.0)

# Axis configuration
ax.set_xlim(0, 100)
ax.set_yscale("log")

ax.set_xlabel("Calidad media de genotipo", fontsize=23)
ax.set_ylabel("Nº variantes (log)", fontsize=22)

ax.tick_params(axis="both", labelsize=22)

# Panel label
ax.text(
    -0.24, 1.10, "(b)",
    transform=ax.transAxes,
    fontsize=24,
    fontweight="bold",
    va="top"
)

# ======================
# (c) Mapping quality per variant (MQ)
# ======================

ax = ax_c

# Define histogram bins
bins_mq = np.arange(35, 66)

# Pre-filter distribution
ax.hist(
    mq_before,
    bins=bins_mq,
    histtype="stepfilled",
    color=COLOR_BEFORE,
    alpha=ALPHA_BEFORE,
    edgecolor=COLOR_EDGE,
    linewidth=0.8,
    label="Prefiltrado"
)

# Post-filter distribution
ax.hist(
    mq_after,
    bins=bins_mq,
    histtype="stepfilled",
    color=COLOR_AFTER,
    alpha=ALPHA_AFTER,
    edgecolor=COLOR_EDGE,
    linewidth=0.8,
    label="Posfiltrado"
)

# Threshold line
ax.axvline(40, color=COLOR_THRESH, linestyle="--", linewidth=3.0)

# Axis configuration
ax.set_xlim(35, 65)

ax.set_xlabel("Calidad de mapeo", fontsize=23)
ax.set_ylabel("Nº variantes", fontsize=22)

ax.tick_params(axis="both", labelsize=22)

# Panel label
ax.text(
    -0.60, 1.10, "(c)",
    transform=ax.transAxes,
    fontsize=24,
    fontweight="bold",
    va="top"
)

# Format y-axis labels
ax.yaxis.set_major_formatter(FuncFormatter(comma_formatter))

# ======================
# Shared legend
# ======================

# Retrieve legend handles
handles, labels = ax.get_legend_handles_labels()

# Compute centered legend position
bbox = ax_b.get_position()
x_center = bbox.x0 + bbox.width / 2

# Add shared legend
fig.legend(
    handles,
    labels,
    loc="lower center",
    bbox_to_anchor=(x_center, -0.02),
    ncol=2,
    frameon=False,
    fontsize=23
)

# ======================
# Final layout adjustments
# ======================

fig.subplots_adjust(
    left=0.06,
    right=0.98,
    top=0.92,
    bottom=0.22
)

# ======================
# Save output figures
# ======================

fig.savefig(
    OUTPUT_DIR / "QC_variants_DP_GQ_MQ_hist_before_after.png",
    dpi=DPI
)

fig.savefig(
    OUTPUT_DIR / "QC_variants_DP_GQ_MQ_hist_before_after.pdf",
    dpi=DPI
)

# Close figure
plt.close(fig)

# ======================
# Print QC summary statistics
# ======================

print("\nResumen de métricas de calidad por variante")
print("==========================================")

print("DP medio por variante:")
print(
    f"  Antes  -> Media: {dp_b[0]:.2f}, "
    f"Mediana: {dp_b[1]:.2f}, "
    f"≥20X: {dp_b[2]}/{dp_b[3]} "
    f"({100*dp_b[2]/dp_b[3]:.2f}%)".replace(".", ",")
)

print(
    f"  Después-> Media: {dp_a[0]:.2f}, "
    f"Mediana: {dp_a[1]:.2f}, "
    f"≥20X: {dp_a[2]}/{dp_a[3]} "
    f"({100*dp_a[2]/dp_a[3]:.2f}%)".replace(".", ",")
)

print()

print("GQ medio por variante:")

print(
    f"  Antes  -> Media: {gq_b[0]:.2f}, "
    f"Mediana: {gq_b[1]:.2f}, "
    f"≥30: {gq_b[2]}/{gq_b[3]} "
    f"({100*gq_b[2]/gq_b[3]:.2f}%)".replace(".", ",")
)

print(
    f"  Después-> Media: {gq_a[0]:.2f}, "
    f"Mediana: {gq_a[1]:.2f}, "
    f"≥30: {gq_a[2]}/{gq_a[3]} "
    f"({100*gq_a[2]/gq_a[3]:.2f}%)".replace(".", ",")
)

print()

print("MQ por variante:")

print(
    f"  Antes  -> Media: {mq_b[0]:.2f}, "
    f"Mediana: {mq_b[1]:.2f}, "
    f"≥40: {mq_b[2]}/{mq_b[3]} "
    f"({100*mq_b[2]/mq_b[3]:.2f}%)".replace(".", ",")
)

print(
    f"  Después-> Media: {mq_a[0]:.2f}, "
    f"Mediana: {mq_a[1]:.2f}, "
    f"≥40: {mq_a[2]}/{mq_a[3]} "
    f"({100*mq_a[2]/mq_a[3]:.2f}%)".replace(".", ",")
)

print()
