#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
Variant filtering and quality control visualization pipeline
===============================================================================

Description
-----------
This script generates multiple publication-quality figures summarizing
variant filtering, genomic distribution, normalization metrics, and
transition/transversion (Ts/Tv) quality statistics from a sequencing dataset.

The generated figures are intended for genomic quality control reporting,
whole-exome sequencing studies, and variant filtering assessment workflows.

The script includes:

    Figure 1:
        Distribution of variant types before and after filtering.

    Figure 2:
        Chromosomal distribution of variants before and after filtering.

    Figure 3:
        Variant quality category distributions.

    Figure 4:
        Number of normalized variants per megabase of captured exonic region.

    Figure 5:
        Combined visualization of chromosomal distribution and normalized
        variant density.

    Figure 6:
        Combined visualization of variant types and chromosomal distribution.

    Figure 7:
        Sample-wise Ts/Tv ratio distribution.

    Figure 8:
        Combined visualization of normalized variants and Ts/Tv ratios.

All figures are exported in both PNG and PDF formats at 600 dpi.

Applications
------------
This workflow is commonly used in:
    - Whole-exome sequencing (WES) studies
    - Variant quality control pipelines
    - Rare variant association studies
    - Population genomics
    - Genomic data preprocessing

Output directory
----------------
All figures are saved inside:

    figuras_filtros/

Dependencies
------------
    numpy
    pandas
    matplotlib

Usage
-----
Run directly from the command line:

    python variant_filtering_visualization.py

Author
------
Prepared for genomic variant quality control and sequencing analysis.

===============================================================================
"""

#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.ticker import FuncFormatter

# ===================== GENERAL CONFIGURATION ===================== #

# Output directory for generated figures
OUTPUT_DIR = Path("figuras_filtros")
OUTPUT_DIR.mkdir(exist_ok=True)

# Figure resolution
DPI = 600

# Default bar width
BAR_WIDTH = 0.65

# Border color for bars
EDGE_COLOR = "black"

# Color palette
COLOR_BEFORE = "#B0BEC5"
COLOR_AFTER  = "#008080"

# Global matplotlib style configuration
plt.rcParams.update({
    "font.family": "sans-serif",
    "axes.linewidth": 1.0
})

# ===================== PERCENTAGE FORMATTER ===================== #

def percent_formatter(x, pos):
    """
    Format axis labels as percentages.
    """
    return f"{round(x)}%"

# Formatter object for y-axis
y_formatter = FuncFormatter(percent_formatter)

# ===================== AUXILIARY FUNCTIONS ===================== #

def normalizar(conteos, total):
    """
    Normalize counts to percentages.

    Parameters
    ----------
    conteos : array-like
        Raw counts.
    total : int
        Total number of variants.

    Returns
    -------
    numpy.ndarray
        Percentages relative to total.
    """
    return (np.array(conteos) / total) * 100


def estilo_ejes(ax, *, ymax=None):
    """
    Apply common axis styling.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object.
    ymax : float, optional
        Upper y-axis limit.
    """

    if ymax is None:
        ax.set_ylim(0, 100)
    else:
        ax.set_ylim(0, ymax * 1.05)

    ax.yaxis.set_major_formatter(y_formatter)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def ajustar_fuentes(ax, *, label_size, tick_size, legend_size):
    """
    Adjust font sizes for axis labels, ticks, and legends.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object.
    label_size : int
        Font size for axis labels.
    tick_size : int
        Font size for tick labels.
    legend_size : int
        Font size for legend text.
    """

    ax.xaxis.label.set_size(label_size)
    ax.yaxis.label.set_size(label_size)

    ax.tick_params(axis="both", labelsize=tick_size)

    legend = ax.get_legend()

    if legend is not None:
        for text in legend.get_texts():
            text.set_fontsize(legend_size)


def guardar_figura(fig, nombre):
    """
    Save figure in PNG and PDF formats.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object.
    nombre : str
        Output filename without extension.
    """

    fig.savefig(OUTPUT_DIR / f"{nombre}.png", dpi=DPI)
    fig.savefig(OUTPUT_DIR / f"{nombre}.pdf", dpi=DPI)

    plt.close(fig)

# ===================== FIGURE 1: VARIANT TYPES ===================== #

# Total number of variants
total = 1_977_270

# Variant categories before filtering
tipos = ["SNP", "Deleción", "Inserción", "STAR"]

antes = [1_752_472, 213_044, 128_145, 52_812]

# Variant categories after filtering
despues_tipos = ["SNP", "Deleción", "Inserción"]

despues = [655_040, 9_669, 6_630]

# Convert counts to percentages
antes_n = normalizar(antes, total)
despues_n = normalizar(despues, total)

# X-axis positions
x = np.arange(len(tipos))

# Create figure
fig, ax = plt.subplots(figsize=(8, 5))

# Plot pre-filtering variants
b1 = ax.bar(
    x,
    antes_n,
    BAR_WIDTH,
    color=COLOR_BEFORE,
    edgecolor=EDGE_COLOR,
    label="Variantes antes del filtrado"
)

# Plot post-filtering variants
idx_after = [tipos.index(t) for t in despues_tipos]

b2 = ax.bar(
    idx_after,
    despues_n,
    BAR_WIDTH,
    color=COLOR_AFTER,
    edgecolor=EDGE_COLOR,
    label="Variantes tras el filtrado"
)

# Axis labels and ticks
ax.set_ylabel("Porcentaje sobre el total de variantes")

ax.set_xticks(x)
ax.set_xticklabels(tipos)

# Legend
ax.legend(frameon=False)

# Apply axis style
estilo_ejes(ax)

# Adjust font sizes
ajustar_fuentes(
    ax,
    label_size=13,
    tick_size=12,
    legend_size=12
)

# Optimize layout
fig.tight_layout()

# Save figure
guardar_figura(fig, "figura_1_tipos_variante")

# ===================== FIGURE 2: CHROMOSOMAL DISTRIBUTION ===================== #

# Chromosome labels
cromosomas = [str(i) for i in range(1, 23)] + ["X", "Y"]

# Variant counts before filtering
antes = [
    201023,145364,110772,77395,86159,84154,97580,68271,87516,83765,
    116043,107015,35781,61384,67632,91204,112143,31302,138285,
    51038,23202,45916,53344,982
]

# Variant counts after filtering
despues = [
    69587,49435,39463,26694,29795,29539,33425,23909,30770,27309,
    42590,35769,11693,19975,21263,34300,40703,9849,51651,
    18933,7816,16871,0,0
]

# Normalize counts to percentages
antes_n = normalizar(antes, total)
despues_n = normalizar(despues, total)

# X-axis positions
x = np.arange(len(cromosomas))

# Create figure
fig, ax = plt.subplots(figsize=(15, 6))

# Plot pre-filtering distribution
b1 = ax.bar(
    x,
    antes_n,
    BAR_WIDTH,
    color=COLOR_BEFORE,
    edgecolor=EDGE_COLOR,
    label="Variantes antes del filtrado"
)

# Plot post-filtering distribution
b2 = ax.bar(
    x,
    despues_n,
    BAR_WIDTH,
    color=COLOR_AFTER,
    edgecolor=EDGE_COLOR,
    label="Variantes tras el filtrado"
)

# Axis configuration
ax.set_xlim(-0.5, len(cromosomas) - 0.5)

ax.set_ylabel("Porcentaje sobre el total de variantes")

ax.set_xticks(x)
ax.set_xticklabels(cromosomas)

# Force vertical legend ordering
handles, labels = ax.get_legend_handles_labels()

ax.legend(
    handles,
    labels,
    frameon=False,
    ncol=1
)

# Apply common axis style
estilo_ejes(ax, ymax=antes_n.max())

# Adjust font sizes
ajustar_fuentes(
    ax,
    label_size=20,
    tick_size=20,
    legend_size=20
)

# Optimize layout
fig.tight_layout()

# Save figure
guardar_figura(fig, "figura_2_distribucion_cromosomas")

# ===================== FIGURE 3: QUALITY CATEGORIES ===================== #

# Quality categories
categorias = [
    "PASS",
    "INDEL\n99,00–99,90",
    "INDEL\n99,90–100,00",
    "SNP\n99,00–99,90",
    "SNP\n99,90–100,00"
]

# Variant counts before filtering
antes = [1_802_502, 17_028, 7_861, 129_074, 20_805]

# Categories retained after filtering
despues_cat = ["PASS", "INDEL\n99,00–99,90"]

# Variant counts after filtering
despues = [1_619_548, 5_874]

# Normalize counts
antes_n = normalizar(antes, total)
despues_n = normalizar(despues, total)

# X-axis positions
x = np.arange(len(categorias))

# Create figure
fig, ax = plt.subplots(figsize=(10, 5))

# Plot pre-filtering categories
b1 = ax.bar(
    x,
    antes_n,
    BAR_WIDTH,
    color=COLOR_BEFORE,
    edgecolor=EDGE_COLOR,
    label="Variantes antes del filtrado"
)

# Plot post-filtering categories
idx_after = [categorias.index(c) for c in despues_cat]

b2 = ax.bar(
    idx_after,
    despues_n,
    BAR_WIDTH,
    color=COLOR_AFTER,
    edgecolor=EDGE_COLOR,
    label="Variantes tras el filtrado"
)

# Axis labels and ticks
ax.set_ylabel("Porcentaje sobre el total de variantes")

ax.set_xticks(x)
ax.set_xticklabels(categorias, ha="center")

# Legend
ax.legend(frameon=False)

# Apply common axis style
estilo_ejes(ax)

# Adjust font sizes
ajustar_fuentes(
    ax,
    label_size=14,
    tick_size=14,
    legend_size=14
)

# Optimize layout
fig.tight_layout()

# Save figure
guardar_figura(fig, "figura_3_categorias_calidad")

# ===================== FIGURE 4: NORMALIZED VARIANTS PER CHROMOSOME ===================== #

# Autosomal chromosome labels
cromosomas = [str(i) for i in range(1, 23)]

# Number of variants normalized by captured exonic megabases
variantes_normalizadas = [
    20047.70863, 19638.26034, 19996.65565, 19471.0276,
    18865.95327, 16973.21772, 20701.08073, 20752.53884,
    22092.50564, 20343.56633, 20921.13924, 19868.0242,
    18417.97284, 18527.62216, 17889.41426, 24029.70436,
    20681.36782, 18128.80675, 23082.07945, 23210.88152,
    22737.05768, 23699.31884
]

# X-axis positions
x = np.arange(len(cromosomas))

# Create figure
fig, ax = plt.subplots(figsize=(15, 6))

# Plot normalized variants
ax.bar(
    x,
    variantes_normalizadas,
    BAR_WIDTH,
    color=COLOR_AFTER,
    edgecolor=EDGE_COLOR
)

# Axis configuration
ax.set_xlim(-0.5, len(cromosomas) - 0.5)

ax.set_ylabel("Variantes por Mb de región exónica capturada")

ax.set_xticks(x)
ax.set_xticklabels(cromosomas)

# Axis styling (non-percentage scale)
ax.set_ylim(0, max(variantes_normalizadas) * 1.05)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Adjust font sizes
ajustar_fuentes(
    ax,
    label_size=20,
    tick_size=20,
    legend_size=20
)

# Optimize layout
fig.tight_layout()

# Save figure
guardar_figura(fig, "figura_4_variantes_normalizadas_cromosoma")

# ===================== FIGURE 5: COMBINED CHROMOSOMAL DISTRIBUTION + NORMALIZED VARIANTS ===================== #

# ---- Figure 2 data ----

# Chromosome labels including sex chromosomes
cromosomas_2 = [str(i) for i in range(1, 23)] + ["X", "Y"]

# Variant counts before filtering
antes_2 = np.array([
    201023,145364,110772,77395,86159,84154,97580,68271,87516,83765,
    116043,107015,35781,61384,67632,91204,112143,31302,138285,
    51038,23202,45916,53344,982
])

# Variant counts after filtering
despues_2 = np.array([
    69587,49435,39463,26694,29795,29539,33425,23909,30770,27309,
    42590,35769,11693,19975,21263,34300,40703,9849,51651,
    18933,7816,16871,0,0
])

# Normalize counts
antes_2_n = normalizar(antes_2, total)
despues_2_n = normalizar(despues_2, total)

# ---- Figure 4 data ----

# Autosomal chromosomes
cromosomas_4 = [str(i) for i in range(1, 23)]

# Normalized variants per chromosome
variantes_norm = np.array([
    20047.70863, 19638.26034, 19996.65565, 19471.0276,
    18865.95327, 16973.21772, 20701.08073, 20752.53884,
    22092.50564, 20343.56633, 20921.13924, 19868.0242,
    18417.97284, 18527.62216, 17889.41426, 24029.70436,
    20681.36782, 18128.80675, 23082.07945, 23210.88152,
    22737.05768, 23699.31884
])

# ---- Create combined figure ----

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 6))

# ===== Left subfigure: Chromosomal distribution =====

x2 = np.arange(len(cromosomas_2))

# Pre-filtering distribution
ax1.bar(
    x2,
    antes_2_n,
    BAR_WIDTH,
    color=COLOR_BEFORE,
    edgecolor=EDGE_COLOR,
    label="Antes del filtrado"
)

# Post-filtering distribution
ax1.bar(
    x2,
    despues_2_n,
    BAR_WIDTH,
    color=COLOR_AFTER,
    edgecolor=EDGE_COLOR,
    label="Tras el filtrado"
)

# Axis configuration
ax1.set_xlim(-0.5, len(cromosomas_2) - 0.5)

ax1.set_ylabel("Porcentaje de variantes")

ax1.set_xticks(x2)
ax1.set_xticklabels(cromosomas_2)

# Legend
ax1.legend(frameon=False)

# Apply axis style
estilo_ejes(ax1, ymax=antes_2_n.max())

# Adjust font sizes
ajustar_fuentes(
    ax1,
    label_size=18,
    tick_size=18,
    legend_size=20
)

# Panel label
ax1.text(
    0.02, 0.95, "(a)",
    transform=ax1.transAxes,
    fontsize=22,
    fontweight="bold",
    va="top",
    ha="left"
)

# ===== Right subfigure: Normalized variants =====

x4 = np.arange(len(cromosomas_4))

# Plot normalized variant density
ax2.bar(
    x4,
    variantes_norm,
    BAR_WIDTH,
    color=COLOR_AFTER,
    edgecolor=EDGE_COLOR
)

# Axis configuration
ax2.set_xlim(-0.5, len(cromosomas_4) - 0.5)

ax2.set_ylabel("Número de variantes")

ax2.set_xticks(x4)
ax2.set_xticklabels(cromosomas_4)

# Axis styling
ax2.set_ylim(0, variantes_norm.max() * 1.05)

ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# Adjust font sizes
ajustar_fuentes(
    ax2,
    label_size=18,
    tick_size=18,
    legend_size=18
)

# Panel label
ax2.text(
    0.02, 0.95, "(b)",
    transform=ax2.transAxes,
    fontsize=22,
    fontweight="bold",
    va="top",
    ha="left"
)

# ---- Save combined figure ----

fig.tight_layout()

guardar_figura(
    fig,
    "figura_5_distribucion_y_normalizacion_cromosomas"
)

# ===================== FIGURE 6: VARIANT TYPES + CHROMOSOMAL DISTRIBUTION ===================== #

# Create combined figure with custom width ratios
fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(22, 9),
    gridspec_kw={"width_ratios": [1, 2.3]}
)

# ===== Subfigure (a): Variant types (excluding STAR variants) =====

# Variant categories
tipos = ["SNP", "Deleción", "Inserción"]

# Variant counts before filtering
antes_tipos = [1_752_472, 213_044, 128_145]

# Variant counts after filtering
despues_tipos = [655_040, 9_669, 6_630]

# Normalize counts
antes_tipos_n = normalizar(antes_tipos, total)
despues_tipos_n = normalizar(despues_tipos, total)

# X-axis positions
x_tipos = np.arange(len(tipos))

# Plot pre-filtering variants
ax1.bar(
    x_tipos,
    antes_tipos_n,
    BAR_WIDTH,
    color=COLOR_BEFORE,
    edgecolor=EDGE_COLOR,
    label="Prefiltrado"
)

# Plot post-filtering variants
ax1.bar(
    x_tipos,
    despues_tipos_n,
    BAR_WIDTH,
    color=COLOR_AFTER,
    edgecolor=EDGE_COLOR,
    label="Posfiltrado"
)

# Axis configuration
ax1.set_ylabel("Variantes (%)")

ax1.set_xticks(x_tipos)
ax1.set_xticklabels(tipos)

# Apply axis styling
estilo_ejes(ax1)

# Format y-axis without percentage symbols
ax1.yaxis.set_major_formatter(
    FuncFormatter(lambda x, pos: f"{round(x)}")
)

# Adjust font sizes
ajustar_fuentes(
    ax1,
    label_size=30,
    tick_size=26,
    legend_size=30
)

# Panel label
ax1.text(
    -0.22, 1.15, "(a)",
    transform=ax1.transAxes,
    fontsize=27,
    fontweight="bold",
    va="top",
    ha="left"
)

# ===== Subfigure (b): Chromosomal distribution =====

# Chromosome labels
cromosomas = [str(i) for i in range(1, 23)] + ["X", "Y"]

# Variant counts before filtering
antes_chr = np.array([
    201023,145364,110772,77395,86159,84154,97580,68271,87516,83765,
    116043,107015,35781,61384,67632,91204,112143,31302,138285,
    51038,23202,45916,53344,982
])

# Variant counts after filtering
despues_chr = np.array([
    69587,49435,39463,26694,29795,29539,33425,23909,30770,27309,
    42590,35769,11693,19975,21263,34300,40703,9849,51651,
    18933,7816,16871,0,0
])

# Normalize counts
antes_chr_n = normalizar(antes_chr, total)
despues_chr_n = normalizar(despues_chr, total)

# X-axis positions
x_chr = np.arange(len(cromosomas))

# Plot pre-filtering chromosomal distribution
ax2.bar(
    x_chr,
    antes_chr_n,
    BAR_WIDTH,
    color=COLOR_BEFORE,
    edgecolor=EDGE_COLOR,
    label="Prefiltrado"
)

# Plot post-filtering chromosomal distribution
ax2.bar(
    x_chr,
    despues_chr_n,
    BAR_WIDTH,
    color=COLOR_AFTER,
    edgecolor=EDGE_COLOR,
    label="Posfiltrado"
)

# Axis configuration
ax2.set_xlim(-0.5, len(cromosomas) - 0.5)

ax2.set_ylabel("Variantes (%)")

ax2.set_xticks(x_chr)
ax2.set_xticklabels(cromosomas)

# Apply axis style
estilo_ejes(ax2, ymax=antes_chr_n.max())

# Format y-axis without percentage symbols
ax2.yaxis.set_major_formatter(
    FuncFormatter(lambda x, pos: f"{round(x)}")
)

# Adjust font sizes
ajustar_fuentes(
    ax2,
    label_size=30,
    tick_size=25,
    legend_size=30
)

# Panel label
ax2.text(
    -0.08, 1.15, "(b)",
    transform=ax2.transAxes,
    fontsize=27,
    fontweight="bold",
    va="top",
    ha="left"
)

# ===== Shared legend =====

# Retrieve legend elements
handles, labels = ax1.get_legend_handles_labels()

# Create shared legend
fig.legend(
    handles,
    labels,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.02),
    ncol=2,
    frameon=False,
    fontsize=30
)

# ---- Final layout adjustments ----

fig.tight_layout()

fig.subplots_adjust(
    wspace=0.20,
    bottom=0.16
)

# Save figure
guardar_figura(
    fig,
    "figura_6_tipos_variante_y_cromosomas"
)

# ===================== FIGURE 7: Ts/Tv RATIO DISTRIBUTION ===================== #

import pandas as pd
from matplotlib.ticker import FuncFormatter

# ---- Load Ts/Tv ratio data ----

df = pd.read_csv(
    "ts_tv_ratio.txt",
    sep="\t",
    decimal=","
)

# Extract Ts/Tv ratios
ratios = df["Ratio Ts/Tv"].values

# X-axis positions (samples)
x = np.arange(len(ratios))

# Mean Ts/Tv ratio
media = ratios.mean()

# ---- Create figure ----

fig, ax = plt.subplots(figsize=(10, 4))

# Scatter plot of Ts/Tv ratios across samples
ax.scatter(
    x,
    ratios,
    s=30,
    color="#696969",
    alpha=0.5
)

# ---- Mean reference line ----

ax.axhline(
    media,
    linestyle="-",
    linewidth=1.5,
    color="#8B0000",
    label="Media"
)

# ---- Axis labels ----

ax.set_xlabel("Muestras", fontsize=14)
ax.set_ylabel("Ratio Ts/Tv", fontsize=14)

# Hide sample IDs on x-axis
ax.set_xticks([])

# Dynamic y-axis limits
ax.set_ylim(
    max(0, ratios.min() - 0.2),
    ratios.max() + 0.2
)

# Legend
ax.legend(frameon=False, fontsize=14)

# Add horizontal grid
ax.grid(True, axis="y", linestyle=":", alpha=0.4)

# ---- Decimal comma formatter for y-axis ----

def comma_formatter(x, pos):
    """
    Format decimal numbers using comma separator.
    """
    return f"{x:.2f}".replace(".", ",")

# Apply formatter
ax.yaxis.set_major_formatter(FuncFormatter(comma_formatter))

# Tick label size
ax.tick_params(axis="y", labelsize=14)

# ---- Axis styling ----

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ---- Save figure ----

fig.tight_layout()

guardar_figura(fig, "figura_7_ratio_ts_tv")

# ===================== FIGURE 8: NORMALIZED VARIANTS + Ts/Tv RATIO ===================== #

# ---- Figure 4 data ----

# Autosomal chromosome labels
cromosomas = [str(i) for i in range(1, 23)]

# Normalized variants per chromosome
variantes_norm = np.array([
    20047.70863, 19638.26034, 19996.65565, 19471.0276,
    18865.95327, 16973.21772, 20701.08073, 20752.53884,
    22092.50564, 20343.56633, 20921.13924, 19868.0242,
    18417.97284, 18527.62216, 17889.41426, 24029.70436,
    20681.36782, 18128.80675, 23082.07945, 23210.88152,
    22737.05768, 23699.31884
])

# ---- Figure 7 data ----

# Reload Ts/Tv ratio data
df = pd.read_csv(
    "ts_tv_ratio.txt",
    sep="\t",
    decimal=","
)

# Extract ratios
ratios = df["Ratio Ts/Tv"].values

# X-axis positions for samples
x_ratios = np.arange(len(ratios))

# Mean Ts/Tv ratio
media = ratios.mean()

# ---- Create combined figure ----

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 6))

# ===== Subfigure (a): Normalized variants per chromosome =====

# X-axis positions for chromosomes
x_chr = np.arange(len(cromosomas))

# Plot normalized variants
ax1.bar(
    x_chr,
    variantes_norm,
    BAR_WIDTH,
    color=COLOR_AFTER,
    edgecolor=EDGE_COLOR
)

# Axis configuration
ax1.set_xlim(-0.5, len(cromosomas) - 0.5)

ax1.set_ylabel("Variantes Normalizadas")

ax1.set_xticks(x_chr)
ax1.set_xticklabels(cromosomas)

# Axis styling
ax1.set_ylim(0, variantes_norm.max() * 1.05)

ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Adjust font sizes
ajustar_fuentes(
    ax1,
    label_size=21,
    tick_size=18,
    legend_size=21
)

# Panel label
ax1.text(
    -0.13, 1.11, "(a)",
    transform=ax1.transAxes,
    fontsize=22,
    fontweight="bold",
    va="top",
    ha="left"
)

# ===== Subfigure (b): Ts/Tv ratio distribution =====

# Scatter plot of sample Ts/Tv ratios
ax2.scatter(
    x_ratios,
    ratios,
    s=30,
    color="#696969",
    alpha=0.5
)

# Mean Ts/Tv ratio line
ax2.axhline(
    media,
    linestyle="-",
    linewidth=1.5,
    color="#8B0000",
    label="Media"
)

# Axis labels
ax2.set_xlabel("Muestras", fontsize=21)
ax2.set_ylabel("Ratio Ts/Tv", fontsize=21)

# Hide x-axis sample labels
ax2.set_xticks([])

# Dynamic y-axis limits
ax2.set_ylim(
    max(0, ratios.min() - 0.2),
    ratios.max() + 0.2
)

# Legend
ax2.legend(frameon=False, fontsize=21)

# Horizontal grid
ax2.grid(True, axis="y", linestyle=":", alpha=0.4)

# Apply decimal comma formatter
ax2.yaxis.set_major_formatter(
    FuncFormatter(lambda x, pos: f"{x:.2f}".replace(".", ","))
)

# Tick label size
ax2.tick_params(axis="y", labelsize=18)

# Axis styling
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# Panel label
ax2.text(
    -0.12, 1.11, "(b)",
    transform=ax2.transAxes,
    fontsize=22,
    fontweight="bold",
    va="top",
    ha="left"
)

# ---- Final layout adjustments ----

fig.tight_layout()

fig.subplots_adjust(wspace=0.17)

# Save figure
guardar_figura(
    fig,
    "figura_8_normalizacion_y_ratio_ts_tv"
)
