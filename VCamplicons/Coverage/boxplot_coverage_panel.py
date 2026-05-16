#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Global mean coverage boxplot per cohort.

Description:
This script generates publication-quality boxplots comparing
the mean sequencing coverage per sample between cohorts.
Outliers are displayed individually using jittered scatter points.
"""

# =========================================================
# Imports
# =========================================================

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import FuncFormatter


# =========================================================
# Number formatting utilities
# =========================================================

def miles_con_espacio_fino(x):
    """
    Format integer values using thin spaces
    as thousands separators.
    """

    return f"{int(x):,}".replace(",", "\u202F")


def formatter_miles(x, pos):
    """
    Formatter wrapper for matplotlib axes.
    """

    return miles_con_espacio_fino(x)


# =========================================================
# Data loading
# =========================================================

df = pd.read_csv(
    "cobertura_media_por_muestra.tsv",
    sep="\t"
)

# Cohort display order

orden = [
    "Di@bet.es",
    "Hortega"
]

# Extract mean coverage values per cohort

datos = [
    df[df["poblacion"] == p]["cobertura_media"].dropna()
    for p in orden
]

# Number of samples per cohort

n_samples = [
    len(y)
    for y in datos
]


# =========================================================
# Cohort color palette
# =========================================================

colores = {
    "Di@bet.es": "#4682B4",
    "Hortega":   "#B0BEC5"
}


# =========================================================
# Global plotting style (Nature-like)
# =========================================================

plt.rcParams.update({

    "font.family": "sans-serif",

    "font.sans-serif": ["Arial"],

    "font.size": 10,

    "axes.linewidth": 1.0,
    "axes.labelsize": 11,

    "xtick.labelsize": 10,
    "ytick.labelsize": 10,

    "xtick.major.size": 4,
    "ytick.major.size": 4,

    "xtick.major.width": 1,
    "ytick.major.width": 1
})


# =========================================================
# Figure initialization
# =========================================================

fig, ax = plt.subplots(
    figsize=(5.2, 4.8)
)


# =========================================================
# Y-axis formatting
# =========================================================

ax.yaxis.set_major_formatter(
    FuncFormatter(formatter_miles)
)


# =========================================================
# Boxplot generation
# =========================================================

bp = ax.boxplot(
    datos,

    widths=0.55,

    patch_artist=True,

    showfliers=False,

    medianprops=dict(
        color="black",
        linewidth=2.0
    ),

    boxprops=dict(
        linewidth=1.2,
        edgecolor="black"
    ),

    whiskerprops=dict(
        linewidth=1.2
    ),

    capprops=dict(
        linewidth=1.2
    )
)

# Apply cohort-specific colors

for box, poblacion in zip(
    bp["boxes"],
    orden
):

    box.set(
        facecolor=colores[poblacion]
    )


# =========================================================
# Outlier visualization
# =========================================================

rng = np.random.default_rng(42)

jitter_sigma = 0.02

for i, y in enumerate(datos, start=1):

    # Interquartile range calculation

    q1, q3 = np.percentile(
        y,
        [25, 75]
    )

    iqr = q3 - q1

    # Outlier detection

    outliers = y[
        (y < q1 - 1.5 * iqr) |
        (y > q3 + 1.5 * iqr)
    ]

    # Jittered x positions

    x = rng.normal(
        i,
        jitter_sigma,
        size=len(outliers)
    )

    # Scatter plot of outliers

    ax.scatter(
        x,
        outliers,

        s=18,

        color=colores[orden[i - 1]],

        alpha=0.9,

        edgecolors="black",

        linewidths=0.3
    )


# =========================================================
# Axis labels and ticks
# =========================================================

ax.set_ylabel(
    "Mean coverage per base (X)"
)

ax.set_xlabel("")

ax.set_xticks([1, 2])

ax.set_xticklabels(orden)


# =========================================================
# Sample size annotations
# =========================================================

ax.set_ylim(600, 4800)

for i, n in enumerate(n_samples, start=1):

    ax.text(
        i,

        ax.get_ylim()[0] -
        (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.06,

        f"(N = {miles_con_espacio_fino(n)})",

        ha="center",
        va="top",

        fontsize=9
    )


# =========================================================
# Figure aesthetics
# =========================================================

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")


# =========================================================
# Figure export
# =========================================================

plt.tight_layout(
    pad=0.6
)

plt.savefig(
    "nivel1_cobertura_global_por_poblacion.png",
    dpi=600
)

plt.savefig(
    "nivel1_cobertura_global_por_poblacion.pdf"
)

plt.close()
