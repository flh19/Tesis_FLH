#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Combined visualization of sequencing coverage and uniformity metrics.

Description:
This script generates a two-panel figure comparing:

(a) Mean sequencing coverage per gene
(b) Sequencing uniformity at ≥50X coverage

across the Di@bet.es and Hortega cohorts.
"""

# =========================================================
# Imports
# =========================================================

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter


# =========================================================
# Global plotting style (Nature-like)
# =========================================================

plt.rcParams.update({

    "font.family": "sans-serif",

    "font.sans-serif": ["Arial"],

    "font.size": 10,

    "axes.linewidth": 1.0,
    "axes.labelsize": 18,

    "xtick.labelsize": 16,
    "ytick.labelsize": 16,

    "xtick.major.size": 4,
    "ytick.major.size": 4,

    "xtick.major.width": 1,
    "ytick.major.width": 1,

    "pdf.fonttype": 42,
    "ps.fonttype": 42
})


# =========================================================
# Cohort color palette
# =========================================================

colores = {

    "Di@bet.es": "#4682B4",

    "Hortega": "#B0BEC5"
}


# =========================================================
# Figure and panel layout
# =========================================================

fig = plt.figure(
    figsize=(12.8, 5.4)
)

gs = GridSpec(
    1,
    2,

    width_ratios=[1, 1],

    wspace=0.18
)

ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])


# =========================================================
# Panel (a):
# Mean coverage per gene
# =========================================================

df = pd.read_csv(
    "cobertura_media_por_gen_y_muestra.tsv",
    sep="\t"
)

genes = [
    "ACAD10",
    "ACSM3",
    "AFMID",
    "GCKR",
    "NAT1"
]

poblaciones = [
    "Di@bet.es",
    "Hortega"
]

x = np.arange(len(genes))

width = 0.35

for i, pob in enumerate(poblaciones):

    # -----------------------------------------------------
    # Extract coverage values per gene
    # -----------------------------------------------------

    datos = [

        df[
            (df["gen"] == g) &
            (df["poblacion"] == pob)
        ]["cobertura"].dropna()

        for g in genes
    ]

    pos = x - width / 2 + i * width

    # -----------------------------------------------------
    # Boxplot generation
    # -----------------------------------------------------

    bp = ax_a.boxplot(
        datos,

        positions=pos,

        widths=width,

        patch_artist=True,

        showfliers=False,

        medianprops=dict(
            color="black",
            linewidth=2
        ),

        boxprops=dict(
            linewidth=1
        ),

        whiskerprops=dict(
            linewidth=1
        ),

        capprops=dict(
            linewidth=1
        )
    )

    # -----------------------------------------------------
    # Apply cohort-specific colors
    # -----------------------------------------------------

    for box in bp["boxes"]:

        box.set_facecolor(
            colores[pob]
        )

        box.set_alpha(0.75)

    # -----------------------------------------------------
    # Outlier visualization
    # -----------------------------------------------------

    for p, y in zip(pos, datos):

        q1, q3 = np.percentile(
            y,
            [25, 75]
        )

        iqr = q3 - q1

        outliers = y[
            (y < q1 - 1.5 * iqr) |
            (y > q3 + 1.5 * iqr)
        ]

        jitter = np.random.normal(
            0,
            0.02,
            size=len(outliers)
        )

        ax_a.scatter(
            np.full(len(outliers), p) + jitter,
            outliers,

            s=10,

            alpha=0.4,

            color=colores[pob],

            edgecolor="black",

            linewidth=0.3
        )

# ---------------------------------------------------------
# Axis formatting
# ---------------------------------------------------------

ax_a.set_xticks(x)

ax_a.set_xticklabels(
    genes,
    fontstyle="italic"
)

ax_a.set_ylabel(
    "Mean coverage per base (X)"
)

# Thin-space thousands separator formatting

ax_a.yaxis.set_major_formatter(
    FuncFormatter(
        lambda y, _:
        f"{int(y):,}".replace(",", "\u2009")
    )
)


# =========================================================
# Panel (b):
# Sequencing uniformity ≥50X
# =========================================================

df = pd.read_csv(
    "uniformidad_50X_por_gen_y_muestra.tsv",
    sep="\t"
)

for i, pob in enumerate(poblaciones):

    # -----------------------------------------------------
    # Extract uniformity values per gene
    # -----------------------------------------------------

    datos = [

        df[
            (df["gen"] == g) &
            (df["poblacion"] == pob)
        ]["uniformidad_50X"].dropna() * 100

        for g in genes
    ]

    pos = x - width / 2 + i * width

    # -----------------------------------------------------
    # Boxplot generation
    # -----------------------------------------------------

    bp = ax_b.boxplot(
        datos,

        positions=pos,

        widths=width,

        patch_artist=True,

        showfliers=False,

        medianprops=dict(
            color="black",
            linewidth=2
        ),

        boxprops=dict(
            linewidth=1
        ),

        whiskerprops=dict(
            linewidth=1
        ),

        capprops=dict(
            linewidth=1
        )
    )

    # -----------------------------------------------------
    # Apply cohort-specific colors
    # -----------------------------------------------------

    for box in bp["boxes"]:

        box.set_facecolor(
            colores[pob]
        )

        box.set_alpha(0.75)

    # -----------------------------------------------------
    # Outlier visualization
    # -----------------------------------------------------

    for p, y in zip(pos, datos):

        q1, q3 = np.percentile(
            y,
            [25, 75]
        )

        iqr = q3 - q1

        outliers = y[
            (y < q1 - 1.5 * iqr) |
            (y > q3 + 1.5 * iqr)
        ]

        jitter = np.random.normal(
            0,
            0.02,
            size=len(outliers)
        )

        ax_b.scatter(
            np.full(len(outliers), p) + jitter,
            outliers,

            s=10,

            alpha=0.4,

            color=colores[pob],

            edgecolor="black",

            linewidth=0.3
        )

# ---------------------------------------------------------
# Axis formatting
# ---------------------------------------------------------

ax_b.set_xticks(x)

ax_b.set_xticklabels(
    genes,
    fontstyle="italic"
)

ax_b.set_ylabel(
    "Bases ≥50X (%)"
)

ax_b.set_ylim(0, 100)

# Percentage formatting

ax_b.yaxis.set_major_formatter(
    FuncFormatter(
        lambda y, _:
        f"{y:.0f}".replace(".", ",")
    )
)


# =========================================================
# Axis cleanup
# =========================================================

for ax in (ax_a, ax_b):

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# =========================================================
# Panel labels
# =========================================================

ax_a.text(
    -0.18,
    1.08,

    "(a)",

    transform=ax_a.transAxes,

    fontsize=16,

    fontweight="bold",

    va="top"
)

ax_b.text(
    -0.16,
    1.08,

    "(b)",

    transform=ax_b.transAxes,

    fontsize=16,

    fontweight="bold",

    va="top"
)


# =========================================================
# Shared legend
# =========================================================

handles = [

    plt.Line2D(
        [0],
        [0],

        color=colores[p],

        lw=10,

        label=p
    )

    for p in poblaciones
]

fig.legend(
    handles=handles,

    frameon=False,

    loc="lower center",

    ncol=2,

    fontsize=16,

    handlelength=2.4,

    columnspacing=3.0,

    bbox_to_anchor=(0.5, 0.0001)
)


# =========================================================
# Fine margin adjustment
# =========================================================

for ax in (ax_a, ax_b):

    ax.tick_params(
        axis="x",
        pad=2
    )

fig.subplots_adjust(
    left=0.07,
    right=0.99,

    bottom=0.15,
    top=0.96,

    wspace=0.18
)


# =========================================================
# Figure export
# =========================================================

plt.savefig(
    "figura_cobertura_uniformidad.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    "figura_cobertura_uniformidad.pdf",
    bbox_inches="tight"
)

plt.close()
