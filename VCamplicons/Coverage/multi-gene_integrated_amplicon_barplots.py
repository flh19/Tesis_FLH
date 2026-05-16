#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Integrated multi-gene amplicon coverage barplots.

Description:
This script generates a multi-panel publication-quality figure
displaying per-amplicon mean coverage for multiple genes
across sequencing cohorts.

The workflow is fully generalizable by modifying the
GENES configuration list.
"""

# =========================================================
# Imports
# =========================================================

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib.ticker import FuncFormatter


# =========================================================
# Global plotting style
# =========================================================

plt.rcParams.update({

    "font.size": 10,

    "axes.labelsize": 11,
    "axes.titlesize": 11,

    "xtick.labelsize": 11,
    "ytick.labelsize": 11,

    "pdf.fonttype": 42,
    "ps.fonttype": 42
})


# =========================================================
# Global configuration
# =========================================================

# ---------------------------------------------------------
# Cohort color palette
# ---------------------------------------------------------

colores = {

    "Di@bet.es": "#4682B4",

    "Hortega": "#B0BEC5"
}

# ---------------------------------------------------------
# Gene configuration
# ---------------------------------------------------------
# Add or remove genes as needed

GENES = [

    {
        "gen": "ACAD10",
        "coverage_column": "cobertura_media"
    },

    {
        "gen": "ACSM3",
        "coverage_column": "cobertura_media"
    },

    {
        "gen": "AFMID",
        "coverage_column": "cobertura_media"
    },

    {
        "gen": "GCKR",
        "coverage_column": "cobertura_media"
    },

    {
        "gen": "NAT1",
        "coverage_column": "cobertura_media"
    }
]

width = 0.35

poblaciones = list(
    colores.keys()
)


# =========================================================
# Figure initialization
# =========================================================

fig = plt.figure(
    figsize=(16, 10)
)

gs = gridspec.GridSpec(
    3,
    4,

    figure=fig,

    hspace=0.55,
    wspace=0.35
)

# ---------------------------------------------------------
# Panel positions
# ---------------------------------------------------------

posiciones = [

    gs[0, 0:2],
    gs[0, 2:4],

    gs[1, 0:2],
    gs[1, 2:4],

    gs[2, 1:3],
]

axes = []


# =========================================================
# Multi-panel barplot generation
# =========================================================

for i, info in enumerate(GENES):

    ax = fig.add_subplot(
        posiciones[i]
    )

    axes.append(ax)

    # -----------------------------------------------------
    # Gene configuration
    # -----------------------------------------------------

    gen = info["gen"]

    col_cob = info["coverage_column"]

    infile = (
        f"nivel3_cobertura_media_por_amplicon_"
        f"{gen}.tsv"
    )

    # -----------------------------------------------------
    # Data loading
    # -----------------------------------------------------

    df = pd.read_csv(
        infile,
        sep="\t"
    )

    amplicones = sorted(
        df["amplicon"].unique()
    )

    x = np.arange(
        len(amplicones)
    )

    # -----------------------------------------------------
    # Cohort-specific barplots
    # -----------------------------------------------------

    for j, pob in enumerate(poblaciones):

        datos = (

            df[
                df["poblacion"] == pob
            ]

            .set_index("amplicon")

            .loc[amplicones]
        )

        ax.bar(

            x - width / 2 + j * width,

            datos[col_cob],

            width=width,

            color=colores[pob],

            alpha=0.8,

            edgecolor="black",

            linewidth=0.5,

            label=pob
        )

    # -----------------------------------------------------
    # Axis formatting
    # -----------------------------------------------------

    ax.set_xticks(x)

    ax.set_xticklabels(
        "",
        rotation=0,
        ha="center"
    )

    ax.set_ylabel(
        "Mean coverage (X)",
        fontsize=12
    )

    ax.set_xlabel(
        f"{gen} ({len(amplicones)} amplicons)",

        fontstyle="italic",

        fontsize=11
    )

    ax.set_title("")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # -----------------------------------------------------
    # Thousands separator formatting
    # -----------------------------------------------------

    ax.yaxis.set_major_formatter(

        FuncFormatter(
            lambda y, _:
            f"{int(y):,}".replace(",", "\u2009")
        )
    )

    # -----------------------------------------------------
    # Panel label
    # -----------------------------------------------------

    letra = f"({chr(ord('a') + i)})"

    ax.text(

        -0.13,
        1.20,

        letra,

        transform=ax.transAxes,

        fontsize=14,

        fontweight="bold",

        va="top",

        ha="left"
    )


# =========================================================
# Shared legend
# =========================================================

handles = [

    plt.Rectangle(
        (0, 0),
        1,
        1,

        color=colores[p],

        alpha=0.8,

        edgecolor="black",

        linewidth=0.5
    )

    for p in poblaciones
]

axes[-1].legend(

    handles,

    poblaciones,

    loc="center left",

    bbox_to_anchor=(1.02, 0.5),

    ncol=1,

    frameon=False,

    fontsize=14,

    handlelength=2.5,

    handleheight=1.5,

    labelspacing=0.6,

    borderpad=0.6
)


# =========================================================
# Figure export
# =========================================================

OUTPREFIX = (
    "nivel3_barplot_cobertura_multigen"
)

plt.savefig(
    f"{OUTPREFIX}.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    f"{OUTPREFIX}.pdf",
    bbox_inches="tight"
)

plt.close()


# =========================================================
# Final message
# =========================================================

print(
    "Integrated figure generated successfully"
)
