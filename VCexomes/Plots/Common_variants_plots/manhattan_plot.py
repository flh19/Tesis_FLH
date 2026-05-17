"""
============================================================
Multi-Phenotype Manhattan Plot Generator
============================================================

This script generates a 2x2 panel Manhattan-style visualization
for multiple phenotypes using REGENIE association results.

"""

# ============================================================
# IMPORTS
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

import warnings
warnings.filterwarnings("ignore")


# ============================================================
# INPUT FILE CONFIGURATION
# ============================================================

# Mapping between phenotype names and input REGENIE files
files = {
    "Glucosa": "glucosa_10_unificado.regenie",
    "HOMA-IR": "HOMA-IR_10_unificado.regenie",
    "Insulina": "insulina_10_unificado.regenie",
    "DM2": "DM2_10_unificado.regenie"
}


# ============================================================
# PHENOTYPE COLORS
# ============================================================

# Each phenotype uses:
# - A darker color for odd chromosomes
# - A lighter color for even chromosomes

phenotype_colors = {
    "Glucosa": ("#5B7C99", "#B8C7D6"),
    "HOMA-IR": ("#B07D62", "#D9B8A6"),
    "Insulina": ("#7A9E7E", "#BDD3BF"),
    "DM2": ("#8A7CA8", "#C9C1D9")
}


# ============================================================
# SIGNIFICANCE THRESHOLDS
# ============================================================

# P-value threshold used for each phenotype

p_thresholds = {
    "Glucosa": 6.18e-7,
    "HOMA-IR": 6.18e-7,
    "Insulina": 6.18e-7,
    "DM2": 6.18e-7
}


# ============================================================
# GENETIC MODEL → MARKER SHAPE
# ============================================================

model_shape_map = {
    "ADD": "o",   # Additive model
    "DOM": "^",   # Dominant model
    "REC": "s"    # Recessive model
}


# ============================================================
# GLOBAL STYLE CONFIGURATION
# ============================================================

mpl.rcParams.update({

    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial"],

    "axes.linewidth": 1.4,

    "xtick.major.width": 1.4,
    "ytick.major.width": 1.4,

    "xtick.major.size": 5,
    "ytick.major.size": 5,

    "xtick.minor.visible": False,
    "ytick.minor.visible": False,

    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

BACKGROUND = "#FFFFFF"
GRID_COLOR = "#E8E8E8"
SPINE_COLOR = "#333333"
LABEL_COLOR = "#1A1A1A"
THRESH_COLOR = "#C0392B"


# ============================================================
# CUSTOM AXIS FORMATTER
# ============================================================

def comma_formatter(x, pos):
    """
    Format tick labels using commas instead of dots
    for decimal separation.

    Example:
        2.5 -> 2,5
    """

    if float(x).is_integer():
        return f"{int(x)}"

    return f"{x:.1f}".replace(".", ",")


# ============================================================
# DATA PREPARATION
# ============================================================

def prepare_equidistant(df):
    """
    Prepare genomic data using an equidistant chromosome layout.

    Steps:
    - Extract p-values
    - Compute -log10(p)
    - Standardize chromosome/position columns
    - Generate evenly spaced cumulative positions

    Parameters
    ----------
    df : pandas.DataFrame
        Input association results dataframe.

    Returns
    -------
    tuple
        (
            processed_dataframe,
            chromosome_tick_positions,
            chromosome_labels,
            chromosome_order
        )
    """

    df = df.copy()

    # --------------------------------------------------------
    # Use combined p-values
    # --------------------------------------------------------

    df["PVAL"] = pd.to_numeric(df["p_comb"], errors="coerce")

    # Remove invalid p-values
    df = df[df["PVAL"].notna()]
    df = df[df["PVAL"] > 0]

    # Compute -log10(p)
    df["log10P"] = -np.log10(df["PVAL"])

    # --------------------------------------------------------
    # Standardize column names
    # --------------------------------------------------------

    df["CHROM"] = df["Chr"].astype(str)
    df["GENPOS"] = pd.to_numeric(df["Pos"], errors="coerce")

    # Best-performing genetic model
    df["MODEL"] = df["best_model"]

    # Sort chromosomes numerically
    chrom_order = sorted(df["CHROM"].unique(), key=lambda x: int(x))

    # Distance between chromosomes
    spacing = 1_000_000

    pos_list = []
    ticks = []
    labels = []

    current_x = 0

    # --------------------------------------------------------
    # Build equidistant chromosome coordinates
    # --------------------------------------------------------

    for chrom in chrom_order:

        sub = df[df["CHROM"] == chrom].copy().sort_values("GENPOS")

        if sub.empty:
            continue

        # If chromosome contains only one variant
        if len(sub) == 1:

            sub["pos_cum"] = current_x + spacing / 2

        else:

            # Spread variants evenly within chromosome space
            sub["pos_cum"] = np.linspace(
                current_x,
                current_x + spacing,
                len(sub)
            )

        pos_list.append(sub)

        # Tick position in chromosome center
        ticks.append(current_x + spacing / 2)

        labels.append(chrom)

        # Add spacing before next chromosome
        current_x += spacing + spacing * 0.2

    return pd.concat(pos_list), ticks, labels, chrom_order


# ============================================================
# LOAD DATA
# ============================================================

prepared = {}

for pheno, file_name in files.items():

    print(f"Loading {pheno} from {file_name}...")

    df = pd.read_csv(file_name, sep="\t")

    prepared[pheno] = prepare_equidistant(df)


# ============================================================
# CREATE FIGURE
# ============================================================

fig, axes = plt.subplots(

    2, 2,

    figsize=(26, 11),

    sharex=True,

    gridspec_kw={
        "hspace": 0.32,
        "wspace": 0.10
    },

    facecolor=BACKGROUND,
)

fig.patch.set_facecolor(BACKGROUND)

axes = axes.flatten()

# Panel order
order = ["DM2", "Glucosa", "Insulina", "HOMA-IR"]


# ============================================================
# Y-AXIS CONSISTENCY
# ============================================================

max_y = max(prepared[p][0]["log10P"].max() for p in order)

y_top = np.ceil(max_y) + 1.5


# ============================================================
# P-VALUE LABEL FORMATTER
# ============================================================

def format_pval(p):
    """
    Format p-values for plot legends using scientific notation.
    """

    mantissa, exponent = f"{p:.2e}".split("e")

    mantissa_formatted = mantissa.rstrip('0').rstrip('.')
    mantissa_formatted = mantissa_formatted.replace(".", ",")

    return rf"$p = {mantissa_formatted} \times 10^{{{int(exponent)}}}$"


# ============================================================
# PANEL PLOTTING
# ============================================================

for ax_idx, (ax, pheno) in enumerate(zip(axes, order)):

    df, ticks, labels, chrom_order = prepared[pheno]

    color_dark, color_light = phenotype_colors[pheno]

    pheno_threshold = p_thresholds[pheno]

    logp_threshold = -np.log10(pheno_threshold)

    ax.set_facecolor(BACKGROUND)

    # --------------------------------------------------------
    # Horizontal grid lines
    # --------------------------------------------------------

    y_ticks_grid = np.arange(0, y_top + 1, 2)

    for yg in y_ticks_grid:

        ax.axhline(
            yg,
            color=GRID_COLOR,
            linewidth=0.6,
            zorder=0
        )

    # --------------------------------------------------------
    # Alternating chromosome background bands
    # --------------------------------------------------------

    for i, chrom in enumerate(chrom_order):

        sub = df[df["CHROM"] == chrom]

        if not sub.empty and i % 2 == 1:

            x0 = sub["pos_cum"].min()
            x1 = sub["pos_cum"].max()

            ax.axvspan(
                x0,
                x1,
                color="#F5F5F5",
                alpha=0.6,
                zorder=0
            )

    # --------------------------------------------------------
    # Plot points chromosome by chromosome
    # --------------------------------------------------------

    for i, chrom in enumerate(chrom_order):

        sub = df[df["CHROM"] == chrom]

        if sub.empty:
            continue

        color = color_dark if i % 2 == 0 else color_light

        # =====================================================
        # NON-SIGNIFICANT VARIANTS
        # =====================================================

        mask_ns = sub["log10P"] < logp_threshold

        ax.scatter(
            sub.loc[mask_ns, "pos_cum"],
            sub.loc[mask_ns, "log10P"],
            c=color,
            s=14,
            alpha=0.65,
            linewidths=0,
            zorder=2,
            rasterized=True
        )

        # =====================================================
        # SIGNIFICANT VARIANTS
        # =====================================================

        mask_sig = sub["log10P"] >= logp_threshold

        sub_sig = sub.loc[mask_sig]

        # Different marker for each genetic model
        for model_type, marker in model_shape_map.items():

            sub_model = sub_sig[sub_sig["MODEL"] == model_type]

            if not sub_model.empty:

                ax.scatter(
                    sub_model["pos_cum"],
                    sub_model["log10P"],
                    c=color,
                    s=135,
                    alpha=1.0,
                    marker=marker,
                    linewidths=0.8,
                    edgecolors="#333333",
                    zorder=4,
                    rasterized=True
                )

    # --------------------------------------------------------
    # Significance threshold line
    # --------------------------------------------------------

    ax.axhline(
        y=logp_threshold,
        color=THRESH_COLOR,
        linestyle="--",
        linewidth=1.4,
        dashes=(5, 4),
        zorder=3,
        alpha=0.85,
        label=format_pval(pheno_threshold)
    )

    # --------------------------------------------------------
    # Individual panel legend
    # --------------------------------------------------------

    ax.legend(
        loc="upper right",
        bbox_to_anchor=(1.0, 1.05),
        frameon=False,
        fontsize=20,
        handlelength=1.5
    )

    # --------------------------------------------------------
    # Phenotype label
    # --------------------------------------------------------

    ax.text(
        0.02,
        0.94,
        pheno,
        transform=ax.transAxes,
        fontsize=24,
        fontweight="bold",
        color=LABEL_COLOR,
        va="top",
    )

    # --------------------------------------------------------
    # Panel letters: (a), (b), (c), (d)
    # --------------------------------------------------------

    panel_letter = chr(97 + ax_idx)

    x_pos = -0.10 if ax_idx % 2 == 0 else -0.08

    ax.text(
        x_pos,
        1.08,
        f"({panel_letter})",
        transform=ax.transAxes,
        fontsize=22,
        fontweight="bold",
        color=LABEL_COLOR,
        va="bottom"
    )

    # --------------------------------------------------------
    # Axis formatting
    # --------------------------------------------------------

    ax.set_ylim(0, y_top)

    ax.yaxis.set_major_formatter(FuncFormatter(comma_formatter))

    # Left panels only
    if ax_idx % 2 == 0:

        ax.set_ylabel(
            r"$-\log_{10}(p)$",
            fontsize=24,
            color=LABEL_COLOR,
            labelpad=8
        )

    else:

        ax.set_ylabel("")

    # Remove top/right borders
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # Customize remaining borders
    for spine in ["left", "bottom"]:

        ax.spines[spine].set_color(SPINE_COLOR)
        ax.spines[spine].set_linewidth(1.4)

    ax.tick_params(
        axis="both",
        which="major",
        labelsize=22,
        width=1.4,
        color=SPINE_COLOR,
        labelcolor=LABEL_COLOR,
    )

    # Hide upper x-axis labels
    if ax_idx < 2:
        ax.tick_params(axis="x", which="major", labelbottom=False)


# ============================================================
# X-AXIS LABELS (BOTTOM PANELS ONLY)
# ============================================================

for ax in axes[2:]:

    ax.tick_params(
        axis="x",
        which="major",
        labelbottom=True,
        labelsize=18
    )

    ax.set_xticks(ticks)

    ax.set_xticklabels(
        labels,
        rotation=0,
        color=LABEL_COLOR
    )

    ax.set_xlabel(
        "Cromosoma",
        fontsize=24,
        color=LABEL_COLOR,
        labelpad=8
    )


# ============================================================
# GLOBAL LEGEND — GENETIC MODELS
# ============================================================

legend_shapes = [

    Line2D(
        [0], [0],
        marker='o',
        linestyle='None',
        markerfacecolor='#666666',
        markeredgecolor='#333333',
        markersize=17,
        label='ADI'
    ),

    Line2D(
        [0], [0],
        marker='^',
        linestyle='None',
        markerfacecolor='#666666',
        markeredgecolor='#333333',
        markersize=17,
        label='DOM'
    ),

    Line2D(
        [0], [0],
        marker='s',
        linestyle='None',
        markerfacecolor='#666666',
        markeredgecolor='#333333',
        markersize=17,
        label='REC'
    ),
]

leg = fig.legend(
    handles=legend_shapes,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.01),
    ncol=3,
    frameon=True,
    fontsize=22,
)

leg.get_frame().set_edgecolor("black")
leg.get_frame().set_linewidth(0.8)
leg.get_frame().set_facecolor("white")


# ============================================================
# SAVE OUTPUT
# ============================================================

plt.tight_layout(rect=[0, 0.05, 1, 1])

# PNG export
plt.savefig(
    "manhattan_common_variants.png",
    dpi=300,
    bbox_inches="tight",
    facecolor=BACKGROUND
)

# PDF export
plt.savefig(
    "manhattan_common_variants.pdf",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.02,
    facecolor=BACKGROUND
)

print("✓ Saved: manhattan_common_variants.png / .pdf")

plt.show()
