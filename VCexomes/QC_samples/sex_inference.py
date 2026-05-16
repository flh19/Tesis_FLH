#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =========================================================
# EXOME SEX INFERENCE PLOTTING SCRIPT
# =========================================================
#
# Description
# -----------
# This script performs sex inference from exome sequencing
# coverage data using X and Y chromosome depth ratios
# relative to autosomal coverage.
#
# The workflow:
#
#   1. Loads per-sample coverage statistics
#   2. Calculates X/Autosome and Y/Autosome ratios
#   3. Infers biological sex based on empirical thresholds
#   4. Generates a publication-quality JointGrid plot
#      containing:
#
#         - Central scatter plot
#         - Marginal density distributions
#         - External legend with sample counts
#
# The final figure is exported as PNG and PDF.
#
#
# Input file
# ----------
# TSV file containing coverage statistics with at least:
#
#   - sample_id
#   - depth_mean
#   - X_depth_mean
#   - Y_depth_mean
#
#
# Sex inference criteria
# ----------------------
# Mujer:
#   X/Autosome >= 1.0
#   Y/Autosome == 0
#
# Hombre:
#   X/Autosome between 0.5 and 0.8
#   Y/Autosome >= 0.3
#
# Otherwise:
#   Indeterminado
#
#
# Output
# ------
# The script generates:
#
#   - inferencia_sexo_exoma.png
#   - inferencia_sexo_exoma.pdf
#
#
# Formatting
# ----------
# - Decimal separator uses commas
# - Thousands are formatted using thin spaces
#   Example:
#       1000 -> 1 000
#
# =========================================================


# =========================================================
# IMPORTS
# =========================================================
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


# =========================================================
# THOUSANDS FORMATTER
# =========================================================
def fmt_miles(x):
    """
    Format integers using thin-space thousands separators.

    Example:
        1000 -> 1 000
    """
    return f"{int(x):,}".replace(",", "\u2009")


# =========================================================
# 1. LOAD DATA
# =========================================================
# Read TSV file containing coverage metrics
df = pd.read_csv(
    "/content/drive/MyDrive/share_ugdg/bioinfo/scripts/FLH_Scripts/TESIS/EXOMAS/QC/relatedness_output.samples.tsv",
    sep="\t"
)

# Set sample identifier as dataframe index
df.set_index("sample_id", inplace=True)


# =========================================================
# 2. AUTOSOMAL DEPTH
# =========================================================
# Use mean autosomal depth as reference coverage
df["Auto_median"] = df["depth_mean"]


# =========================================================
# 3. CALCULATE COVERAGE RATIOS
# =========================================================
# X chromosome normalized depth
df["XAutoRatio"] = df["X_depth_mean"] / df["Auto_median"]

# Y chromosome normalized depth
df["YAutoRatio"] = df["Y_depth_mean"] / df["Auto_median"]


# =========================================================
# 4. SEX INFERENCE
# =========================================================
# Initialize all samples as undetermined
df["sexo_inferido"] = "Indeterminado"

# Female inference criteria
df.loc[
    (df["XAutoRatio"] >= 1.0) & (df["YAutoRatio"] == 0),
    "sexo_inferido"
] = "Mujer"

# Male inference criteria
df.loc[
    (df["XAutoRatio"].between(0.5, 0.8)) & (df["YAutoRatio"] >= 0.3),
    "sexo_inferido"
] = "Hombre"


# =========================================================
# 5. COLOR PALETTE
# =========================================================
# Define consistent colors for each inferred sex category
sex_palette = {
    "Mujer": "#FF7F50",
    "Hombre": "#008B8B",
    "Indeterminado": "#808080"
}


# =========================================================
# 6. AXIS FORMATTERS
# =========================================================
def coma_decimal(x, pos):
    """
    Format decimal numbers using commas instead of dots.

    Example:
        1.5 -> 1,5
    """
    return f"{x:.1f}".replace(".", ",")


def formatter_y(y, pos):
    """
    Custom formatter for Y-axis values.

    - Hides negative labels
    - Uses comma decimal separator
    """
    if y < 0:
        return ""

    return f"{y:.1f}".replace(".", ",")


# =========================================================
# 7. JOINTGRID STRUCTURE
# =========================================================
# Create JointGrid layout:
#   - Central scatter plot
#   - Marginal density plots
g = sns.JointGrid(
    data=df,
    x="XAutoRatio",
    y="YAutoRatio",
    height=4.8,
    ratio=4,
    space=0
)

# ---------------------------------------------------------
# REAL FIGURE SIZE (WIDTH, HEIGHT)
# ---------------------------------------------------------
g.figure.set_size_inches(7.7, 3.8)


# =========================================================
# 8. CENTRAL SCATTER PLOT
# =========================================================
sns.scatterplot(
    data=df,
    x="XAutoRatio",
    y="YAutoRatio",
    hue="sexo_inferido",
    palette=sex_palette,
    alpha=0.7,
    s=35,
    ax=g.ax_joint,
    legend=False
)


# =========================================================
# 9. MARGINAL DENSITY DISTRIBUTIONS
# =========================================================
# X-axis density distribution
sns.kdeplot(
    data=df,
    x="XAutoRatio",
    hue="sexo_inferido",
    palette=sex_palette,
    fill=True,
    alpha=0.4,
    common_norm=False,
    ax=g.ax_marg_x,
    legend=False
)

# Y-axis density distribution
sns.kdeplot(
    data=df,
    y="YAutoRatio",
    hue="sexo_inferido",
    palette=sex_palette,
    fill=True,
    alpha=0.4,
    common_norm=False,
    ax=g.ax_marg_y,
    legend=False
)


# =========================================================
# 10. Y-AXIS LIMITS
# =========================================================
# Prevent cutting upper points
ymax = df["YAutoRatio"].max()

g.ax_joint.set_ylim(-0.05, ymax * 1.05)


# =========================================================
# 11. LABELS AND AXIS FORMAT
# =========================================================
# Axis labels
g.ax_joint.set_xlabel("ratio X / Autosomas", fontsize=13)
g.ax_joint.set_ylabel("ratio Y / Autosomas", fontsize=13)

# Apply decimal comma formatting
g.ax_joint.xaxis.set_major_formatter(FuncFormatter(coma_decimal))
g.ax_joint.yaxis.set_major_formatter(FuncFormatter(formatter_y))

# Tick label size
g.ax_joint.tick_params(labelsize=13)


# =========================================================
# 12. EXTERNAL LEGEND WITH SAMPLE COUNTS
# =========================================================
# Count samples per inferred category
conteos = df["sexo_inferido"].value_counts()

# Custom legend markers
handles = [
    plt.Line2D(
        [0], [0],
        marker="o",
        color="w",
        markerfacecolor=sex_palette["Hombre"],
        markersize=8
    ),

    plt.Line2D(
        [0], [0],
        marker="o",
        color="w",
        markerfacecolor=sex_palette["Mujer"],
        markersize=8
    ),

    plt.Line2D(
        [0], [0],
        marker="o",
        color="w",
        markerfacecolor=sex_palette["Indeterminado"],
        markersize=8
    )
]

# Legend labels with formatted sample counts
labels = [
    f"Hombre (N = {fmt_miles(conteos.get('Hombre', 0))})",
    f"Mujer (N = {fmt_miles(conteos.get('Mujer', 0))})",
    f"Indeterminado (N = {fmt_miles(conteos.get('Indeterminado', 0))})"
]

# Leave space for external legend
g.figure.subplots_adjust(right=0.78)

# Create external legend
legend = g.figure.legend(
    handles,
    labels,
    loc="center left",
    bbox_to_anchor=(0.70, 0.22),
    frameon=True,
    fontsize=11
)

# Legend frame styling
legend.get_frame().set_edgecolor("black")
legend.get_frame().set_linewidth(1.0)
legend.get_frame().set_alpha(1.0)


# =========================================================
# 13. SAVE FIGURES
# =========================================================
# Export PNG version
g.figure.savefig(
    "/content/drive/MyDrive/share_ugdg/bioinfo/scripts/FLH_Scripts/TESIS/EXOMAS/QC/inferencia_sexo_exoma.png",
    dpi=600,
    bbox_inches="tight"
)

# Export PDF version
g.figure.savefig(
    "/content/drive/MyDrive/share_ugdg/bioinfo/scripts/FLH_Scripts/TESIS/EXOMAS/QC/inferencia_sexo_exoma.pdf",
    dpi=600,
    bbox_inches="tight"
)


# =========================================================
# 14. DISPLAY FIGURE
# =========================================================
plt.show()
