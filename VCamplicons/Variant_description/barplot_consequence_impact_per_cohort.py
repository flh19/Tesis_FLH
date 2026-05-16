#!/usr/bin/env python3

"""
============================================================
VARIANT CONSEQUENCE AND IMPACT COMPARISON PLOT
============================================================

Description:
    This script generates a two-panel horizontal stacked barplot
    comparing the distribution of variant consequences and impact
    categories between two cohorts.

    The plot summarizes:
        - Functional consequence annotations
        - Variant impact classes (HIGH, MODERATE, LOW, MODIFIER)

    Variant annotations are expected to come from VEP-like outputs.

Main features:
    - Consequence label simplification
    - Optional variant deduplication
    - Publication-quality multipanel figure
    - PNG and PDF export

------------------------------------------------------------
INPUT FILE REQUIREMENTS
------------------------------------------------------------

Input files must be tab-separated tables containing at least
the following columns:

    - ID
    - Consecuencia
    - Impacto

Example:

    ID              Consecuencia              Impacto
    rs12345         missense_variant          MODERATE
    rs67890         intron_variant            MODIFIER

------------------------------------------------------------
USAGE
------------------------------------------------------------

Basic execution:

    python script.py cohortA.tsv cohortB.tsv

Custom labels:

    python script.py cohortA.tsv cohortB.tsv \
        --labelA "Cases" \
        --labelB "Controls"

Enable variant deduplication:

    python script.py cohortA.tsv cohortB.tsv --dedup

Custom output directory:

    python script.py cohortA.tsv cohortB.tsv \
        --outdir results/

Custom colors:

    python script.py cohortA.tsv cohortB.tsv \
        --colors "HIGH=#D88B8B,MODERATE=#AED8C8,LOW=#F9CEAE,MODIFIER=#CAD9FB"

------------------------------------------------------------
OUTPUT
------------------------------------------------------------

The script generates:

    comparacion_consecuencia_impacto.png
    comparacion_consecuencia_impacto.pdf

inside the selected output directory.

------------------------------------------------------------
DEPENDENCIES
------------------------------------------------------------

Required Python packages:

    - pandas
    - matplotlib

Install with:

    pip install pandas matplotlib

============================================================
"""

# ============================================================
# IMPORTS
# ============================================================

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# GLOBAL CONFIGURATION
# ============================================================

# Ordered impact categories used in the stacked bars
IMPACT_ORDER = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

# Spanish labels displayed in the legend
IMPACT_LABELS_ES = {
    "HIGH": "ALTO",
    "MODERATE": "MODERADO",
    "LOW": "BAJO",
    "MODIFIER": "MODIFICADOR",
}

# Mapping between raw VEP consequence annotations
# and simplified labels used in the plot
CONSEQUENCE_MAP = {
    "3_prime_UTR_variant": "3' UTR",
    "5_prime_UTR_variant": "5' UTR",
    "frameshift_variant": "Frameshift",
    "frameshift_variant&splice_region_variant": "Frameshift",
    "inframe_deletion": "Deleción inframe",
    "inframe_insertion": "Inserción inframe",
    "intron_variant": "Intrónica",
    "missense_variant": "Missense",
    "missense_variant&splice_region_variant": "Missense",
    "splice_acceptor_variant": "Aceptora splicing",
    "splice_donor_variant": "Donadora splicing",
    "splice_region_variant&intron_variant": "Splicing\nno canónico",
    "splice_region_variant&synonymous_variant": "Splicing\nno canónico",
    "stop_gained": "Ganancia stop",
    "synonymous_variant": "Sinónima",
}


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def parse_colors(colors_str: str | None) -> dict[str, str]:
    """
    Parse a comma-separated string of colors into a dictionary.

    Example input:
        "HIGH=#D88B8B,MODERATE=#AED8C8"

    Returns:
        {
            "HIGH": "#D88B8B",
            "MODERATE": "#AED8C8"
        }
    """

    if not colors_str:
        return {}

    out = {}

    for chunk in colors_str.split(","):
        k, v = chunk.split("=")
        out[k.strip().upper()] = v.strip()

    return out


def italicize_english_words(label: str) -> str:
    """
    Apply matplotlib LaTeX italics formatting to selected
    English biological terms appearing in consequence labels.
    """

    if label is None:
        return ""

    text = str(label)

    substitutions = [
        (r"\bUTR\b", lambda m: r"$\mathit{UTR}$"),
        (r"\binframe\b", lambda m: r"$\mathit{inframe}$"),
        (r"\bFrameshift\b", lambda m: r"$\mathit{Frameshift}$"),
        (r"\bMissense\b", lambda m: r"$\mathit{Missense}$"),
        (r"\bSplicing\b", lambda m: r"$\mathit{Splicing}$"),
        (r"\bstart\b", lambda m: r"$\mathit{start}$"),
        (r"\bstop\b", lambda m: r"$\mathit{stop}$"),
    ]

    for pattern, replacement in substitutions:
        text = re.sub(pattern, replacement, text)

    return text


def map_consequence(raw: str) -> str:
    """
    Convert raw consequence annotations into simplified labels.
    """

    if pd.isna(raw):
        return ""

    return CONSEQUENCE_MAP.get(raw, raw)


# ============================================================
# DATA LOADING AND PREPROCESSING
# ============================================================

def load_table(path: Path, dedup: bool = False) -> pd.DataFrame:
    """
    Load and preprocess the input variant table.

    Steps:
    1. Read tab-separated file
    2. Standardize impact labels
    3. Translate and format consequence annotations
    4. Optionally remove duplicated variant IDs
    """

    df = pd.read_csv(path, sep="\t", dtype=str)

    # Standardize impact labels
    df["Impacto"] = df["Impacto"].str.upper().str.strip()

    # Translate and format consequence annotations
    df["Consecuencia_es"] = (
        df["Consecuencia"]
        .apply(map_consequence)
        .apply(italicize_english_words)
    )

    # Remove duplicated variants if requested
    if dedup:
        df = df.drop_duplicates(subset=["ID"])

    return df


# ============================================================
# COUNT TABLE PREPARATION
# ============================================================

def prepare_counts(df: pd.DataFrame, top_n: int = 14) -> pd.DataFrame:
    """
    Generate a contingency table containing the number of variants
    per consequence category and impact class.
    """

    counts = (
        df.groupby(["Consecuencia_es", "Impacto"])
        .size()
        .unstack(fill_value=0)
    )

    # Ensure all impact categories are present
    for impact in IMPACT_ORDER:
        if impact not in counts.columns:
            counts[impact] = 0

    # Reorder columns consistently
    counts = counts[IMPACT_ORDER]

    # Calculate total variants per consequence category
    counts["TOTAL"] = counts.sum(axis=1)

    # Keep only the top categories
    counts = counts.sort_values("TOTAL", ascending=False).head(top_n)

    return counts


# ============================================================
# MULTIPANEL PLOT
# ============================================================

def plot_multipanel(
    counts_A: pd.DataFrame,
    counts_B: pd.DataFrame,
    label_A: str,
    label_B: str,
    outbase: Path,
    colors: dict[str, str],
    dpi: int,
    xpad_frac: float,
    bar_height: float = 0.95,
    y_gap: float = 2.2,
):
    """
    Generate a two-panel horizontal stacked barplot comparing
    variant consequences and impact categories between cohorts.
    """

    # --------------------------------------------------------
    # GLOBAL TYPOGRAPHY SETTINGS
    # --------------------------------------------------------

    plt.rcParams.update({
        "font.size": 13,
        "axes.titlesize": 18,
        "axes.labelsize": 18,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 18,
        "legend.title_fontsize": 18,
    })

    # --------------------------------------------------------
    # ENSURE CONSISTENT CATEGORY ORDER BETWEEN PANELS
    # --------------------------------------------------------

    order = counts_A.index
    counts_B = counts_B.reindex(order).fillna(0)

    # --------------------------------------------------------
    # FIGURE DIMENSIONS
    # --------------------------------------------------------

    max_x = max(counts_A["TOTAL"].max(), counts_B["TOTAL"].max())
    n = len(order)

    fig_h = max(8.0, 0.65 * n + 2.5)

    fig, axes = plt.subplots(
        ncols=2,
        figsize=(16, fig_h),
        sharey=True,
    )

    # --------------------------------------------------------
    # PANEL GENERATION
    # --------------------------------------------------------

    for ax, counts, title in zip(
        axes,
        [counts_A, counts_B],
        [label_A, label_B],
    ):

        y = [i * y_gap for i in range(n)]

        # Cumulative positions for stacked bars
        cum = pd.Series(0, index=counts.index)

        # ----------------------------------------------------
        # STACKED HORIZONTAL BARS
        # ----------------------------------------------------

        for impact in IMPACT_ORDER:

            ax.barh(
                y,
                counts[impact].values,
                left=cum.values,
                height=bar_height,
                label=IMPACT_LABELS_ES[impact],
                color=colors.get(impact),
                edgecolor="black",
                linewidth=0.8,
            )

            cum += counts[impact]

        # ----------------------------------------------------
        # TOTAL COUNTS AT BAR ENDS
        # ----------------------------------------------------

        for yi, total in zip(y, cum.values):

            ax.text(
                total + max_x * 0.01,
                yi,
                f"{int(total)}",
                va="center",
                ha="left",
                fontsize=17,
            )

        # ----------------------------------------------------
        # PANEL FORMATTING
        # ----------------------------------------------------

        ax.set_title(title, fontweight="bold", pad=10)

        ax.set_xlim(0, max_x * (1 + xpad_frac))

        ax.xaxis.grid(True, linestyle="--", alpha=0.35)

        ax.set_axisbelow(True)

        ax.invert_yaxis()

    # --------------------------------------------------------
    # SHARED Y-AXIS LABELS
    # --------------------------------------------------------

    axes[0].set_yticks([i * y_gap for i in range(n)])
    axes[0].set_yticklabels(order)

    # --------------------------------------------------------
    # X-AXIS LABELS
    # --------------------------------------------------------

    axes[0].set_xlabel("Número de variantes")
    axes[1].set_xlabel("Número de variantes")

    # --------------------------------------------------------
    # SHARED LEGEND
    # --------------------------------------------------------

    handles, labels = axes[0].get_legend_handles_labels()

    legend = fig.legend(
        handles,
        labels,
        title="Impacto",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=4,
        frameon=True,
    )

    legend.get_frame().set_edgecolor("black")
    legend.get_frame().set_linewidth(0.8)
    legend.get_frame().set_facecolor("white")

    # --------------------------------------------------------
    # LAYOUT ADJUSTMENT AND EXPORT
    # --------------------------------------------------------

    fig.tight_layout(rect=[0, 0.08, 1, 1])

    outbase.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        outbase.with_suffix(".png"),
        dpi=dpi,
        bbox_inches="tight"
    )

    fig.savefig(
        outbase.with_suffix(".pdf"),
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():

    # --------------------------------------------------------
    # ARGUMENT PARSER
    # --------------------------------------------------------

    parser = argparse.ArgumentParser()

    parser.add_argument("cohorte_A")
    parser.add_argument("cohorte_B")

    parser.add_argument("--labelA", default="Cohorte A")
    parser.add_argument("--labelB", default="Cohorte B")

    parser.add_argument("-o", "--outdir", default="plots")

    parser.add_argument("--dpi", type=int, default=300)

    parser.add_argument("--dedup", action="store_true")

    parser.add_argument("--xpad", type=float, default=0.15)

    parser.add_argument(
        "--colors",
        default=(
            "HIGH=#D88B8B,"
            "MODERATE=#AED8C8,"
            "LOW=#F9CEAE,"
            "MODIFIER=#CAD9FB"
        ),
    )

    args = parser.parse_args()

    # --------------------------------------------------------
    # COLOR PARSING
    # --------------------------------------------------------

    colors = parse_colors(args.colors)

    # --------------------------------------------------------
    # LOAD INPUT TABLES
    # --------------------------------------------------------

    dfA = load_table(
        Path(args.cohorte_A),
        dedup=args.dedup
    )

    dfB = load_table(
        Path(args.cohorte_B),
        dedup=args.dedup
    )

    # --------------------------------------------------------
    # PREPARE COUNT TABLES
    # --------------------------------------------------------

    counts_A = prepare_counts(dfA)

    counts_B = prepare_counts(dfB)

    # --------------------------------------------------------
    # OUTPUT FILE PREFIX
    # --------------------------------------------------------

    outbase = (
        Path(args.outdir)
        / "comparacion_consecuencia_impacto"
    )

    # --------------------------------------------------------
    # GENERATE FIGURE
    # --------------------------------------------------------

    plot_multipanel(
        counts_A,
        counts_B,
        args.labelA,
        args.labelB,
        outbase,
        colors,
        args.dpi,
        args.xpad,
    )

    # --------------------------------------------------------
    # SUCCESS MESSAGE
    # --------------------------------------------------------

    print(f"[OK] {outbase}.png / .pdf")


# ============================================================
# SCRIPT ENTRY POINT
# ============================================================

if __name__ == "__main__":
    main()
