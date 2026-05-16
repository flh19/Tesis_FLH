#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
Variant consequence and impact visualization pipeline
===============================================================================

Description
-----------
This script generates publication-quality visualizations summarizing the
distribution of variant consequences and predicted functional impacts from
variant annotation datasets (e.g., Ensembl VEP annotations).

The generated figure combines:

    1. Horizontal stacked bar plots:
        - Distribution of the most frequent variant consequence categories.
        - Stratified by predicted functional impact.

    2. Pie chart:
        - Global distribution of impact categories.

Variant consequences are translated, grouped, and standardized to improve
visual interpretation for publication-quality figures.

Applications
------------
This workflow is commonly used in:
    - Whole-exome sequencing studies,
    - Variant annotation summaries,
    - Rare variant analyses,
    - Functional genomics,
    - Population genetics studies,
    - Doctoral theses and scientific manuscripts.

Input
-----
One or more tab-delimited annotation files containing at least:

    - Consecuencia
    - Impacto
    - ID

Expected annotation nomenclature follows Ensembl VEP conventions.

Output
------
For each input file, the script generates:

    - PNG figure
    - PDF figure

Output filename format:

    <input_name>_tipo_impacto.png
    <input_name>_tipo_impacto.pdf

Usage
-----
Run from command line:

    python plot_variant_consequences.py input.tsv

Multiple files:

    python plot_variant_consequences.py file1.tsv file2.tsv

Optional arguments:
    --outdir
    --dpi
    --dedup
    --xpad
    --colors

Dependencies
------------
    pandas
    matplotlib

Author
------
Prepared for genomic variant annotation visualization and sequencing studies.

===============================================================================
"""

#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# ============================================================
# Internal impact order (VEP hierarchy)
# ============================================================

IMPACT_ORDER = [
    "HIGH",
    "MODERATE",
    "LOW",
    "MODIFIER"
]

# ============================================================
# Spanish labels for legend
# ============================================================

IMPACT_LABELS_ES = {
    "HIGH": "ALTO",
    "MODERATE": "MODERADO",
    "LOW": "BAJO",
    "MODIFIER": "MODIFICADOR",
}

# ============================================================
# Consequence translation and grouping
# ============================================================

CONSEQUENCE_MAP = {

    # UTR variants
    "3_prime_UTR_variant": "3' UTR",
    "3_prime_UTR_variant&NMD_transcript_variant": "3' UTR",

    "5_prime_UTR_variant": "5' UTR",
    "5_prime_UTR_variant&NMD_transcript_variant": "5' UTR",

    # Coding sequence
    "coding_sequence_variant": "Secuencia codificante",
    "coding_sequence_variant&NMD_transcript_variant":
        "Secuencia codificante",

    # Downstream
    "downstream_gene_variant": "Downstream",

    # Frameshift
    "frameshift_variant": "Frameshift",
    "frameshift_variant&NMD_transcript_variant": "Frameshift",
    "frameshift_variant&splice_region_variant": "Frameshift",
    "frameshift_variant&splice_region_variant&NMD_transcript_variant":
        "Frameshift",
    "frameshift_variant&start_lost": "Frameshift",
    "frameshift_variant&start_lost&NMD_transcript_variant":
        "Frameshift",
    "frameshift_variant&stop_lost": "Frameshift",

    # Terminal codon
    "incomplete_terminal_codon_variant&coding_sequence_variant":
        "Codón terminal incompleto",

    # In-frame insertion
    "inframe_insertion&stop_retained_variant":
        "Inserción in-frame",

    # Intergenic
    "intergenic_variant": "Intergénica",

    # Intronic
    "intron_variant": "Intrónica",
    "intron_variant&NMD_transcript_variant": "Intrónica",
    "intron_variant&non_coding_transcript_variant":
        "Intrónica",

    # miRNA
    "mature_miRNA_variant": "miRNA maduro",

    # Missense
    "missense_variant": "Missense",
    "missense_variant&NMD_transcript_variant": "Missense",
    "missense_variant&splice_region_variant": "Missense",
    "missense_variant&splice_region_variant&NMD_transcript_variant":
        "Missense",

    # Non-coding transcript exon
    "non_coding_transcript_exon_variant": "Exón ARNnc",

    # Splice acceptor
    "splice_acceptor_variant": "Aceptora splicing",
    "splice_acceptor_variant&NMD_transcript_variant":
        "Aceptora splicing",
    "splice_acceptor_variant&non_coding_transcript_variant":
        "Aceptora splicing",

    # Non-canonical splicing
    "splice_donor_5th_base_variant&intron_variant":
        "Splicing no canónico",

    "splice_donor_5th_base_variant&intron_variant&NMD_transcript_variant":
        "Splicing no canónico",

    "splice_donor_5th_base_variant&intron_variant&non_coding_transcript_variant":
        "Splicing no canónico",

    "splice_donor_region_variant&intron_variant":
        "Splicing no canónico",

    "splice_donor_region_variant&intron_variant&NMD_transcript_variant":
        "Splicing no canónico",

    "splice_donor_region_variant&intron_variant&non_coding_transcript_variant":
        "Splicing no canónico",

    "splice_polypyrimidine_tract_variant&intron_variant":
        "Splicing no canónico",

    "splice_region_variant&intron_variant":
        "Splicing no canónico",

    "splice_region_variant&synonymous_variant":
        "Splicing no canónico",

    # Donor splice
    "splice_donor_variant": "Donadora splicing",
    "splice_donor_variant&NMD_transcript_variant":
        "Donadora splicing",

    # Start/stop variants
    "start_lost": "Pérdida start",
    "start_retained_variant": "Start retenido",

    "stop_gained": "Ganancia stop",
    "stop_lost": "Pérdida stop",

    "stop_retained_variant": "Stop retenido",

    # Synonymous
    "synonymous_variant": "Sinónima",

    # Upstream
    "upstream_gene_variant": "Upstream",
}

# ============================================================
# Utility functions
# ============================================================

def parse_colors(colors_str: str | None) -> dict[str, str]:
    """
    Parse custom color definitions from command line.
    """

    if not colors_str:
        return {}

    out: dict[str, str] = {}

    for chunk in colors_str.split(","):

        chunk = chunk.strip()

        if "=" not in chunk:
            raise ValueError(
                f"Formato inválido en --colors: {chunk}"
            )

        k, v = chunk.split("=", 1)

        out[k.strip().upper()] = v.strip()

    return out


def format_thousands(x):
    """
    Format integers using thin-space thousands separator.
    """
    return f"{int(x):,}".replace(",", "\u2009")


def italicize_english_words(label: str) -> str:
    """
    Italicize selected English biological terms using LaTeX formatting.
    """

    if label is None:
        return ""

    text = str(label)

    subs = [

        (r"\bUTR\b",
         lambda m: r"$\mathit{UTR}$"),

        (r"\bFrameshift\b",
         lambda m: r"$\mathit{Frameshift}$"),

        (r"\bMissense\b",
         lambda m: r"$\mathit{Missense}$"),

        (r"\bSplicing\b",
         lambda m: r"$\mathit{Splicing}$"),

        (r"\bsplicing\b",
         lambda m: r"$\mathit{splicing}$"),

        (r"\bUpstream\b",
         lambda m: r"$\mathit{Upstream}$"),

        (r"\bDownstream\b",
         lambda m: r"$\mathit{Downstream}$"),

        (r"\bstart\b",
         lambda m: r"$\mathit{start}$"),

        (r"\bstop\b",
         lambda m: r"$\mathit{stop}$"),
    ]

    for pat, repl in subs:
        text = re.sub(pat, repl, text)

    return text


def map_consequence(raw: str) -> str:
    """
    Translate and group VEP consequence terms.
    """

    s = "" if pd.isna(raw) else str(raw).strip()

    if s in CONSEQUENCE_MAP:
        return CONSEQUENCE_MAP[s]

    sl = s.lower()

    if "missense_variant" in sl:
        return "missense"

    if "frameshift_variant" in sl:
        return "Frameshift"

    if "splice" in sl:
        return "Splicing no canónico"

    return s


def load_table(path: Path, dedup: bool = False) -> pd.DataFrame:
    """
    Load and preprocess annotation table.
    """

    df = pd.read_csv(path, sep="\t", dtype=str)

    # Standardize impact labels
    df["Impacto"] = (
        df["Impacto"]
        .str.upper()
        .str.strip()
    )

    # Translate and stylize consequences
    df["Consecuencia_es"] = (
        df["Consecuencia"]
        .apply(map_consequence)
        .apply(italicize_english_words)
    )

    # Optional deduplication by variant ID
    if dedup:
        df = df.drop_duplicates(subset=["ID"])

    return df

# ============================================================
# Plotting function
# ============================================================

def make_plot(
    df: pd.DataFrame,
    title: str,
    outbase: Path,
    dpi: int,
    xpad_frac: float,
    colors: dict[str, str],
    bar_height: float = 0.55,
    y_gap: float = 1,
    bar_edge_width: float = 0.8,
) -> None:

    """
    Generate stacked consequence-impact plot with inset pie chart.
    """

    # Font configuration
    plt.rcParams.update({

        "font.size": 28,
        "axes.titlesize": 34,
        "axes.labelsize": 34,
        "xtick.labelsize": 28,
        "ytick.labelsize": 28,
        "legend.fontsize": 28,
        "legend.title_fontsize": 30,
    })

    # Count variants by consequence and impact
    counts = (
        df.groupby(["Consecuencia_es", "Impacto"])
        .size()
        .unstack(fill_value=0)
    )

    # Ensure all impact categories exist
    for imp in IMPACT_ORDER:

        if imp not in counts.columns:
            counts[imp] = 0

    counts = counts[IMPACT_ORDER]

    # Total variants per consequence
    counts["TOTAL"] = counts.sum(axis=1)

    # Keep top 13 categories
    counts = (
        counts
        .sort_values("TOTAL", ascending=True)
        .tail(13)
    )

    # ========================================================
    # Global impact proportions
    # ========================================================

    impact_totals = df["Impacto"].value_counts()

    impact_totals = (
        impact_totals
        .reindex(IMPACT_ORDER)
        .fillna(0)
    )

    impact_perc = (
        impact_totals /
        impact_totals.sum() * 100
    )

    # ========================================================
    # Figure size configuration
    # ========================================================

    n = len(counts)

    FIG_W = 25.0
    ROW_H = 0.75

    FIG_H = max(7.0, n * ROW_H + 3.0)

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

    # Y positions
    y = [i * y_gap for i in range(n - 1, -1, -1)]

    cum = pd.Series(0, index=counts.index)

    # ========================================================
    # Stacked horizontal bars
    # ========================================================

    for imp in IMPACT_ORDER:

        kwargs = {}

        if imp in colors:
            kwargs["color"] = colors[imp]

        ax.barh(
            y,
            counts[imp].values,
            left=cum.values,
            height=bar_height,
            label=IMPACT_LABELS_ES[imp],
            edgecolor="black",
            linewidth=bar_edge_width,
            **kwargs,
        )

        cum += counts[imp]

    # Axis labels and ticks
    ax.set_yticks(y)

    ax.set_yticklabels(
        counts.index,
        fontsize=30
    )

    ax.tick_params(axis="y", pad=10)

    ax.set_xlabel(
        "Número de variantes",
        fontsize=34
    )

    # Optional title
    if title:
        ax.set_title(
            title,
            fontweight="bold",
            pad=18
        )

    # Grid
    ax.xaxis.grid(
        True,
        linestyle="--",
        linewidth=0.9,
        alpha=0.35
    )

    ax.set_axisbelow(True)

    # Remove upper/right borders
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(top=False, right=False)

    ax.spines["left"].set_linewidth(1.2)

    # Thousands formatter
    ax.xaxis.set_major_formatter(
        FuncFormatter(
            lambda x, pos: format_thousands(x)
        )
    )

    # X-axis limits
    max_total = counts["TOTAL"].max()

    ax.set_xlim(
        0,
        max_total * (1 + xpad_frac) + 1
    )

    # ========================================================
    # Labels at bar ends
    # ========================================================

    for yi, total in zip(y, counts["TOTAL"]):

        ax.text(
            total + max_total * 0.012,
            yi,
            format_thousands(total),
            va="center",
            ha="left",
            fontsize=26,
        )

    # ========================================================
    # Legend
    # ========================================================

    leg = ax.legend(
        title="Impacto",
        loc="upper right",
        frameon=True,
        fontsize=28,
        title_fontsize=30,
    )

    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.8)

    # ========================================================
    # Inset pie chart
    # ========================================================

    ax_pie = ax.inset_axes([0.28, 0.34, 0.50, 0.62])

    pie_colors = [
        colors.get(imp, None)
        for imp in IMPACT_ORDER
    ]

    def autopct_es(pct):
        """
        Percentage formatter using decimal comma.
        """
        return f"{pct:.2f}%".replace(".", ",")

    wedges, texts, autotexts = ax_pie.pie(
        impact_totals,
        labels=None,
        colors=pie_colors,
        startangle=90,
        counterclock=False,
        autopct=autopct_es,
        pctdistance=0.68,
        textprops={"fontsize": 30},
        wedgeprops={
            "edgecolor": "black",
            "linewidth": 1
        },
    )

    # Move HIGH-impact percentage label outside
    idx_alto = IMPACT_ORDER.index("HIGH")

    autotexts[idx_alto].set_position((0, 1.08))
    autotexts[idx_alto].set_ha("center")
    autotexts[idx_alto].set_va("bottom")

    ax_pie.set_title(
        "Impacto",
        fontsize=30,
        pad=18
    )

    # ========================================================
    # Save figure
    # ========================================================

    fig.tight_layout(pad=1.5)

    outbase.parent.mkdir(
        parents=True,
        exist_ok=True
    )

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
# Main execution
# ============================================================

def main():

    # Command-line parser
    p = argparse.ArgumentParser()

    p.add_argument("inputs", nargs="+")

    p.add_argument(
        "-o",
        "--outdir",
        default="plots"
    )

    p.add_argument(
        "--dpi",
        type=int,
        default=300
    )

    p.add_argument(
        "--dedup",
        action="store_true"
    )

    p.add_argument(
        "--xpad",
        type=float,
        default=0.15
    )

    p.add_argument(
        "--colors",
        default="",
        help=(
            'Ej: '
            '"HIGH=#D88B8B,'
            'MODERATE=#AED8C8,'
            'LOW=#F9CEAE,'
            'MODIFIER=#CAD9FB"'
        ),
    )

    args = p.parse_args()

    # Parse color configuration
    colors = parse_colors(args.colors)

    # Output directory
    outdir = Path(args.outdir)

    # Process each input file
    for f in args.inputs:

        path = Path(f)

        # Load annotation table
        df = load_table(
            path,
            dedup=args.dedup
        )

        # Output base name
        outbase = (
            outdir /
            f"{path.stem}_tipo_impacto"
        )

        # Generate figure
        make_plot(
            df=df,
            title="",
            outbase=outbase,
            dpi=args.dpi,
            xpad_frac=args.xpad,
            colors=colors,
        )

        print(f"[OK] {outbase}.png / .pdf")

# ============================================================
# Script entry point
# ============================================================

if __name__ == "__main__":
    main()
