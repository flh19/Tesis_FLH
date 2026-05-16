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

#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# ===================== CONFIGURACIÓN ===================== #

IMPACT_ORDER = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

IMPACT_LABELS_ES = {
    "HIGH": "ALTO",
    "MODERATE": "MODERADO",
    "LOW": "BAJO",
    "MODIFIER": "MODIFICADOR",
}

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
    "splice_region_variant&intron_variant": "Splicing no canónico",
    "splice_region_variant&synonymous_variant": "Splicing no canónico",
    "stop_gained": "Ganancia stop",
    "synonymous_variant": "Sinónima",
}


# ===================== UTILIDADES ===================== #

def parse_colors(colors_str: str | None) -> dict[str, str]:
    if not colors_str:
        return {}
    out = {}
    for chunk in colors_str.split(","):
        k, v = chunk.split("=")
        out[k.strip().upper()] = v.strip()
    return out


def italicize_english_words(label: str) -> str:
    if label is None:
        return ""

    text = str(label)

    subs = [
        (r"\bUTR\b", lambda m: r"$\mathit{UTR}$"),
        (r"\binframe\b", lambda m: r"$\mathit{inframe}$"),
        (r"\bFrameshift\b", lambda m: r"$\mathit{Frameshift}$"),
        (r"\bMissense\b", lambda m: r"$\mathit{Missense}$"),
        (r"\bSplicing\b", lambda m: r"$\mathit{Splicing}$"),
        (r"\bstart\b", lambda m: r"$\mathit{start}$"),
        (r"\bstop\b", lambda m: r"$\mathit{stop}$"),
    ]

    for pat, repl in subs:
        text = re.sub(pat, repl, text)

    return text


def map_consequence(raw: str) -> str:
    if pd.isna(raw):
        return ""
    return CONSEQUENCE_MAP.get(raw, raw)


def load_table(path: Path, dedup: bool = False) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)

    df["Impacto"] = df["Impacto"].str.upper().str.strip()
    df["Consecuencia_es"] = (
        df["Consecuencia"]
        .apply(map_consequence)
        .apply(italicize_english_words)
    )

    if dedup:
        df = df.drop_duplicates(subset=["ID"])

    return df


# ===================== PREPARAR CONTEOS ===================== #

def prepare_counts(df: pd.DataFrame, top_n: int = 14) -> pd.DataFrame:
    counts = (
        df.groupby(["Consecuencia_es", "Impacto"])
        .size()
        .unstack(fill_value=0)
    )

    for imp in IMPACT_ORDER:
        if imp not in counts.columns:
            counts[imp] = 0

    counts = counts[IMPACT_ORDER]
    counts["TOTAL"] = counts.sum(axis=1)
    counts = counts.sort_values("TOTAL", ascending=False).head(top_n)

    return counts


# ===================== PLOT MULTIPANEL ===================== #

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

    # ---------- TIPOGRAFÍA GLOBAL ----------
    plt.rcParams.update({
        "font.size": 25,
        "axes.titlesize": 25,
        "axes.labelsize": 25,
        "xtick.labelsize": 25,
        "ytick.labelsize": 25,
        "legend.fontsize": 22,
        "legend.title_fontsize": 22,
    })

    # ---------- ORDEN CONSISTENTE ----------
    order = counts_A.index
    counts_B = counts_B.reindex(order).fillna(0)

    max_x = max(counts_A["TOTAL"].max(), counts_B["TOTAL"].max())
    n = len(order)

    fig_h = max(8.0, 0.65 * n + 2.5)
    fig, axes = plt.subplots(
        ncols=2,
        figsize=(20, fig_h),
        sharey=True,
    )

    for ax, counts, title in zip(
        axes,
        [counts_A, counts_B],
        [label_A, label_B],
    ):
        y = [i * y_gap for i in range(n)]
        cum = pd.Series(0, index=counts.index)

        for imp in IMPACT_ORDER:
            ax.barh(
                y,
                counts[imp].values,
                left=cum.values,
                height=bar_height,
                label=IMPACT_LABELS_ES[imp],
                color=colors.get(imp),
                edgecolor="black",
                linewidth=0.8,
            )
            cum += counts[imp]

        # ---------- NÚMEROS AL FINAL ----------
        for yi, total in zip(y, cum.values):
            ax.text(
                total + max_x * 0.01,
                yi,
                f"{int(total)}",
                va="center",
                ha="left",
                fontsize=22,
            )

        ax.set_title(title, fontweight="bold", pad=10)
        ax.set_xlim(0, max_x * (1 + xpad_frac))
        ax.xaxis.grid(True, linestyle="--", alpha=0.35)
        ax.set_axisbelow(True)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Mantener solo ejes X e Y
        ax.spines["left"].set_linewidth(1.5)
        ax.spines["bottom"].set_linewidth(1.5) 
        ax.invert_yaxis()

    axes[0].set_yticks([i * y_gap for i in range(n)])
    axes[0].set_yticklabels(order)

    axes[0].set_xlabel("Número de variantes")
    axes[1].set_xlabel("Número de variantes")

    # ---------- LEYENDA ÚNICA ----------
    handles, labels = axes[0].get_legend_handles_labels()
    leg = fig.legend(
        handles,
        labels,
        title="Impacto",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=4,
        frameon=True,
    )
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.8)
    leg.get_frame().set_facecolor("white")

    fig.tight_layout(rect=[0, 0.08, 1, 1])
    outbase.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outbase.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


# ===================== MAIN ===================== #

def main():
    p = argparse.ArgumentParser()
    p.add_argument("cohorte_A")
    p.add_argument("cohorte_B")
    p.add_argument("--labelA", default="Cohorte A")
    p.add_argument("--labelB", default="Cohorte B")
    p.add_argument("-o", "--outdir", default="plots")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--dedup", action="store_true")
    p.add_argument("--xpad", type=float, default=0.15)
    p.add_argument(
        "--colors",
        default="HIGH=#D88B8B,MODERATE=#AED8C8,LOW=#F9CEAE,MODIFIER=#CAD9FB",
    )

    args = p.parse_args()
    colors = parse_colors(args.colors)

    dfA = load_table(Path(args.cohorte_A), dedup=args.dedup)
    dfB = load_table(Path(args.cohorte_B), dedup=args.dedup)

    counts_A = prepare_counts(dfA)
    counts_B = prepare_counts(dfB)

    outbase = Path(args.outdir) / "comparacion_consecuencia_impacto"

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

    print(f"[OK] {outbase}.png / .pdf")


if __name__ == "__main__":
    main()
