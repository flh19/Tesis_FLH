#!/usr/bin/env python3
"""
FastQC_post_hortega.py

Boxplots "Per base sequence quality" ENTRE MUESTRAS
AGRUPADOS EN VENTANAS DE 5 bp
"""

import argparse
import zipfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


ALLOWED_VALUE_COLS = {"mean", "median", "q1", "q3", "p10", "p90"}

COLOR_BOX  = "grey"
COLOR_LINE = "#A9A9A9"
COLOR_THRESH_ORANGE = "#f37a48"
COLOR_THRESH_GREEN  = "#478d48"
GRID_ALPHA = 0.35
WINDOW_BP = 5
MAX_BASE = 150


@dataclass
class AggRow:
    start: int
    end: int
    mean: float
    median: float
    q1: float
    q3: float
    p10: float
    p90: float
    n: int


def parse_base_range(s: str) -> Tuple[int, int]:
    if "-" in s:
        a, b = s.split("-", 1)
        return int(a), int(b)
    return int(s), int(s)


def read_first_matching_from_zip(zpath: str, suffix: str) -> Optional[str]:
    with zipfile.ZipFile(zpath, "r") as z:
        for name in z.namelist():
            if name.endswith(suffix):
                with z.open(name) as f:
                    return f.read().decode("utf-8", errors="replace")
    return None


def extract_per_base_quality_values(txt: str, value_col: str) -> Dict[int, float]:
    col_idx = {"mean": 1, "median": 2, "q1": 3, "q3": 4, "p10": 5, "p90": 6}[value_col]
    out: Dict[int, float] = {}
    in_mod = False

    for line in txt.splitlines():
        if line.startswith(">>Per base sequence quality"):
            in_mod = True
            continue
        if not in_mod:
            continue
        if line.startswith(">>END_MODULE"):
            break
        if not line or line.startswith("#"):
            continue

        parts = line.split("\t")
        if len(parts) < 7:
            continue

        try:
            v = float(parts[col_idx])
        except ValueError:
            continue
        if not np.isfinite(v):
            continue

        a, b = parse_base_range(parts[0])
        for pos in range(a, b + 1):
            if pos <= MAX_BASE:
                out[pos] = v

    if not out:
        raise ValueError("No se pudieron extraer valores")
    return out


def process_zip(zpath: str, value_col: str) -> Dict[int, float]:
    txt = read_first_matching_from_zip(zpath, "fastqc_data.txt")
    if txt is None:
        raise ValueError(f"fastqc_data.txt no encontrado en {zpath}")

    d = extract_per_base_quality_values(txt, value_col)

    for pos in range(1, MAX_BASE + 1):
        d.setdefault(pos, np.nan)

    return d


def aggregate_across_samples(dicts: List[Dict[int, float]], min_n: int) -> List[AggRow]:
    window_vals: Dict[Tuple[int, int], List[float]] = {}

    for d in dicts:
        for pos, val in d.items():
            if pos > MAX_BASE:
                continue

            w_start = ((pos - 1) // WINDOW_BP) * WINDOW_BP + 1
            w_end = w_start + WINDOW_BP - 1
            window_vals.setdefault((w_start, w_end), []).append(val)

    rows: List[AggRow] = []
    for (a, b) in sorted(window_vals):
        vals = np.array(window_vals[(a, b)], dtype=float)
        vals = vals[np.isfinite(vals)]

        if vals.size < min_n:
            continue

        rows.append(
            AggRow(
                start=a,
                end=b,
                mean=float(vals.mean()),
                median=float(np.percentile(vals, 50)),
                q1=float(np.percentile(vals, 25)),
                q3=float(np.percentile(vals, 75)),
                p10=float(np.percentile(vals, 10)),
                p90=float(np.percentile(vals, 90)),
                n=int(vals.size),
            )
        )

    if not rows:
        raise ValueError("No quedan ventanas tras min-n")
    return rows


def plot_merged(
    rows: List[AggRow],
    out_prefix: str,
    title: str,
    y_min: float = 18,
    y_max: float = 40,
    dpi: int = 200,
):

    n_bins = len(rows)
    fig, ax = plt.subplots(figsize=(18, 6))

    fig.patch.set_alpha(0)
    ax.set_facecolor("none")

    ax.axhline(20, color=COLOR_THRESH_ORANGE, linestyle="--", linewidth=3)
    ax.axhline(28, color=COLOR_THRESH_GREEN, linestyle="--", linewidth=3)

    box_w = 0.75
    cap_w = 0.25

    for i, r in enumerate(rows, start=1):
        x0 = i - box_w / 2

        ax.add_patch(Rectangle(
            (x0, r.q1),
            box_w,
            r.q3 - r.q1,
            facecolor=COLOR_BOX,
            edgecolor=COLOR_BOX,
            alpha=0.22
        ))

        ax.vlines(i, r.p10, r.p90, color=COLOR_BOX)
        ax.hlines([r.p10, r.p90], i - cap_w, i + cap_w, color=COLOR_BOX)
        ax.hlines(r.median, x0, x0 + box_w, color=COLOR_BOX, linewidth=1.2)

    ax.plot(range(1, n_bins + 1),
            [r.mean for r in rows],
            color=COLOR_LINE, linewidth=1.3)

    ax.set_xlim(0.5, n_bins + 0.5)
    ax.set_ylim(y_min, y_max)

    ax.set_yticks([20, 25, 30, 35, 40])
    ax.set_yticklabels(["20", "25", "30", "35", "40"], fontsize=20)

    ax.set_ylabel("Calidad Phred (Q)", fontsize=22)
    ax.set_xlabel("Intervalos de posición (pb)", fontsize=22)
    ax.set_title(title, fontsize=18)

    ticks = list(range(1, n_bins + 1))
    labels = [f"{r.start}-{r.end}" for r in rows]

    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=20)

    ax.grid(True, alpha=GRID_ALPHA)
    fig.tight_layout()

    fig.savefig(f"{out_prefix}_merged.png", dpi=dpi, transparent=True)
    fig.savefig(f"{out_prefix}_merged.pdf", dpi=dpi, transparent=True)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file-list", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--value-col", default="mean", choices=sorted(ALLOWED_VALUE_COLS))
    ap.add_argument("--min-n", type=int, default=1500)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--title", default="")

    args = ap.parse_args()

    paths = [p.strip() for p in open(args.file_list) if p.strip()]
    if not paths:
        raise ValueError("file-list vacío")

    all_dicts: List[Dict[int, float]] = []

    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = [ex.submit(process_zip, p, args.value_col) for p in paths]
        for f in as_completed(futs):
            all_dicts.append(f.result())

    rows = aggregate_across_samples(all_dicts, args.min_n)

    plot_merged(
        rows=rows,
        out_prefix=args.out_prefix,
        title=args.title,
        dpi=args.dpi
    )


if __name__ == "__main__":
    main()
