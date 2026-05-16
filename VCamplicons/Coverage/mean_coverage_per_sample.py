#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mean coverage calculation per sample from .cov files.

Description:
This script processes coverage files generated for two sequencing cohorts
and calculates the mean sequencing depth per sample. Results from both
cohorts are merged into a single output table.
"""

# =========================================================
# Imports
# =========================================================

import os
import pandas as pd

# =========================================================
# Coverage processing functions
# =========================================================

def procesar_poblacion(ruta, poblacion):
    """
    Process all .cov files from a cohort directory and
    calculate the mean coverage per sample.

    Parameters
    ----------
    ruta : str
        Path to the directory containing .cov files.

    poblacion : str
        Cohort name.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing:
        - sample name
        - cohort name
        - mean coverage
    """

    resultados = []

    # -----------------------------------------------------
    # Iterate through coverage files
    # -----------------------------------------------------

    for fichero in os.listdir(ruta):

        if not fichero.endswith(".cov"):
            continue

        path = os.path.join(ruta, fichero)

        # -------------------------------------------------
        # Read coverage file
        # -------------------------------------------------

        df = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=[
                "muestra",
                "amplicon",
                "posicion",
                "cobertura"
            ]
        )

        # -------------------------------------------------
        # Convert coverage column to numeric values
        # -------------------------------------------------

        df["cobertura"] = pd.to_numeric(
            df["cobertura"],
            errors="coerce"
        )

        # -------------------------------------------------
        # Calculate mean coverage
        # -------------------------------------------------

        cobertura_media = df["cobertura"].mean()

        # -------------------------------------------------
        # Store results
        # -------------------------------------------------

        resultados.append({
            "muestra": fichero.replace(
                ".primerclipped.cov",
                ""
            ),

            "poblacion": poblacion,

            "cobertura_media": cobertura_media
        })

    return pd.DataFrame(resultados)


# =========================================================
# Input cohort processing
# =========================================================

# Cohort 1:
# Path to coverage files from Di@bet.es cohort

pob1 = procesar_poblacion(
    "../Diabet_study/cov_files/",
    "Di@bet.es"
)

# Cohort 2:
# Path to coverage files from Hortega cohort

pob2 = procesar_poblacion(
    "../Hortega_study/cov_files/",
    "Hortega"
)


# =========================================================
# Merge cohort results
# =========================================================

df_final = pd.concat(
    [pob1, pob2],
    ignore_index=True
)


# =========================================================
# Export results
# =========================================================

df_final.to_csv(
    "mean_coverage_per_sample.tsv",
    sep="\t",
    index=False
)


# =========================================================
# Final message
# =========================================================

print(
    "Output file generated: mean_coverage_per_sample.tsv"
)
