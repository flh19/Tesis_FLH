#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Per-gene sequencing uniformity calculation at 50X coverage threshold.

Description:
This script processes sequencing coverage files from two cohorts
and calculates per-gene coverage uniformity for each sample.
Uniformity is defined as the proportion of positions with
coverage greater than or equal to 50X.
"""

# =========================================================
# Imports
# =========================================================

import os

import pandas as pd


# =========================================================
# Global configuration
# =========================================================

UMBRAL = 50  # 50X coverage threshold


# =========================================================
# Uniformity processing functions
# =========================================================

def procesar_poblacion(ruta, poblacion):
    """
    Process all .cov files from a cohort directory and
    calculate per-gene sequencing uniformity.

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
        - gene name
        - uniformity at 50X
        - sample name
        - cohort name
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
        # Ensure numeric coverage values
        # -------------------------------------------------

        df["cobertura"] = pd.to_numeric(
            df["cobertura"],
            errors="coerce"
        ).fillna(0)

        # -------------------------------------------------
        # Extract gene name from amplicon identifier
        # -------------------------------------------------

        df["gen"] = (
            df["amplicon"]
            .str.split("-")
            .str[0]
        )

        # -------------------------------------------------
        # Calculate per-gene uniformity
        # -------------------------------------------------

        resumen = (
            df
            .groupby("gen")
            .apply(
                lambda x:
                (
                    x["cobertura"] >= UMBRAL
                ).sum() / len(x)
            )
            .reset_index(
                name="uniformidad_50X"
            )
        )

        # -------------------------------------------------
        # Add sample and cohort metadata
        # -------------------------------------------------

        resumen["muestra"] = fichero.replace(
            ".primerclipped.cov",
            ""
        )

        resumen["poblacion"] = poblacion

        resultados.append(resumen)

    return pd.concat(
        resultados,
        ignore_index=True
    )


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

df_uniformidad = pd.concat(
    [pob1, pob2],
    ignore_index=True
)


# =========================================================
# Export results
# =========================================================

df_uniformidad.to_csv(
    "uniformidad_50X_por_gen_y_muestra.tsv",
    sep="\t",
    index=False
)


# =========================================================
# Final message
# =========================================================

print(
    "Output file generated: "
    "uniformidad_50X_por_gen_y_muestra.tsv"
)
