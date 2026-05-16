#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Per-gene mean coverage calculation from sequencing coverage files.

Description:
This script processes sequencing coverage files from two cohorts
and calculates the mean coverage per gene for each sample.
Gene names are extracted from amplicon identifiers.
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
    calculate the mean coverage per gene.

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
        - mean coverage
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
        )

        # -------------------------------------------------
        # Extract gene name from amplicon identifier
        # -------------------------------------------------

        df["gen"] = (
            df["amplicon"]
            .str.split("-")
            .str[0]
        )

        # -------------------------------------------------
        # Calculate mean coverage per gene
        # -------------------------------------------------

        resumen = (
            df
            .groupby("gen")["cobertura"]
            .mean()
            .reset_index()
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

df_gen = pd.concat(
    [pob1, pob2],
    ignore_index=True
)


# =========================================================
# Export results
# =========================================================

df_gen.to_csv(
    "cobertura_media_por_gen_y_muestra.tsv",
    sep="\t",
    index=False
)


# =========================================================
# Final message
# =========================================================

print(
    "Output file generated: "
    "cobertura_media_por_gen_y_muestra.tsv"
)
