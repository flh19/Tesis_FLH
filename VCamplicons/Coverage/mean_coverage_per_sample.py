#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd

def procesar_poblacion(ruta, poblacion):
    resultados = []

    for fichero in os.listdir(ruta):
        if not fichero.endswith(".cov"):
            continue

        path = os.path.join(ruta, fichero)

        df = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=["muestra", "amplicon", "posicion", "cobertura"]
        )

        df["cobertura"] = pd.to_numeric(df["cobertura"], errors="coerce")
        cobertura_media = df["cobertura"].mean()

        resultados.append({
            "muestra": fichero.replace(".primerclipped.cov", ""),
            "poblacion": poblacion,
            "cobertura_media": cobertura_media
        })

    return pd.DataFrame(resultados)


# ====== RUTAS ======
pob1 = procesar_poblacion("../Diabet_study/cov_files/", "Di@bet.es") #path to the cov_files cohort 1
pob2 = procesar_poblacion("../Hortega_study/cov_files/", "Hortega") #path to the cov_files cohort 2

# Combinar
df_final = pd.concat([pob1, pob2], ignore_index=True)

# Guardar
df_final.to_csv("mean_coverage_per_sample.tsv", sep="\t", index=False)

print("Archivo generado: mean_coverage_per_sample.tsv")
