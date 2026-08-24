"""
06b_robustness_checks.py
Verificaciones de robustez usando archivos ya existentes (Script 06),
sin recomputar WOT:
  1. Pseudorreplicacion: la caida de entropia, se sostiene a nivel de
     animal (n=2-6), no solo a nivel celula (n=1563-4481)?
  2. Efecto composicional vs intrinseco: dentro del Estado 1 (el unico
     con n razonable en ambos intervalos), cae la entropia por si
     sola, o es solo redistribucion de estados?
"""
import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from config import PROCESSED, RESULTS_IRREV, PROJECT_ROOT

print("Chequeo 1: entropia a nivel animal (no solo celula)\n")

entropy = pd.read_csv(os.path.join(RESULTS_IRREV, "entropy_por_celula.csv"))
meta = pd.read_csv(os.path.join(PROCESSED, "cell_metadata.csv"))
entropy = entropy.merge(meta[["id", "animal"]], on="id", how="left")

animal_summary = entropy.groupby(["intervalo", "animal"])["entropy_norm"].agg(
    ["mean", "count"]
).reset_index()
print(animal_summary.to_string(index=False))

n_animales_0_15 = animal_summary[animal_summary.intervalo == "0_15"]["animal"].nunique()
n_animales_15_30 = animal_summary[animal_summary.intervalo == "15_30"]["animal"].nunique()
print(f"\nAnimales distintos en 0->15sem: {n_animales_0_15}")
print(f"Animales distintos en 15->30sem: {n_animales_15_30}")
print("\nCon este n de animales no se puede calcular un p-valor formal")
print("robusto a nivel animal por si solo (potencia insuficiente); lo que")
print("si se puede reportar es si la direccion del efecto se mantiene al")
print("promediar por animal (impreso arriba). Documentar como limitacion.\n")

animal_summary.to_csv(os.path.join(RESULTS_IRREV, "entropy_por_animal.csv"), index=False)

# "Animal1" en 0_15 y "Animal1" en 15_30 son animales distintos (diseno
# transversal, cohortes independientes por timepoint) - la comparacion
# correcta es entre grupos independientes, no pareada.
medias_0_15 = animal_summary.loc[animal_summary.intervalo == "0_15", "mean"].values
medias_15_30 = animal_summary.loc[animal_summary.intervalo == "15_30", "mean"].values

print(f"Medias por animal, 0->15sem (n={len(medias_0_15)}): {sorted(medias_0_15)}")
print(f"Medias por animal, 15->30sem (n={len(medias_15_30)}): {sorted(medias_15_30)}")

u_stat, p_mw = mannwhitneyu(medias_0_15, medias_15_30, alternative="two-sided")
print(f"\nMann-Whitney U = {u_stat}, p = {p_mw:.4f}\n")

print("Chequeo 2: efecto composicional vs intrinseco (Estado 1)\n")

por_estado = pd.read_csv(os.path.join(RESULTS_IRREV, "entropy_por_estado.csv"))
estado1 = por_estado[por_estado.seurat_clusters == 1]
print(estado1.to_string(index=False))
print("\nSi la entropia del Estado 1 cae entre intervalos, la caida no es")
print("puramente composicional (hay reduccion dentro del mismo estado).")
print("Nota: n desigual entre intervalos - interpretar con cautela.\n")

print("Chequeos completos. Archivo nuevo: entropy_por_animal.csv")
