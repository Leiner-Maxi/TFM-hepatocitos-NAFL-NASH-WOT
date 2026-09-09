"""
08_generate_figures.py
Genera las figuras de resultados a partir de los CSVs producidos por
los Scripts 06, 06b y 07 (fig10-fig12), y combina las figuras UMAP de
los Scripts 02-03 en paneles multi-panel (a,b,c). No recalcula ningun
resultado cientifico nuevo, salvo el test de Mann-Whitney (recomputado
desde el resumen ya guardado en entropy_por_animal.csv - sin
aleatoriedad, resultado identico al ya reportado por Script 06b, que
no lo guarda en archivo).

Todo se guarda en PDF (vectorial): formato requerido por Nature y la
mayoria de revistas Q1 para figuras de lineas/barras/texto - ver
nature.com/nature/for-authors/final-submission. La conversion a PNG
para manuscrito/presentacion queda centralizada en Script 09, que debe
correrse DESPUES de este script.

Input:  results/irreversibility/{bootstrap_ci,entropy_por_animal,
        driver_genes_alta_confianza,permutation_test}.{csv,txt}
        results/figures/main/{umap_por_cluster,umap_por_condicion,
        umap_por_fraccion,umap_hepatocitos_por_estado,
        umap_hepatocitos_por_condicion}.pdf (Scripts 02-03, R)
Output: results/figures/main/{fig10_entropia_bootstrap,
        fig11_validacion_animal,fig12_driver_genes,
        panel_umap_global,panel_umap_hepatocitos}.pdf
"""
import os
import sys
import re
import numpy as np
import pandas as pd
import pymupdf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(__file__))
from config import RESULTS_IRREV, PROJECT_ROOT

FIGURES_MAIN = os.path.join(PROJECT_ROOT, "results", "figures", "main")
os.makedirs(FIGURES_MAIN, exist_ok=True)

plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})

# ============================================================
# Figura 10: entropia normalizada por intervalo, con bootstrap CI
# ============================================================
bootstrap = pd.read_csv(os.path.join(RESULTS_IRREV, "bootstrap_ci.csv")).set_index("intervalo")

with open(os.path.join(RESULTS_IRREV, "permutation_test.txt")) as f:
    perm_text = f.read()
p_match = re.search(r"p-valor.*?:\s*([\d.]+)", perm_text)
p_value = float(p_match.group(1)) if p_match else None
p_label = "p < 0.001" if (p_value is not None and p_value < 0.001) else f"p = {p_value:.4f}"

intervalos = ["Chow -> NAFLD\n(0-15 sem)", "NAFLD -> NASH\n(15-30 sem)"]
medias = [bootstrap.loc["0_15", "media"], bootstrap.loc["15_30", "media"]]
ci_low = [bootstrap.loc["0_15", "ci_low"], bootstrap.loc["15_30", "ci_low"]]
ci_high = [bootstrap.loc["0_15", "ci_high"], bootstrap.loc["15_30", "ci_high"]]
errores = [[m - lo for m, lo in zip(medias, ci_low)],
           [hi - m for m, hi in zip(medias, ci_high)]]
colores = ["#4C72B0", "#C44E52"]

fig, ax = plt.subplots(figsize=(6, 5))
ax.bar(intervalos, medias, yerr=errores, capsize=8, color=colores,
       edgecolor="black", linewidth=0.8, width=0.55)
ax.set_ylabel("Entropía normalizada de destino (H / log N)")
ax.set_ylim(0, 0.9)
ax.set_title("Entropía normalizada por intervalo temporal\n(IC 95%, bootstrap n=1000)")
for i, m in enumerate(medias):
    ax.text(i, m + 0.03, f"{m:.3f}", ha="center", fontweight="bold")
ax.text(0.5, 0.85, f"{p_label} (test de permutación, n=1000)",
        ha="center", transform=ax.transAxes, fontsize=9, style="italic")
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_MAIN, "fig10_entropia_bootstrap.pdf"))
plt.close()
print("Figura 10 guardada (entropia normalizada por intervalo).")

# ============================================================
# Figura 11: validacion a nivel de replica biologica
# ============================================================
animal_summary = pd.read_csv(os.path.join(RESULTS_IRREV, "entropy_por_animal.csv"))
medias_0_15 = animal_summary.loc[animal_summary.intervalo == "0_15", "mean"].values
medias_15_30 = animal_summary.loc[animal_summary.intervalo == "15_30", "mean"].values

u_stat, p_mw = mannwhitneyu(medias_0_15, medias_15_30, alternative="two-sided")

fig, ax = plt.subplots(figsize=(6, 5))
# Dos generadores independientes (semillas 42 y 1), no uno secuencial:
# reproduce exactamente la disposicion horizontal de los puntos del
# script original.
x0 = np.random.RandomState(42).normal(0, 0.04, len(medias_0_15))
x1 = np.random.RandomState(1).normal(1, 0.04, len(medias_15_30))

ax.scatter(x0, medias_0_15, s=90, color="#4C72B0", edgecolor="black",
           zorder=3, label=f"0→15 sem (n={len(medias_0_15)} animales)")
ax.scatter(x1, medias_15_30, s=90, color="#C44E52", edgecolor="black",
           zorder=3, label=f"15→30 sem (n={len(medias_15_30)} animales)")
ax.hlines(np.mean(medias_0_15), -0.2, 0.2, color="#4C72B0", linewidth=2)
ax.hlines(np.mean(medias_15_30), 0.8, 1.2, color="#C44E52", linewidth=2)
ax.set_xticks([0, 1])
ax.set_xticklabels(["Chow→NAFLD", "NAFLD→NASH"])
ax.set_ylabel("Entropía normalizada media por animal")
# CORREGIDO: titulo principal sin estadistico (evita peso visual
# igual al p<0.001 de Figura 10); el estadistico pasa a anotacion
# secundaria, mas pequenia y con nota de potencia limitada, dejando
# claro que esta es evidencia CONFIRMATORIA, no la evidencia primaria.
ax.set_title("Validación a nivel de réplica biológica", fontsize=12, fontweight="bold", pad=22)
ax.text(0.5, 1.02, f"Mann-Whitney U={u_stat:.0f}, p={p_mw:.3f} "
        f"(confirmatorio; potencia limitada por n=6 vs n=3)",
        ha="center", transform=ax.transAxes, fontsize=8.5, style="italic", color="#444444")
ax.set_xlim(-0.5, 1.5)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_MAIN, "fig11_validacion_animal.pdf"))
plt.close()
print(f"Figura 11 guardada (Mann-Whitney U={u_stat:.0f}, p={p_mw:.4f}).")

# ============================================================
# Figura 12: genes con mayor asociacion
# ============================================================
driver = pd.read_csv(os.path.join(RESULTS_IRREV, "driver_genes_alta_confianza.csv"))
neg = driver[driver.spearman_rho < 0].sort_values("spearman_rho").head(10)
pos = driver[driver.spearman_rho > 0].sort_values("spearman_rho", ascending=False).head(10)

fig, axes = plt.subplots(1, 2, figsize=(11, 6), sharey=False)
axes[0].barh(neg["gene"][::-1], neg["spearman_rho"][::-1], color="#C44E52", edgecolor="black")
axes[0].set_xlabel("ρ de Spearman")
axes[0].set_title("Asociados a mayor compromiso\n(menor entropía)")
axes[0].axvline(0, color="black", linewidth=0.8)

axes[1].barh(pos["gene"][::-1], pos["spearman_rho"][::-1], color="#4C72B0", edgecolor="black")
axes[1].set_xlabel("ρ de Spearman")
axes[1].set_title("Asociados a mayor plasticidad\n(mayor entropía)")
axes[1].axvline(0, color="black", linewidth=0.8)

fig.suptitle("Genes con mayor asociación a la entropía normalizada de destino\n"
             "(intervalo 15→30 semanas, FDR<0.05, |ρ|>0.15)", y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_MAIN, "fig12_driver_genes.pdf"), bbox_inches="tight")
plt.close()
print("Figura 12 guardada.")

# ============================================================
# Paneles combinados: UMAP global (a,b,c) y UMAP hepatocitario (a,b)
# Combina figuras PDF ya generadas por Scripts 02-03 (R); no genera
# datos nuevos, solo compone imagenes ya producidas.
# ============================================================
PANELES = [
    {
        "nombre_salida": "panel_umap_global",
        "archivos": ["umap_por_cluster", "umap_por_condicion", "umap_por_fraccion"],
        "cols": 3,
    },
    {
        "nombre_salida": "panel_umap_hepatocitos",
        "archivos": ["umap_hepatocitos_por_estado", "umap_hepatocitos_por_condicion"],
        "cols": 2,
    },
]


def pdf_a_array(pdf_path, dpi=200):
    zoom = dpi / 72
    doc = pymupdf.open(pdf_path)
    pix = doc[0].get_pixmap(matrix=pymupdf.Matrix(zoom, zoom))
    doc.close()
    return np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3)


for panel in PANELES:
    rutas = [os.path.join(FIGURES_MAIN, f"{n}.pdf") for n in panel["archivos"]]
    faltantes = [r for r in rutas if not os.path.exists(r)]
    if faltantes:
        print(f"AVISO: faltan {faltantes} (correr Scripts 02/03 primero) - se omite {panel['nombre_salida']}")
        continue

    imagenes = [pdf_a_array(r) for r in rutas]
    cols = panel["cols"]
    fig, axes = plt.subplots(1, cols, figsize=(6.2 * cols, 6.2 / imagenes[0].shape[1] * imagenes[0].shape[0]))
    if cols == 1:
        axes = [axes]

    for ax, img, letra in zip(axes, imagenes, "abcdefgh"):
        ax.imshow(img)
        ax.axis("off")
        ax.text(0.0, 1.02, letra, transform=ax.transAxes, fontsize=20,
                fontweight="bold", va="bottom", ha="left")

    plt.tight_layout()
    out_path = os.path.join(FIGURES_MAIN, f"{panel['nombre_salida']}.pdf")
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"{panel['nombre_salida']}: {len(imagenes)} paneles ({list('abcdefgh'[:len(imagenes)])}) -> {out_path}")

print("\nScript 08 completo. Correr Script 09 a continuacion para generar los PNG.")
