"""
07b_gene_robustness_loo.py
Validacion leave-one-animal-out de los 20 genes que sustentan la
Figura 12, la Tabla 4 y la discusion biologica (apartado 5.5): para
cada gen, recalcula la correlacion de Spearman excluyendo un animal a
la vez (3 animales en el intervalo 15->30sem, ver Script 06b), con
correccion FDR (Benjamini-Hochberg, mismo metodo que Script 07)
recalculada dentro de cada subconjunto de 20 genes.

Clasificacion de robustez (sin umbral numerico nuevo, reutiliza el
FDR<0.05 ya establecido en la metodologia del TFM):
  Robusto              - signo igual al original en las 3 exclusiones
                         Y significativo (FDR<0.05) en las 3.
  Parcialmente robusto - signo igual en las 3, pero pierde
                         significancia en al menos una.
  No robusto           - el signo se invierte en al menos una
                         exclusion.

Input:  results/irreversibility/{driver_genes_alta_confianza,
        entropy_por_celula}.csv
        data/processed/{cell_metadata,expression_matrix}.{csv,txt}
Output: results/irreversibility/gene_loo_robustness.csv
        results/irreversibility/gene_loo_summary.md
        results/figures/main/gene_loo_heatmap.pdf
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from config import PROCESSED, RESULTS_IRREV, PROJECT_ROOT

FIGURES_MAIN = os.path.join(PROJECT_ROOT, "results", "figures", "main")
os.makedirs(FIGURES_MAIN, exist_ok=True)

N_TOP_POR_DIRECCION = 10  # mismo criterio que arma la Tabla 4 (10+10)
FDR_THRESHOLD = 0.05      # mismo umbral ya establecido en Script 07


def bh_fdr(p_vals):
    """Benjamini-Hochberg, identico al usado en Script 07."""
    n = len(p_vals)
    order = np.argsort(p_vals)
    ranked_p = p_vals[order] * n / (np.arange(n) + 1)
    ranked_p = np.minimum.accumulate(ranked_p[::-1])[::-1]
    q_vals = np.empty(n)
    q_vals[order] = np.clip(ranked_p, 0, 1)
    return q_vals


# --- Cargar y seleccionar los 20 genes dinamicamente (sin lista manual) ---
driver = pd.read_csv(os.path.join(RESULTS_IRREV, "driver_genes_alta_confianza.csv"))
neg20 = driver[driver.spearman_rho < 0].sort_values("spearman_rho").head(N_TOP_POR_DIRECCION)
pos20 = driver[driver.spearman_rho > 0].sort_values("spearman_rho", ascending=False).head(N_TOP_POR_DIRECCION)
genes_20 = pd.concat([neg20, pos20], ignore_index=True)
print(f"Genes seleccionados dinamicamente desde driver_genes_alta_confianza.csv: {len(genes_20)}")

# --- Cargar datos reales del pipeline ---
entropy = pd.read_csv(os.path.join(RESULTS_IRREV, "entropy_por_celula.csv"))
meta = pd.read_csv(os.path.join(PROCESSED, "cell_metadata.csv"))
expr = pd.read_csv(os.path.join(PROCESSED, "expression_matrix.txt")).set_index("id")

df = entropy[entropy.intervalo == "15_30"].merge(meta[["id", "animal"]], on="id", how="left")
df = df[df["id"].isin(expr.index)].reset_index(drop=True)
expr_aligned = expr.loc[df["id"]]

animales = sorted(df["animal"].dropna().unique())
print(f"Animales en intervalo 15->30sem: {animales}")
if len(animales) != 3:
    print(f"AVISO: se esperaban 3 animales, se encontraron {len(animales)}.")

# --- Recalcular Spearman + FDR en cada subconjunto leave-one-out ---
rho_loo = {a: [] for a in animales}
p_loo = {a: [] for a in animales}

for _, row in genes_20.iterrows():
    gene = row["gene"]
    expr_gene = expr_aligned[gene].values
    for animal_excluido in animales:
        mask = (df["animal"] != animal_excluido).values
        rho, p = spearmanr(expr_gene[mask], df.loc[mask, "entropy_norm"].values)
        rho_loo[animal_excluido].append(rho)
        p_loo[animal_excluido].append(p)

# FDR recalculado DENTRO de cada subconjunto (20 genes), no global
q_loo = {a: bh_fdr(np.array(p_loo[a])) for a in animales}

# --- Construir tabla de resultados ---
filas = []
for i, (_, row) in enumerate(genes_20.iterrows()):
    gene = row["gene"]
    rho_original = row["spearman_rho"]
    signo_original = np.sign(rho_original)

    rhos_i = [rho_loo[a][i] for a in animales]
    qs_i = [q_loo[a][i] for a in animales]

    mantiene_signo = all(np.sign(r) == signo_original for r in rhos_i)
    sig_en_las_3 = all(q < FDR_THRESHOLD for q in qs_i)

    if not mantiene_signo:
        clasificacion = "No robusto"
    elif sig_en_las_3:
        clasificacion = "Robusto"
    else:
        clasificacion = "Parcialmente robusto"

    fila = {"gene": gene, "spearman_rho_original": rho_original}
    for a, r, q in zip(animales, rhos_i, qs_i):
        fila[f"rho_sin_{a}"] = r
        fila[f"fdr_sin_{a}"] = q
    fila["rango_variacion"] = max(rhos_i) - min(rhos_i)
    fila["cambio_maximo_absoluto"] = max(abs(rho_original - r) for r in rhos_i)
    fila["mantiene_signo"] = mantiene_signo
    fila["clasificacion"] = clasificacion
    filas.append(fila)

resultados = pd.DataFrame(filas)
csv_path = os.path.join(RESULTS_IRREV, "gene_loo_robustness.csv")
resultados.to_csv(csv_path, index=False)
print(f"\nGuardado: {csv_path}")

# --- Resumen ---
n_total = len(resultados)
n_robusto = (resultados.clasificacion == "Robusto").sum()
n_parcial = (resultados.clasificacion == "Parcialmente robusto").sum()
n_no_robusto = (resultados.clasificacion == "No robusto").sum()

print(f"\nRobusto: {n_robusto} | Parcialmente robusto: {n_parcial} | No robusto: {n_no_robusto} (de {n_total})")

# --- gene_loo_summary.md: texto generado deterministicamente por el
# script a partir de los numeros calculados arriba, no por un modelo
# de lenguaje - mismo principio que gene_loo_robustness.csv. ---
genes_robustos = resultados[resultados.clasificacion == "Robusto"]["gene"].tolist()
genes_parciales = resultados[resultados.clasificacion == "Parcialmente robusto"]["gene"].tolist()
genes_no_robustos = resultados[resultados.clasificacion == "No robusto"]["gene"].tolist()

md_lines = [
    "# Validación leave-one-animal-out de genes asociados al compromiso de destino",
    "",
    "## Metodología",
    "",
    f"Se evaluaron los {n_total} genes que sustentan la Figura 12, la Tabla 4 y la "
    f"discusión biológica del apartado 5.5 (los {N_TOP_POR_DIRECCION} de mayor |ρ| en "
    "cada dirección, seleccionados dinámicamente desde `driver_genes_alta_confianza.csv`). "
    f"Para cada gen se recalculó la correlación de Spearman entre expresión génica y "
    f"entropía normalizada de destino excluyendo un animal a la vez "
    f"({', '.join(animales)}, intervalo 15→30 semanas), con corrección FDR "
    "(Benjamini-Hochberg) recalculada dentro de cada subconjunto de "
    f"{n_total} genes.",
    "",
    "## Criterios de clasificación",
    "",
    f"- **Robusto**: signo idéntico al original en las 3 exclusiones y significativo "
    f"(FDR<{FDR_THRESHOLD}) en las 3.",
    f"- **Parcialmente robusto**: signo idéntico en las 3 exclusiones, pierde "
    f"significancia (FDR≥{FDR_THRESHOLD}) en al menos una.",
    "- **No robusto**: el signo se invierte en al menos una exclusión.",
    "",
    "No se introdujo ningún umbral numérico nuevo: el criterio de significancia "
    f"(FDR<{FDR_THRESHOLD}) es el mismo ya establecido y utilizado en el Script 07 "
    "para la identificación original de estos genes.",
    "",
    "## Resultados",
    "",
    f"- Genes evaluados: **{n_total}**",
    f"- Robustos: **{n_robusto}** de {n_total} ({100*n_robusto/n_total:.1f}%)",
    f"- Parcialmente robustos: **{n_parcial}** de {n_total} ({100*n_parcial/n_total:.1f}%)",
    f"- No robustos: **{n_no_robusto}** de {n_total} ({100*n_no_robusto/n_total:.1f}%)",
    "",
    f"**Genes robustos**: {', '.join(genes_robustos) if genes_robustos else 'ninguno'}.",
    "",
    f"**Genes parcialmente robustos**: {', '.join(genes_parciales) if genes_parciales else 'ninguno'}.",
    "",
    f"**Genes no robustos**: {', '.join(genes_no_robustos) if genes_no_robustos else 'ninguno'}.",
    "",
    "## Interpretación metodológica",
    "",
    "Un gen clasificado como robusto muestra el mismo patrón de asociación con el "
    "compromiso de destino independientemente de qué animal se excluya del análisis, "
    "lo que reduce la probabilidad de que el hallazgo sea un artefacto de un único "
    "individuo, en línea con la preocupación de pseudorreplicación documentada por "
    "Squair et al. (2021) y ya aplicada a la métrica de entropía en el apartado 5.4. "
    "Un gen no robusto no debe descartarse como falso positivo sin más evidencia, pero "
    "su asociación debe interpretarse con mayor cautela que la de un gen robusto.",
    "",
]

md_path = os.path.join(RESULTS_IRREV, "gene_loo_summary.md")
with open(md_path, "w", encoding="utf-8") as f:
    f.write("\n".join(md_lines))
print(f"Guardado: {md_path}")

# --- Heatmap ---
plot_df = resultados.set_index("gene")[
    ["spearman_rho_original"] + [f"rho_sin_{a}" for a in animales]
]
orden_clasificacion = {"Robusto": 0, "Parcialmente robusto": 1, "No robusto": 2}
orden_genes = resultados.assign(
    _orden=resultados.clasificacion.map(orden_clasificacion)
).sort_values(["_orden", "spearman_rho_original"], ascending=[True, False])["gene"]
plot_df = plot_df.loc[orden_genes]

fig, ax = plt.subplots(figsize=(6, 0.35 * len(plot_df) + 1.5))
im = ax.imshow(plot_df.values, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
ax.set_xticks(range(len(plot_df.columns)))
ax.set_xticklabels(["ρ original"] + [f"sin {a}" for a in animales], rotation=30, ha="right")
ax.set_yticks(range(len(plot_df)))
ax.set_yticklabels(plot_df.index, style="italic")
for i, gene in enumerate(plot_df.index):
    clasif = resultados.set_index("gene").loc[gene, "clasificacion"]
    color_borde = {"Robusto": "#2ca02c", "Parcialmente robusto": "#ff7f0e", "No robusto": "#d62728"}[clasif]
    ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(plot_df.columns), 1,
                                fill=False, edgecolor=color_borde, linewidth=2))
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("ρ de Spearman")
ax.set_title("Estabilidad leave-one-animal-out de los 20 genes principales\n"
              "(borde verde=robusto, naranja=parcial, rojo=no robusto)", fontsize=10)
plt.tight_layout()
heatmap_path = os.path.join(FIGURES_MAIN, "gene_loo_heatmap.pdf")
plt.savefig(heatmap_path, bbox_inches="tight")
plt.close()
print(f"Guardado: {heatmap_path}")

print("\nScript 07b completo.")
