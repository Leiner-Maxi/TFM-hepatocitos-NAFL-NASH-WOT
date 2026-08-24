"""
07_driver_gene_analysis.py
Correlacion de Spearman entre expresion genica (celulas origen, 15
semanas) y entropia normalizada de destino (transport map 15->30sem).

Se usa Spearman, no GLM binomial negativo, porque expression_matrix.txt
contiene datos normalizados/log-transformados, no conteos crudos.
Correccion FDR (Benjamini-Hochberg) sobre las pruebas.

Input:  results/irreversibility/entropy_por_celula.csv (Script 06)
        data/processed/expression_matrix.txt
Output: results/irreversibility/driver_genes_full.csv
        results/irreversibility/driver_genes_top.csv (FDR<0.05)
        results/irreversibility/driver_genes_alta_confianza.csv
        (FDR<0.05 y |rho| por encima del umbral de tamano de efecto)
"""
import os
import sys
import datetime
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from config import PROCESSED, RESULTS_IRREV, PROJECT_ROOT, DRIVER_GENE_EFFECT_SIZE_THRESHOLD

LOGS = os.path.join(PROJECT_ROOT, "logs")
ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = os.path.join(LOGS, f"07_driver_gene_analysis_{ts}.log")
log_lines = []


def log(msg):
    print(msg)
    log_lines.append(str(msg))


log("Script 07 - Identificacion de genes asociados al compromiso de destino")
log(f"Timestamp: {ts}\n")

entropy_path = os.path.join(RESULTS_IRREV, "entropy_por_celula.csv")
if not os.path.exists(entropy_path):
    log(f"ERROR CRITICO: no existe {entropy_path}. Correr Script 06 primero.")
    sys.exit(1)

df_entropy = pd.read_csv(entropy_path)
df_15_30 = df_entropy[df_entropy["intervalo"] == "15_30"].copy()
log(f"Celulas origen (15 semanas) con entropia calculada: {len(df_15_30)}")

n_estado7 = (df_15_30["seurat_clusters"] == 7).sum()
log(f"De estas, {n_estado7} pertenecen al Estado 7 (entropia media mas baja).\n")

expr = pd.read_csv(os.path.join(PROCESSED, "expression_matrix.txt"))
expr = expr.set_index("id")

common_ids = df_15_30["id"].isin(expr.index)
log(f"Celulas del intervalo 15->30 encontradas en expression_matrix.txt: "
    f"{common_ids.sum()} de {len(df_15_30)}")
df_15_30 = df_15_30[common_ids.values].reset_index(drop=True)

expr_subset = expr.loc[df_15_30["id"]]
assert list(expr_subset.index) == list(df_15_30["id"]), \
    "Los IDs de expr_subset y df_15_30 deben coincidir en orden exacto"

gene_names = expr_subset.columns.tolist()
expr_matrix = expr_subset.values
entropy_values = df_15_30["entropy_norm"].values

log(f"Matriz de expresion alineada: {expr_matrix.shape[0]} celulas x "
    f"{expr_matrix.shape[1]} genes\n")

log("Calculando correlacion de Spearman gen-por-gen...")
n_genes = expr_matrix.shape[1]
rho_values = np.zeros(n_genes)
p_values = np.zeros(n_genes)

for i in range(n_genes):
    gene_expr = expr_matrix[:, i]
    if np.std(gene_expr) == 0:
        rho_values[i] = np.nan
        p_values[i] = np.nan
        continue
    rho, p = spearmanr(gene_expr, entropy_values)
    rho_values[i] = rho
    p_values[i] = p

results = pd.DataFrame({
    "gene": gene_names,
    "spearman_rho": rho_values,
    "p_value": p_values
}).dropna()

log(f"Genes con varianza cero (excluidos): {n_genes - len(results)}")
log(f"Genes evaluados: {len(results)}\n")


def bh_fdr(p_vals):
    """Benjamini-Hochberg manual, evita dependencia de statsmodels."""
    n = len(p_vals)
    order = np.argsort(p_vals)
    ranked_p = p_vals[order]
    q_ranked = ranked_p * n / (np.arange(n) + 1)
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    q_ranked = np.clip(q_ranked, 0, 1)
    q_vals = np.empty(n)
    q_vals[order] = q_ranked
    return q_vals


results["q_value"] = bh_fdr(results["p_value"].values)
results = results.sort_values("p_value").reset_index(drop=True)

n_sig = (results["q_value"] < 0.05).sum()
log(f"Genes significativos a FDR < 0.05: {n_sig} de {len(results)}\n")

results.to_csv(os.path.join(RESULTS_IRREV, "driver_genes_full.csv"), index=False)

sig_results = results[results["q_value"] < 0.05].copy()
sig_results.to_csv(os.path.join(RESULTS_IRREV, "driver_genes_top.csv"), index=False)

top_negative = sig_results.sort_values("spearman_rho").head(15)
top_positive = sig_results.sort_values("spearman_rho", ascending=False).head(15)

# rho negativo: mayor expresion -> menor entropia -> mayor compromiso.
# rho positivo: mayor expresion -> mayor entropia -> mayor plasticidad.
log("Top 15 genes asociados a compromiso de destino:")
log(top_negative[["gene", "spearman_rho", "q_value"]].to_string(index=False))
log("")

log("Top 15 genes asociados a mantener plasticidad:")
log(top_positive[["gene", "spearman_rho", "q_value"]].to_string(index=False))
log("")

human_risk_genes_mouse = ["Pnpla3", "Tm6sf2", "Hsd17b13", "Mboat7"]
log("Cruce con genes de riesgo NAFLD/NASH humanos conocidos:")
for g in human_risk_genes_mouse:
    match = results[results["gene"] == g]
    if len(match) > 0:
        row = match.iloc[0]
        log(f"  {g}: rho={row.spearman_rho:.4f}, q={row.q_value:.4f}"
            f"{' (significativo)' if row.q_value < 0.05 else ' (no significativo)'}")
    else:
        log(f"  {g}: no presente en los {len(results)} HVG evaluados")
log("")

# Filtro de tamano de efecto: separa significancia estadistica de
# relevancia biologica. Con n~4481 celulas, FDR<0.05 por si solo
# retiene demasiados genes para ser interpretable; esta lista final
# es la que se reporta como hallazgo principal en Resultados.
driver_genes_alta_confianza = results[
    (results["q_value"] < 0.05) &
    (results["spearman_rho"].abs() > DRIVER_GENE_EFFECT_SIZE_THRESHOLD)
].copy()
driver_genes_alta_confianza = driver_genes_alta_confianza.sort_values(
    "spearman_rho", key=lambda x: x.abs(), ascending=False
)
driver_genes_alta_confianza.to_csv(
    os.path.join(RESULTS_IRREV, "driver_genes_alta_confianza.csv"), index=False
)

log(f"Filtro de tamano de efecto (|rho| > {DRIVER_GENE_EFFECT_SIZE_THRESHOLD}):")
log(f"Genes con FDR<0.05 solamente: {len(sig_results)} de {len(results)}")
log(f"Genes con FDR<0.05 Y |rho|>{DRIVER_GENE_EFFECT_SIZE_THRESHOLD}: "
    f"{len(driver_genes_alta_confianza)} de {len(results)}\n")

with open(log_file, "w") as f:
    f.write("\n".join(log_lines))

log(f"Script 07 completo. Log: {log_file}")
