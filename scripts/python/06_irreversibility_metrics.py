"""
06_irreversibility_metrics.py
Metrica de irreversibilidad (fate commitment): entropia de Shannon
normalizada de la distribucion de destino de cada celula, bootstrap,
test de permutacion y control ortogonal.

Input:  results/wot/tmaps/hepatocytes_*.h5ad, data/processed/cell_metadata.csv,
        data/processed/expression_matrix.txt
Output: results/irreversibility/{entropy_por_celula,entropy_resumen,
        bootstrap_ci,permutation_test,raw_expression_entropy}.csv

Nota de lenguaje: esto mide "fate commitment" (determinismo del destino
hacia adelante), no reversibilidad literal en sentido de transporte
inverso independiente.
"""
import os
import sys
import datetime
import numpy as np
import pandas as pd
import anndata

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from config import PROCESSED, RESULTS_WOT, RESULTS_IRREV, PROJECT_ROOT, N_BOOTSTRAP, N_PERMUTATIONS, SEED_GLOBAL

LOGS = os.path.join(PROJECT_ROOT, "logs")
os.makedirs(RESULTS_IRREV, exist_ok=True)

ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = os.path.join(LOGS, f"06_irreversibility_metrics_{ts}.log")
log_lines = []


def log(msg):
    print(msg)
    log_lines.append(str(msg))


np.random.seed(SEED_GLOBAL)

log("Script 06 - Metricas de irreversibilidad (fate commitment)")
log(f"Timestamp: {ts}\n")

# --- Cargar transport maps ---
tmap_dir = os.path.join(RESULTS_WOT, "tmaps")
path_0_15 = os.path.join(tmap_dir, "hepatocytes_0.0_15.0.h5ad")
path_15_30 = os.path.join(tmap_dir, "hepatocytes_15.0_30.0.h5ad")

for p in [path_0_15, path_15_30]:
    if not os.path.exists(p):
        log(f"ERROR CRITICO: no existe {p}")
        log("Verificar que Script 05 se corrio completo.")
        sys.exit(1)

tmap_0_15 = anndata.read_h5ad(path_0_15)
tmap_15_30 = anndata.read_h5ad(path_15_30)

log(f"Transport map 0->15sem: {tmap_0_15.shape[0]} celulas origen x "
    f"{tmap_0_15.shape[1]} celulas destino")
log(f"Transport map 15->30sem: {tmap_15_30.shape[0]} celulas origen x "
    f"{tmap_15_30.shape[1]} celulas destino\n")


def row_entropy(X):
    """Entropia de Shannon por fila, normalizando cada fila a una
    distribucion de probabilidad. Filas con suma cero -> NaN."""
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    row_sums = X.sum(axis=1, keepdims=True)
    entropies = np.full(X.shape[0], np.nan)
    valid = row_sums.flatten() > 0
    P = np.zeros_like(X)
    P[valid] = X[valid] / row_sums[valid]
    with np.errstate(divide="ignore", invalid="ignore"):
        logP = np.where(P > 0, np.log(P), 0)
    entropies[valid] = -(P[valid] * logP[valid]).sum(axis=1)
    return entropies


log("Calculando entropia por celula (intervalo 0->15sem)...")
entropy_0_15 = row_entropy(tmap_0_15.X)
log("Calculando entropia por celula (intervalo 15->30sem)...")
entropy_15_30 = row_entropy(tmap_15_30.X)

df_entropy_0_15 = pd.DataFrame({
    "id": tmap_0_15.obs_names, "entropy_raw": entropy_0_15, "intervalo": "0_15"
})
df_entropy_15_30 = pd.DataFrame({
    "id": tmap_15_30.obs_names, "entropy_raw": entropy_15_30, "intervalo": "15_30"
})

# Normalizacion critica: los dos transport maps tienen distinto numero
# de celulas destino (techo teorico log(N) distinto para cada uno).
# Sin normalizar, parte de la diferencia de entropia seria artefacto
# del tamano del espacio de destino, no senal biologica.
n_destino_0_15 = tmap_0_15.shape[1]
n_destino_15_30 = tmap_15_30.shape[1]
log_max_0_15 = np.log(n_destino_0_15)
log_max_15_30 = np.log(n_destino_15_30)

df_entropy_0_15["entropy_norm"] = df_entropy_0_15["entropy_raw"] / log_max_0_15
df_entropy_15_30["entropy_norm"] = df_entropy_15_30["entropy_raw"] / log_max_15_30

log(f"Techo teorico de entropia: log({n_destino_0_15})={log_max_0_15:.4f} "
    f"(0->15sem) vs log({n_destino_15_30})={log_max_15_30:.4f} (15->30sem)\n")

df_entropy_all = pd.concat([df_entropy_0_15, df_entropy_15_30], ignore_index=True)

n_nan_0_15 = df_entropy_0_15["entropy_raw"].isna().sum()
n_nan_15_30 = df_entropy_15_30["entropy_raw"].isna().sum()
log(f"Celulas con masa transportada cero (excluidas): "
    f"{n_nan_0_15} en 0->15sem, {n_nan_15_30} en 15->30sem\n")

# --- Unir con estados transcriptomicos (Script 03) ---
cell_meta = pd.read_csv(os.path.join(PROCESSED, "cell_metadata.csv"))
df_entropy_all = df_entropy_all.merge(
    cell_meta[["id", "week", "seurat_clusters"]], on="id", how="left"
)
df_entropy_all.to_csv(os.path.join(RESULTS_IRREV, "entropy_por_celula.csv"), index=False)

# --- Resumen: entropia normalizada es la metrica principal ---
resumen_intervalo = df_entropy_all.groupby("intervalo")[["entropy_norm", "entropy_raw"]].agg(
    ["mean", "median", "std", "count"]
).reset_index()
log("Entropia normalizada por intervalo:")
log(resumen_intervalo.to_string(index=False))
log("")

resumen_estado = df_entropy_all.groupby(["intervalo", "seurat_clusters"])["entropy_norm"].agg(
    ["mean", "count"]
).reset_index()
log("Entropia normalizada por estado y por intervalo:")
log(resumen_estado.to_string(index=False))
log("")

resumen_intervalo.to_csv(os.path.join(RESULTS_IRREV, "entropy_resumen.csv"))
resumen_estado.to_csv(os.path.join(RESULTS_IRREV, "entropy_por_estado.csv"), index=False)

mean_0_15 = df_entropy_all.loc[df_entropy_all.intervalo == "0_15", "entropy_norm"].mean()
mean_15_30 = df_entropy_all.loc[df_entropy_all.intervalo == "15_30", "entropy_norm"].mean()
diff_observada = mean_15_30 - mean_0_15
log(f"Entropia normalizada media 0->15sem: {mean_0_15:.4f}")
log(f"Entropia normalizada media 15->30sem: {mean_15_30:.4f}")
log(f"Diferencia (15_30 - 0_15): {diff_observada:.4f}\n")

# --- Bootstrap ---
def bootstrap_mean_ci(values, n_boot=N_BOOTSTRAP, alpha=0.05, seed=SEED_GLOBAL):
    values = values.dropna().values
    rng = np.random.default_rng(seed)
    boot_means = np.array([
        rng.choice(values, size=len(values), replace=True).mean()
        for _ in range(n_boot)
    ])
    lo, hi = np.percentile(boot_means, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return values.mean(), lo, hi


log(f"Bootstrap ({N_BOOTSTRAP} iteraciones), IC 95%:")
bootstrap_results = []
for interval_name, group in df_entropy_all.groupby("intervalo"):
    m, lo, hi = bootstrap_mean_ci(group["entropy_norm"])
    bootstrap_results.append({"intervalo": interval_name, "media": m, "ci_low": lo, "ci_high": hi})
    log(f"  {interval_name}: media={m:.4f}, IC95%=[{lo:.4f}, {hi:.4f}]")
pd.DataFrame(bootstrap_results).to_csv(os.path.join(RESULTS_IRREV, "bootstrap_ci.csv"), index=False)
log("")

# --- Test de permutacion ---
log(f"Test de permutacion ({N_PERMUTATIONS} iteraciones):")
combined = df_entropy_all.dropna(subset=["entropy_norm"]).copy()
n_0_15 = (combined.intervalo == "0_15").sum()

rng = np.random.default_rng(SEED_GLOBAL)
values = combined["entropy_norm"].values
perm_diffs = np.zeros(N_PERMUTATIONS)
for i in range(N_PERMUTATIONS):
    shuffled = rng.permutation(values)
    perm_diffs[i] = shuffled[n_0_15:].mean() - shuffled[:n_0_15].mean()

p_value = (np.abs(perm_diffs) >= np.abs(diff_observada)).mean()
log(f"Diferencia observada: {diff_observada:.4f}")
log(f"p-valor (permutacion, dos colas): {p_value:.4f}")

with open(os.path.join(RESULTS_IRREV, "permutation_test.txt"), "w") as pf:
    pf.write(f"Diferencia observada: {diff_observada:.6f}\n")
    pf.write(f"p-valor (permutacion, {N_PERMUTATIONS} iteraciones): {p_value:.6f}\n")
log("")

# --- Control ortogonal: entropia del transcriptoma crudo, independiente de WOT ---
log("Control ortogonal: entropia per-celula del transcriptoma crudo:")
expr = pd.read_csv(os.path.join(PROCESSED, "expression_matrix.txt"))
expr_ids = expr["id"]
expr_vals = expr.drop(columns=["id"]).values
expr_vals = np.clip(expr_vals, 0, None)

row_sums = expr_vals.sum(axis=1, keepdims=True)
valid = row_sums.flatten() > 0
P = np.zeros_like(expr_vals)
P[valid] = expr_vals[valid] / row_sums[valid]
with np.errstate(divide="ignore", invalid="ignore"):
    logP = np.where(P > 0, np.log(P), 0)
raw_entropy = np.full(expr_vals.shape[0], np.nan)
raw_entropy[valid] = -(P[valid] * logP[valid]).sum(axis=1)

df_raw_entropy = pd.DataFrame({"id": expr_ids, "raw_entropy": raw_entropy})
df_raw_entropy = df_raw_entropy.merge(cell_meta[["id", "week"]], on="id", how="left")

resumen_raw = df_raw_entropy.groupby("week")["raw_entropy"].agg(["mean", "count"]).reset_index()
log(resumen_raw.to_string(index=False))
log("")

df_raw_entropy.to_csv(os.path.join(RESULTS_IRREV, "raw_expression_entropy.csv"), index=False)

# --- Growth rate aprendido por WOT vs estimado ---
g_path = os.path.join(tmap_dir, "hepatocytes_g.txt")
if os.path.exists(g_path):
    learned_g = pd.read_csv(g_path, sep=None, engine="python")
    log("Growth rate aprendido por WOT (tras growth_iters) vs estimado:")
    log(f"Columnas: {list(learned_g.columns)}")
    log(learned_g.describe().to_string())
    learned_g.to_csv(os.path.join(RESULTS_IRREV, "growth_rate_aprendido.csv"), index=False)
    log("")

with open(log_file, "w") as f:
    f.write("\n".join(log_lines))

log(f"\nScript 06 completo. Log: {log_file}")
