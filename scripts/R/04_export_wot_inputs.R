#####################################################################
# 04_export_wot_inputs.R
#
# Exportar archivos de entrada para WOT (Python).
#
# INPUT:  CHECKPOINTS/checkpoint_03_hepatocytes.rds (7.682 hepatocitos,
#         11 estados, 3 timepoints reales: 0/15/30 semanas)
# OUTPUT: PROCESSED/{cell_days,growth_rates,expression_matrix}.txt
#         PROCESSED/{cell_metadata,matrix,umap_coords}.csv
#         LOGS/04_export_wot_inputs_<timestamp>.log
#
# Growth rates: transformacion sigmoide de scores de proliferacion/
# apoptosis a tasas de nacimiento/muerte, combinadas como factor
# multiplicativo exp(beta-delta) (Schiebinger et al. 2019).
#####################################################################

source(here::here("scripts", "R", "config.R"))

library(Seurat)
library(dplyr)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(LOGS, paste0("04_export_wot_inputs_", ts, ".log"))
sink(log_file, split = TRUE)
cat("Script 04 - Exportacion de inputs para WOT\n")
cat("Timestamp:", ts, "\n\n")

hep_obj <- readRDS(file.path(CHECKPOINTS, "checkpoint_03_hepatocytes.rds"))
cat("Objeto cargado:", ncol(hep_obj), "hepatocitos x", nrow(hep_obj), "genes\n")
stopifnot(ncol(hep_obj) == 7682)
stopifnot(identical(sort(unique(hep_obj$week)), as.numeric(HEPATOCYTE_TIMEPOINTS_REAL)))
cat("Timepoints confirmados: 0, 15, 30 semanas\n\n")

# --- Growth rates ---
PROLIFERATION_MARKERS <- c("Mki67", "Pcna", "Top2a", "Mcm2", "Ccnb1", "Cdk1")
APOPTOSIS_MARKERS <- c("Casp3", "Casp8", "Casp9", "Bax", "Bak1", "Bad")
BETA_MAX <- 1.7; BETA_MIN <- 0.3
DELTA_MAX <- 1.7; DELTA_MIN <- 0.3

prolif_presentes <- intersect(PROLIFERATION_MARKERS, rownames(hep_obj))
apop_presentes <- intersect(APOPTOSIS_MARKERS, rownames(hep_obj))
cat("Marcadores de proliferacion presentes:", paste(prolif_presentes, collapse = ", "), "\n")
cat("Marcadores de apoptosis presentes:", paste(apop_presentes, collapse = ", "), "\n\n")

if (length(prolif_presentes) < 3 || length(apop_presentes) < 3) {
  warning("Pocos marcadores presentes (<3) - el score de growth rate ",
          "puede ser poco robusto. Revisar antes de interpretar resultados de WOT.")
}

hep_obj <- AddModuleScore(hep_obj, features = list(prolif_presentes),
                           name = "ProlifScore", seed = SEED_GLOBAL)
hep_obj <- AddModuleScore(hep_obj, features = list(apop_presentes),
                           name = "ApopScore", seed = SEED_GLOBAL)

sigmoid_to_rate <- function(score, rate_min, rate_max) {
  rate_min + (rate_max - rate_min) / (1 + exp(-score))
}

hep_obj$birth_rate <- sigmoid_to_rate(hep_obj$ProlifScore1, BETA_MIN, BETA_MAX)
hep_obj$death_rate <- sigmoid_to_rate(hep_obj$ApopScore1, DELTA_MIN, DELTA_MAX)
# Factor multiplicativo positivo, no diferencia lineal: una diferencia
# puede dar negativo y rompe la matematica interna de wot (potencias
# fraccionarias sobre este valor). tau = ln(2)/(beta-delta) implica
# que el factor de crecimiento real es exp(beta-delta).
hep_obj$growth_rate <- exp(hep_obj$birth_rate - hep_obj$death_rate)

cat("Resumen de growth rates calculados:\n")
print(summary(hep_obj$growth_rate))
cat("\nGrowth rate medio por condicion:\n")
print(hep_obj@meta.data %>% group_by(week) %>%
        summarise(growth_medio = mean(growth_rate), .groups = "drop"))
cat("\n")

# --- Exportacion ---
# Extension .txt, no .csv: wot solo reconoce mtx/txt/h5ad/loom.
cell_days <- data.frame(id = colnames(hep_obj), day = hep_obj$week)
write.csv(cell_days, file.path(PROCESSED, "cell_days.txt"), row.names = FALSE)
cat("cell_days.txt exportado:", nrow(cell_days), "celulas\n")

cell_metadata <- hep_obj@meta.data %>%
  mutate(id = colnames(hep_obj)) %>%
  select(id, week, seurat_clusters, fraction, animal, capture,
         cell_type_official, growth_rate, birth_rate, death_rate,
         nCount_RNA, nFeature_RNA, percent.mt)
write.csv(cell_metadata, file.path(PROCESSED, "cell_metadata.csv"), row.names = FALSE)
cat("cell_metadata.csv exportado:", nrow(cell_metadata), "celulas,",
    ncol(cell_metadata), "columnas\n")

# Columnas "id"/"cell_growth_rate" exactas requeridas por wot.
growth_rates_export <- data.frame(id = colnames(hep_obj),
                                   cell_growth_rate = hep_obj$growth_rate)
write.csv(growth_rates_export, file.path(PROCESSED, "growth_rates.txt"), row.names = FALSE)
cat("growth_rates.txt exportado\n")

pca_coords <- as.data.frame(Embeddings(hep_obj, "pca"))
pca_coords <- cbind(id = rownames(pca_coords), pca_coords)
write.csv(pca_coords, file.path(PROCESSED, "matrix.csv"), row.names = FALSE)
cat("matrix.csv exportado (PCA", N_PCS_HEPATOCYTES, "D)\n")

umap_coords <- as.data.frame(Embeddings(hep_obj, "umap"))
umap_coords <- cbind(id = rownames(umap_coords), umap_coords)
write.csv(umap_coords, file.path(PROCESSED, "umap_coords.csv"), row.names = FALSE)
cat("umap_coords.csv exportado\n")

# Genes MT excluidos: su %MT varia sistematicamente con la condicion
# (Script 01), lo que podria confundir senal tecnica con senal
# genuina de estado celular en el PCA local que WOT calcula.
hvg <- setdiff(VariableFeatures(hep_obj), MT_GENES_MOUSE)
cat("HVG totales:", length(VariableFeatures(hep_obj)),
    "| excluyendo", length(VariableFeatures(hep_obj)) - length(hvg),
    "genes MT:", length(hvg), "\n")
expr_mat <- as.data.frame(t(as.matrix(GetAssayData(hep_obj, layer = "data")[hvg, ])))
expr_mat <- cbind(id = rownames(expr_mat), expr_mat)
write.csv(expr_mat, file.path(PROCESSED, "expression_matrix.txt"), row.names = FALSE)
cat("expression_matrix.txt exportado:", length(hvg), "HVG x", ncol(hep_obj), "celulas\n\n")

n_ids <- length(unique(cell_days$id))
stopifnot(n_ids == nrow(cell_metadata))
stopifnot(n_ids == nrow(growth_rates_export))
stopifnot(n_ids == nrow(pca_coords))
stopifnot(n_ids == nrow(umap_coords))
stopifnot(n_ids == nrow(expr_mat))
cat("Verificacion de consistencia: los", n_ids, "IDs de celula coinciden\n")
cat("exactamente en los 6 archivos exportados.\n\n")

cat("Archivos guardados en:", PROCESSED, "\n\n")
cat("---- sessionInfo() ----\n")
print(sessionInfo())
sink()
cat("Script 04 completo. Log:", log_file, "\n")
