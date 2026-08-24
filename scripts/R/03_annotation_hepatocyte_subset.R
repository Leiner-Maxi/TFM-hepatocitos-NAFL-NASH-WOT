#####################################################################
# 03_annotation_hepatocyte_subset.R
#
# Subsetting hepatocitario, re-clustering especifico y verificacion
# de zonacion.
#
# INPUT:  CHECKPOINTS/checkpoint_02_clustered.rds
# OUTPUT: CHECKPOINTS/checkpoint_03_hepatocytes.rds
#         RESULTS_HEP/composicion_por_condicion.csv, elbow_plot_hepatocitos.pdf
#         FIGURES_MAIN/umap_hepatocitos_por_estado.pdf,
#                      umap_hepatocitos_por_condicion.pdf, zonacion_markers.pdf
#         LOGS/03_annotation_hepatocyte_subset_<timestamp>.log
#
# Criterio de seleccion: fraction=="Hepatocyte" & cell_type_official==
# "Hepatocytes". Por construccion excluye la contaminacion NPC_34weeks
# (todas tienen fraction=="NPC"); se agrega ademas un stopifnot como
# doble verificacion.
#####################################################################

source(here::here("scripts", "R", "config.R"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(LOGS, paste0("03_annotation_hepatocyte_subset_", ts, ".log"))
sink(log_file, split = TRUE)
cat("Script 03 - Subsetting hepatocitario\n")
cat("Timestamp:", ts, "\n\n")

seurat_obj <- readRDS(file.path(CHECKPOINTS, "checkpoint_02_clustered.rds"))
cat("Objeto cargado:", ncol(seurat_obj), "celulas x", nrow(seurat_obj), "genes\n")
stopifnot(ncol(seurat_obj) == 65158)

# --- Subsetting hepatocitario ---
n_por_fraction_solo <- sum(seurat_obj$fraction == "Hepatocyte")
n_por_ambos <- sum(seurat_obj$fraction == "Hepatocyte" &
                    seurat_obj$cell_type_official == "Hepatocytes")
cat("Celulas con fraction=='Hepatocyte':", n_por_fraction_solo, "\n")
cat("Celulas con fraction=='Hepatocyte' Y cell_type_official=='Hepatocytes':",
    n_por_ambos, "\n\n")

hep_obj <- subset(seurat_obj, subset = fraction == "Hepatocyte" &
                                        cell_type_official == "Hepatocytes")

n_npc34w_en_subset <- sum(hep_obj$is_npc34w_hep_contamination)
cat("Celulas de contaminacion NPC_34weeks en el subset final:",
    n_npc34w_en_subset, "(debe ser 0)\n")
stopifnot(n_npc34w_en_subset == 0)

weeks_presentes <- sort(unique(hep_obj$week))
cat("Timepoints presentes:", paste(weeks_presentes, collapse = ", "), "\n")
stopifnot(identical(weeks_presentes, as.numeric(HEPATOCYTE_TIMEPOINTS_REAL)))

cat("Composicion final del subset hepatocitario:\n")
print(table(hep_obj$week))
cat("\n")

# --- Re-procesamiento especifico ---
# HVGs, PCA y clustering recalculados solo sobre hepatocitos: los HVG
# globales estarian dominados por marcadores no-hepatocitarios
# (inmune, endotelial), poco informativos aqui.
hep_obj <- FindVariableFeatures(hep_obj, selection.method = "vst",
                                 nfeatures = N_HVG, verbose = FALSE)
hep_obj <- ScaleData(hep_obj, features = VariableFeatures(hep_obj), verbose = FALSE)
hep_obj <- RunPCA(hep_obj, features = VariableFeatures(hep_obj),
                   npcs = N_PCS_HEPATOCYTES, seed.use = SEED_PCA, verbose = FALSE)

total_variance_real <- hep_obj[["pca"]]@misc$total.variance
pct_var <- (hep_obj[["pca"]]@stdev^2) / total_variance_real * 100
cum_var <- cumsum(pct_var)
cat("Varianza acumulada REAL en", N_PCS_HEPATOCYTES, "PCs (hepatocitos):",
    round(cum_var[N_PCS_HEPATOCYTES], 1), "%\n\n")

p_elbow_hep <- ElbowPlot(hep_obj, ndims = N_PCS_HEPATOCYTES) +
  ggtitle(sprintf("Elbow plot hepatocitos - %.1f%% varianza acumulada",
                   cum_var[N_PCS_HEPATOCYTES]))
ggsave(file.path(RESULTS_HEP, "elbow_plot_hepatocitos.pdf"), p_elbow_hep, width = 7, height = 5)

hep_obj <- RunUMAP(hep_obj, dims = 1:N_PCS_HEPATOCYTES,
                    n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST,
                    seed.use = SEED_UMAP, verbose = FALSE)

hep_obj <- FindNeighbors(hep_obj, dims = 1:N_PCS_HEPATOCYTES, verbose = FALSE)
hep_obj <- FindClusters(hep_obj, resolution = CLUSTERING_RESOLUTION_HEP,
                         random.seed = SEED_CLUSTERING, verbose = FALSE)

n_estados <- length(unique(Idents(hep_obj)))
cat("Estados transcriptomicos hepatocitarios encontrados (resolucion =",
    CLUSTERING_RESOLUTION_HEP, "):", n_estados, "\n\n")
cat("Tamano de cada estado:\n")
print(table(Idents(hep_obj)))
cat("\n")

# --- Composicion por condicion (input clave para Scripts 04-06) ---
composicion <- hep_obj@meta.data %>%
  group_by(week, seurat_clusters) %>%
  summarise(n_celulas = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = week, values_from = n_celulas, values_fill = 0)

cat("Composicion de estados por condicion temporal:\n")
print(composicion)
write.csv(composicion, file.path(RESULTS_HEP, "composicion_por_condicion.csv"),
          row.names = FALSE)
cat("\n")

# --- Zonacion ---
genes_perivenosos <- c("Cyp2e1", "Glul", "Oat")
genes_periportales <- c("Cps1", "Ass1", "Cyp2f2")
genes_perivenosos_presentes <- intersect(genes_perivenosos, rownames(hep_obj))
genes_periportales_presentes <- intersect(genes_periportales, rownames(hep_obj))

cat("Marcadores perivenosos presentes:", paste(genes_perivenosos_presentes, collapse = ", "), "\n")
cat("Marcadores periportales presentes:", paste(genes_periportales_presentes, collapse = ", "), "\n\n")

if (length(genes_perivenosos_presentes) > 0 && length(genes_periportales_presentes) > 0) {
  hep_obj <- AddModuleScore(hep_obj, features = list(genes_perivenosos_presentes),
                             name = "PerivenousScore", seed = SEED_GLOBAL)
  hep_obj <- AddModuleScore(hep_obj, features = list(genes_periportales_presentes),
                             name = "PeriportalScore", seed = SEED_GLOBAL)

  zonacion_por_condicion <- hep_obj@meta.data %>%
    group_by(week) %>%
    summarise(perivenoso_medio = mean(PerivenousScore1),
              periportal_medio = mean(PeriportalScore1), .groups = "drop")
  cat("Scores de zonacion por condicion temporal:\n")
  print(zonacion_por_condicion)
  write.csv(zonacion_por_condicion, file.path(RESULTS_HEP, "zonacion_por_condicion.csv"),
            row.names = FALSE)

  # CORREGIDO: combine=FALSE + titulos manuales, en vez de dejar que
  # FeaturePlot use "PerivenousScore1"/"PeriportalScore1" (nombres de
  # columna internos de AddModuleScore) como titulo de panel.
  p_zonacion_list <- FeaturePlot(hep_obj, features = c("PerivenousScore1", "PeriportalScore1"),
                                  reduction = "umap", combine = FALSE)
  p_zonacion_list[[1]] <- p_zonacion_list[[1]] + ggtitle("Puntuación pericentral")
  p_zonacion_list[[2]] <- p_zonacion_list[[2]] + ggtitle("Puntuación periportal")
  p_zonacion <- wrap_plots(p_zonacion_list)
  ggsave(file.path(FIGURES_MAIN, "zonacion_markers.pdf"), p_zonacion, width = 10, height = 5)
  cat("\n")
}

# --- Visualizaciones ---
p_estados <- DimPlot(hep_obj, reduction = "umap", label = TRUE) +
  ggtitle(sprintf("UMAP hepatocitos por estado (%d estados)", n_estados))
ggsave(file.path(FIGURES_MAIN, "umap_hepatocitos_por_estado.pdf"), p_estados, width = 8, height = 6)

# CORREGIDO: "condición" con tilde (antes "condicion").
p_condicion_hep <- DimPlot(hep_obj, reduction = "umap", group.by = "week") +
  ggtitle("UMAP hepatocitos por condición temporal")
ggsave(file.path(FIGURES_MAIN, "umap_hepatocitos_por_condicion.pdf"), p_condicion_hep,
       width = 8, height = 6)

cat("Figuras guardadas en", FIGURES_MAIN, "y", RESULTS_HEP, "\n\n")

saveRDS(hep_obj, file.path(CHECKPOINTS, "checkpoint_03_hepatocytes.rds"))
cat("Checkpoint guardado: checkpoint_03_hepatocytes.rds\n")
cat("Dimensiones finales:", ncol(hep_obj), "hepatocitos x", nrow(hep_obj), "genes,",
    n_estados, "estados transcriptomicos, 3 timepoints reales (0/15/30 semanas)\n\n")
cat("---- sessionInfo() ----\n")
print(sessionInfo())
sink()
cat("Script 03 completo. Log:", log_file, "\n")
