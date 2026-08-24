#####################################################################
# 02_normalization_clustering.R
#
# Normalizacion, HVGs, PCA, UMAP y clustering global.
#
# INPUT:  CHECKPOINTS/checkpoint_01_qc_filtered.rds (65.158 celulas,
#         25.127 genes)
# OUTPUT: CHECKPOINTS/checkpoint_02_clustered.rds
#         RESULTS_CLUST/elbow_plot.pdf, cluster_composition.csv
#         FIGURES_MAIN/umap_por_condicion.pdf, umap_por_cluster.pdf,
#                      umap_por_fraccion.pdf
#         LOGS/02_normalization_clustering_<timestamp>.log
#####################################################################

source(here::here("scripts", "R", "config.R"))

library(Seurat)
library(dplyr)
library(ggplot2)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(LOGS, paste0("02_normalization_clustering_", ts, ".log"))
sink(log_file, split = TRUE)
cat("Script 02 - Normalizacion y clustering global\n")
cat("Timestamp:", ts, "\n\n")

seurat_obj <- readRDS(file.path(CHECKPOINTS, "checkpoint_01_qc_filtered.rds"))
cat("Objeto cargado:", ncol(seurat_obj), "celulas x", nrow(seurat_obj), "genes\n")
stopifnot(ncol(seurat_obj) == 65158)

cat("\nDistribucion por condicion y fraccion:\n")
print(table(seurat_obj$week, seurat_obj$fraction))
cat("\n")

# --- Normalizacion ---
# LogNormalize, no SCTransform: preserva la senal temporal que WOT
# necesita mas adelante (SCTransform introduce una regresion de
# varianza que puede atenuarla).
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize",
                             scale.factor = 10000, verbose = FALSE)
cat("Normalizacion (LogNormalize) completa.\n")

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                    nfeatures = N_HVG, verbose = FALSE)
cat("HVGs seleccionados:", length(VariableFeatures(seurat_obj)), "(objetivo:", N_HVG, ")\n\n")

# --- PCA ---
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj),
                      npcs = N_PCS_FULL, seed.use = SEED_PCA, verbose = FALSE)

# Varianza real: el denominador correcto es la varianza total de los
# HVG (misc$total.variance, calculado por Seurat internamente), no la
# suma de los propios PCs retenidos (eso siempre da 100% por
# definicion y no dice nada sobre la varianza real capturada).
var_explicada <- seurat_obj[["pca"]]@stdev^2
total_variance_real <- seurat_obj[["pca"]]@misc$total.variance
pct_var <- var_explicada / total_variance_real * 100
cum_var <- cumsum(pct_var)
cat("Varianza TOTAL en el espacio de", N_HVG, "HVG:", round(total_variance_real, 1), "\n")
cat("Varianza acumulada REAL en PC", N_PCS_FULL, ":", round(cum_var[N_PCS_FULL], 1), "%\n\n")

p_elbow <- ElbowPlot(seurat_obj, ndims = N_PCS_FULL) +
  ggtitle(sprintf("Elbow plot - %.1f%% varianza acumulada en %d PCs",
                   cum_var[N_PCS_FULL], N_PCS_FULL))
ggsave(file.path(RESULTS_CLUST, "elbow_plot.pdf"), p_elbow, width = 7, height = 5)

# --- UMAP ---
seurat_obj <- RunUMAP(seurat_obj, dims = 1:N_PCS_FULL,
                       n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST,
                       seed.use = SEED_UMAP, verbose = FALSE)
cat("UMAP calculado (n.neighbors =", UMAP_N_NEIGHBORS, ", min.dist =", UMAP_MIN_DIST, ")\n\n")

# --- Clustering ---
# Sin correccion de batch: el batch coincide con el tiempo biologico
# de interes; corregirlo eliminaria la senal que WOT necesita modelar.
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:N_PCS_FULL, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = CLUSTERING_RESOLUTION_FULL,
                            random.seed = SEED_CLUSTERING, verbose = FALSE)

n_clusters <- length(unique(Idents(seurat_obj)))
cat("Clusters encontrados (resolucion =", CLUSTERING_RESOLUTION_FULL, "):", n_clusters, "\n\n")
cat("Tamano de cada cluster:\n")
print(table(Idents(seurat_obj)))
cat("\n")

# --- Validacion cruzada contra anotacion oficial de Su et al. ---
cross_tab <- table(cluster = Idents(seurat_obj), cell_type = seurat_obj$cell_type_official)
cat("Composicion de cada cluster por CellType oficial (Su et al.):\n")
print(cross_tab)

pureza_por_cluster <- apply(cross_tab, 1, function(x) max(x) / sum(x) * 100)
cat("\nPureza por cluster (%, tipo oficial dominante / total del cluster):\n")
print(round(pureza_por_cluster, 1))
cat("\n")

write.csv(as.data.frame.matrix(cross_tab),
          file.path(RESULTS_CLUST, "cluster_composition.csv"))

# --- Visualizaciones ---
# CORREGIDO: "condición" y "fracción" con tilde (antes "condicion" y
# "fraccion", inconsistente con la ortografia del resto del TFM).
p_cluster <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters",
                      label = TRUE) + ggtitle("UMAP por cluster")
ggsave(file.path(FIGURES_MAIN, "umap_por_cluster.pdf"), p_cluster, width = 8, height = 6)

p_condicion <- DimPlot(seurat_obj, reduction = "umap", group.by = "week") +
  ggtitle("UMAP por condición temporal (semanas)")
ggsave(file.path(FIGURES_MAIN, "umap_por_condicion.pdf"), p_condicion, width = 8, height = 6)

p_fraccion <- DimPlot(seurat_obj, reduction = "umap", group.by = "fraction") +
  ggtitle("UMAP por fracción (Hepatocyte vs NPC)")
ggsave(file.path(FIGURES_MAIN, "umap_por_fraccion.pdf"), p_fraccion, width = 8, height = 6)

cat("Figuras guardadas en", FIGURES_MAIN, "y", RESULTS_CLUST, "\n\n")

saveRDS(seurat_obj, file.path(CHECKPOINTS, "checkpoint_02_clustered.rds"))
cat("Checkpoint guardado: checkpoint_02_clustered.rds\n\n---- sessionInfo() ----\n")
print(sessionInfo())
sink()
cat("Script 02 completo. Log:", log_file, "\n")
