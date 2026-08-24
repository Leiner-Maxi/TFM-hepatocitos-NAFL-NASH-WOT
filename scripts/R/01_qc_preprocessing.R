#####################################################################
# 01_qc_preprocessing.R
#
# Carga de datos crudos GSE166504 y control de calidad.
#
# INPUT:  RAW_DATA/GSE166504_cell_metadata_20220204_tsv.gz
#         RAW_DATA/GSE166504_cell_raw_counts_20220204_txt.gz
# OUTPUT: CHECKPOINTS/checkpoint_01_qc_filtered.rds
#         RESULTS_QC/qc_summary.csv
#         FIGURES_SUPP/qc_violin_prefiltro.pdf, qc_violin_postfiltro.pdf
#         LOGS/01_qc_preprocessing_<timestamp>.log
#
# Verificado directamente sobre los datos crudos: matriz 25.127 genes
# x 82.168 celulas, 38 muestras (coincide con la lista GSM de GEO). No
# existe captura "Hepatocyte_34weeks_*"; las 6.048 celulas a 34 semanas
# provienen exclusivamente de NPC_34weeks_Animal1 (Capture1: 2.990,
# Capture2: 3.058). De esas, 441 estan anotadas CellType="Hepatocytes"
# por Su et al. (216 Capture1, 225 Capture2) - contaminacion
# documentada, se marca aqui y se excluye en Script 03.
#####################################################################

source(here::here("scripts", "R", "config.R"))

library(data.table)
library(Matrix)
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(LOGS, paste0("01_qc_preprocessing_", ts, ".log"))
sink(log_file, split = TRUE)
cat("Script 01 - QC preprocessing\n")
cat("Timestamp:", ts, "\n")
cat("Semilla global:", SEED_GLOBAL, "\n\n")

# --- 1. Carga de la matriz cruda ---
# 25.127 x 82.168 ~ 2.06 mil millones de celdas; como matriz densa
# double serian ~16.5 GB, riesgoso en 32GB RAM. Se lee con fread() y
# se convierte de inmediato a dgCMatrix (dispersa), liberando el
# intermedio denso con rm()+gc().
#
# Localizacion por patron, no por nombre exacto: los nombres de GEO
# han variado al pasar por Windows (guion bajo vs punto, mayusculas).

find_raw_file <- function(pattern, label) {
  candidates <- list.files(RAW_DATA, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(candidates) == 0) {
    stop("No se encontro ningun archivo de '", label, "' (patron '", pattern,
         "') en ", RAW_DATA, ".\nArchivos presentes en esa carpeta: ",
         paste(list.files(RAW_DATA), collapse = ", "))
  }
  if (length(candidates) > 1) {
    warning("Mas de un archivo coincide con '", label, "': ",
            paste(basename(candidates), collapse = ", "),
            " -- usando el primero: ", basename(candidates[1]))
  }
  cat("Archivo de", label, "localizado:", basename(candidates[1]), "\n")
  candidates[1]
}

raw_path <- find_raw_file("raw.*counts.*\\.gz$", "raw counts")

cat("Leyendo matriz cruda (archivo ~4GB descomprimido)...\n")
t0 <- Sys.time()

# El header se lee por separado: fread() con 82K columnas puede
# asignar como nombres de columna los VALORES de la primera fila de
# datos en vez de los barcodes reales del encabezado (falla observada
# en Windows). Leer el header aparte elimina esa ambiguedad.
con_header <- gzfile(raw_path, "r")
header_line <- readLines(con_header, n = 1)
close(con_header)
cell_barcodes <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
cat("Barcodes leidos directamente del header:", length(cell_barcodes), "\n")
stopifnot(length(cell_barcodes) == 82168)
stopifnot(!grepl("^[0-9.]+$", cell_barcodes[1]))

dt_raw <- fread(raw_path, header = FALSE, skip = 1, sep = "\t",
                 quote = "", showProgress = TRUE)
cat("Dimensiones leidas:", nrow(dt_raw), "genes x", ncol(dt_raw) - 1, "celulas\n")
stopifnot(ncol(dt_raw) - 1 == length(cell_barcodes))

# Tolerancia documentada: diferencia marginal (<0.01%) en conteo de
# genes entre plataformas se acepta si ningun gen MT falta (verificado
# abajo); una diferencia >10 genes detiene el pipeline.
n_genes_read <- nrow(dt_raw)
cat("Genes leidos:", n_genes_read, "(referencia esperada: 25127)\n")
if (abs(n_genes_read - 25127) > 10) {
  stop("Diferencia mayor a 10 genes respecto a lo esperado (", n_genes_read,
       " vs 25127) - investigar antes de continuar.")
}
cat("Tiempo de lectura:", round(difftime(Sys.time(), t0, units = "mins"), 1), "minutos\n")

gene_names <- dt_raw[[1]]

# Genes mitocondriales de este dataset en mayusculas (no el estandar
# Nd1/Cox1 de raton); comparacion case-insensitive por seguridad.
gene_names_upper <- toupper(gene_names)
mt_genes_upper <- toupper(MT_GENES_MOUSE)
mt_found_early <- MT_GENES_MOUSE[mt_genes_upper %in% gene_names_upper]
cat("Genes mitocondriales presentes:", length(mt_found_early), "de 13\n")
if (length(mt_found_early) < 13) {
  stop("Faltan genes mitocondriales criticos: ",
       paste(setdiff(MT_GENES_MOUSE, mt_found_early), collapse = ", "),
       " - esto SI afecta el calculo de %MT, resolver antes de continuar.")
} else {
  MT_GENES_ACTUAL_CASE <- gene_names[match(mt_genes_upper, gene_names_upper)]
  cat("Los 13 genes MT estan completos (case real: ",
      paste(head(MT_GENES_ACTUAL_CASE, 3), collapse = ", "), "...).\n\n")
}

dt_raw[, 1 := NULL]
gc()
mat_sparse <- Matrix(as.matrix(dt_raw), sparse = TRUE)
rownames(mat_sparse) <- gene_names
colnames(mat_sparse) <- cell_barcodes
rm(dt_raw); gc()

cat("Matriz dispersa construida:", nrow(mat_sparse), "x", ncol(mat_sparse), "\n")
cat("Tamano en memoria (dispersa):", round(object.size(mat_sparse) / 1e9, 2), "GB\n\n")

# --- 2. Metadata desde los barcodes ---
# Match contra la lista de las 38 muestras reales (mas robusto que
# asumir un sufijo de longitud fija de barcode).

KNOWN_FILENAMES <- c(
  "Hepatocyte_15weeks_Animal1_Capture1","Hepatocyte_15weeks_Animal1_Capture2",
  "Hepatocyte_15weeks_Animal2_Capture1","Hepatocyte_15weeks_Animal2_Capture2",
  "Hepatocyte_15weeks_Animal3_Capture1","Hepatocyte_15weeks_Animal3_Capture2",
  "Hepatocyte_30weeks_Animal1_Capture1","Hepatocyte_30weeks_Animal1_Capture2",
  "Hepatocyte_30weeks_Animal2_Capture1","Hepatocyte_30weeks_Animal2_Capture2",
  "Hepatocyte_Chow_Animal1_Capture1","Hepatocyte_Chow_Animal2_Capture1",
  "Hepatocyte_Chow_Animal3_Capture1","Hepatocyte_Chow_Animal4_Capture1",
  "Hepatocyte_Chow_Animal5_Capture1","Hepatocyte_Chow_Animal6_Capture1",
  "NPC_15weeks_Animal1_Capture1","NPC_15weeks_Animal1_Capture2",
  "NPC_15weeks_Animal2_Capture1","NPC_15weeks_Animal2_Capture2",
  "NPC_15weeks_Animal3_Capture1","NPC_15weeks_Animal3_Capture2",
  "NPC_30weeks_Animal1_Capture1","NPC_30weeks_Animal1_Capture2",
  "NPC_30weeks_Animal2_Capture1","NPC_30weeks_Animal2_Capture2",
  "NPC_30weeks_Animal3_Capture1","NPC_30weeks_Animal3_Capture2",
  "NPC_30weeks_Animal4_Capture1","NPC_30weeks_Animal4_Capture2",
  "NPC_34weeks_Animal1_Capture1","NPC_34weeks_Animal1_Capture2",
  "NPC_Chow_Animal1_Capture1","NPC_Chow_Animal1_Capture2",
  "NPC_Chow_Animal2_Capture1","NPC_Chow_Animal2_Capture2",
  "NPC_Chow_Animal3_Capture1","NPC_Chow_Animal3_Capture2"
)
stopifnot(length(KNOWN_FILENAMES) == 38)
# Orden por longitud descendente: evita que un nombre corto capture
# por error un prefijo de un nombre mas largo que lo contiene.
KNOWN_FILENAMES_SORTED <- KNOWN_FILENAMES[order(-nchar(KNOWN_FILENAMES))]

match_filename <- function(x) {
  hit <- KNOWN_FILENAMES_SORTED[startsWith(x, KNOWN_FILENAMES_SORTED)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

parse_meta <- data.frame(full_id = cell_barcodes, stringsAsFactors = FALSE) %>%
  mutate(
    file_name  = vapply(full_id, match_filename, character(1), USE.NAMES = FALSE),
    barcode_nt = ifelse(is.na(file_name), NA_character_,
                         substr(full_id, nchar(file_name) + 2, nchar(full_id))),
    fraction   = ifelse(grepl("^Hepatocyte", file_name), "Hepatocyte", "NPC"),
    week_raw   = case_when(
      grepl("Chow", file_name)    ~ "0",
      grepl("15weeks", file_name) ~ "15",
      grepl("30weeks", file_name) ~ "30",
      grepl("34weeks", file_name) ~ "34",
      TRUE ~ NA_character_
    ),
    week = as.numeric(week_raw),
    animal = str_extract(file_name, "Animal\\d+"),
    capture = str_extract(file_name, "Capture\\d+")
  )

n_na_week <- sum(is.na(parse_meta$week))
cat("Celulas sin match a ninguna de las 38 muestras conocidas:", n_na_week,
    sprintf("(%.4f%% del total)\n", 100 * n_na_week / nrow(parse_meta)))

if (n_na_week > 0) {
  cat("Ejemplos de full_id sin match:\n")
  print(head(parse_meta$full_id[is.na(parse_meta$week)], 10))
}

if (n_na_week / nrow(parse_meta) > 0.001) {
  stop("Mas del 0.1% de celulas sin match a muestra conocida - ",
       "investigar antes de continuar (posible corrupcion del archivo).")
} else if (n_na_week > 0) {
  warning(n_na_week, " celulas sin match a muestra conocida - se ",
          "excluiran mas adelante. Documentar en Limitaciones del TFM.")
}

cat("Distribucion de muestras (fraction x week):\n")
print(table(parse_meta$fraction, parse_meta$week))
cat("\n")

# --- 3. Cruce con la anotacion oficial de Su et al. ---

meta_path <- find_raw_file("metadata.*\\.gz$", "metadata")
official_meta <- fread(meta_path, header = TRUE, sep = "\t", quote = "")
colnames(official_meta) <- c("file_name_official", "cell_type_official", "cell_id")
official_meta$full_id <- paste0(official_meta$file_name_official, "_", official_meta$cell_id)

parse_meta <- parse_meta %>%
  left_join(official_meta %>% select(full_id, cell_type_official), by = "full_id")

cat("Celulas sin anotacion oficial de CellType (deberia ser 0):",
    sum(is.na(parse_meta$cell_type_official)), "\n\n")

write.csv(parse_meta, file.path(PROCESSED_META, "cell_metadata_cruzada.csv"),
          row.names = FALSE)

# --- 4. Contaminacion NPC_34weeks etiquetada "Hepatocytes" ---
# Confirmado: 441 celulas (216 Capture1 + 225 Capture2), todas de
# NPC_34weeks_Animal1.

parse_meta <- parse_meta %>%
  mutate(is_npc34w_hep_contamination =
           (!is.na(week) & fraction == "NPC" & week == 34 & cell_type_official == "Hepatocytes"))

n_contam <- sum(parse_meta$is_npc34w_hep_contamination)
cat("Celulas marcadas como contaminacion hepatocitaria NPC_34weeks:",
    n_contam, "(esperado: 441)\n")
if (n_contam != 441) {
  warning("El numero de celulas contaminantes no coincide con lo verificado ",
          "manualmente (441). Revisar el parseo antes de continuar.")
}
cat("\n")

rownames(parse_meta) <- parse_meta$full_id

# --- 5. Objeto Seurat ---

seurat_obj <- CreateSeuratObject(
  counts = mat_sparse,
  project = "GSE166504_LM",
  meta.data = parse_meta[colnames(mat_sparse), ]
)
rm(mat_sparse); gc()

# --- 6. Porcentaje mitocondrial ---

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj, features = intersect(MT_GENES_ACTUAL_CASE, rownames(seurat_obj))
)
genes_mt_encontrados <- intersect(MT_GENES_ACTUAL_CASE, rownames(seurat_obj))
cat("Genes mitocondriales encontrados:", length(genes_mt_encontrados), "de 13\n\n")

# --- 7. QC exploratorio pre-filtro ---

p_prefiltro <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        group.by = "fraction", pt.size = 0, ncol = 3)
ggsave(file.path(FIGURES_SUPP, "qc_violin_prefiltro.pdf"), p_prefiltro, width = 12, height = 5)

cat("Resumen %MT por condicion (pre-filtro):\n")
print(seurat_obj@meta.data %>% group_by(week) %>%
        summarise(pct_mt_mediana = median(percent.mt), n = n()))
cat("\n")

# --- 8. Filtrado de calidad (miQC preferido, fallback MAD-based) ---

USE_MIQC <- requireNamespace("miQC", quietly = TRUE) &&
            requireNamespace("flexmix", quietly = TRUE)

if (USE_MIQC) {
  cat("Usando miQC (Hippen, Falco et al. 2021, PLoS Comp Biol).\n")
  library(miQC); library(SingleCellExperiment); library(scater)

  sce <- as.SingleCellExperiment(seurat_obj)
  model <- mixtureModel(sce)
  saveRDS(model, file.path(CHECKPOINTS, "miqc_model.rds"))

  p_miqc <- plotModel(sce, model)
  ggsave(file.path(FIGURES_SUPP, "miqc_model.pdf"), p_miqc, width = 7, height = 5)

  sce_filtered <- filterCells(sce, model)
  keep_cells <- colnames(sce_filtered)
  cat("miQC retuvo", length(keep_cells), "de", ncol(seurat_obj), "celulas\n")

} else {
  cat("miQC no disponible (BiocManager::install(c('miQC','flexmix'))).\n")
  cat("Usando umbrales MAD-based como fallback:\n")

  med_feat <- median(seurat_obj$nFeature_RNA); mad_feat <- mad(seurat_obj$nFeature_RNA)
  min_features_threshold <- max(200, med_feat - 3 * mad_feat)

  med_mt <- median(seurat_obj$percent.mt); mad_mt <- mad(seurat_obj$percent.mt)
  max_mt_threshold <- min(30, med_mt + 3 * mad_mt)

  cat("Umbral nFeature_RNA (MAD-based):", round(min_features_threshold), "\n")
  cat("Umbral percent.mt (MAD-based, tope 30%):", round(max_mt_threshold, 1), "\n")

  keep_cells <- colnames(seurat_obj)[
    seurat_obj$nFeature_RNA > min_features_threshold &
    seurat_obj$percent.mt < max_mt_threshold
  ]
}

n_antes <- ncol(seurat_obj)
seurat_obj <- subset(seurat_obj, cells = keep_cells)
n_despues <- ncol(seurat_obj)

cat("\nCelulas antes:", n_antes, "| despues:", n_despues,
    sprintf("| eliminadas: %d (%.1f%%)\n", n_antes - n_despues,
            100 * (n_antes - n_despues) / n_antes))

n_contam_post <- sum(seurat_obj$is_npc34w_hep_contamination)
cat("De las 441 celulas contaminantes NPC_34weeks,", n_contam_post,
    "sobrevivieron el QC (se excluyen en Script 03 por origen de muestra).\n\n")

# --- 9. QC post-filtro y tabla resumen ---

p_postfiltro <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                         group.by = "fraction", pt.size = 0, ncol = 3)
ggsave(file.path(FIGURES_SUPP, "qc_violin_postfiltro.pdf"), p_postfiltro, width = 12, height = 5)

qc_summary <- seurat_obj@meta.data %>%
  group_by(week, fraction) %>%
  summarise(n_celulas = n(), nFeature_mediana = median(nFeature_RNA),
            nCount_mediana = median(nCount_RNA), pct_mt_mediana = median(percent.mt),
            .groups = "drop")
write.csv(qc_summary, file.path(RESULTS_QC, "qc_summary.csv"), row.names = FALSE)
print(qc_summary)

# --- 10. Checkpoint ---

saveRDS(seurat_obj, file.path(CHECKPOINTS, "checkpoint_01_qc_filtered.rds"))
cat("\nCheckpoint guardado.\n\n---- sessionInfo() ----\n")
print(sessionInfo())
sink()
cat("Script 01 completo. Log:", log_file, "\n")
