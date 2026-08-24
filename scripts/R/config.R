#####################################################################
# config.R - fuente unica de rutas, semillas y parametros.
# Cargar al inicio de cada script:
#   source(here::here("scripts", "R", "config.R"))
#####################################################################
library(here)

# --- Rutas ---
RAW_DATA        <- here("data", "raw")
PROCESSED       <- here("data", "processed")
PROCESSED_META  <- here("data", "processed", "metadata")
EXTERNAL        <- here("data", "external")
GENE_SETS       <- here("data", "external", "gene_sets")
CHECKPOINTS     <- here("checkpoints")
LOGS            <- here("logs")
RESULTS         <- here("results")
RESULTS_QC      <- here("results", "qc")
RESULTS_CLUST   <- here("results", "clustering")
RESULTS_HEP     <- here("results", "hepatocytes")
RESULTS_WOT     <- here("results", "wot")
RESULTS_IRREV   <- here("results", "irreversibility")
FIGURES_MAIN    <- here("results", "figures", "main")
FIGURES_SUPP    <- here("results", "figures", "supplementary")
TABLES          <- here("results", "tables")
REPORTS         <- here("reports")

# --- Semillas ---
SEED_GLOBAL     <- 42
SEED_PCA        <- SEED_GLOBAL
SEED_UMAP       <- SEED_GLOBAL
SEED_CLUSTERING <- SEED_GLOBAL
SEED_BOOTSTRAP  <- SEED_GLOBAL
set.seed(SEED_GLOBAL)

# --- Script 01 (QC) ---
QC_MIN_FEATURES <- 200
QC_MAX_MT_PCT   <- 30
# GSE166504 anota los genes mitocondriales en mayusculas (COX1, ND1,
# CYTB...) en vez del prefijo estandar de raton/Ensembl (^mt-).
MT_GENES_MOUSE <- c("ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
                     "COX1","COX2","COX3","ATP6","ATP8","CYTB")

# --- Script 02 (clustering global) ---
N_HVG                      <- 2000
N_PCS_FULL                 <- 20
CLUSTERING_RESOLUTION_FULL <- 0.5
UMAP_N_NEIGHBORS           <- 30
UMAP_MIN_DIST              <- 0.3

# --- Script 03 (subset hepatocitario) ---
N_PCS_HEPATOCYTES         <- 15
CLUSTERING_RESOLUTION_HEP <- 0.3
# 441 celulas "Hepatocytes" dentro de NPC_34weeks_Animal1 son
# contaminacion documentada por los propios autores; excluir sin
# excepcion (ver docs/methodology_decisions.md).
EXCLUDE_NPC_34W <- TRUE
HEPATOCYTE_TIMEPOINTS_REAL <- c(0, 15, 30)

# --- WOT (uso real en scripts/python/config.py; aqui solo para
#     trazabilidad de un unico valor de referencia) ---
WOT_EPSILON <- 0.05
WOT_LAMBDA1 <- 1
WOT_LAMBDA2 <- 50
# Defaults de Schiebinger et al. 2019, robustos en rango amplio
# (epsilon 5e-5 a 0.5, lambda 0.1 a 32).

# --- Creacion de carpetas ---
# Ausente en la version anterior de este archivo; sin esto, cualquier
# script falla en su primer intento de escritura si las carpetas no
# existen ya en disco.
for (d in c(RAW_DATA, PROCESSED, PROCESSED_META, EXTERNAL, GENE_SETS,
            CHECKPOINTS, LOGS, RESULTS_QC, RESULTS_CLUST, RESULTS_HEP,
            RESULTS_WOT, RESULTS_IRREV, FIGURES_MAIN, FIGURES_SUPP,
            TABLES, REPORTS)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
