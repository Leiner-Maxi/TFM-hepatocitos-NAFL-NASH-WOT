#####################################################################
# config.py - fuente unica de verdad: rutas, semillas, parametros (WOT)
#####################################################################
import os
import numpy as np

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PROCESSED    = os.path.join(PROJECT_ROOT, "data", "processed")
RESULTS_WOT  = os.path.join(PROJECT_ROOT, "results", "wot")
RESULTS_IRREV = os.path.join(PROJECT_ROOT, "results", "irreversibility")
CHECKPOINTS  = os.path.join(PROJECT_ROOT, "checkpoints")

SEED_GLOBAL = 42
np.random.seed(SEED_GLOBAL)

WOT_EPSILON = 0.05
WOT_LAMBDA1 = 1
WOT_LAMBDA2 = 50
EPSILON_GRID = [0.03, 0.05, 0.10]
LAMBDA2_GRID = [30, 50, 80]

HEPATOCYTE_TIMEPOINTS_REAL = [0, 15, 30]

PROLIFERATION_MARKERS = ["Mki67", "Pcna", "Top2a", "Mcm2", "Ccnb1", "Cdk1"]
APOPTOSIS_MARKERS      = ["Casp3", "Casp8", "Casp9", "Bax", "Bak1", "Bad"]

BETA_MAX = 1.7
BETA_MIN = 0.3
DELTA_MAX = 1.7
DELTA_MIN = 0.3
GROWTH_ITERS = 3

N_BOOTSTRAP    = 1000
N_PERMUTATIONS = 1000

# --- Genes driver (Script 07) ---
# Con n~4481 celulas, correlaciones de Spearman muy debiles ya cruzan
# significancia estadistica (error estandar ~1/sqrt(n) ~ 0.015). El
# umbral de tamano de efecto separa significancia estadistica de
# relevancia biologica - un gen debe cumplir AMBOS criterios para
# reportarse como "driver" en Resultados, no solo FDR<0.05.
DRIVER_GENE_EFFECT_SIZE_THRESHOLD = 0.15
