# Irreversibilidad transcriptómica en hepatocitos durante la progresión NAFL-NASH

Análisis de trayectorias celulares mediante Waddington Optimal Transport (WOT)
sobre datos de secuenciación de célula única (GSE166504, Su et al. 2021),
cuantificando el compromiso de destino transcripcional de hepatocitos durante
la progresión de esteatosis hepática no alcohólica (NAFL) a esteatohepatitis
(NASH).

## Resultados principales

- **65.158 células** retenidas tras control de calidad (de 82.168 iniciales).
- **7.682 hepatocitos** identificados mediante criterio doble de selección,
  reagrupados en **11 estados transcripcionales**, tras verificar y excluir
  una fracción de contaminación no parenquimal a 34 semanas (441 células,
  confirmado contra la anotación oficial de los autores originales).
- Entropía normalizada de destino (compromiso de destino transcripcional):
  **0,767 → 0,627** entre los intervalos control-NAFLD y NAFLD-NASH
  (bootstrap IC95% sin solapamiento; test de permutación p<0,001).
- Validado de forma independiente a nivel de réplica biológica (no solo
  celular): separación completa entre 6 y 3 animales (Mann-Whitney U=18,
  p=0,024), descartando pseudorreplicación.
- **294 genes** asociados significativamente al compromiso de destino
  (FDR<0,05 y |ρ|>0,15), con *Gnmt* como asociación más fuerte (ρ=0,410),
  y los genes de riesgo genético humano *Pnpla3* e *Hsd17b13* también
  significativos.
- Los **20 genes principales** (10 por dirección, los que sustentan la
  Tabla 4 y la discusión biológica) se validaron mediante leave-one-animal-out:
  **20 de 20 mantienen signo y significancia** al excluir cualquiera de los
  3 animales del intervalo 15→30 semanas.
- Las matrices de transporte agregadas por estado (masa y probabilidad
  normalizada) se visualizan directamente como evidencia complementaria al
  resultado de entropía.

## Estructura del repositorio

```
scripts/
  R/            Scripts 01-04: control de calidad, clustering, extracción
                hepatocitaria, exportación de datos para WOT
  python/       Scripts 05-09: transporte óptimo, métricas de
                irreversibilidad, validación de robustez, visualización
                de matrices de transporte, genes asociados, figuras finales
run_analysis.R  Script maestro que corre 01-04 en secuencia
data/
  raw/          Datos crudos de GSE166504 (no incluidos, ver abajo)
  processed/    Archivos de entrada para WOT (generados por Script 04)
results/        Figuras, tablas, matrices de transporte y resultados numéricos
```

### Scripts completos (15 archivos)

| # | Archivo | Función |
|---|---|---|
| — | `scripts/R/config.R` | Rutas, semillas, parámetros (R) |
| 01 | `scripts/R/01_qc_preprocessing.R` | Control de calidad y filtrado |
| 02 | `scripts/R/02_normalization_clustering.R` | Normalización, PCA, UMAP, clustering global |
| 03 | `scripts/R/03_annotation_hepatocyte_subset.R` | Extracción hepatocitaria, zonación |
| 04 | `scripts/R/04_export_wot_inputs.R` | Exportación de inputs para WOT |
| — | `scripts/python/config.py` | Rutas, semillas, parámetros (Python) |
| 05 | `scripts/python/05_wot_transport_maps.py` | Cálculo de mapas de transporte (WOT) |
| 06 | `scripts/python/06_irreversibility_metrics.py` | Entropía normalizada, bootstrap, permutación |
| 06b | `scripts/python/06b_robustness_checks.py` | Validación a nivel de réplica biológica |
| 06c | `scripts/python/06c_transport_maps.py` | Matrices de transporte agregadas por estado |
| 07 | `scripts/python/07_driver_gene_analysis.py` | Identificación de genes asociados |
| 07b | `scripts/python/07b_gene_robustness_loo.py` | Validación leave-one-animal-out de genes |
| 08 | `scripts/python/08_generate_figures.py` | Figuras finales (entropía, animal, genes, paneles UMAP) |
| 09 | `scripts/python/09_pdf_to_png_conversion.py` | Conversión de figuras a PNG para manuscrito |
| — | `run_analysis.R` | Corre 01-04 en secuencia, sesión limpia |

## Reproducir el análisis

### 1. Datos

Este repositorio no incluye los datos crudos (~4GB descomprimidos). Descargar
manualmente de GEO ([GSE166504](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166504))
y colocar en `data/raw/`:

- `GSE166504_cell_raw_counts.20220204.txt.gz`
- `GSE166504_cell_metadata.20220204.tsv.gz`

### 2. Entorno R

```r
# R 4.5.2. Paquetes principales: Seurat 5.5.0, data.table 1.18.2,
# dplyr 1.2.1, ggplot2, stringr, here, patchwork
source("scripts/R/config.R")   # crea la estructura de carpetas automáticamente
source("scripts/R/01_qc_preprocessing.R")
source("scripts/R/02_normalization_clustering.R")
source("scripts/R/03_annotation_hepatocyte_subset.R")
source("scripts/R/04_export_wot_inputs.R")
```

O de una sola vez, en sesión limpia (`Session > Restart R` antes):

```r
source("run_analysis.R")
```

### 3. Entorno Python (WOT)

```bash
conda create -n wot_env python=3.8
conda activate wot_env
pip install wot==1.0.8 anndata pymupdf

python scripts/python/05_wot_transport_maps.py
python scripts/python/06_irreversibility_metrics.py
python scripts/python/06b_robustness_checks.py
python scripts/python/06c_transport_maps.py
python scripts/python/07_driver_gene_analysis.py
python scripts/python/07b_gene_robustness_loo.py
python scripts/python/08_generate_figures.py
python scripts/python/09_pdf_to_png_conversion.py
```

## Reproducibilidad

- Semilla fija (42) en todo paso estocástico, centralizada en `config.R` /
  `config.py`.
- Cada script genera un log con timestamp y `sessionInfo()` (R) o
  parámetros exactos (Python).
- **Nota sobre control de calidad**: el filtrado del Script 01 usa miQC
  (Hippen et al. 2021) si está disponible en el entorno, o un umbral
  MAD-based como fallback en caso contrario. Los resultados reportados
  en este trabajo (65.158 células retenidas) corresponden al fallback
  MAD-based, dado que miQC no estaba instalado en el entorno de
  desarrollo. Si se corre con miQC instalado, el número de células
  retenidas puede diferir.
- Verificado de extremo a extremo en sesión limpia (`Session > Restart R`
  seguido de los 15 scripts en orden), reproduciendo cada resultado
  numérico exacto.

## Metodología

Ver el documento completo del TFM (`manuscript/`) para la justificación
detallada de cada decisión metodológica, resultados completos, discusión
y limitaciones.

## Datos y código disponibles

- Datos: [GSE166504](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166504) (Su et al., 2021, *iScience*)
- WOT: [broadinstitute/wot](https://github.com/broadinstitute/wot) (Schiebinger et al., 2019, *Cell*)

## Declaración de uso de IA

Se utilizó Claude (Anthropic) como herramienta de apoyo extensivo en el
desarrollo de este pipeline, incluyendo generación asistida de código,
depuración activa de errores técnicos y correcciones metodológicas
identificadas durante el desarrollo. Ver Anexo I del TFM para la
declaración completa.
