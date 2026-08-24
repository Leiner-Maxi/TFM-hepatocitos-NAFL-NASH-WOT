#####################################################################
# run_analysis.R
# Script maestro: corre el pipeline de R completo en una sesion limpia,
# con log de sessionInfo() por paso (Session > Restart R antes de correr).
#####################################################################

source(here::here("scripts", "R", "config.R"))

log_step <- function(step_name, expr) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  message(sprintf("[%s] Iniciando: %s", ts, step_name))
  result <- eval(expr, envir = parent.frame())
  log_file <- file.path(LOGS, sprintf("%s_%s.log", ts, step_name))
  writeLines(c(sprintf("Paso: %s", step_name), sprintf("Timestamp: %s", ts),
               "---- sessionInfo() ----", capture.output(sessionInfo())), log_file)
  message(sprintf("  -> log guardado en %s", log_file))
  result
}

log_step("01_qc_preprocessing",             quote(source(here::here("scripts","R","01_qc_preprocessing.R"))))
log_step("02_normalization_clustering",     quote(source(here::here("scripts","R","02_normalization_clustering.R"))))
log_step("03_annotation_hepatocyte_subset", quote(source(here::here("scripts","R","03_annotation_hepatocyte_subset.R"))))
log_step("04_export_wot_inputs",            quote(source(here::here("scripts","R","04_export_wot_inputs.R"))))

message("\nLado R completo. Continuar con Python (conda activate wot_env):")
message("  scripts/python/05_wot_transport_maps.py")
message("  scripts/python/06_irreversibility_metrics.py")
message("  scripts/python/06b_robustness_checks.py")
message("  scripts/python/06c_transport_maps.py")
message("  scripts/python/07_driver_gene_analysis.py")
message("  scripts/python/07b_gene_robustness_loo.py")
message("  scripts/python/08_generate_figures.py")
message("  scripts/python/09_pdf_to_png_conversion.py")
