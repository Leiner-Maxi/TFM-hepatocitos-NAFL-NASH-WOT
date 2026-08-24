"""
05_wot_transport_maps.py

Calculo de matrices de transporte optimo (WOT), via interfaz de linea
de comandos (mas confiable que la API interna de Python para wot 1.0.8,
cuya documentacion publica es limitada).

Input:  data/processed/{expression_matrix,cell_days,growth_rates}.txt
Output: results/wot/tmaps/hepatocytes_*.h5ad (0->15, 15->30 semanas)
Entorno: conda activate wot_env (Python 3.8, WOT 1.0.8, POT 0.9.3)
"""
import subprocess
import sys
import os
import datetime

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from config import (WOT_EPSILON, WOT_LAMBDA1, WOT_LAMBDA2, GROWTH_ITERS,
                     PROCESSED, RESULTS_WOT, PROJECT_ROOT)

LOGS = os.path.join(PROJECT_ROOT, "logs")
os.makedirs(LOGS, exist_ok=True)
os.makedirs(RESULTS_WOT, exist_ok=True)
os.makedirs(os.path.join(RESULTS_WOT, "tmaps"), exist_ok=True)

ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = os.path.join(LOGS, f"05_wot_transport_maps_{ts}.log")


def log(msg, f):
    print(msg)
    f.write(msg + "\n")
    f.flush()


with open(log_file, "w") as f:
    log(f"Script 05 - WOT transport maps", f)
    log(f"Timestamp: {ts}", f)
    log(f"Parametros: epsilon={WOT_EPSILON}, lambda1={WOT_LAMBDA1}, "
        f"lambda2={WOT_LAMBDA2}, growth_iters={GROWTH_ITERS}\n", f)

    # Verificacion previa de flags disponibles en la instalacion real,
    # antes de comprometerse a la corrida completa.
    log("Verificando instalacion de wot...", f)
    try:
        help_result = subprocess.run(
            ["wot", "optimal_transport", "-h"],
            capture_output=True, text=True, timeout=180
        )
        log(help_result.stdout, f)
        if help_result.returncode != 0:
            log("ADVERTENCIA: 'wot optimal_transport -h' devolvio error:", f)
            log(help_result.stderr, f)
    except subprocess.TimeoutExpired:
        log("ERROR: 'wot optimal_transport -h' no respondio en 180 segundos.", f)
        log("Verificar manualmente: where wot / wot optimal_transport -h", f)
        sys.exit(1)
    except FileNotFoundError:
        log("ERROR CRITICO: comando 'wot' no encontrado. Verificar que", f)
        log("el entorno conda wot_env este activado y wot instalado.", f)
        sys.exit(1)

    matrix_path = os.path.join(PROCESSED, "expression_matrix.txt")
    days_path = os.path.join(PROCESSED, "cell_days.txt")
    growth_path = os.path.join(PROCESSED, "growth_rates.txt")
    out_prefix = os.path.join(RESULTS_WOT, "tmaps", "hepatocytes")

    for p in [matrix_path, days_path, growth_path]:
        if not os.path.exists(p):
            log(f"ERROR CRITICO: no existe el archivo {p}", f)
            log("Verificar que Script 04 (R) se corrio completo.", f)
            sys.exit(1)
    log(f"\nArchivos de entrada verificados:\n  {matrix_path}\n  "
        f"{days_path}\n  {growth_path}\n", f)

    # Con 3 timepoints reales (0, 15, 30 semanas): 2 transport maps.
    cmd = [
        "wot", "optimal_transport",
        "--matrix", matrix_path,
        "--cell_days", days_path,
        "--cell_growth_rates", growth_path,
        "--growth_iters", str(GROWTH_ITERS),
        "--epsilon", str(WOT_EPSILON),
        "--lambda1", str(WOT_LAMBDA1),
        "--lambda2", str(WOT_LAMBDA2),
        "--out", out_prefix,
        "--verbose"
    ]

    log("Comando: " + " ".join(cmd), f)
    log("\nPuede tardar minutos a horas (Sinkhorn OT sin GPU). La salida", f)
    log("se muestra en vivo linea por linea a medida que wot progresa.\n", f)

    t0 = datetime.datetime.now()
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1, universal_newlines=True
    )
    for line in process.stdout:
        log(line.rstrip(), f)
    process.wait()
    elapsed = (datetime.datetime.now() - t0).total_seconds() / 60
    returncode = process.returncode

    log(f"\nTiempo de ejecucion: {elapsed:.1f} minutos", f)
    log(f"Codigo de salida: {returncode}", f)

    if returncode != 0:
        log("\nERROR: el comando wot optimal_transport termino con error.", f)
        sys.exit(1)

    log("\nVerificando archivos de salida...", f)
    tmap_dir = os.path.join(RESULTS_WOT, "tmaps")
    output_files = [x for x in os.listdir(tmap_dir) if x.startswith("hepatocytes")]
    log(f"Archivos generados en {tmap_dir}:", f)
    for of in output_files:
        log(f"  {of}", f)

    n_expected_maps = len([0, 15, 30]) - 1
    log(f"\nTransport maps esperados: {n_expected_maps} (0->15, 15->30)", f)
    log(f"Archivos .h5ad encontrados: "
        f"{len([x for x in output_files if x.endswith('.h5ad')])}", f)

    log("\nScript 05 completo.", f)

print(f"\nLog completo guardado en: {log_file}")
