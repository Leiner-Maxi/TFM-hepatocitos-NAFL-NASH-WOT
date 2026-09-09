"""
09_pdf_to_png_conversion.py
Convierte a PNG las figuras PDF generadas por los Scripts 01-03 (R) y
08 (Python), para insertarlas en el manuscrito o presentaciones que no
manejan bien archivos vectoriales (Word, PowerPoint). Estas son las 13
figuras completas del TFM/paper (11 individuales + 2 paneles
combinados generados por Script 10).

Los PDF originales (en results/figures/main y .../supplementary) siguen
siendo el archivo maestro, vectorial, apto para publicacion. Los PNG
generados aqui son una copia derivada para manuscrito/presentacion, no
una fuente de verdad nueva.

Requiere: pip install PyMuPDF

Input:  results/figures/supplementary/{qc_violin_prefiltro,
        qc_violin_postfiltro}.pdf
        results/figures/main/{umap_por_cluster,umap_por_condicion,
        umap_por_fraccion,umap_hepatocitos_por_estado,
        umap_hepatocitos_por_condicion,zonacion_markers,
        fig10_entropia_bootstrap,fig11_validacion_animal,
        fig12_driver_genes}.pdf
Output: results/figures/png_para_manuscrito/*.png (200 dpi)
"""
import os
import sys
import pymupdf

sys.path.insert(0, os.path.dirname(__file__))
from config import PROJECT_ROOT

DPI = 200
FIGURES_SUPP = os.path.join(PROJECT_ROOT, "results", "figures", "supplementary")
FIGURES_MAIN = os.path.join(PROJECT_ROOT, "results", "figures", "main")
OUT_DIR = os.path.join(PROJECT_ROOT, "results", "figures", "png_para_manuscrito")
os.makedirs(OUT_DIR, exist_ok=True)

# (nombre base sin extension, carpeta de origen)
ARCHIVOS = [
    # Scripts 01-03 (R) - control de calidad y clustering
    ("qc_violin_prefiltro", FIGURES_SUPP),
    ("qc_violin_postfiltro", FIGURES_SUPP),
    ("umap_por_cluster", FIGURES_MAIN),
    ("umap_por_condicion", FIGURES_MAIN),
    ("umap_por_fraccion", FIGURES_MAIN),
    ("umap_hepatocitos_por_estado", FIGURES_MAIN),
    ("umap_hepatocitos_por_condicion", FIGURES_MAIN),
    ("zonacion_markers", FIGURES_MAIN),
    # Script 08 (Python) - resultados de irreversibilidad
    ("fig10_entropia_bootstrap", FIGURES_MAIN),
    ("fig11_validacion_animal", FIGURES_MAIN),
    ("fig12_driver_genes", FIGURES_MAIN),
    # Script 10 (Python) - paneles combinados
    ("panel_umap_global", FIGURES_MAIN),
    ("panel_umap_hepatocitos", FIGURES_MAIN),
]

# PyMuPDF renderiza a 72 dpi por defecto (1 punto PDF = 1 pixel);
# el factor de zoom escala eso al DPI deseado.
zoom = DPI / 72
matrix = pymupdf.Matrix(zoom, zoom)

convertidos = 0
for nombre, carpeta in ARCHIVOS:
    pdf_path = os.path.join(carpeta, f"{nombre}.pdf")
    if not os.path.exists(pdf_path):
        print(f"AVISO: no existe {pdf_path} (correr Scripts 01-03 primero) - se omite.")
        continue

    doc = pymupdf.open(pdf_path)
    pix = doc[0].get_pixmap(matrix=matrix)
    out_path = os.path.join(OUT_DIR, f"{nombre}.png")
    pix.save(out_path)
    doc.close()

    print(f"{nombre}.pdf -> {nombre}.png ({pix.width}x{pix.height}px, {DPI} dpi)")
    convertidos += 1

print(f"\n{convertidos} de {len(ARCHIVOS)} figuras convertidas.")
print(f"Salida: {OUT_DIR}")
