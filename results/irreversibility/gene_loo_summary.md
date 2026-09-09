# Validación leave-one-animal-out de genes asociados al compromiso de destino

## Metodología

Se evaluaron los 20 genes que sustentan la Figura 12, la Tabla 4 y la discusión biológica del apartado 5.5 (los 10 de mayor |ρ| en cada dirección, seleccionados dinámicamente desde `driver_genes_alta_confianza.csv`). Para cada gen se recalculó la correlación de Spearman entre expresión génica y entropía normalizada de destino excluyendo un animal a la vez (Animal1, Animal2, Animal3, intervalo 15→30 semanas), con corrección FDR (Benjamini-Hochberg) recalculada dentro de cada subconjunto de 20 genes.

## Criterios de clasificación

- **Robusto**: signo idéntico al original en las 3 exclusiones y significativo (FDR<0.05) en las 3.
- **Parcialmente robusto**: signo idéntico en las 3 exclusiones, pierde significancia (FDR≥0.05) en al menos una.
- **No robusto**: el signo se invierte en al menos una exclusión.

No se introdujo ningún umbral numérico nuevo: el criterio de significancia (FDR<0.05) es el mismo ya establecido y utilizado en el Script 07 para la identificación original de estos genes.

## Resultados

- Genes evaluados: **20**
- Robustos: **20** de 20 (100.0%)
- Parcialmente robustos: **0** de 20 (0.0%)
- No robustos: **0** de 20 (0.0%)

**Genes robustos**: Hspa5, Ces3a, Slco1b2, Neat1, Cp, Malat1, Pzp, Mug1, H2-K1, Fgg, Gnmt, Rps2, Wfdc21, Gstp2, Gstp1, Fabp2, Car3, Acaa1b, Serpina12, Cth.

**Genes parcialmente robustos**: ninguno.

**Genes no robustos**: ninguno.

## Interpretación metodológica

Un gen clasificado como robusto muestra el mismo patrón de asociación con el compromiso de destino independientemente de qué animal se excluya del análisis, lo que reduce la probabilidad de que el hallazgo sea un artefacto de un único individuo, en línea con la preocupación de pseudorreplicación documentada por Squair et al. (2021) y ya aplicada a la métrica de entropía en el apartado 5.4. Un gen no robusto no debe descartarse como falso positivo sin más evidencia, pero su asociación debe interpretarse con mayor cautela que la de un gen robusto.
