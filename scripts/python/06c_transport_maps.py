"""
06c_transport_maps.py

Aggregates cell-level WOT transport maps into state-to-state summary
matrices and renders them as a single multi-panel figure.

Methodology: state-to-state aggregation by summation, mathematically
equivalent to `wot transition_table`; the row-normalized matrix is
equivalent to the WOT fate matrix. Both follow Schiebinger et al.,
"Optimal-Transport Analysis of Single-Cell Gene Expression Identifies
Developmental Trajectories in Reprogramming," Cell 176(4):928-943,
2019, and the official documentation at
https://broadinstitute.github.io/wot/tutorial/ (Notebooks 5-6).

These matrices and the entropy metric in Script 06 are computed from
the same transport-map objects (Script 05 output). Script 06 performs
the quantitative analysis, using Shannon entropy of the destination
distribution. This script only provides a visual, interpretable
representation of that same mathematical object; the patterns shown
here are compatible with the entropy result but do not themselves
constitute quantitative evidence.

Read-only with respect to Scripts 05/06/06b/07: consumes their output
files and never modifies them. Panel dimensions and color scales are
computed from the input data, not assumed in advance.

Input
-----
results/wot/tmaps/hepatocytes_0.0_15.0.h5ad
results/wot/tmaps/hepatocytes_15.0_30.0.h5ad
data/processed/cell_metadata.csv

Output
------
results/wot_visualization/transport_states_0_15_mass.csv
results/wot_visualization/transport_states_15_30_mass.csv
results/wot_visualization/transport_states_15_30_normalized.csv
results/wot_visualization/heatmap_transport_0_15.pdf
results/wot_visualization/heatmap_transport_0_15.png
results/wot_visualization/heatmap_transport_15_30.pdf
results/wot_visualization/heatmap_transport_15_30.png
results/wot_visualization/heatmap_transport_15_30_normalized.pdf
results/wot_visualization/heatmap_transport_15_30_normalized.png
results/wot_visualization/transport_summary.pdf
results/wot_visualization/transport_summary.png
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import NamedTuple

import anndata
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

sys.path.insert(0, str(Path(__file__).parent))
from config import PROCESSED, RESULTS_WOT, PROJECT_ROOT

STATE_COLUMN = "seurat_clusters"
OUT_DIR = Path(PROJECT_ROOT) / "results" / "wot_visualization"
COLORMAP = "viridis"
FIGURE_DPI = 600
CELL_SIZE_INCHES = 0.42
MIN_PANEL_INCHES = 2.3
MAX_PANEL_INCHES = 5.5
PANEL_LABEL_STYLE = dict(fontsize=13, fontweight="bold", ha="left", va="bottom")


class StateAggregation(NamedTuple):
    mass: pd.DataFrame
    n_origin_cells: pd.Series


def read_transport_map(h5ad_path: Path) -> anndata.AnnData:
    """Load a WOT transport map, failing clearly if it is absent."""
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Missing transport map: {h5ad_path}. Run Script 05 first.")
    return anndata.read_h5ad(h5ad_path)


def assign_states(cell_ids: pd.Index, metadata: pd.DataFrame) -> pd.Series:
    """Map transport-map cell identifiers to their transcriptional state.

    Parameters
    ----------
    cell_ids : pd.Index
        Cell identifiers from a transport map's obs_names or var_names.
    metadata : pd.DataFrame
        cell_metadata.csv contents, with an "id" column and STATE_COLUMN.

    Returns
    -------
    pd.Series
        State label per cell, aligned with cell_ids.
    """
    indexed = metadata.set_index("id")
    missing = cell_ids.difference(indexed.index)
    if len(missing) > 0:
        raise ValueError(f"{len(missing)} cells absent from cell_metadata.csv (e.g. {missing[0]})")
    return indexed.loc[cell_ids, STATE_COLUMN].astype(str)


def aggregate_transport_by_state(
    transport: np.ndarray,
    origin_states: pd.Series,
    dest_states: pd.Series,
) -> StateAggregation:
    """Aggregate a cell-level transport matrix into state-to-state mass.

    Implements sum_{x in A_i} sum_{y in B_j} transport(x, y), the
    formulation used by `wot transition_table` (Schiebinger et al.,
    2019; WOT documentation, Notebook 6). States with zero origin
    cells are excluded from the result rather than filled with zeros,
    since absent data and zero transported mass are not equivalent.

    Parameters
    ----------
    transport : np.ndarray
        Dense cell-by-cell transport matrix (rows=origin, cols=dest).
    origin_states : pd.Series
        State label per origin cell, aligned with transport's rows.
    dest_states : pd.Series
        State label per destination cell, aligned with transport's
        columns.

    Returns
    -------
    StateAggregation
        Aggregated mass matrix and origin cell counts per state.
    """
    origin_labels = sorted(origin_states.unique(), key=int)
    dest_labels = sorted(dest_states.unique(), key=int)

    origin_indicator = pd.get_dummies(origin_states)[origin_labels].to_numpy(dtype=float)
    dest_indicator = pd.get_dummies(dest_states)[dest_labels].to_numpy(dtype=float)

    mass = origin_indicator.T @ transport @ dest_indicator
    mass_df = pd.DataFrame(mass, index=origin_labels, columns=dest_labels)
    mass_df.index.name, mass_df.columns.name = "origin_state", "destination_state"

    n_origin_cells = origin_states.value_counts().reindex(origin_labels)
    return StateAggregation(mass_df, n_origin_cells)


def normalize_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    """Row-normalize a mass matrix into a fate probability matrix.

    Equivalent to the WOT fate matrix (Schiebinger et al., 2019;
    Notebook 5): each row becomes a probability distribution over
    destination states, P(destination | origin).
    """
    return matrix.div(matrix.sum(axis=1), axis=0)


def _panel_extent(matrix: pd.DataFrame) -> tuple[float, float]:
    """Axes width/height in inches so every panel's cells share a
    consistent physical size, regardless of how many states a given
    dataset happens to produce. Bounded to MIN/MAX_PANEL_INCHES so an
    unusually small or large number of states cannot degenerate the
    layout, the same way a published figure keeps a sensible panel
    size regardless of table length."""
    width = float(np.clip(matrix.shape[1] * CELL_SIZE_INCHES, MIN_PANEL_INCHES, MAX_PANEL_INCHES))
    height = float(np.clip(matrix.shape[0] * CELL_SIZE_INCHES, MIN_PANEL_INCHES, MAX_PANEL_INCHES))
    return width, height


def _draw_panel(
    ax: plt.Axes,
    matrix: pd.DataFrame,
    n_origin_cells: pd.Series,
    title: str,
    norm: Normalize,
    panel_label: str = "",
) -> plt.cm.ScalarMappable:
    mesh = ax.pcolormesh(
        matrix.to_numpy(), cmap=COLORMAP, norm=norm, edgecolors="white", linewidth=0.6
    )
    ax.set_xticks(np.arange(matrix.shape[1]) + 0.5)
    ax.set_xticklabels(matrix.columns, fontsize=8)
    ax.set_yticks(np.arange(matrix.shape[0]) + 0.5)
    ax.set_yticklabels([f"{s} (n={n_origin_cells[s]})" for s in matrix.index], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Destination state", fontsize=9)
    ax.set_ylabel("Origin state", fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.set_aspect("auto")
    if panel_label:
        ax.text(-0.05, 1.12, panel_label, transform=ax.transAxes, **PANEL_LABEL_STYLE)
    return mesh


def plot_transport_summary(
    mass_0_15: StateAggregation,
    mass_15_30: StateAggregation,
    normalized_15_30: pd.DataFrame,
    mass_norm: Normalize,
    prob_norm: Normalize,
    out_path_stem: Path,
) -> None:
    """Render the A/B/C transport summary figure (A and B on top, C below).

    Panels A and B share a colorbar and value range (raw transported
    mass), making their magnitudes directly comparable. Panel C uses
    an independent [0, 1] scale, since it shows row-normalized fate
    probabilities rather than mass. Panel size follows the number of
    states in each matrix, bounded to a fixed range so the layout
    stays stable regardless of how many states a given dataset
    produces.
    """
    width_a, height_a = _panel_extent(mass_0_15.mass)
    width_b, height_b = _panel_extent(mass_15_30.mass)
    _, height_c = _panel_extent(normalized_15_30)
    top_row_height = max(height_a, height_b)
    cbar_width = 0.35

    fig_width = width_a + width_b + cbar_width + 1.6
    fig_height = top_row_height + height_c + 1.6
    fig = plt.figure(figsize=(fig_width, fig_height))
    grid = fig.add_gridspec(
        2, 3,
        width_ratios=[width_a, width_b, cbar_width],
        height_ratios=[top_row_height, height_c],
        hspace=0.6, wspace=0.35,
        left=0.08, right=0.90, top=0.90, bottom=0.08,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    cax_mass = fig.add_subplot(grid[0, 2])
    ax_c = fig.add_subplot(grid[1, 0:2])
    cax_prob = fig.add_subplot(grid[1, 2])

    mesh_mass = _draw_panel(
        ax_a, mass_0_15.mass, mass_0_15.n_origin_cells,
        "Transport mass (0\u219215 weeks)", mass_norm, panel_label="A",
    )
    _draw_panel(
        ax_b, mass_15_30.mass, mass_15_30.n_origin_cells,
        "Transport mass (15\u219230 weeks)", mass_norm, panel_label="B",
    )
    mesh_prob = _draw_panel(
        ax_c, normalized_15_30, mass_15_30.n_origin_cells,
        "Normalized transport probability (15\u219230 weeks)", prob_norm, panel_label="C",
    )

    cbar_mass = fig.colorbar(mesh_mass, cax=cax_mass)
    cbar_mass.set_label("Transported mass", fontsize=9)
    cbar_prob = fig.colorbar(mesh_prob, cax=cax_prob)
    cbar_prob.set_label("P(destination | origin)", fontsize=9)

    fig.savefig(out_path_stem.with_suffix(".pdf"))
    fig.savefig(out_path_stem.with_suffix(".png"), dpi=FIGURE_DPI)
    plt.close(fig)


def render_individual_heatmap(
    matrix: pd.DataFrame,
    n_origin_cells: pd.Series,
    title: str,
    colorbar_label: str,
    norm: Normalize,
    out_path_stem: Path,
) -> None:
    """Render a single heatmap using the same drawing routine as the
    combined figure, operating on an already-computed matrix."""
    width, height = _panel_extent(matrix)
    fig, ax = plt.subplots(figsize=(width + 1.4, height + 1.1))
    mesh = _draw_panel(ax, matrix, n_origin_cells, title, norm)
    cbar = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label(colorbar_label, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path_stem.with_suffix(".pdf"))
    fig.savefig(out_path_stem.with_suffix(".png"), dpi=FIGURE_DPI)
    plt.close(fig)


def process_interval(h5ad_path: Path, metadata: pd.DataFrame) -> StateAggregation:
    tmap = read_transport_map(h5ad_path)
    transport = np.asarray(tmap.X.todense()) if hasattr(tmap.X, "todense") else np.asarray(tmap.X)
    origin_states = assign_states(tmap.obs_names, metadata)
    dest_states = assign_states(tmap.var_names, metadata)
    return aggregate_transport_by_state(transport, origin_states, dest_states)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    metadata_path = Path(PROCESSED) / "cell_metadata.csv"
    if not metadata_path.exists():
        raise FileNotFoundError(f"Missing {metadata_path}. Run Script 04 first.")
    metadata = pd.read_csv(metadata_path)

    tmap_dir = Path(RESULTS_WOT) / "tmaps"
    mass_0_15 = process_interval(tmap_dir / "hepatocytes_0.0_15.0.h5ad", metadata)
    mass_15_30 = process_interval(tmap_dir / "hepatocytes_15.0_30.0.h5ad", metadata)
    normalized_15_30 = normalize_rows(mass_15_30.mass)

    mass_0_15.mass.to_csv(OUT_DIR / "transport_states_0_15_mass.csv")
    mass_15_30.mass.to_csv(OUT_DIR / "transport_states_15_30_mass.csv")
    normalized_15_30.to_csv(OUT_DIR / "transport_states_15_30_normalized.csv")

    shared_vmax = max(mass_0_15.mass.to_numpy().max(), mass_15_30.mass.to_numpy().max())
    mass_norm = Normalize(vmin=0.0, vmax=shared_vmax)
    prob_norm = Normalize(vmin=0.0, vmax=1.0)

    render_individual_heatmap(
        mass_0_15.mass, mass_0_15.n_origin_cells,
        "Transport mass (0\u219215 weeks)", "Transported mass", mass_norm,
        OUT_DIR / "heatmap_transport_0_15",
    )
    render_individual_heatmap(
        mass_15_30.mass, mass_15_30.n_origin_cells,
        "Transport mass (15\u219230 weeks)", "Transported mass", mass_norm,
        OUT_DIR / "heatmap_transport_15_30",
    )
    render_individual_heatmap(
        normalized_15_30, mass_15_30.n_origin_cells,
        "Normalized transport probability (15\u219230 weeks)", "P(destination | origin)", prob_norm,
        OUT_DIR / "heatmap_transport_15_30_normalized",
    )

    plot_transport_summary(
        mass_0_15, mass_15_30, normalized_15_30, mass_norm, prob_norm, OUT_DIR / "transport_summary"
    )

    print(f"Aggregated matrices: {mass_0_15.mass.shape}, {mass_15_30.mass.shape}")
    print(f"Output written to {OUT_DIR}")


if __name__ == "__main__":
    main()
