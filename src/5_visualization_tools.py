"""
MSI Delocalization Visualisation Tool
─────────────────────────────────────
Reads a scored h5ad file (produced by the scoring pipeline) and provides
interactive visualisation of individual m/z features:

  • Heatmap view   – spatial intensity map with tissue contour overlay
  • Distance map   – arrows from nonzero BG pixels to nearest tissue pixel

Folder structure (relative to the working directory):
  <working_dir>/
  └── results/           ← input  (scored h5ad from the scoring pipeline)
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
from skimage.measure import find_contours

matplotlib.use("Qt5Agg")  # Required for VSCode interactivity

# ═══════════════════════════════════════════════════════════════════════════════
# RESOLVE PATHS RELATIVE TO THE WORKING DIRECTORY
# ═══════════════════════════════════════════════════════════════════════════════
WORKING_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
RESULTS_DIR = os.path.join(WORKING_DIR, "delocalization_score_result_h5ad")

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

# Scored h5ad to visualise
INPUT_FILENAME = "sample_1_scored.h5ad"
INPUT_FILE = os.path.join(RESULTS_DIR, INPUT_FILENAME)

# m/z value to visualise (set to None to use the top-ranked feature)
TARGET_MZ = None

# Distance-map arrow count
N_ARROWS = 100

# Custom colourmap (dark → navy → blue → purple → red → orange → yellow → white)
CUSTOM_CMAP = LinearSegmentedColormap.from_list("msi_heatmap", [
    (0.00,       "#454545"),
    (0.00000001, "#000000"),
    (0.10,       "#000080"),
    (0.15,       "#0000FF"),
    (0.30,       "#8000FF"),
    (0.45,       "#FF0000"),
    (0.60,       "#FF8000"),
    (0.75,       "#FFFF00"),
    (1.00,       "#FFFFFF"),
])


# ═══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def _to_dense_col(X, col_idx):
    """Extract a single column from X as a dense 1-D array."""
    if sp.issparse(X):
        return X[:, col_idx].toarray().ravel()
    return np.asarray(X[:, col_idx]).ravel()


def _resolve_mz_index(adata, mz_value):
    """Return the column index of the closest m/z in adata.var['mzs']."""
    mz_axis = adata.var["mzs"].values.astype(float)
    idx = np.argmin(np.abs(mz_axis - mz_value))
    return idx, mz_axis[idx]


def tissue_contours(adata):
    """
    Compute tissue boundary contours from the binary tissue mask.

    Returns
    -------
    list of (contour_x, contour_y) arrays
    """
    x = adata.obs["x"].values
    y = adata.obs["y"].values
    tissue = adata.obs["tissue"].values == 1

    grid_x = np.unique(x)
    grid_y = np.unique(y)
    x_map = {v: i for i, v in enumerate(grid_x)}
    y_map = {v: i for i, v in enumerate(grid_y)}

    mask_2d = np.zeros((len(grid_y), len(grid_x)), dtype=bool)
    for xi, yi, t in zip(x, y, tissue):
        if t:
            mask_2d[y_map[yi], x_map[xi]] = True

    contours = find_contours(mask_2d.astype(float), level=0.5)
    return [
        (grid_x[c[:, 1].astype(int)], grid_y[c[:, 0].astype(int)])
        for c in contours
    ]


def compute_bg_to_tissue_pairs(adata, mz_index):
    """
    For a single m/z column, find every nonzero background pixel and its
    nearest tissue pixel (coordinates + distance).

    Returns
    -------
    pd.DataFrame with columns: bg_x, bg_y, tissue_x, tissue_y, distance
    (empty DataFrame if no valid BG pixels)
    """
    tissue_mask = adata.obs["tissue"].values == 1
    bg_mask = ~tissue_mask

    intensities = _to_dense_col(adata.X, mz_index)
    bg_int = intensities[bg_mask]

    valid = bg_int > 0
    if valid.sum() == 0:
        return pd.DataFrame(columns=["bg_x", "bg_y", "tissue_x", "tissue_y", "distance"])

    bg_x = adata.obs.loc[bg_mask, "x"].values[valid]
    bg_y = adata.obs.loc[bg_mask, "y"].values[valid]

    tissue_coords = np.column_stack([
        adata.obs.loc[tissue_mask, "x"].values,
        adata.obs.loc[tissue_mask, "y"].values,
    ])

    tree = cKDTree(tissue_coords)
    distances, indices = tree.query(np.column_stack([bg_x, bg_y]))
    nearest = tissue_coords[indices]

    return pd.DataFrame({
        "bg_x": bg_x,
        "bg_y": bg_y,
        "tissue_x": nearest[:, 0],
        "tissue_y": nearest[:, 1],
        "distance": distances,
    })


# ═══════════════════════════════════════════════════════════════════════════════
# HEATMAP VISUALISATION
# ═══════════════════════════════════════════════════════════════════════════════

def visualise_heatmap(adata, mz_value, contours, cmap=None, save_path=None):
    """
    Plot spatial intensity heatmap for a single m/z with tissue contour overlay.

    Parameters
    ----------
    adata : AnnData       – scored h5ad
    mz_value : float      – target m/z
    contours : list       – from tissue_contours()
    cmap : colormap       – default: CUSTOM_CMAP
    save_path : str|None  – if provided, save figure and close
    """
    if cmap is None:
        cmap = CUSTOM_CMAP

    mz_idx, matched_mz = _resolve_mz_index(adata, mz_value)

    x = adata.obs["x"].values
    y = adata.obs["y"].values
    intensities = _to_dense_col(adata.X, mz_idx)

    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(x, y, c=intensities, s=12, marker="s", cmap=cmap)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(f"Intensity @ m/z {matched_mz:.4f}", fontsize=18)
    cbar.ax.tick_params(labelsize=16)

    for cx, cy in contours:
        ax.plot(cx, cy, color="red", linewidth=2)

    ax.set_title(f"m/z {matched_mz:.4f}", fontsize=20)
    ax.set_xlabel("x", fontsize=18)
    ax.set_ylabel("y", fontsize=18)
    ax.tick_params(labelsize=16)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300)
        plt.close(fig)
    else:
        plt.show()


# ═══════════════════════════════════════════════════════════════════════════════
# DISTANCE MAP VISUALISATION
# ═══════════════════════════════════════════════════════════════════════════════

def visualise_distance_map(adata, mz_value, n_arrows=100, save_path=None):
    """
    Plot nonzero background pixels and arrows to their nearest tissue pixel.

    Parameters
    ----------
    adata : AnnData       – scored h5ad
    mz_value : float      – target m/z
    n_arrows : int        – max arrows to draw (sampled randomly)
    save_path : str|None  – if provided, save figure and close
    """
    mz_idx, matched_mz = _resolve_mz_index(adata, mz_value)
    df = compute_bg_to_tissue_pairs(adata, mz_idx)

    if df.empty:
        print(f"  No nonzero BG pixels for m/z {matched_mz:.4f}")
        return

    fig, ax = plt.subplots(figsize=(10, 8))

    ax.scatter(df["bg_x"], df["bg_y"], c="blue", s=10)
    ax.scatter(df["tissue_x"], df["tissue_y"], c="red", s=10)

    # Draw a random sample of arrows
    sample = df.sample(n=min(n_arrows, len(df)), random_state=42)
    for _, row in sample.iterrows():
        ax.arrow(
            row["bg_x"], row["bg_y"],
            row["tissue_x"] - row["bg_x"],
            row["tissue_y"] - row["bg_y"],
            color="gray", alpha=0.8, head_width=0.5,
        )

    # Legend
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="blue",
               markersize=10, label="Nonzero BG pixels"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="red",
               markersize=10, label="Nearest tissue pixel"),
        Line2D([0], [0], color="gray", lw=2, label="BG → Tissue direction"),
    ]

    ax.set_title(f"BG → Tissue Distance Map  ·  m/z {matched_mz:.4f}", fontsize=18)
    ax.set_xlabel("x", fontsize=16)
    ax.set_ylabel("y", fontsize=16)
    ax.tick_params(labelsize=14)
    ax.legend(handles=legend_handles, fontsize=13)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300)
        plt.close(fig)
    else:
        plt.show()


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    # --- Load scored h5ad ---
    print(f"Loading: {INPUT_FILE}")
    adata = sc.read_h5ad(INPUT_FILE)
    print(f"  Shape: {adata.shape}")

    # --- Determine target m/z ---
    if TARGET_MZ is not None:
        mz = TARGET_MZ
    elif "delocalization_score" in adata.var.columns:
        # Use the top-ranked (most delocalized) m/z
        mz = float(adata.var.sort_values("delocalization_score", ascending=False).index[0])
        print(f"  Using top-ranked m/z: {mz}")
    else:
        mz = float(adata.var_names[0])
        print(f"  No score column found — using first m/z: {mz}")

    # --- Compute contours once ---
    contours = tissue_contours(adata)

    # --- Show heatmap ---
    print(f"  Plotting heatmap for m/z {mz:.4f} …")
    visualise_heatmap(adata, mz, contours)

    # --- Show distance map ---
    print(f"  Plotting distance map for m/z {mz:.4f} …")
    visualise_distance_map(adata, mz, n_arrows=N_ARROWS)


if __name__ == "__main__":
    main()
