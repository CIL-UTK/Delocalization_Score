"""
MSI Delocalization Scoring Pipeline
────────────────────────────────────
Quantifies lipid delocalization from tissue into the background region
using spatial metrics (mean distance, nonzero area).

Workflow:
  1. Load overlaid h5ad (with tissue mask from the annotation tool)
  2. TIC-normalise intensities
  3. Remove m/z features that are not enriched in tissue (fold-change filter)
  4. Mask weak background signals (below a fraction of tissue max)
  5. Compute delocalization metrics per m/z (mean BG-to-tissue distance, area)
  6. Score and rank m/z features
  7. Store metrics in adata.var and save scored h5ad

Folder structure (relative to the working directory):
  <working_dir>/
  ├── overlaid_h5ad/   ← input  (h5ad with tissue column)
  └── results/         ← output (scored h5ad)
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
from scipy.spatial import cKDTree
from tqdm import tqdm

# ═══════════════════════════════════════════════════════════════════════════════
# RESOLVE PATHS RELATIVE TO THE WORKING DIRECTORY
# ═══════════════════════════════════════════════════════════════════════════════
WORKING_DIR = "/home/Delocalization_Score_V2"
OVERLAID_DIR = os.path.join(WORKING_DIR, "overlaid_h5ad")
RESULTS_DIR = os.path.join(WORKING_DIR, "delocalization_score_result_h5ad")

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

# File names
INPUT_FILENAME = "sample_1_overlaid.h5ad"
OUTPUT_FILENAME = "sample_1_scored.h5ad"

# Full paths (assembled automatically)
INPUT_FILE = os.path.join(OVERLAID_DIR, INPUT_FILENAME)
OUTPUT_FILE = os.path.join(RESULTS_DIR, OUTPUT_FILENAME)

# Tissue enrichment filter
TISSUE_FOLD_CHANGE = 1  # keep m/z where mean_tissue >= FC × mean_background

# Background intensity masking
# Set BG pixel to 0 if intensity < this fraction of tissue max for that m/z
MIN_BG_INTENSITY_FRAC = 0.1

# Delocalization score weighting
AREA_WEIGHT = 0.95  # weight for normalised area (1 - weight → mean distance)


# ═══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def _to_dense_col(X, col_idx):
    """Extract a single column from X as a dense 1-D array."""
    if sp.issparse(X):
        return X[:, col_idx].toarray().ravel()
    return np.asarray(X[:, col_idx]).ravel()


def tic_normalise(adata):
    """
    Total Ion Current normalisation (in-place).

    Stores original intensities in adata.layers["raw"] and per-pixel
    TIC values in adata.obs["tic"].
    """
    tic = np.asarray(adata.X.sum(axis=1)).ravel()
    tic[tic == 0] = 1.0

    adata.layers["raw"] = adata.X.copy()

    if sp.issparse(adata.X):
        adata.X = adata.X.multiply(1.0 / tic[:, np.newaxis])
    else:
        adata.X = adata.X / tic[:, np.newaxis]

    adata.obs["tic"] = tic
    print(f"  TIC normalisation applied ({adata.n_obs} pixels)")


def filter_mz_by_tissue(adata, fold_change=2.0):
    """
    Remove m/z features whose mean tissue intensity is below
    `fold_change × mean_background` intensity.
    """
    tissue = adata.obs["tissue"] == 1
    bg = adata.obs["tissue"] == 0

    mean_tissue = np.asarray(adata.X[tissue].mean(axis=0)).ravel()
    mean_bg = np.asarray(adata.X[bg].mean(axis=0)).ravel()

    keep = mean_tissue >= fold_change * mean_bg
    n_drop = (~keep).sum()
    print(f"  Dropping {n_drop} m/z features (tissue < {fold_change}× background)")
    return adata[:, keep].copy()


def mask_low_bg_intensities(adata, frac=0.1):
    """
    Zero out background pixels whose intensity is below
    `frac × max_tissue_intensity` for each m/z.
    """
    X = adata.X
    tissue_mask = adata.obs["tissue"].values == 1
    bg_indices = np.where(~tissue_mask)[0]

    if sp.issparse(X):
        X = X.tocsr()
        max_tissue = X[tissue_mask].max(axis=0).toarray().ravel()
    else:
        max_tissue = X[tissue_mask].max(axis=0)

    thresholds = frac * max_tissue

    for i in bg_indices:
        if sp.issparse(X):
            row = X[i].tocoo()
            row.data[row.data <= thresholds[row.col]] = 0
            X[i] = row.tocsr()
        else:
            below = X[i] <= thresholds
            X[i, below] = 0

    adata.X = X
    print(f"  Masked BG intensities < {frac}× tissue max ({len(bg_indices)} BG pixels)")


# ═══════════════════════════════════════════════════════════════════════════════
# DELOCALIZATION METRICS
# ═══════════════════════════════════════════════════════════════════════════════

def compute_delocalization_metrics(adata):
    """
    For each m/z, compute:
      - bg_dist_mean      : mean distance from nonzero BG pixels to nearest tissue pixel
      - bg_nonzero_area   : number of nonzero BG pixels

    Returns
    -------
    pd.DataFrame with one row per m/z
    """
    mz_values = adata.var["mzs"].values
    n_mz = len(mz_values)

    # Pre-compute tissue / background coordinate arrays
    tissue_mask = adata.obs["tissue"].values == 1
    bg_mask = ~tissue_mask

    tissue_coords = np.column_stack([
        adata.obs.loc[tissue_mask, "x"].values,
        adata.obs.loc[tissue_mask, "y"].values,
    ])
    bg_x = adata.obs.loc[bg_mask, "x"].values
    bg_y = adata.obs.loc[bg_mask, "y"].values

    # Build KD-tree for tissue pixels (shared across all m/z)
    tree = cKDTree(tissue_coords) if tissue_coords.shape[0] > 0 else None

    rows = []
    for j in tqdm(range(n_mz), desc="Computing metrics"):
        intensities = _to_dense_col(adata.X, j)
        bg_int = intensities[bg_mask]

        row = {"mz": mz_values[j]}

        # Nonzero background pixels
        valid = bg_int > 0
        n_valid = valid.sum()
        row["bg_nonzero_area"] = n_valid

        if n_valid > 0 and tree is not None:
            valid_coords = np.column_stack([bg_x[valid], bg_y[valid]])
            distances, _ = tree.query(valid_coords)
            row["bg_dist_mean"] = distances.mean()
        else:
            row["bg_dist_mean"] = np.nan

        rows.append(row)

    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════════════════
# DELOCALIZATION SCORING
# ═══════════════════════════════════════════════════════════════════════════════

def score_delocalization(df, area_weight=0.95):
    """
    Add normalised columns and a weighted delocalization score.

    score = area_weight × norm_area + (1 − area_weight) × norm_mean_dist
    """
    def _minmax(s):
        lo, hi = s.min(), s.max()
        return (s - lo) / (hi - lo) if hi > lo else pd.Series(0.0, index=s.index)

    df["norm_bg_nonzero_area"] = _minmax(df["bg_nonzero_area"])
    df["norm_bg_dist_mean"] = _minmax(df["bg_dist_mean"])
    df["delocalization_score"] = (
        area_weight * df["norm_bg_nonzero_area"]
        + (1 - area_weight) * df["norm_bg_dist_mean"]
    )
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    # --- Load ---
    print(f"Loading: {INPUT_FILE}")
    adata = sc.read_h5ad(INPUT_FILE)

    # --- TIC normalise ---
    tic_normalise(adata)

    # --- Filter m/z by tissue enrichment ---
    adata = filter_mz_by_tissue(adata, fold_change=TISSUE_FOLD_CHANGE)

    # --- Mask weak background signals ---
    mask_low_bg_intensities(adata, frac=MIN_BG_INTENSITY_FRAC)

    # --- Compute delocalization metrics ---
    df_metrics = compute_delocalization_metrics(adata)

    # --- Score and rank ---
    df_metrics = score_delocalization(df_metrics, area_weight=AREA_WEIGHT)

    # --- Store metrics in adata.var ---
    df_metrics = df_metrics.set_index("mz")
    for col in df_metrics.columns:
        adata.var[col] = df_metrics[col].values

    print(f"  Stored metrics in adata.var: {list(df_metrics.columns)}")

    # --- Save scored AnnData ---
    os.makedirs(RESULTS_DIR, exist_ok=True)
    adata.write(OUTPUT_FILE)
    print(f"  Saved: {OUTPUT_FILE}  ({adata.shape})")


if __name__ == "__main__":
    main()
