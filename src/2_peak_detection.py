"""
MSI Peak Processing Pipeline
─────────────────────────────
Processes raw h5ad mass spectrometry imaging data through:
  1. Global spectrum construction
  2. Baseline correction (ALS) + Savitzky-Golay smoothing
  3. Peak detection with minimum distance filtering
  4. Per-pixel intensity summation around selected peaks
  5. Lock-mass polynomial drift correction

Folder structure (relative to this working directory):
  <working_dir>/
  ├── raw_h5ad/                ← input  h5ad files
  └── processed_h5ad/          ← output h5ad files (created automatically)
"""

import os
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import scipy.sparse as sp
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.signal import savgol_filter, find_peaks
from numpy.polynomial import Polynomial
from tqdm import tqdm

# ═══════════════════════════════════════════════════════════════════════════════
# RESOLVE PATHS RELATIVE TO THE WORKING DIRECTORY
# ═══════════════════════════════════════════════════════════════════════════════
WORKING_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
RAW_DIR = os.path.join(WORKING_DIR, "raw_h5ad")
PROCESSED_DIR = os.path.join(WORKING_DIR, "processed_h5ad")

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

# File names (inside the folders above)
INPUT_FILENAME = "sample_1.h5ad"
OUTPUT_FILENAME = "sample_1_peaks.h5ad"

# Full paths (assembled automatically)
INPUT_FILE = os.path.join(RAW_DIR, INPUT_FILENAME)
OUTPUT_FILE = os.path.join(PROCESSED_DIR, OUTPUT_FILENAME)

# Baseline correction (Asymmetric Least Squares)
BASELINE_LAMBDA = 1e5
BASELINE_P = 0.01
BASELINE_NITER = 10

# Peak detection
TOP_MZ_COUNT = 1000
MIN_PROMINENCE = 0.00001
MIN_PEAK_DISTANCE_DA = 0.024

# Savitzky-Golay smoothing
SG_WINDOW_LENGTH = 9
SG_POLYORDER = 3

# Per-pixel intensity summation
NEIGHBOUR_POINTS = 2
ADD_THRESHOLD_DA = 0.012

# Lock-mass alignment
# [LPC 16:0, PC 32:0, PC 34:1, PC 38:6]
REF_MZ = [496.339, 734.569, 760.585, 806.569]
PPM_TOLERANCE = 30

# Decimal precision for corrected m/z values
MZ_DECIMALS = 4


# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: format m/z values to exactly MZ_DECIMALS decimal places
# ═══════════════════════════════════════════════════════════════════════════════

def _fmt_mz(mz):
    """Format a single m/z value as a string with exactly MZ_DECIMALS places."""
    return f"{mz:.{MZ_DECIMALS}f}"


def _fmt_mz_array(mzs):
    """Format an array of m/z values as a list of fixed-decimal strings."""
    return [f"{mz:.{MZ_DECIMALS}f}" for mz in mzs]


# ═══════════════════════════════════════════════════════════════════════════════
# CORE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def create_global_spectrum(X, mz_axis):
    """Sum all pixel spectra to produce a single global (mean) spectrum."""
    if sp.issparse(X):
        X = X.toarray()
    return mz_axis, X.sum(axis=0)


def als_baseline(y, lam=1e5, p=0.01, niter=10):
    """
    Asymmetric Least Squares baseline estimation.

    Parameters
    ----------
    y : array-like   – 1-D spectrum
    lam : float      – smoothness penalty (larger → smoother baseline)
    p : float        – asymmetry weight (0 < p < 1)
    niter : int      – number of re-weighted iterations

    Returns
    -------
    z : ndarray – estimated baseline
    """
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, 1, 2], shape=(L - 2, L))
    DTD = lam * D.T @ D
    w = np.ones(L)
    for _ in range(niter):
        W = sparse.diags(w, 0)
        z = spsolve(W + DTD, w * y)
        w = np.where(y > z, p, 1 - p)
    return z


def detect_top_peaks(
    mz_axis, global_spec, *,
    window_length=9, polyorder=3,
    min_prominence=0.01, top_n=100,
    min_distance_da=0.024,
    baseline_lam=1e5, baseline_p=0.01, baseline_niter=10,
):
    """
    Detect top-N peaks from the global summed spectrum.

    Steps: baseline correction → smoothing → peak finding →
           minimum-distance filtering → top-N selection.
    """
    # Baseline correction
    baseline = als_baseline(global_spec, lam=baseline_lam, p=baseline_p, niter=baseline_niter)
    corrected = np.clip(global_spec - baseline, 0, None)

    # Savitzky-Golay smoothing
    smoothed = savgol_filter(corrected, window_length=window_length, polyorder=polyorder)

    # Initial peak detection
    peaks, _ = find_peaks(smoothed, prominence=min_prominence)
    print(f"  Initial detected peaks : {len(peaks)}")

    mz_peaks = mz_axis[peaks]
    intensities = smoothed[peaks]

    # Sort by descending intensity
    order = np.argsort(intensities)[::-1]
    mz_sorted = mz_peaks[order]
    int_sorted = intensities[order]

    # Greedy filter: keep a peak only if no stronger peak exists within ±da
    kept_mz, kept_int = [], []
    accepted = np.zeros(len(mz_sorted), dtype=bool)

    for i in range(len(mz_sorted)):
        if accepted[i]:
            continue
        kept_mz.append(mz_sorted[i])
        kept_int.append(int_sorted[i])
        accepted |= np.abs(mz_sorted - mz_sorted[i]) <= min_distance_da

    kept_mz = np.array(kept_mz)
    kept_int = np.array(kept_int)
    print(f"  After ±{min_distance_da} Da filter : {len(kept_mz)}")

    # Top-N selection
    if top_n is not None and len(kept_mz) > top_n:
        top_idx = np.argsort(kept_int)[-top_n:][::-1]
        kept_mz = kept_mz[top_idx]
        kept_int = kept_int[top_idx]

    # Final sort by descending intensity
    final_order = np.argsort(kept_int)[::-1]
    return kept_mz[final_order], kept_int[final_order]


def sum_peaks_per_pixel(adata, target_mzs, flank=2, mz_threshold=0.02):
    """
    For each pixel, sum intensities in a window around each target m/z.

    Parameters
    ----------
    adata : AnnData        – full-resolution data
    target_mzs : array     – selected peak m/z values
    flank : int            – number of index neighbours on each side
    mz_threshold : float   – max Dalton distance for a neighbour to count

    Returns
    -------
    AnnData with shape (n_spots, n_targets)
    """
    mz_axis = adata.var_names.astype(float).values
    X = adata.X
    is_sparse = sp.issparse(X)

    n_spots = adata.n_obs
    n_targets = len(target_mzs)
    result = np.zeros((n_spots, n_targets), dtype=np.float32)

    # Pre-compute valid index sets per target
    nearest_idx = np.clip(np.searchsorted(mz_axis, target_mzs), 0, len(mz_axis) - 1)
    valid_indices = []
    for i, tmz in enumerate(target_mzs):
        lo = max(nearest_idx[i] - flank, 0)
        hi = min(nearest_idx[i] + flank + 1, len(mz_axis))
        candidates = np.arange(lo, hi)
        mask = np.abs(mz_axis[candidates] - tmz) <= mz_threshold
        valid_indices.append(candidates[mask])

    # Pixel-wise summation
    for i in tqdm(range(n_spots), desc="Summing intensities per pixel"):
        spectrum = X[i, :].toarray().ravel() if is_sparse else X[i, :]
        for j, idx in enumerate(valid_indices):
            result[i, j] = spectrum[idx].sum()

    new_adata = anndata.AnnData(X=result, obs=adata.obs.copy())
    new_adata.var = pd.DataFrame(index=_fmt_mz_array(target_mzs))
    return new_adata


def lock_mass_correction(adata, lock_mzs, ppm_tol=20, max_deg=2):
    """
    Correct m/z drift using lock-mass references and best-fit polynomial.

    Selects the polynomial degree (1 or 2) that minimises MSE on the
    lock-mass residuals, then applies the correction across all m/z values.
    """
    if "mzs" not in adata.var.columns:
        raise ValueError("adata.var must contain a 'mzs' column with numeric m/z values.")

    mzs = adata.var["mzs"].values.astype(float)
    measured, deltas = [], []

    for ref in lock_mzs:
        window = ref * ppm_tol * 1e-6
        candidates = np.where((mzs >= ref - window) & (mzs <= ref + window))[0]
        if len(candidates) == 0:
            raise ValueError(f"No m/z found within {ppm_tol} ppm of lock-mass {ref:.3f}")
        closest = candidates[np.argmin(np.abs(mzs[candidates] - ref))]
        measured.append(mzs[closest])
        deltas.append(mzs[closest] - ref)

    measured = np.array(measured)
    deltas = np.array(deltas)

    # Select best polynomial degree by MSE
    best_poly, best_deg, best_mse = None, 1, np.inf
    for deg in range(1, max_deg + 1):
        poly = Polynomial.fit(measured, deltas, deg=deg)
        mse = np.mean((poly(measured) - deltas) ** 2)
        if mse < best_mse:
            best_mse, best_deg, best_poly = mse, deg, poly

    # Apply correction
    corrected = mzs - best_poly(mzs)
    adata.var["mz_corrected"] = np.round(corrected, MZ_DECIMALS)

    print(f"  Lock-mass correction : degree-{best_deg} polynomial, "
          f"lock-masses {lock_mzs}")
    return adata


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    # --- Load ---
    print(f"Reading : {INPUT_FILE}")
    adata = sc.read_h5ad(INPUT_FILE)

    # --- Global spectrum ---
    mz_axis, global_spec = create_global_spectrum(
        adata.X, adata.var_names.astype(float).values
    )

    # --- Peak detection ---
    print("Peak detection …")
    selected_mz, selected_int = detect_top_peaks(
        mz_axis, global_spec,
        window_length=SG_WINDOW_LENGTH,
        polyorder=SG_POLYORDER,
        min_prominence=MIN_PROMINENCE,
        top_n=TOP_MZ_COUNT,
        min_distance_da=MIN_PEAK_DISTANCE_DA,
        baseline_lam=BASELINE_LAMBDA,
        baseline_p=BASELINE_P,
        baseline_niter=BASELINE_NITER,
    )

    # --- Per-pixel intensity summation ---
    print("Per-pixel summation …")
    adata_peaks = sum_peaks_per_pixel(
        adata, selected_mz,
        flank=NEIGHBOUR_POINTS,
        mz_threshold=ADD_THRESHOLD_DA,
    )

    # --- Lock-mass alignment ---
    adata_peaks.var["mzs"] = adata_peaks.var.index.astype(float).values
    adata_peaks = lock_mass_correction(adata_peaks, REF_MZ, ppm_tol=PPM_TOLERANCE)

    # Rename columns & index with exactly MZ_DECIMALS decimal places
    adata_peaks.var["mzs_original"] = adata_peaks.var["mzs"]
    adata_peaks.var["mzs"] = adata_peaks.var["mz_corrected"]
    adata_peaks.var_names = pd.Index(_fmt_mz_array(adata_peaks.var["mz_corrected"]))
    adata_peaks.var_names.name = None
    adata_peaks.var.drop(columns=["mz_corrected"], inplace=True)

    # Sort by corrected m/z
    adata_peaks = adata_peaks[:, adata_peaks.var.sort_values("mzs").index]

    # --- Save ---
    os.makedirs(PROCESSED_DIR, exist_ok=True)
    adata_peaks.write(OUTPUT_FILE)
    print(f"Saved: {OUTPUT_FILE}  ({adata_peaks.shape})")


if __name__ == "__main__":
    main()
