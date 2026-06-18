"""
imzML → AnnData Conversion Pipeline
====================================
Converts mass spectrometry imaging (imzML) data into an AnnData (.h5ad) object
with spatial metadata, TIC normalization, and experiment annotations.

Folder structure (relative to WORKING_DIR):
  WORKING_DIR/
  ├── imzml_data/          ← input  .imzML + .ibd files
  └── raw_h5ad/            ← output .h5ad files (created automatically)
"""

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from pyimzml.ImzMLParser import ImzMLParser
from scipy.sparse import csr_matrix


# ═══════════════════════════════════════════════════════════════════════════════
# RESOLVE PATHS RELATIVE TO THE WORKING DIRECTORY
# ═══════════════════════════════════════════════════════════════════════════════
WORKING_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
IMZML_DIR   = os.path.join(WORKING_DIR, "imzml_data")
RAW_DIR     = os.path.join(WORKING_DIR, "raw_h5ad")

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class PipelineConfig:
    """All user-facing parameters in one place."""

    # File names (just the stem — paths are assembled automatically)
    input_filename:  str = "sample_1.imzml"
    output_filename: str = "sample_1.h5ad"

    # Experiment metadata
    sample:         str = "1-1"
    batch:          str = "Slide_1"
    age_group:      str = "Aged"
    disease_status: str = "AD"

    # Spatial / processing
    image_size_cm:    float = 20.0
    image_resolution: float = 0.001   # cm per pixel
    batch_size:       int   = 5_000
    mz_decimals:      int   = 6       # rounding precision for m/z alignment

    # --- Derived full paths (read-only) ---
    @property
    def input_path(self) -> str:
        return os.path.join(IMZML_DIR, self.input_filename)

    @property
    def output_path(self) -> str:
        return os.path.join(RAW_DIR, self.output_filename)


# ──────────────────────────────────────────────
# imzML Reader
# ──────────────────────────────────────────────

class ImzMLReader:
    """Thin wrapper around pyimzml that exposes coordinates and spectra."""

    def __init__(self, path: str):
        self.parser = ImzMLParser(path)
        self.coordinates = np.array(
            [(x, y) for x, y, *_ in self.parser.coordinates]
        )

    @property
    def n_spectra(self) -> int:
        return len(self.coordinates)

    def get_spectrum(self, index: int):
        mzs, ints = self.parser.getspectrum(index)
        return np.asarray(mzs, dtype=np.float64), np.asarray(ints, dtype=np.float64)


# ──────────────────────────────────────────────
# Core Conversion
# ──────────────────────────────────────────────

def _collect_global_mzs(reader: ImzMLReader, decimals: int) -> np.ndarray:
    """Single pass to collect and deduplicate all m/z values."""
    all_mzs = []
    for i in range(reader.n_spectra):
        mzs, _ = reader.get_spectrum(i)
        all_mzs.append(np.round(mzs, decimals))
    return np.unique(np.concatenate(all_mzs))


def _build_sparse_matrix(
    reader: ImzMLReader,
    global_mzs: np.ndarray,
    decimals: int,
    batch_size: int,
) -> csr_matrix:
    """
    Build a (spectra × m/z) sparse matrix using vectorised per-spectrum
    operations instead of a Python-level inner loop.
    """
    mz_to_col = {mz: i for i, mz in enumerate(global_mzs)}
    n_spectra, n_mzs = reader.n_spectra, len(global_mzs)

    rows, cols, vals = [], [], []

    for start in range(0, n_spectra, batch_size):
        end = min(start + batch_size, n_spectra)
        for i in range(start, end):
            mzs, intensities = reader.get_spectrum(i)
            if mzs.size == 0:
                continue

            # Drop zeros & round in one shot
            nonzero = intensities != 0
            mzs = np.round(mzs[nonzero], decimals)
            intensities = intensities[nonzero]

            # Vectorised column lookup
            col_indices = np.array(
                [mz_to_col[mz] for mz in mzs if mz in mz_to_col]
            )
            mask = np.isin(mzs, global_mzs)
            rows.extend([i] * col_indices.size)
            cols.extend(col_indices)
            vals.extend(intensities[mask])

        print(f"  spectra {start:>7,d} – {end:>7,d}  /  {n_spectra:,d}")

    return csr_matrix(
        (np.array(vals, dtype=np.float64),
         (np.array(rows, dtype=np.int64),
          np.array(cols, dtype=np.int64))),
        shape=(n_spectra, n_mzs),
    )


def _attach_spatial(adata: sc.AnnData, size_cm: float, resolution: float):
    """
    Create a minimal Squidpy/Scanpy-compatible spatial dict.
    Uses a small placeholder image to avoid multi-GB allocations.
    """
    PLACEHOLDER_PX = 512
    blank = np.zeros((PLACEHOLDER_PX, PLACEHOLDER_PX, 3), dtype=np.uint8)

    scale_factor = PLACEHOLDER_PX / (size_cm / resolution)

    adata.uns["spatial"] = {
        "spatial": {
            "scalefactors": {
                "tissue_hires_scalef": scale_factor,
                "spot_diameter_fullres": 0.05,
            },
            "images": {
                "hires": blank,
                "lowres": blank,
            },
        }
    }


def build_anndata(reader: ImzMLReader, cfg: PipelineConfig) -> sc.AnnData:
    """End-to-end construction of an annotated AnnData from an imzML reader."""

    print(f"Total spectra: {reader.n_spectra:,d}")

    # 1 ── Global m/z index
    print("Collecting global m/z values …")
    global_mzs = _collect_global_mzs(reader, cfg.mz_decimals)
    print(f"  → {len(global_mzs):,d} unique m/z values")

    # 2 ── Sparse intensity matrix
    print("Building sparse matrix …")
    X = _build_sparse_matrix(reader, global_mzs, cfg.mz_decimals, cfg.batch_size)

    # 3 ── Assemble AnnData
    adata = sc.AnnData(X=X)

    # m/z as variable metadata & index
    adata.var["mz"] = global_mzs
    adata.var_names = [str(mz) for mz in global_mzs]

    # Normalised coordinates
    coords = reader.coordinates.astype(np.float64)
    coords -= coords.min(axis=0)
    adata.obs["x"] = coords[:, 0]
    adata.obs["y"] = coords[:, 1]
    adata.obsm["spatial"] = coords  # Squidpy convention

    # TIC & TIC-normalised layer
    tic = np.asarray(X.sum(axis=1)).ravel()
    adata.obs["TIC"] = tic
    tic_safe = np.where(tic == 0, 1.0, tic)
    adata.layers["tic_normalized"] = X.multiply(1.0 / tic_safe[:, None]).tocsr()

    # Experiment metadata (categorical to save memory)
    for key, value in {
        "sample": cfg.sample,
        "batch": cfg.batch,
        "age_group": cfg.age_group,
        "disease_status": cfg.disease_status,
    }.items():
        adata.obs[key] = pd.Categorical([value] * adata.n_obs)

    # Store raw counts before any downstream filtering
    adata.raw = adata.copy()

    # Spatial visualisation scaffold
    _attach_spatial(adata, cfg.image_size_cm, cfg.image_resolution)

    print(f"AnnData ready — shape: {adata.shape}")
    return adata


# ──────────────────────────────────────────────
# Entry Point
# ──────────────────────────────────────────────

if __name__ == "__main__":
    cfg = PipelineConfig()

    # Ensure output directory exists
    Path(cfg.output_path).parent.mkdir(parents=True, exist_ok=True)

    reader = ImzMLReader(cfg.input_path)
    adata = build_anndata(reader, cfg)

    adata.write(cfg.output_path)
    print(f"Saved → {cfg.output_path}")
