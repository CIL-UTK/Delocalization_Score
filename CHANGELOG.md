# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [2.0.0] — 2026-06-18

### Changed
- **Complete rewrite** from Jupyter notebooks (`.ipynb`) to standalone Python scripts (`.py`) for better version control, modularity, and reproducibility.
- All scripts now use a clean configuration block at the top — no inline parameter changes needed.
- Replaced notebook-specific globals with `dataclass` and module-level constants.
- Each script is self-contained with its own `if __name__ == "__main__"` entry point.

### Added
- `requirements.txt` — pinned dependency list for reproducible environments.
- `CITATION.cff` — machine-readable citation metadata (renders on GitHub).
- `.gitignore` — excludes large data files, virtual environments, and IDE artifacts.
- `CHANGELOG.md` — this file.
- Comprehensive `README.md` with pipeline architecture diagram, per-step usage instructions, parameter reference table, and mathematical framework section.
- Lock-mass polynomial drift correction in `2_peak_detection.py`.
- Asymmetric Least Squares (ALS) baseline correction in `2_peak_detection.py`.
- KD-tree accelerated nearest-tissue-pixel distance computation in `4_delocalization_score_calculator.py`.
- Interactive H&E alignment GUI with real-time sliders in `3_mz_overlay.py`.
- Background-to-tissue distance arrow visualization in `5_visualization_tools.py`.
- Tissue contour overlay on spatial heatmaps in `5_visualization_tools.py`.

### Removed
- Jupyter notebook files (`*.ipynb`). The old notebooks are preserved in the `v1.0.0` tag for reference.

## [1.0.0] — 2025-11-23

### Added
- Initial release as Jupyter notebooks.
- imzML parsing and AnnData conversion.
- Peak detection (Savitzky–Golay + SciPy `find_peaks`).
- H&E overlay and polygon ROI annotation.
- Composite delocalization score (area + mean distance).
- Basic visualization of spatial heatmaps.
