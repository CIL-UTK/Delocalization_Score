# Delocalization Score for MALDI-MSI

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facsmeasuresciau.5c00148-blue)](https://doi.org/10.1021/acsmeasuresciau.5c00148)

A Python-based pipeline for **quantifying analyte delocalization** in matrix-assisted laser desorption/ionization mass spectrometry imaging (MALDI-MSI) data. This toolkit implements a composite scoring approach that integrates spatial metrics — nonzero off-tissue area and mean background-to-tissue distance — into a single interpretable delocalization score.

> **Associated Publication:**  
> Jarrahi, A., Jones, A., Tang, W., Qi, H., & Crouch, A. C. (2026). *Mathematical Framework for Quantifying Delocalization in MALDI-MSI via a Composite Scoring Approach*. ACS Measurement Science Au, 6(1), 134–149.  
> [https://doi.org/10.1021/acsmeasuresciau.5c00148](https://doi.org/10.1021/acsmeasuresciau.5c00148)

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Installation](#installation)
- [Directory Structure](#directory-structure)
- [Usage](#usage)
  - [Step 1 — imzML to AnnData Conversion](#step-1--imzml-to-anndata-conversion)
  - [Step 2 — Peak Detection & Lock-Mass Correction](#step-2--peak-detection--lock-mass-correction)
  - [Step 3 — H&E Overlay & Tissue Annotation](#step-3--he-overlay--tissue-annotation)
  - [Step 4 — Delocalization Score Calculation](#step-4--delocalization-score-calculation)
  - [Step 5 — Visualization](#step-5--visualization)
- [Configuration Reference](#configuration-reference)
- [Mathematical Framework](#mathematical-framework)
- [Data Availability](#data-availability)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Overview

Delocalization — the spatial spread of analytes beyond tissue boundaries — is a common artifact in MALDI-MSI that degrades spatial resolution and sensitivity. This pipeline provides a **reproducible, data-driven metric** to quantify delocalization across all detected m/z features in an MSI dataset.

**Key capabilities:**

| Feature | Description |
|---|---|
| **imzML ingestion** | Parses raw imzML/ibd files into sparse AnnData objects with full spatial metadata |
| **Peak detection** | Baseline correction (ALS), Savitzky–Golay smoothing, and proximity-filtered peak picking |
| **Lock-mass correction** | Polynomial drift correction using user-specified reference masses |
| **Interactive registration** | GUI-based alignment of H&E histology images with MSI ion heatmaps |
| **Tissue annotation** | Polygon-based ROI drawing to define tissue vs. background regions |
| **Delocalization scoring** | Composite metric combining normalized off-tissue area and mean BG-to-tissue distance |
| **Visualization** | Spatial heatmaps with tissue contours and BG → tissue distance arrow maps |

---

## Pipeline Architecture

```
┌─────────────────────┐
│   imzML + ibd files  │
└─────────┬───────────┘
          ▼
┌─────────────────────┐     ┌──────────────────────────────────────────────┐
│  1_imzml_reader.py   │────▶│  raw_h5ad/sample_1.h5ad                     │
│  imzML → AnnData     │     │  (sparse matrix, spatial coords, TIC norm)  │
└─────────────────────┘     └──────────────────┬───────────────────────────┘
                                               ▼
┌─────────────────────┐     ┌──────────────────────────────────────────────┐
│  2_peak_detection.py │────▶│  processed_h5ad/sample_1_peaks.h5ad         │
│  ALS + SG + peaks    │     │  (top-N peaks, lock-mass corrected m/z)     │
└─────────────────────┘     └──────────────────┬───────────────────────────┘
                                               ▼
┌─────────────────────┐     ┌──────────────────────────────────────────────┐
│  3_mz_overlay.py     │────▶│  overlaid_h5ad/sample_1_overlaid.h5ad       │
│  H&E alignment + ROI │     │  (binary tissue mask in obs["tissue"])      │
└─────────────────────┘     └──────────────────┬───────────────────────────┘
                                               ▼
┌──────────────────────────┐ ┌─────────────────────────────────────────────┐
│ 4_delocalization_score   │─▶│ results/sample_1_scored.h5ad               │
│ _calculator.py            │ │ (delocalization scores in var columns)      │
└──────────────────────────┘ └──────────────────┬──────────────────────────┘
                                                ▼
┌──────────────────────────┐ ┌─────────────────────────────────────────────┐
│ 5_visualization_tools.py │─▶│ Heatmaps, distance maps, score rankings    │
└──────────────────────────┘ └─────────────────────────────────────────────┘
```

---

## Installation

### Prerequisites

- Python ≥ 3.9
- A Qt5-compatible backend for matplotlib (required for interactive GUI tools in Steps 3 and 5)

### Setup

```bash
# Clone the repository
git clone https://github.com/CIL-UTK/Delocalization_Score.git
cd Delocalization_Score

# Create and activate a virtual environment (recommended)
python -m venv venv
source venv/bin/activate        # Linux / macOS
# venv\Scripts\activate         # Windows

# Install dependencies
pip install -r requirements.txt
```

### Dependencies

| Package | Purpose |
|---|---|
| `numpy` | Numerical computation |
| `scipy` | Sparse matrices, signal processing, spatial trees |
| `pandas` | Tabular data handling |
| `anndata` | AnnData storage format |
| `scanpy` | Single-cell / spatial data utilities |
| `pyimzml` | imzML file parsing |
| `matplotlib` | Plotting and interactive GUI |
| `scikit-image` | Image transformation and contour detection |
| `tqdm` | Progress bars |
| `PyQt5` | Qt5 backend for matplotlib (interactive overlay/annotation) |

---

## Directory Structure

Before running the pipeline, organize your working directory as follows:

```
Delocalization_Score/
├── src/                              # Pipeline scripts
│   ├── 1_imzml_reader.py
│   ├── 2_peak_detection.py
│   ├── 3_mz_overlay.py
│   ├── 4_delocalization_score_calculator.py
│   └── 5_visualization_tools.py
│
├── data/                             # ← Create this; place your data here
│   ├── imzml_data/                   # Input imzML + ibd files
│   ├── H&E_images/                   # H&E stained tissue images
│   ├── raw_h5ad/                     # Auto-created by Step 1
│   ├── processed_h5ad/               # Auto-created by Step 2
│   ├── overlaid_h5ad/                # Auto-created by Step 3
│   └── delocalization_score_result_h5ad/  # Auto-created by Step 4
│
├── notebooks/                        # Jupyter notebooks (examples / tutorials)
│   └── example_workflow.ipynb
│
├── requirements.txt
├── LICENSE
├── CITATION.cff
├── .gitignore
└── README.md
```

> **Note:** Update the `WORKING_DIR` variable at the top of each script to point to your local data directory. All sub-folder paths are derived from it automatically.

---

## Data

### 📥 Download the Dataset

The dataset (22 MALDI-MSI samples + H&E images) is hosted externally due to file size.

**[Download Dataset (SharePoint)](https://liveutk-my.sharepoint.com/:f:/g/personal/ajarrah1_vols_utk_edu/EqUysi_9h2BOnOSajIEcq64BhPRXaRExNfXPcdQ12BxlRg?e=wpZNMZ)**

Password:
data@CIL-UTK

After downloading, place the files in the corresponding subfolders under `data/`. See [`data/README.md`](data/README.md) for full details.

### File Naming Conventions

**imzML files** follow this pattern:

```
YYYYMMDD_<model><group>_<age><sex>_<slide>_pos_<animal#>_<section>
```

Examples:
- `20230828_3xTgControl_youngfemale_5_pos_1_1` → Control, young female, slide 1, #5, section 1
- `20230914_3xTgAD_youngfemale_1_pos_3_1` → AD, young female, slide 3, #1, section 1

**H&E image files** use a shorthand:

| Code | Meaning | Example |
|---|---|---|
| `y` / `o` | Young / Old (Aged) | `yfm1_c` → Young female 1, Control |
| `f` / `m` | Female / Male | `ofm1_ad` → Old female 1, AD |
| `_c` / `_ad` | Control / AD | `ofm2_ad` → Old female 2, AD |

---

## Usage

Each script is designed to run independently. Edit the **CONFIGURATION** block at the top of each file to set filenames, parameters, and metadata before running.

### Step 1 — imzML to AnnData Conversion

**Script:** `1_imzml_reader.py`

Converts raw imzML/ibd mass spectrometry imaging files into AnnData (`.h5ad`) objects with:
- A sparse (CSR) intensity matrix (pixels × m/z bins)
- Spatial coordinates in `obs["x"]`, `obs["y"]`, and `obsm["spatial"]`
- TIC-normalized intensities in `layers["tic_normalized"]`
- Experiment metadata as categorical columns (sample, batch, age group, disease status)

**Configuration:**
```python
# In PipelineConfig:
input_filename  = "sample_1.imzml"     # imzML file in imzml_data/
output_filename = "sample_1.h5ad"       # output in raw_h5ad/
sample          = "1-1"
batch           = "Slide_1"
age_group       = "Aged"
disease_status  = "AD"
mz_decimals     = 6                     # rounding precision for m/z alignment
```

**Run:**
```bash
python src/1_imzml_reader.py
```

---

### Step 2 — Peak Detection & Lock-Mass Correction

**Script:** `2_peak_detection.py`

Processes the raw h5ad through:
1. **Global spectrum construction** — summing all pixel spectra
2. **Baseline correction** — Asymmetric Least Squares (ALS)
3. **Smoothing** — Savitzky–Golay filter
4. **Peak detection** — `scipy.signal.find_peaks` with prominence filtering
5. **Proximity filtering** — greedy minimum-distance deduplication (default ±0.024 Da)
6. **Top-N selection** — retain the 1000 most intense peaks
7. **Per-pixel intensity summation** — integrate signal around each peak (±2 indices, ±0.012 Da)
8. **Lock-mass polynomial correction** — drift correction using reference masses

**Configuration:**
```python
INPUT_FILENAME  = "sample_1.h5ad"
OUTPUT_FILENAME = "sample_1_peaks.h5ad"

# Baseline correction (ALS)
BASELINE_LAMBDA = 1e5
BASELINE_P      = 0.01

# Peak detection
TOP_MZ_COUNT          = 1000
MIN_PEAK_DISTANCE_DA  = 0.024

# Savitzky-Golay smoothing
SG_WINDOW_LENGTH = 9
SG_POLYORDER     = 3

# Lock-mass references [LPC 16:0, PC 32:0, PC 34:1, PC 38:6]
REF_MZ        = [496.339, 734.569, 760.585, 806.569]
PPM_TOLERANCE = 30
```

**Run:**
```bash
python src/2_peak_detection.py
```

---

### Step 3 — H&E Overlay & Tissue Annotation

**Script:** `3_mz_overlay.py`

> ⚠️ **Requires a display** — this step launches interactive matplotlib windows.

Interactive GUI workflow:
1. Generates a heatmap for a target m/z (default: 800.55)
2. Opens an alignment window with sliders for **rotation**, **scale X/Y**, **translation X/Y**, and **alpha blending** to register the H&E image onto the MSI heatmap
3. After alignment (close the window), opens a polygon annotation tool for drawing a tissue ROI
4. Saves the binary tissue mask into `adata.obs["tissue"]`

**Configuration:**
```python
INPUT_FILENAME    = "sample_1_peaks.h5ad"
HE_IMAGE_FILENAME = "sample_1_HE.jpg"
OUTPUT_FILENAME   = "sample_1_overlaid.h5ad"
TARGET_MZ         = 800.55
MZ_TOLERANCE      = 0.1
```

**Run:**
```bash
python src/3_mz_overlay.py
```

---

### Step 4 — Delocalization Score Calculation

**Script:** `4_delocalization_score_calculator.py`

Core analysis step:
1. **TIC normalization** — corrects for pixel-level signal variation
2. **Fold-change filtering** — removes m/z features not enriched in tissue (configurable threshold)
3. **Background masking** — zeros out BG pixels with intensity < 10% of tissue max (per m/z)
4. **Metric computation** — for each m/z:
   - `bg_nonzero_area`: count of nonzero off-tissue pixels
   - `bg_dist_mean`: mean Euclidean distance from nonzero BG pixels to nearest tissue pixel (via KD-tree)
5. **Composite scoring** — min-max normalization followed by weighted linear combination:

$$
s_i = \alpha \cdot \tilde{a}_i + (1 - \alpha) \cdot \tilde{d}_i
$$

where $\alpha = 0.95$ by default (emphasizing area over distance).

**Configuration:**
```python
INPUT_FILENAME         = "sample_1_overlaid.h5ad"
OUTPUT_FILENAME        = "sample_1_scored.h5ad"
TISSUE_FOLD_CHANGE     = 1       # fold-change threshold for tissue enrichment
MIN_BG_INTENSITY_FRAC  = 0.1     # BG masking threshold (fraction of tissue max)
AREA_WEIGHT            = 0.95    # α in the composite score
```

**Run:**
```bash
python src/4_delocalization_score_calculator.py
```

**Output columns added to `adata.var`:**

| Column | Description |
|---|---|
| `bg_nonzero_area` | Number of nonzero off-tissue pixels |
| `bg_dist_mean` | Mean distance from nonzero BG pixels to tissue |
| `norm_bg_nonzero_area` | Min-max normalized area |
| `norm_bg_dist_mean` | Min-max normalized mean distance |
| `delocalization_score` | Final composite score ∈ [0, 1] |

---

### Step 5 — Visualization

**Script:** `5_visualization_tools.py`

> ⚠️ **Requires a display** — launches interactive matplotlib windows.

Two visualization modes:

1. **Heatmap view** — spatial intensity map for any m/z feature with tissue boundary contours overlaid in red
2. **Distance map** — arrows connecting nonzero background pixels to their nearest tissue pixel, showing the spatial direction and magnitude of delocalization

**Configuration:**
```python
INPUT_FILENAME = "sample_1_scored.h5ad"
TARGET_MZ      = None      # None → auto-select the most delocalized m/z
N_ARROWS       = 100       # max arrows on the distance map
```

**Run:**
```bash
python src/5_visualization_tools.py
```

---

## Configuration Reference

### Key Parameters Across Scripts

| Parameter | Script | Default | Description |
|---|---|---|---|
| `mz_decimals` | 1 | 6 | Decimal places for m/z rounding during alignment |
| `batch_size` | 1 | 5000 | Spectra processed per batch (memory control) |
| `BASELINE_LAMBDA` | 2 | 1e5 | ALS smoothness penalty |
| `BASELINE_P` | 2 | 0.01 | ALS asymmetry weight |
| `SG_WINDOW_LENGTH` | 2 | 9 | Savitzky–Golay window |
| `SG_POLYORDER` | 2 | 3 | Savitzky–Golay polynomial order |
| `TOP_MZ_COUNT` | 2 | 1000 | Number of peaks retained |
| `MIN_PEAK_DISTANCE_DA` | 2 | 0.024 | Minimum m/z separation between peaks (Da) |
| `REF_MZ` | 2 | `[496.339, ...]` | Lock-mass reference m/z values |
| `PPM_TOLERANCE` | 2 | 30 | Lock-mass matching window (ppm) |
| `TARGET_MZ` | 3 | 800.55 | m/z for heatmap generation during overlay |
| `TISSUE_FOLD_CHANGE` | 4 | 1 | Minimum tissue/BG fold-change for feature retention |
| `MIN_BG_INTENSITY_FRAC` | 4 | 0.1 | BG intensity threshold (fraction of tissue max) |
| `AREA_WEIGHT` | 4 | 0.95 | α weight for area in composite score |

---

## Mathematical Framework

### Composite Delocalization Score

For each of $M$ analyte features, two spatial metrics are computed:

- **$a_i$** — the number of nonzero-intensity off-tissue pixels (background area)
- **$d_i$** — the mean Euclidean distance from nonzero off-tissue pixels to the nearest tissue pixel

Both metrics are min-max normalized to $[0,1]$:

$$
\tilde{a}_i = \frac{a_i - \min_j a_j}{\max_j a_j - \min_j a_j}, \quad \tilde{d}_i = \frac{d_i - \min_j d_j}{\max_j d_j - \min_j d_j}
$$

The final delocalization score is:

$$
s_i = \alpha \cdot \tilde{a}_i + (1 - \alpha) \cdot \tilde{d}_i, \quad s_i \in [0, 1]
$$

A higher $s_i$ indicates greater spatial delocalization. The default $\alpha = 0.95$ reflects the finding that off-tissue area is a stronger indicator of delocalization than mean distance, as validated through visual inspection across multiple tissue samples.

---

## Data Availability

The MALDI-MSI datasets used in the associated publication (murine brain sections from 3xTg-AD and control mice) are available via the link provided in the [GitHub repository](https://github.com/CIL-UTK/Delocalization_Score). Data were acquired on a Waters Synapt G2-Si mass spectrometer in positive ion mode (m/z 50–2000) at 60 μm pixel size using DHB matrix (40 mg/mL in 70% methanol).

---

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit your changes (`git commit -m "Add your feature"`)
4. Push to the branch (`git push origin feature/your-feature`)
5. Open a Pull Request

Please ensure code follows the existing style (PEP 8, type hints where practical) and include docstrings for new functions.

---

## Citation

If you use this software in your research, please cite:

```bibtex
@article{jarrahi2026delocalization,
  title     = {Mathematical Framework for Quantifying Delocalization in {MALDI-MSI} via a Composite Scoring Approach},
  author    = {Jarrahi, Amin and Jones, Allison and Tang, Weisheng and Qi, Hairong and Crouch, Anna Colleen},
  journal   = {ACS Measurement Science Au},
  volume    = {6},
  number    = {1},
  pages     = {134--149},
  year      = {2026},
  doi       = {10.1021/acsmeasuresciau.5c00148},
  publisher = {American Chemical Society}
}
```

---

## License

This project is licensed under the [MIT License](LICENSE).

---

## Contact

**Crouch Imaging Lab**  
Department of Biomedical Engineering  
Tickle College of Engineering, University of Tennessee, Knoxville  

- 🌐 [crouchimaging.utk.edu](https://crouchimaging.utk.edu)  
- 📧 acrouch5@utk.edu [Dr. Colleen Crouch]
- 📧 alexajarrahi@gmail.com [Alex Jarrahi]  
- 📷 Instagram: [@crouch.imaging.lab](https://instagram.com/crouch.imaging.lab)  
- 💻 GitHub: [CIL-UTK](https://github.com/CIL-UTK)
