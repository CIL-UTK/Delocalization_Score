"""
MSI H&E Overlay & Tissue Annotation Tool
─────────────────────────────────────────
Interactive tool for overlaying MSI heatmaps onto H&E images,
drawing polygon ROIs, and saving tissue masks into AnnData objects.

Workflow:
  1. Load processed h5ad and generate an m/z heatmap
  2. Load the matching H&E image
  3. Interactively align H&E to MSI (rotation, scale, translation, alpha)
  4. Draw a polygon ROI on the blended overlay
  5. Save the AnnData with a binary "tissue" column to overlaid_h5ad/

Folder structure (relative to the working directory):
  <working_dir>/
  ├── processed_h5ad/             ← input h5ad files
  ├── H&E_images/                 ← input H&E image files
  └── overlaid_h5ad/              ← output h5ad files (created automatically)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scanpy as sc
from matplotlib.widgets import PolygonSelector, Slider, CheckButtons
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap
from skimage.transform import resize, rotate

matplotlib.use("Qt5Agg")  # Required for VSCode interactivity

# ═══════════════════════════════════════════════════════════════════════════════
# RESOLVE PATHS RELATIVE TO The WORKING DIRECTORY
# ═══════════════════════════════════════════════════════════════════════════════
WORKING_DIR = "/home/Delocalization_Score_V2"
PROCESSED_DIR = os.path.join(WORKING_DIR, "processed_h5ad")
HE_IMAGE_DIR = os.path.join(WORKING_DIR, "H&E_images")
OVERLAID_DIR = os.path.join(WORKING_DIR, "overlaid_h5ad")

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

# File names (inside the folders above)
INPUT_FILENAME = "sample_1_peaks.h5ad"
HE_IMAGE_FILENAME = "sample_1_HE.jpg"
OUTPUT_FILENAME = "sample_1_overlaid.h5ad"

# Full paths (assembled automatically)
INPUT_FILE = os.path.join(PROCESSED_DIR, INPUT_FILENAME)
HE_IMAGE_FILE = os.path.join(HE_IMAGE_DIR, HE_IMAGE_FILENAME)
OUTPUT_FILE = os.path.join(OVERLAID_DIR, OUTPUT_FILENAME)

# Target m/z for the heatmap visualisation
TARGET_MZ = 800.55
MZ_TOLERANCE = 0.1

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

def sort_adata_by_mz(adata, mz_col="mzs"):
    """Sort AnnData variables by numeric m/z value."""
    adata.var[mz_col] = pd.to_numeric(adata.var_names, errors="coerce")
    n_nan = adata.var[mz_col].isna().sum()
    if n_nan > 0:
        print(f"  ⚠ {n_nan} var_names could not be converted to float.")
    return adata[:, adata.var[mz_col].sort_values().index]


def get_mz_heatmap(adata, target_mz, tol=0.1):
    """
    Extract a 2-D intensity image for the closest m/z within tolerance.

    Returns
    -------
    image : ndarray  – (height, width) intensity array
    matched_mz : float
    """
    mz_axis = adata.var_names.astype(float).values
    mz_diff = np.abs(mz_axis - target_mz)
    if mz_diff.min() > tol:
        raise ValueError(f"No m/z found within ±{tol} of {target_mz}")

    idx = np.argmin(mz_diff)
    matched_mz = mz_axis[idx]

    intensities = (
        adata.X[:, idx].toarray().ravel()
        if hasattr(adata.X, "toarray")
        else adata.X[:, idx]
    )

    x = adata.obs["x"].values.astype(int)
    y = adata.obs["y"].values.astype(int)

    image = np.zeros((y.max() + 1, x.max() + 1))
    image[y, x] = intensities
    return image, matched_mz


def load_he_image(path):
    """
    Read an H&E image file and normalise to [0, 1] RGB.

    Parameters
    ----------
    path : str – path to the image file (JPG, PNG, etc.)

    Returns
    -------
    ndarray – (H, W, 3) float array in [0, 1]
    """
    img = mpimg.imread(path)
    if img.ndim == 2:
        img = np.stack([img] * 3, axis=-1)
    return np.clip(img / 255.0 if img.max() > 1 else img, 0, 1)


def normalise_image(image):
    """Min-max normalise a 2-D array to [0, 1]."""
    lo, hi = image.min(), image.max()
    return (image - lo) / (hi - lo) if hi > lo else np.zeros_like(image)


# ═══════════════════════════════════════════════════════════════════════════════
# INTERACTIVE OVERLAY ALIGNMENT
# ═══════════════════════════════════════════════════════════════════════════════

def interactive_overlay(heatmap_rgb, he_image):
    """
    Open a matplotlib window with sliders for rotation, scale, translation,
    and alpha blending between the H&E and MSI heatmap.

    Returns
    -------
    img_ax : AxesImage – handle to the displayed overlay (used later for ROI)
    """
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.45)

    img_ax = ax.imshow(np.zeros_like(heatmap_rgb), cmap="gray")
    ax.set_title("Adjust sliders to align H&E → MSI, then close window")

    # ── slider factory ──
    def _slider(rect, label, vmin, vmax, vinit):
        return Slider(plt.axes(rect), label, vmin, vmax, valinit=vinit)

    s_angle   = _slider([0.25, 0.37, 0.65, 0.03], "Angle",   -180, 180, 0)
    s_scale_x = _slider([0.25, 0.32, 0.65, 0.03], "Scale X",  0.1, 3.0, 1.0)
    s_scale_y = _slider([0.25, 0.27, 0.65, 0.03], "Scale Y",  0.1, 3.0, 1.0)
    s_tx      = _slider([0.25, 0.22, 0.65, 0.03], "Move X",  -200, 200, 0)
    s_ty      = _slider([0.25, 0.17, 0.65, 0.03], "Move Y",  -200, 200, 0)
    s_alpha   = _slider([0.25, 0.12, 0.65, 0.03], "Alpha",    0.0, 1.0, 0.5)

    cb_ax = plt.axes([0.025, 0.12, 0.15, 0.1])
    cb_uniform = CheckButtons(cb_ax, ["Uniform Scale"], [False])

    def _update(_=None):
        angle = s_angle.val
        sx, sy = s_scale_x.val, s_scale_y.val
        tx, ty = s_tx.val, s_ty.val
        alpha = s_alpha.val

        if cb_uniform.get_status()[0]:
            sy = sx

        # Transform H&E
        rotated = rotate(he_image, angle, resize=False, preserve_range=True)
        target_h = int(heatmap_rgb.shape[0] * sy)
        target_w = int(heatmap_rgb.shape[1] * sx)
        resized = resize(rotated, (target_h, target_w), preserve_range=True)

        # Place transformed H&E onto a canvas matching the heatmap size
        h, w = heatmap_rgb.shape[:2]
        canvas = np.zeros_like(heatmap_rgb)

        y_start = int((h - target_h) / 2 + ty)
        x_start = int((w - target_w) / 2 + tx)
        y_end = min(y_start + target_h, h)
        x_end = min(x_start + target_w, w)

        src_y = max(-y_start, 0)
        src_x = max(-x_start, 0)

        if y_start < h and x_start < w:
            dy = y_end - max(y_start, 0)
            dx = x_end - max(x_start, 0)
            canvas[max(y_start, 0):y_end,
                   max(x_start, 0):x_end] = resized[src_y:src_y + dy,
                                                     src_x:src_x + dx]

        overlay = np.clip((1 - alpha) * canvas + alpha * heatmap_rgb, 0, 1)
        img_ax.set_data(overlay)
        fig.canvas.draw_idle()

    for s in (s_angle, s_scale_x, s_scale_y, s_tx, s_ty, s_alpha):
        s.on_changed(_update)
    cb_uniform.on_clicked(lambda _: _update())

    _update()
    plt.show(block=True)
    return img_ax


# ═══════════════════════════════════════════════════════════════════════════════
# POLYGON ROI ANNOTATOR
# ═══════════════════════════════════════════════════════════════════════════════

class PolygonAnnotator:
    """
    Display an image and let the user draw a single closed polygon.
    Double-click (or press Enter) to finish. The window closes automatically.

    Attributes
    ----------
    coords : list of (x, y) or None
    """

    def __init__(self, image, cmap="hot"):
        self.image = image
        self.coords = None

        self.fig, self.ax = plt.subplots()
        self.ax.imshow(image, cmap=cmap)
        self.ax.set_title("Click to draw polygon · Double-click to finish")

        self.selector = PolygonSelector(
            self.ax, self._on_select, useblit=True,
            props=dict(color="cyan", linewidth=2),
            handle_props=dict(marker="o", markersize=5, mec="cyan", mfc="cyan"),
        )
        plt.show(block=True)

    def _on_select(self, verts):
        self.coords = verts
        print(f"  Polygon drawn with {len(verts)} vertices.")
        plt.close(self.fig)

    def get_mask(self):
        """Return a boolean mask with the same (H, W) as the displayed image."""
        if self.coords is None:
            print("  No polygon was drawn.")
            return None
        h, w = self.image.shape[:2]
        yy, xx = np.mgrid[:h, :w]
        points = np.column_stack((xx.ravel(), yy.ravel()))
        return Path(self.coords).contains_points(points).reshape(h, w)


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    # ── Load data ──
    adata = sc.read_h5ad(INPUT_FILE)
    adata = sort_adata_by_mz(adata)
    print(f"Loaded AnnData : {INPUT_FILE}  {adata.shape}")

    # ── Build MSI heatmap ──
    image, matched_mz = get_mz_heatmap(adata, TARGET_MZ, tol=MZ_TOLERANCE)
    print(f"Heatmap m/z    : {TARGET_MZ} → matched {matched_mz}")

    heatmap_norm = normalise_image(image)
    heatmap_rgb = CUSTOM_CMAP(heatmap_norm)[..., :3]

    # ── Load H&E ──
    he_image = load_he_image(HE_IMAGE_FILE)
    print(f"H&E image      : {HE_IMAGE_FILE}")

    # ── Interactive alignment ──
    img_ax = interactive_overlay(heatmap_rgb, he_image)

    # ── Polygon ROI on the final overlay ──
    final_overlay = img_ax.get_array()
    annotator = PolygonAnnotator(final_overlay, cmap=None)
    mask = annotator.get_mask()

    if mask is None:
        print("No tissue mask created — aborting save.")
        return

    # Preview
    plt.imshow(np.where(mask, image, 0), cmap=CUSTOM_CMAP)
    plt.title(f"Masked m/z = {matched_mz:.4f}")
    plt.show()

    # ── Write tissue mask into AnnData ──
    x = adata.obs["x"].astype(int).values
    y = adata.obs["y"].astype(int).values
    adata.obs["tissue"] = mask[y, x].astype(int)
    print(f"Tissue mask counts:\n{adata.obs['tissue'].value_counts().to_string()}")

    # ── Save ──
    os.makedirs(OVERLAID_DIR, exist_ok=True)
    adata.write(OUTPUT_FILE)
    print(f"Saved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
