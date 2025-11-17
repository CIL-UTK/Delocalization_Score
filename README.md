MALDI-MSI Data Processing Pipeline.\
\
This repository contains a complete MALDI-MSI (Matrix-Assisted Laser Desorption/Ionization Mass Spectrometry Imaging) data processing pipeline. The pipeline is designed to process raw MSI data (.imzML + .ibd) acquired from Waters Synapt Mass Spectrometers, extract m/z metrics, organize the data in AnnData objects, perform signal processing, remove matrix related analytes, normalize intensities, overlay tissue optical images, detect tissue margins, group m/z features, and compute a delocalization score for downstream analysis.\
The pipeline is implemented in Python, using pandas, numpy, anndata, scipy, plotly, matplotlib, skimage, scanpy, tqdm, os,and seaborn tools.\
\
Overview\
\
This pipeline provides a robust, reproducible workflow for multi-sample MALDI-MSI studies, allowing users to:
1. Load raw MSI datasets (.imzML + .ibd).
2. Extract intensity metrics per m/z and store in AnnData format.
3. Remove matrix signals.
4. Normalize tissue and background intensities, including log-transformed values.
5. Overlay optical images for spatial context manually.
6. Annotate tissue margins manually.
7. Group m/z features based on m/z distance.
8. Compute delocalization scores to prioritize relevant features.
