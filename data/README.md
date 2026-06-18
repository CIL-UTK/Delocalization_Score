# Data Directory

This folder contains the input and intermediate data files for the Delocalization Score pipeline. **Data files are not included in this repository** due to their size — they are hosted externally.

## 📥 Download the Dataset

Download all data files from the link below and place them in the corresponding subfolders:

**[Download Dataset (SharePoint)](https://liveutk-my.sharepoint.com/:f:/g/personal/ajarrah1_vols_utk_edu/EqUysi_9h2BOnOSajIEcq64BhPRXaRExNfXPcdQ12BxlRg?e=wpZNMZ)**

Password:
data@CIL-UTK

## 📁 Folder Structure

After downloading, organize the files as follows:

```
data/
├── imzml_data/              ← Raw imzML + ibd files
│   ├── sample_1.imzml
│   └── sample_1.ibd
│
├── H&E_images/              ← H&E stained tissue images
│   └── sample_1_HE.jpg
│
├── raw_h5ad/                ← Auto-created by Step 1
│   └── (generated .h5ad files)
│
├── processed_h5ad/          ← Auto-created by Step 2
│   └── (generated .h5ad files)
│
├── overlaid_h5ad/           ← Auto-created by Step 3
│   └── (generated .h5ad files)
│
└── delocalization_score_result_h5ad/  ← Auto-created by Step 4
    └── (generated .h5ad files)
```

## Notes

- Only `imzml_data/` and `H&E_images/` need to be populated manually. The other subfolders are created automatically by the pipeline scripts.
- All `.h5ad`, `.imzml`, `.ibd`, and image files are excluded from version control via `.gitignore`.
- Update the `WORKING_DIR` variable in each script under `src/` to point to this `data/` directory.
