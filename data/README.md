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
├── imzml_data/              ← Raw imzML + ibd files (22 samples)
│   ├── 20230828_3xTgControl_youngfemale_5_pos_1_1.imzml
│   ├── 20230828_3xTgControl_youngfemale_5_pos_1_1.ibd
│   ├── 20230914_3xTgAD_youngfemale_1_pos_3_1.imzml
│   ├── 20230914_3xTgAD_youngfemale_1_pos_3_1.ibd
│   └── ...
│
├── H&E_images/              ← H&E stained tissue images
│   ├── yfm1_c.jpg
│   ├── ofm1_c.jpg
│   ├── ofm1_ad.jpg
│   ├── ofm2_ad.jpg
│   └── ...
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

## File Naming Conventions

### imzML Files

The imzML filenames follow this pattern:

```
YYYYMMDD_<model><group>_<age><sex>_<slide>_pos_<animal#>_<section>
```

| Component | Example | Meaning |
|---|---|---|
| `20230828` | Date | Acquisition date |
| `3xTgControl` | Model + Group | 3xTg mouse model, Control group |
| `3xTgAD` | Model + Group | 3xTg mouse model, AD group |
| `youngfemale` | Age + Sex | Young female animal |
| `5` | Slide | Slide 5 |
| `pos` | Polarity | Positive ion mode |
| `1_1` | Animal number_Section | Animal 1, Section 1 |

### H&E Image Files

The H&E image filenames follow a shorthand pattern:

```
<age><sex><animal#>_<group>.jpg
```

| Code | Meaning | Example |
|---|---|---|
| `y` | Young | `yfm1_c` → Young female 1, Control |
| `o` | Old (Aged) | `ofm1_c` → Old/Aged female 1, Control |
| `fm` | Female | |
| `m` | Male | |
| `_c` | Control | `yfm1_c` → Young female 1, Control |
| `_ad` | AD (Alzheimer's Disease) | `ofm2_ad` → Old/Aged female 2, AD |

## Notes

- Only `imzml_data/` and `H&E_images/` need to be populated manually. The other subfolders are created automatically by the pipeline scripts.
- All `.h5ad`, `.imzml`, `.ibd`, and image files are excluded from version control via `.gitignore`.
- The `WORKING_DIR` variable in each script under `src/` auto-resolves to this `data/` directory — no manual path editing needed.
