# HPA Downloader – Usage Protocol

## 1. Objective
To extract structured immunohistochemistry (IHC) cancer image datasets from the Human Protein Atlas (HPA) using a reproducible and selective workflow.

## 2. Requirements
- Windows 10+ (for executable)
- Internet connection
- Optional: Python 3.9+ (for development mode)

## 3. Procedure

### Step 1 — Launch
Open `ImageDownloaderHPA.exe` (or run `python ImageDownloaderHPA.py`).

### Step 2 — Provide URL
Paste a valid HPA cancer/IHC page URL. Example:
- `https://www.proteinatlas.org/ENSG.../cancer/...#img`

### Step 3 — Preview
Click **Preview**.  
The tool will detect:
- Gene name
- Antibody IDs
- Cancer subtypes
- Number of images available per group

### Step 4 — Select Targets
Select:
- Entire antibodies, and/or
- Specific cancer subtypes per antibody

If nothing is selected, the tool downloads **everything** detected in the preview.

### Step 5 — Choose Export Format
Select output image format:
- `.jpg`
- `.tif`

### Step 6 — Download
Click **Download** and wait until completion.

## 4. Output and Validation
Outputs are stored under:
`~/Downloads/HPA Images/GENE_NAME/`

Validation checklist:
- Confirm downloaded counts match preview counts
- Confirm CSV summaries exist
- Spot-check file naming and folder structure

## 5. Reproducibility Notes
- Folder hierarchy preserves provenance (gene, antibody, cancer subtype)
- CSV metadata enables downstream analysis and dataset auditing

## 6. Data Use and Ethics
Users must comply with the Human Protein Atlas data usage policy and licensing terms.