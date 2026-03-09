# HPA Cancer/IHC Image Downloader

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18926273.svg)](https://doi.org/10.5281/zenodo.18926273)

A lightweight graphical application for structured and selective downloading of immunohistochemistry (IHC) cancer images from the Human Protein Atlas (HPA).

![HPA Downloader Concept](assets/screenshots/HPADownloader_concept.png)


This tool supports reproducible digital pathology workflows by organizing downloads by **Gene → Antibody → Cancer subtype** and exporting **metadata CSV summaries**.

---

## Features

- URL-based automatic parsing of HPA cancer/IHC pages  
- Preview tree: **Gene → Antibody → Cancer subtype → image count**  
- Selective download (choose specific antibodies and/or cancer subtypes)  
- Export formats: **.tif** or **.jpg** 
- Configurable output directory
- Automatic metadata extraction (patient/staining fields when available)  
- Retry-safe downloads  
- Standalone Windows executable available in **Releases**

---

## Quick Start (Windows Executable)

1. Download the latest `ImageDownloaderHPA.exe` from **Releases**
2. Run the executable
3. Paste a valid HPA cancer/IHC URL
4. Click **Preview**
5. Select what you need
6. Click **Download**

> Output is saved to: `~/Downloads/HPA Images/`

---

## Usage (Preview → Select → Download)

1. Paste a Human Protein Atlas cancer/IHC URL, for example:  
   `https://www.proteinatlas.org/ENSG00000077782-FGFR1/cancer/lung+cancer#img`

2. Click **Preview** to build the tree view and selection panel.

3. Choose:
   - Full antibodies
   - Specific cancer subtypes within each antibody

4. Choose output format: `.tif` (default) or `.jpg`.

5. Click **Download**.

![HPA Downloader Concept](assets/screenshots/gui_main.png)


---

## Installation (Python)

Requires Python ≥3.9

Clone the repository:

``` bash
git clone https://github.com/Juaco2r/HPADownloader.git
cd HPADownloader
```

Install dependencies:

``` bash
pip install -r requirements.txt
```

Run the application:

``` bash
python src/ImageDownloaderHPA.py
```

---

## Dependencies

The application requires only three external Python packages:

- requests
- beautifulsoup4
- lxml

All dependencies can be installed automatically using:

```bash
pip install -r requirements.txt

---

## Example Output

Example dataset generated from:

https://www.proteinatlas.org/ENSG00000077782-FGFR1/cancer/lung+cancer

The software identified:

- 2 antibodies
- 4 different cancer subtypes (2 for CAP033614 and 3 for HPA056402)
- 61 total images

Images will be organized automatically and a metadata CSV file will be generated.

---

## Output Folder Structure

```text
HPA Images/
└── GENE_NAME/
    ├── HPAxxxxxx/
    │   ├── CancerSubtypeA/
    │   │   ├── ID_patient_1.jpg
    │   │   ├── ID_patient_2.jpg
    │   │   └── CancerSubtypeA_summary.csv
    │   └── HPAxxxxxx_summary_all_images.csv
    └── GENE_all_antibodies_summary.csv

```
---

## Citation

If you use this software in research, please cite the archived version:

Rodriguez-Rojas J. (2026).  
**HPA Cancer/IHC Image Downloader (v1.2)**.  
Zenodo.  
https://doi.org/10.5281/zenodo.18926273

Source code: https://github.com/Juaco2r/HPADownloader