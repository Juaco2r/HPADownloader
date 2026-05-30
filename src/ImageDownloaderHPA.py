"""
HPA IHC Image Downloader

GUI-based tool for previewing, selecting, and downloading immunohistochemistry images from cancer and normal tissue pages from the Human Protein Atlas (HPA).
The software organizes downloaded images by gene, antibody, HPA cancer/tissue category,
and HPA diagnostic subtype when available, and exports metadata summaries in CSV format.
"""

import re
import csv
import json
import time
import threading
import queue
import requests
import sys
import argparse
import hashlib
import platform
from bs4 import BeautifulSoup
from html import unescape
from pathlib import Path
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import webbrowser
import xml.etree.ElementTree as ET
from urllib.parse import urlparse, urljoin, unquote
from datetime import datetime

__version__ = "1.3"
__author__ = "José Rodríguez-Rojas"
__github_url__ = "https://github.com/Juaco2r/HPADownloader"
__doi__ = "https://doi.org/10.5281/zenodo.20465365"

HEADERS = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome Safari/537.36"}

# Reuse HTTP connections across HTML, XML and image requests.
# This keeps the downloader robust while reducing network overhead during bulk jobs.
SESSION = requests.Session()
SESSION.headers.update(HEADERS)

# Cache per-gene XML metadata so overview pages do not fetch the same ENSG XML
# repeatedly for every cancer/tissue category.
_XML_METADATA_CACHE = {}

BASE_IMG_HOST = "https://images.proteinatlas.org"
USER_DOWNLOAD_DIR = Path.home() / "Downloads"
ROOT_DIR = USER_DOWNLOAD_DIR / "HPA Images"


# ---------------------------------------------------------------------
# Core parsing and download utilities
# ---------------------------------------------------------------------

def resource_path(relative_path):
    """Get absolute path to resource (works for PyInstaller)."""
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = Path(__file__).resolve().parent.parent
    return Path(base_path) / relative_path


def download_html(url):
    """Download and return the HTML content of an HPA page."""
    r = SESSION.get(url, timeout=30)
    r.raise_for_status()
    return r.text


def is_valid_hpa_url(url):
    """Return True if the URL appears to belong to the Human Protein Atlas."""
    return bool(url) and "proteinatlas.org" in url.lower()


def extract_gene_name(soup):
    """Extract the gene name from the parsed HPA HTML page."""
    gn = soup.find("div", class_="gene_name")
    if gn:
        if gn.has_attr("data-gene_name"):
            gene = gn["data-gene_name"].strip()
            if gene:
                return gene
        text_gene = gn.get_text(strip=True)
        if text_gene:
            return text_gene
    return "UnknownGene"


def extract_antibody_ids(soup):
    """Extract unique antibody identifiers from the parsed HPA page."""
    text = soup.get_text(" ", strip=True)
    ids = sorted(set(re.findall(r"\b(HPA\d{6}|CAB\d{6})\b", text)))
    return ids


def clean_html_title(html_title):
    """Convert HTML-formatted title text into plain text with line breaks."""
    txt = unescape(html_title)
    txt = re.sub(r"</?b>", "", txt)
    txt = re.sub(r"<br\s*/?>", "\n", txt)
    return txt.strip()


def infer_hpa_section_from_url(url):
    """Infer the HPA section from the URL so output labels are not cancer-only."""
    u = (url or "").lower()
    if "/cancer" in u:
        return "cancer"
    if "/tissue" in u:
        return "tissue"
    return "ihc"


def is_image_href(href):
    """Return True for links that look like HPA image links."""
    if not href:
        return False
    href_l = href.lower()
    return (
        "images.proteinatlas.org" in href_l
        or href_l.startswith("/images/")
        or href_l.startswith("/images_static/")
        or re.search(r"/\d+/[^/]+\.(jpg|tif|tiff)$", href_l) is not None
    )


def safe_folder_name(value, fallback="Unknown"):
    """Create a safe folder name while preserving readable labels."""
    value = (value or "").strip()
    if not value:
        value = fallback
    value = value.replace("NOS", "").strip()
    value = re.sub(r"[<>:\"/\\|?*]+", "_", value)
    value = re.sub(r"\s+", " ", value).strip(" ._")
    return value or fallback


def clean_text_value(value, max_len=1200):
    """Return a compact one-line string safe for CSV/Excel cells."""
    if value is None:
        return ""
    value = str(value)
    value = unescape(value)
    value = re.sub(r"[\r\n\t]+", " ", value)
    value = re.sub(r"\s+", " ", value).strip()
    if max_len and len(value) > max_len:
        value = value[:max_len - 3].rstrip() + "..."
    return value


def is_plausible_annotation_label(value):
    """Return True for HPA cell-type/compartment labels, not full page sections."""
    value = clean_text_value(value, max_len=500)
    if not value or len(value) > 90:
        return False
    bad = [
        "protein expression", "rna expression", "ntpm", "gtex", "fantom",
        "consensus", "sample id", "average", "max subtype", "show all",
        "show less", "antibody ", "information about", "scaled tags"
    ]
    low = value.lower()
    if any(b in low for b in bad):
        return False
    if re.search(r"\b(?:GTEX|FF:)\b", value):
        return False
    # Avoid long numeric expression rows.
    if len(re.findall(r"\d", value)) > 8:
        return False
    return True


def compact_annotation_summary(annotations, max_len=1200):
    """Build a short annotation summary from parsed HPA annotation dictionaries."""
    parts = []
    for ann in annotations or []:
        name = clean_text_value(ann.get("AnnotationType") or "Annotation", max_len=120)
        if not is_plausible_annotation_label(name):
            continue
        vals = []
        for out_key, in_key in [
            ("Staining", "Staining"),
            ("Intensity", "Intensity"),
            ("Quantity", "Quantity"),
            ("Location", "Location"),
        ]:
            val = clean_text_value(ann.get(in_key), max_len=180)
            if val:
                vals.append(f"{out_key}={val}")
        if vals:
            parts.append(f"{name}: " + "; ".join(vals))
    return clean_text_value(" | ".join(parts), max_len=max_len)


def sanitize_csv_row(row):
    """Clean all row values before writing to CSV."""
    cleaned = {}
    for k, v in row.items():
        max_len = 1200 if k == "AnnotationSummary" else 400
        cleaned[k] = clean_text_value(v, max_len=max_len)
    return cleaned


def _split_name_code(line):
    """Split strings such as 'Colon (T-67000)' into label and code."""
    m = re.match(r"^(.*?)\s*\(([^()]+)\)\s*$", line.strip())
    if not m:
        return line.strip(), None
    return m.group(1).strip(), m.group(2).strip()


def parse_metadata_from_title(title_text):
    """
    Parse metadata embedded in HPA image titles.

    This supports both cancer pages and normal tissue pages. Cancer pages
    usually contain one annotation set, while tissue pages can contain
    repeated annotation blocks by cell type, e.g. Endothelial cells,
    Glandular cells, Peripheral nerve/ganglion.
    """
    meta = {
        "Gender": None,
        "Age": None,
        "PatientID": None,
        "LocationsRaw": [],
        "Tissue": None,
        "TissueCode": None,
        "Diagnosis": None,
        "DiagnosisCode": None,
        "CancerType": None,
        "CancerCode": None,
        "AntibodyStaining": None,
        "Intensity": None,
        "Quantity": None,
        "Location": None,
        "AnnotationType": None,
        "AnnotationSummary": "",
        "AllAnnotations": [],
    }

    lines = [l.strip() for l in title_text.split("\n") if l.strip()]

    # Typical first lines can be: Tissue / Antibody / Female, age 55 / Colon (T-67000) ...
    for line in lines:
        m = re.match(r"^(Male|Female)\s*,\s*age\s*(\d+)", line, re.IGNORECASE)
        if m:
            meta["Gender"] = m.group(1).capitalize()
            meta["Age"] = m.group(2)
            continue

        m = re.match(r"Patient id:\s*(\d+)", line, re.IGNORECASE)
        if m:
            meta["PatientID"] = m.group(1)
            continue

        label, code = _split_name_code(line)
        if code:
            if code.startswith("T-"):
                meta["Tissue"] = label
                meta["TissueCode"] = code
                meta["LocationsRaw"].append(line)
            elif code.startswith("M-"):
                meta["Diagnosis"] = label
                meta["DiagnosisCode"] = code
                # In HPA cancer pages this is the cancer type; in normal tissue
                # pages this can be 'Normal tissue'.
                if "normal tissue" not in label.lower():
                    meta["CancerType"] = label
                    meta["CancerCode"] = code
                else:
                    meta["LocationsRaw"].append(line)

    # Parse repeated annotation blocks. A block starts with a line without ':'
    # followed by Staining/Intensity/Quantity/Location lines.
    annotations = []
    current = None

    ignored_plain = set()
    for line in lines:
        if re.fullmatch(r"(HPA|CAB)\d{6}", line):
            ignored_plain.add(line)
        if re.match(r"^(Male|Female)\s*,\s*age\s*\d+", line, re.IGNORECASE):
            ignored_plain.add(line)
        if re.match(r"Patient id:", line, re.IGNORECASE):
            ignored_plain.add(line)
        if _split_name_code(line)[1]:
            ignored_plain.add(line)

    for line in lines:
        if line in ignored_plain:
            continue

        m = re.match(r"^(Antibody staining|Staining|Intensity|Quantity|Location):\s*(.*)$", line, re.IGNORECASE)
        if m:
            key = m.group(1).strip().lower()
            value = m.group(2).strip()
            if current is None:
                current = {"AnnotationType": None, "Staining": None, "Intensity": None, "Quantity": None, "Location": None}
            if key in ("antibody staining", "staining"):
                current["Staining"] = value
            elif key == "intensity":
                current["Intensity"] = value
            elif key == "quantity":
                current["Quantity"] = value
            elif key == "location":
                current["Location"] = value
            continue

        # A plain non-metadata line begins a new annotation type only if it
        # looks like a real cell type/compartment label. This prevents full
        # page text (RNA/GTEx/FANTOM tables) from becoming CSV metadata.
        if not is_plausible_annotation_label(line):
            continue
        if current and any(v for k, v in current.items() if k != "AnnotationType"):
            annotations.append(current)
        current = {"AnnotationType": line, "Staining": None, "Intensity": None, "Quantity": None, "Location": None}

    if current and any(v for k, v in current.items() if k != "AnnotationType"):
        annotations.append(current)

    meta["AllAnnotations"] = annotations
    if annotations:
        first = annotations[0]
        meta["AnnotationType"] = first.get("AnnotationType")
        meta["AntibodyStaining"] = first.get("Staining")
        meta["Intensity"] = first.get("Intensity")
        meta["Quantity"] = first.get("Quantity")
        meta["Location"] = first.get("Location")
        meta["AnnotationSummary"] = compact_annotation_summary(annotations)

    return meta


def normalize_category_folder(category):
    """Generic folder normalizer for cancer subtype, tissue, or other HPA groups."""
    return safe_folder_name(category, fallback="Unknown")




def valid_hpa_diagnostic_subtype(value, section_category=""):
    """
    Return a safe diagnostic/histology subtype only when it was explicitly
    parsed from HPA metadata/table labels.

    We intentionally do NOT invent labels such as Unknown subtype. If HPA does
    not expose a subtype, the image is stored directly in the cancer category
    folder and the DiagnosticCategory fields stay empty.
    """
    value = clean_text_value(value, max_len=180)
    section_category = clean_text_value(section_category, max_len=180)
    if not value:
        return ""
    low = value.lower().strip()
    if low in {"unknown", "none", "n/a", "na", "not available", "not detected"}:
        return ""
    if section_category and low == section_category.lower().strip():
        return ""
    if low in {"cancer", "tissue", "ihc"}:
        return ""
    if len(value) > 90 or any(b in low for b in ["protein expression", "rna expression", "gtex", "fantom", "sample id", "show all"]):
        return ""
    return value

def normalize_cancer_folder(cancer_type):
    """Backward-compatible wrapper for older cancer-only naming."""
    return normalize_category_folder(cancer_type)


def normalize_antibody_numeric(antibody_id):
    """
    Extract normalized numeric antibody ID from strings like HPA000007 or CAB014263.
    HPA000007 -> '7'
    CAB014263 -> '14263'
    """
    m = re.fullmatch(r"(HPA|CAB)(\d+)", antibody_id.strip())
    if not m:
        return None
    return str(int(m.group(2)))


def extract_antibody_numeric_from_href(href):
    """
    Extract antibody numeric folder from an HPA image URL.

    Examples
    --------
    https://images.proteinatlas.org/7/1339_B_1_1.tif -> '7'
    https://images.proteinatlas.org/14263/32472_B_1_1.jpg -> '14263'
    //images.proteinatlas.org/35305/84031_A_8_3_medium.jpg -> '35305'
    """
    if not href:
        return None
    m = re.search(r"(?:images\.proteinatlas\.org|^)/(\d+)/", href)
    if not m:
        return None
    return str(int(m.group(1)))


def normalize_image_href(href):
    """Convert relative/partial HPA image hrefs into absolute URLs."""
    if not href:
        return href
    href = href.strip()
    if href.startswith("//"):
        return "https:" + href
    if href.startswith("/images/"):
        return BASE_IMG_HOST + href[len("/images"):]
    if href.startswith("/images_static/"):
        return BASE_IMG_HOST + href[len("/images_static"):]
    return href


def collect_title_candidate(a_tag):
    """Collect image metadata text from several possible HPA attributes."""
    candidates = []
    img = a_tag.find("img") if a_tag else None
    for tag in (img, a_tag):
        if not tag:
            continue
        for attr in ["title", "alt", "data-title", "data-original-title", "aria-label"]:
            value = tag.get(attr)
            if value:
                candidates.append(str(value))
    # Some HPA layouts place useful title text in a nearby wrapper, but very
    # large parent attributes often contain entire page sections/RNA tables.
    parent = a_tag.parent if a_tag else None
    if parent is not None:
        for attr in ["title", "data-title", "data-original-title", "aria-label"]:
            value = parent.get(attr)
            if value and len(str(value)) < 2500:
                candidates.append(str(value))
    candidates = [c for c in candidates if c and len(c) < 2500]
    best = max(candidates, key=len) if candidates else ""
    return clean_html_title(best) if best else ""


def image_link_for_format(href, img_ext):
    """Return the requested HPA image URL while avoiding invalid _medium.tif links."""
    link = href
    if img_ext == ".tif":
        # Older HPA IHC links often expose .jpg/.tif alternatives.  Thumbnail links
        # such as *_medium.jpg should not become *_medium.tif; instead request .tif
        # from the same base image when possible.
        link = re.sub(r"_medium\.jpg($|\?)", r".tif\1", link, flags=re.IGNORECASE)
        link = re.sub(r"\.jpg($|\?)", r".tif\1", link, flags=re.IGNORECASE)
    return link


def infer_tissue_name_from_url(url):
    """Infer tissue label from URLs like /tissue/colon."""
    m = re.search(r"/tissue/([^/?#]+)", url or "", re.IGNORECASE)
    if not m:
        return None
    return m.group(1).replace("_", " ").replace("-", " ").title()




def local_xml_tag(tag):
    """Return XML tag name without namespace."""
    return str(tag).split("}", 1)[-1].lower()


def clean_image_key(value):
    """Normalize an HPA image URL/path to a stable key such as 84031_A_8_3."""
    if not value:
        return None
    value = str(value).strip()
    value = value.split("?", 1)[0].split("#", 1)[0]
    name = value.rsplit("/", 1)[-1]
    name = re.sub(r"_medium(?=\.(jpg|jpeg|tif|tiff)$)", "", name, flags=re.IGNORECASE)
    name = re.sub(r"\.(jpg|jpeg|tif|tiff)$", "", name, flags=re.IGNORECASE)
    return name.lower() if name else None


def extract_image_links_from_text(text):
    """Extract HPA image links from plain text or XML attributes."""
    if not text:
        return []
    return re.findall(r"https?://images\.proteinatlas\.org/[^\s<>\"']+", str(text))


def extract_ensg_from_url(url):
    """Extract ENSG identifier from an HPA gene URL."""
    m = re.search(r"/(ENSG\d+)", url or "")
    return m.group(1) if m else None


def merge_missing_metadata(base_meta, extra_meta):
    """Fill empty metadata values in base_meta using extra_meta."""
    if not extra_meta:
        return base_meta
    aliases = {
        "PatientID": ["PatientID", "patientId", "patientid", "PatientId", "patient", "sampleId", "sample"],
        "Gender": ["Gender", "Sex", "sex", "gender"],
        "Age": ["Age", "age"],
        "Tissue": ["Tissue", "tissue", "organ", "tissueName"],
        "TissueCode": ["TissueCode", "tissueCode", "snomedCode", "snomed"],
        "Diagnosis": ["Diagnosis", "diagnosis", "sampleDescription", "pathology"],
        "DiagnosisCode": ["DiagnosisCode", "diagnosisCode"],
        "CancerType": ["CancerType", "cancerType", "tumorType"],
        "CancerCode": ["CancerCode", "cancerCode"],
        "Intensity": ["Intensity", "intensity"],
        "Quantity": ["Quantity", "quantity", "fraction"],
        "Location": ["Location", "location", "mainLocation"],
        "AntibodyStaining": ["AntibodyStaining", "staining", "level", "score"],
        "AnnotationType": ["AnnotationType", "cellType", "celltype"],
    }
    for target, keys in aliases.items():
        if base_meta.get(target):
            continue
        for k in keys:
            v = extra_meta.get(k)
            if v not in (None, ""):
                if target == "Gender" and str(v).lower() in ("male", "female"):
                    v = str(v).capitalize()
                base_meta[target] = str(v).strip()
                break
    return base_meta


def collect_xml_context_metadata(node, parent_map, max_up=7):
    """
    Collect metadata around an XML image node.

    HPA XML layouts have changed over releases, so this is intentionally
    defensive: it climbs through parent sample/assay nodes and collects likely
    metadata leaves such as patientId, sex/gender, age, tissue, diagnosis,
    intensity, quantity and location.
    """
    useful_tags = {
        "patientid", "patient", "sampleid", "sample", "sex", "gender", "age",
        "tissue", "tissuename", "tissuedescription", "tissuecode", "snomedcode", "snomed",
        "diagnosis", "diagnosiscode", "pathology", "cancertype", "cancercode", "tumortype",
        "staining", "level", "score", "intensity", "quantity", "fraction", "location",
        "mainlocation", "celltype", "celltypeid", "antibody", "antibodyid", "imageurl"
    }
    metadata = {}
    cur = node
    visited = []
    for _ in range(max_up):
        if cur is None:
            break
        visited.append(cur)
        cur = parent_map.get(cur)

    for ctx in visited:
        # Text visible on the context node can contain compact labels.
        ctx_text = " ".join(t.strip() for t in ctx.itertext() if t and t.strip())
        m = re.search(r"\b(Male|Female)\b", ctx_text, re.IGNORECASE)
        if m and not metadata.get("Gender"):
            metadata["Gender"] = m.group(1).capitalize()
        m = re.search(r"\bage\s*(\d{1,3})\b", ctx_text, re.IGNORECASE)
        if m and not metadata.get("Age"):
            metadata["Age"] = m.group(1)
        m = re.search(r"Patient\s*id\s*:?\s*(\d+)", ctx_text, re.IGNORECASE)
        if m and not metadata.get("PatientID"):
            metadata["PatientID"] = m.group(1)

        for el in ctx.iter():
            tag = local_xml_tag(el.tag)
            if tag not in useful_tags:
                continue
            txt = (el.text or "").strip()
            if not txt or len(txt) > 250:
                continue
            # Preserve the first useful value for each raw XML tag.
            metadata.setdefault(tag, txt)

            if tag in ("patientid", "patient", "sampleid", "sample"):
                mm = re.search(r"\d+", txt)
                if mm:
                    metadata.setdefault("PatientID", mm.group(0))
            elif tag in ("sex", "gender"):
                if txt.lower() in ("male", "female"):
                    metadata.setdefault("Gender", txt.capitalize())
            elif tag == "age":
                mm = re.search(r"\d{1,3}", txt)
                if mm:
                    metadata.setdefault("Age", mm.group(0))
            elif tag in ("tissue", "tissuename", "tissuedescription"):
                metadata.setdefault("Tissue", txt)
            elif tag in ("tissuecode", "snomedcode", "snomed"):
                metadata.setdefault("TissueCode", txt)
            elif tag in ("diagnosis", "pathology"):
                metadata.setdefault("Diagnosis", txt)
            elif tag == "diagnosiscode":
                metadata.setdefault("DiagnosisCode", txt)
            elif tag in ("cancertype", "tumortype"):
                metadata.setdefault("CancerType", txt)
            elif tag == "cancercode":
                metadata.setdefault("CancerCode", txt)
            elif tag in ("staining", "level", "score"):
                metadata.setdefault("AntibodyStaining", txt)
            elif tag == "intensity":
                metadata.setdefault("Intensity", txt)
            elif tag in ("quantity", "fraction"):
                metadata.setdefault("Quantity", txt)
            elif tag in ("location", "mainlocation"):
                metadata.setdefault("Location", txt)
            elif tag in ("celltype", "celltypeid"):
                metadata.setdefault("AnnotationType", txt)
    return metadata


def build_xml_image_metadata_map(url, log=None):
    """
    Fetch the per-gene HPA XML and build image-key -> metadata mapping.

    This is mainly needed for normal tissue pages, because the visible page may
    expose thumbnails as direct JPG links without patient metadata in the title.
    The HPA help page states that single-gene XML can be fetched by appending
    .xml to the ENSG entry URL.

    Results are cached by ENSG identifier. This avoids downloading and parsing
    the same XML repeatedly when an overview page is expanded into many child
    cancer/tissue category pages.
    """
    ensg = extract_ensg_from_url(url)
    if not ensg:
        return {}

    if ensg in _XML_METADATA_CACHE:
        return _XML_METADATA_CACHE[ensg]

    xml_url = f"https://www.proteinatlas.org/{ensg}.xml"
    try:
        r = SESSION.get(xml_url, timeout=30)
        r.raise_for_status()
        root = ET.fromstring(r.content)
    except Exception as e:
        if log:
            log(f"XML metadata unavailable: {e}")
        _XML_METADATA_CACHE[ensg] = {}
        return {}

    parent_map = {child: parent for parent in root.iter() for child in parent}
    image_meta = {}

    for el in root.iter():
        chunks = []
        if el.text:
            chunks.append(el.text)
        for v in el.attrib.values():
            chunks.append(v)
        joined = " ".join(chunks)
        links = extract_image_links_from_text(joined)

        # Also catch raw filenames or relative paths in XML text/attributes.
        if not links and re.search(r"\.(jpg|jpeg|tif|tiff)\b", joined, re.IGNORECASE):
            links = [joined]

        for link in links:
            key = clean_image_key(link)
            if not key:
                continue
            meta = collect_xml_context_metadata(el, parent_map)
            if not meta:
                continue
            if key not in image_meta:
                image_meta[key] = meta
            else:
                for k, v in meta.items():
                    image_meta[key].setdefault(k, v)

    _XML_METADATA_CACHE[ensg] = image_meta
    return image_meta

def extract_tissue_staining_by_antibody(soup, antibody_ids):
    """
    Extract a compact tissue-level antibody staining summary from real HTML tables only.

    This deliberately avoids parsing the whole page text as a fallback, because
    normal tissue pages also contain large RNA/GTEx/FANTOM tables that can pollute
    CSV rows with thousands of unrelated values.
    """
    result = {ab: [] for ab in antibody_ids}
    valid_values = re.compile(r"^(not detected|low|medium|high)$", re.IGNORECASE)

    for table in soup.find_all("table"):
        rows = []
        for tr in table.find_all("tr"):
            cells = [clean_text_value(c.get_text(" ", strip=True), max_len=250) for c in tr.find_all(["th", "td"])]
            cells = [c for c in cells if c]
            if cells:
                rows.append(cells)
        if not rows:
            continue

        header_idx = None
        header_abs = []
        for i, cells in enumerate(rows):
            joined = " ".join(cells)
            found = [ab for ab in antibody_ids if ab in joined]
            if found:
                header_idx = i
                header_abs = found
                break
        if header_idx is None or not header_abs:
            continue

        for cells in rows[header_idx + 1:]:
            if len(cells) < 2:
                continue
            cell_type = cells[0]
            if not is_plausible_annotation_label(cell_type):
                continue
            values = cells[1:1 + len(header_abs)]
            for ab, val in zip(header_abs, values):
                val = clean_text_value(val, max_len=80)
                if valid_values.match(val):
                    val = "Not detected" if val.lower() == "not detected" else val.capitalize()
                    pair = (cell_type, val)
                    if pair not in result.setdefault(ab, []):
                        result[ab].append(pair)

    summary = {}
    for ab, pairs in result.items():
        if not pairs:
            continue
        annotations = []
        for cell_type, staining in pairs:
            annotations.append({
                "AnnotationType": cell_type,
                "Staining": staining,
                "Intensity": None,
                "Quantity": None,
                "Location": None,
            })
        summary[ab] = {
            "AllAnnotations": annotations,
            "AnnotationType": annotations[0]["AnnotationType"],
            "AntibodyStaining": annotations[0]["Staining"],
            "AnnotationSummary": compact_annotation_summary(annotations),
        }
    return summary

def download_image_with_retry(image_link, image_path, max_retries=3, log=None):
    """Download a single image file with retry support and streaming writes."""
    for attempt in range(max_retries):
        try:
            if log:
                log(f"Downloading: {image_path.name}")

            # Stream the response to disk instead of loading the complete TIFF/JPG
            # into memory. This is safer for large HPA TIFF files.
            with SESSION.get(image_link, timeout=30, stream=True) as r:
                r.raise_for_status()
                with open(image_path, "wb") as out:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            out.write(chunk)

            if log:
                log(f"✓ OK: {image_path.name}")
            return True
        except Exception as e:
            # Remove incomplete files so a failed partial download is not treated
            # as an already existing valid image in a future run.
            try:
                if image_path.exists():
                    image_path.unlink()
            except Exception:
                pass

            if log:
                log(f"✗ Error on attempt {attempt + 1}: {e}")
            if attempt < max_retries - 1:
                time.sleep(1)
            else:
                if log:
                    log(f"✗ Final download failure: {image_link}")
                return False


# ---------------------------------------------------------------------
# Preview inventory construction
# ---------------------------------------------------------------------

def _build_preview_inventory_single(url, img_ext=".tif"):
    """
    Build an in-memory preview inventory without downloading files.

    This version keeps the old cancer-table logic when possible, but adds a
    generic HPA image scanner for tissue pages and newer/alternative layouts.
    Items are organized as: gene -> antibody -> category. Category means cancer
    subtype for cancer pages or tissue name for tissue pages.
    """
    html = download_html(url)
    soup = BeautifulSoup(html, "lxml")

    gene_name = extract_gene_name(soup)
    antibody_ids = extract_antibody_ids(soup)
    if not antibody_ids:
        raise RuntimeError("No antibody IDs were found on the page.")

    antibody_by_number = {}
    for ab in antibody_ids:
        num = normalize_antibody_numeric(ab)
        if num:
            antibody_by_number[num] = ab

    section = infer_hpa_section_from_url(url)
    page_category = infer_page_category_from_url(url, section)
    inferred_tissue = infer_tissue_name_from_url(url)
    xml_image_metadata = build_xml_image_metadata_map(url)
    tissue_staining_by_ab = extract_tissue_staining_by_antibody(soup, antibody_ids) if section == "tissue" else {}

    inv = {}
    total = 0

    def add_item(antibody_id, category, title_text, href, fallback_category_type):
        nonlocal total

        meta = parse_metadata_from_title(title_text or "")

        # Enrich title-derived metadata with XML metadata when HPA thumbnails
        # do not expose patient/sample data in HTML attributes.
        image_key = clean_image_key(href)
        if image_key and image_key in xml_image_metadata:
            meta = merge_missing_metadata(meta, xml_image_metadata[image_key])

        # Add tissue-level antibody staining table summary when the individual
        # thumbnail title does not contain the annotation blocks.
        if section == "tissue" and antibody_id in tissue_staining_by_ab:
            stain_meta = tissue_staining_by_ab[antibody_id]
            for k, v in stain_meta.items():
                if k == "AllAnnotations" and not meta.get("AllAnnotations"):
                    meta[k] = v
                elif not meta.get(k):
                    meta[k] = v

        # Rebuild a compact final summary after all enrichments.
        if meta.get("AllAnnotations"):
            meta["AnnotationSummary"] = compact_annotation_summary(meta.get("AllAnnotations"))

        diagnostic_category = ""
        diagnostic_folder = ""

        if section == "cancer":
            # Top folder = explicit HPA cancer category from the URL when present
            # (e.g. colorectal cancer, breast cancer, lung cancer).
            # Subfolder = HPA diagnostic/histology label only when it is actually
            # present in the HPA table/title/XML (e.g. Adenocarcinoma, Duct
            # carcinoma, Squamous cell carcinoma). No synthetic subtype is made.
            raw_diagnostic = meta.get("CancerType") or meta.get("Diagnosis") or category or ""
            diagnostic_category = valid_hpa_diagnostic_subtype(raw_diagnostic, section_category=page_category)
            final_category = page_category or diagnostic_category or "Cancer"
            category_type = "CancerCategory" if page_category else "CancerType"
            if diagnostic_category:
                diagnostic_folder = normalize_category_folder(diagnostic_category)
        elif section == "tissue":
            final_category = meta.get("Tissue") or page_category or category or inferred_tissue or "Tissue"
            category_type = "Tissue"
            if not meta.get("Tissue") and (page_category or inferred_tissue):
                meta["Tissue"] = page_category or inferred_tissue
        else:
            diagnostic_category = meta.get("CancerType") or meta.get("Diagnosis") or ""
            final_category = page_category or meta.get("CancerType") or meta.get("Tissue") or meta.get("Diagnosis") or category or "IHC"
            category_type = fallback_category_type or "Category"

        category_folder = normalize_category_folder(final_category)
        image_link = image_link_for_format(href, img_ext)

        inv.setdefault(antibody_id, {})
        bucket = inv[antibody_id].setdefault(category_folder, {"count": 0, "items": [], "_seen": set(), "_counters": {}})

        dedupe_key = (antibody_id, category_folder, image_link)
        if dedupe_key in bucket["_seen"]:
            return
        bucket["_seen"].add(dedupe_key)

        patient_id = meta["PatientID"] if meta.get("PatientID") else "unknown"
        counter_key = (category_folder, patient_id)
        idx = bucket["_counters"].get(counter_key, 0) + 1
        bucket["_counters"][counter_key] = idx

        image_name = f"ID_{patient_id}_{idx}{img_ext}"
        bucket["count"] += 1
        row = {
            "Gene": gene_name,
            "HPASection": section,
            "AntibodyID": antibody_id,
            "CategoryFolder": category_folder,
            "CategoryType": category_type,
            "Category": final_category,
            "CancerFolder": category_folder,
            "SectionCategory": final_category,
            "SectionCategoryFolder": category_folder,
            "DiagnosticCategory": diagnostic_category,
            "DiagnosticCategoryFolder": diagnostic_folder,
            "DownloadSubfolder": diagnostic_folder,
            "ImageName": image_name,
            "ImageLink": image_link,
            "ImageKey": image_key or "",
            "PatientID": patient_id,
            "Gender": meta.get("Gender"),
            "Age": meta.get("Age"),
            "Tissue": meta.get("Tissue"),
            "TissueCode": meta.get("TissueCode"),
            "Diagnosis": meta.get("Diagnosis"),
            "DiagnosisCode": meta.get("DiagnosisCode"),
            "CancerType": meta.get("CancerType") or (diagnostic_category if section == "cancer" else None),
            "CancerCode": meta.get("CancerCode"),
            "LocationCodes": "; ".join(meta.get("LocationsRaw", [])) if meta.get("LocationsRaw") else "",
            "AnnotationType": meta.get("AnnotationType"),
            "AntibodyStaining": meta.get("AntibodyStaining"),
            "Intensity": meta.get("Intensity"),
            "Quantity": meta.get("Quantity"),
            "Location": meta.get("Location"),
            "AnnotationSummary": meta.get("AnnotationSummary", ""),
        }
        bucket["items"].append(sanitize_csv_row(row))
        total += 1

    # ------------------------------------------------------------------
    # 1) Original cancer-specific parser. This preserves the behavior that
    #    was already working on HPA cancer IHC pages.
    # ------------------------------------------------------------------
    if section == "cancer":
        rows = soup.find_all("tr")
        for antibody_id in antibody_ids:
            antibody_number = normalize_antibody_numeric(antibody_id)
            if antibody_number is None:
                continue
            for tr in rows:
                ths = tr.find_all("th")
                tds = tr.find_all("td")
                if len(ths) < 1 or len(tds) < 1:
                    continue
                for col_idx, th in enumerate(ths):
                    cancer_label = th.get_text(" ", strip=True)
                    if not cancer_label or col_idx >= len(tds):
                        continue
                    td = tds[col_idx]
                    cancer_divs = td.find_all("div", class_="cancerAnnoations")
                    for div in cancer_divs:
                        for a in div.find_all("a", href=True):
                            href = normalize_image_href(a["href"])
                            if not is_image_href(href):
                                continue
                            href_antibody_number = extract_antibody_numeric_from_href(href)
                            if href_antibody_number != antibody_number:
                                continue
                            title_text = collect_title_candidate(a)
                            add_item(antibody_id, cancer_label, title_text, href, "CancerType")

    # ------------------------------------------------------------------
    # 2) Generic parser. Important for tissue pages, where the links may be
    #    simple anchors to images.proteinatlas.org and not cancerAnnoations
    #    divs. Also acts as fallback if HPA changes the cancer markup.
    # ------------------------------------------------------------------
    for a in soup.find_all("a", href=True):
        href = normalize_image_href(a["href"])
        if not is_image_href(href):
            continue

        href_antibody_number = extract_antibody_numeric_from_href(href)
        if not href_antibody_number:
            continue

        antibody_id = antibody_by_number.get(href_antibody_number)
        if not antibody_id:
            continue

        title_text = collect_title_candidate(a)

        # If a title contains an explicit different antibody, skip it.
        ids_in_title = set(re.findall(r"\b(HPA\d{6}|CAB\d{6})\b", title_text or ""))
        if ids_in_title and antibody_id not in ids_in_title:
            continue

        fallback_category = inferred_tissue if section == "tissue" else None
        add_item(antibody_id, fallback_category, title_text, href, "Category")

    # Clean helper fields and remove empty antibodies/categories.
    clean_inv = {}
    for ab_id, categories in inv.items():
        clean_categories = {}
        for category, payload in categories.items():
            payload.pop("_seen", None)
            payload.pop("_counters", None)
            if payload.get("items"):
                clean_categories[category] = payload
        if clean_categories:
            clean_inv[ab_id] = clean_categories

    if not clean_inv:
        raise RuntimeError(
            "No downloadable IHC image records were found. "
            "The parser found antibody IDs, but no matching image links. "
            "Try opening the page in a browser and confirm it contains IHC thumbnails."
        )

    return gene_name, clean_inv, total



def get_gene_url_segment(url):
    """Return the HPA gene URL segment, e.g. ENSG00000146648-EGFR."""
    try:
        parts = [p for p in urlparse(url).path.split("/") if p]
    except Exception:
        return ""
    for p in parts:
        if p.startswith("ENSG") and "-" in p:
            return p
    return parts[0] if parts else ""


def infer_page_category_from_url(url, section=None):
    """
    Return the explicit HPA page category from URLs such as:
      /ENSG...-GENE/cancer/colorectal+cancer -> colorectal cancer
      /ENSG...-GENE/tissue/colon -> colon

    This is intentionally different from diagnosis/cell-type metadata. For
    cancer overview downloads it lets us keep folders like colorectal cancer,
    breast cancer, lung cancer, etc., instead of mixing all images under
    histological labels such as Adenocarcinoma.
    """
    section = (section or infer_hpa_section_from_url(url) or "").lower()
    try:
        parts = [unquote(p) for p in urlparse(url).path.split("/") if p]
    except Exception:
        return ""

    for i, part in enumerate(parts):
        if part.lower() == section and i + 1 < len(parts):
            category = parts[i + 1].split("#", 1)[0].split("?", 1)[0].strip()
            if category:
                return clean_text_value(category.replace("+", " ").replace("-", " "), max_len=120)
    return ""


def is_hpa_section_overview_url(url, section=None):
    """True for overview pages like /ENSG...-GENE/cancer or /tissue."""
    try:
        parts = [unquote(p) for p in urlparse(url).path.split("/") if p]
    except Exception:
        return False
    if len(parts) != 2:
        return False
    if not (parts[0].startswith("ENSG") and "-" in parts[0]):
        return False
    if section:
        return parts[1].lower() == section.lower()
    return parts[1].lower() in {"cancer", "tissue"}


def discover_available_category_urls(url, soup=None, section=None):
    """
    Discover child category pages from an HPA overview page.

    Example:
      /ENSG00000146648-EGFR/cancer
    becomes:
      /ENSG00000146648-EGFR/cancer/colorectal+cancer
      /ENSG00000146648-EGFR/cancer/lung+cancer
      etc.
    """
    section = section or infer_hpa_section_from_url(url)
    if section not in {"cancer", "tissue"}:
        return []

    if soup is None:
        soup = BeautifulSoup(download_html(url), "lxml")

    gene_seg = get_gene_url_segment(url)
    if not gene_seg:
        return []

    category_urls = []
    seen = set()
    wanted_prefix = f"/{gene_seg}/{section}/"

    for a in soup.find_all("a", href=True):
        href_abs = urljoin(url, a["href"])
        parsed = urlparse(href_abs)
        path = unquote(parsed.path)

        if not path.startswith(wanted_prefix):
            continue

        category_part = path.split(wanted_prefix, 1)[1].strip("/")
        if not category_part:
            continue

        # Only one category level here. Avoid unrelated deeper pages.
        if "/" in category_part:
            continue

        # Ignore non-category anchors that may be injected by scripts/search.
        category_part = category_part.split("#", 1)[0].split("?", 1)[0].strip()
        if not category_part:
            continue

        normalized_url = f"https://www.proteinatlas.org/{gene_seg}/{section}/{category_part}"
        key = normalized_url.lower()
        if key not in seen:
            seen.add(key)
            category_urls.append(normalized_url)

    return category_urls


def merge_single_gene_inventory(target_inv, source_inv):
    """Merge antibody -> category inventory dictionaries."""
    for ab_id, categories in source_inv.items():
        ab_bucket = target_inv.setdefault(ab_id, {})
        for cat, payload in categories.items():
            target = ab_bucket.setdefault(cat, {"count": 0, "items": []})
            existing = set((r.get("AntibodyID"), r.get("CategoryFolder"), r.get("ImageLink")) for r in target.get("items", []))
            for item in payload.get("items", []):
                key = (item.get("AntibodyID"), item.get("CategoryFolder"), item.get("ImageLink"))
                if key in existing:
                    continue
                existing.add(key)
                target["items"].append(item)
                target["count"] = target.get("count", 0) + 1


def build_preview_inventory(url, img_ext=".tif"):
    """
    Build preview inventory from either a direct category page or an overview page.

    Direct pages such as /cancer/colorectal+cancer are parsed directly.
    Overview pages such as /cancer are first expanded into their available
    cancer category pages, then each child page is parsed and merged.
    """
    # Overview pages usually do not contain all IHC image links directly.
    # HPA exposes cancer/tissue categories as child pages, so we expand them.
    section = infer_hpa_section_from_url(url)
    if is_hpa_section_overview_url(url, section) and section in {"cancer", "tissue"}:
        html = download_html(url)
        soup = BeautifulSoup(html, "lxml")
        gene_name = extract_gene_name(soup)
        category_urls = discover_available_category_urls(url, soup=soup, section=section)

        if category_urls:
            combined_inv = {}
            total = 0
            failures = []
            for child_url in category_urls:
                try:
                    child_gene, child_inv, child_total = _build_preview_inventory_single(child_url, img_ext=img_ext)
                    if child_gene and child_gene != "UnknownGene":
                        gene_name = child_gene
                    merge_single_gene_inventory(combined_inv, child_inv)
                    total += child_total
                except Exception as e:
                    failures.append(f"{child_url}: {e}")

            if combined_inv:
                return gene_name, combined_inv, total

            if failures:
                raise RuntimeError(
                    "No image records were found after expanding the overview page. "
                    "First category error: " + failures[0]
                )

    # Direct page or fallback behavior.
    return _build_preview_inventory_single(url, img_ext=img_ext)


# ---------------------------------------------------------------------
# Bulk preview and download execution
# ---------------------------------------------------------------------

def parse_bulk_urls(text):
    """Parse a pasted text block into a clean list of unique HPA URLs."""
    urls = []
    seen = set()
    for raw in re.split(r"[\r\n]+", text or ""):
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        # Allow lines like "FGFR2 https://..." by extracting the first HPA URL.
        m = re.search(r"https?://(?:www\.)?proteinatlas\.org/\S+", raw)
        url = m.group(0).strip() if m else raw
        url = url.rstrip(",;)")
        if is_valid_hpa_url(url) and url not in seen:
            seen.add(url)
            urls.append(url)
    return urls


def merge_inventory_into_bulk(bulk_inv, gene_name, inv):
    """Merge a single-gene preview inventory into a gene -> antibody -> category inventory."""
    gene_bucket = bulk_inv.setdefault(gene_name, {})
    for ab_id, categories in inv.items():
        ab_bucket = gene_bucket.setdefault(ab_id, {})
        for cat, payload in categories.items():
            target = ab_bucket.setdefault(cat, {"count": 0, "items": [], "_seen": set()})
            for item in payload.get("items", []):
                key = (
                    item.get("Gene"),
                    item.get("HPASection"),
                    item.get("AntibodyID"),
                    item.get("CategoryFolder") or item.get("CancerFolder"),
                    item.get("ImageLink"),
                    item.get("ImageKey"),
                )
                if key in target["_seen"]:
                    continue
                target["_seen"].add(key)
                target["items"].append(item)
                target["count"] += 1


def clean_bulk_inventory(bulk_inv):
    """Remove temporary helper fields and empty branches from bulk inventory."""
    clean = {}
    total = 0
    for gene, inv in bulk_inv.items():
        gene_clean = {}
        for ab_id, categories in inv.items():
            ab_clean = {}
            for cat, payload in categories.items():
                payload.pop("_seen", None)
                if payload.get("items"):
                    ab_clean[cat] = payload
                    total += payload.get("count", len(payload.get("items", [])))
            if ab_clean:
                gene_clean[ab_id] = ab_clean
        if gene_clean:
            clean[gene] = gene_clean
    return clean, total


def build_bulk_preview_inventory(urls, img_ext=".tif", log=None):
    """Build a preview inventory for multiple HPA URLs."""
    bulk_inv = {}
    for i, url in enumerate(urls, start=1):
        if log:
            log(f"[{i}/{len(urls)}] Previewing: {url}")
        gene, inv, total = build_preview_inventory(url, img_ext=img_ext)
        if log:
            log(f"    ✓ {gene}: {total} items")
        merge_inventory_into_bulk(bulk_inv, gene, inv)
    clean, total = clean_bulk_inventory(bulk_inv)
    if not clean:
        raise RuntimeError("No downloadable IHC image records were found in the provided URLs.")
    return clean, total


NO_SUBTYPE_KEY = "__NO_HPA_DIAGNOSTIC_SUBTYPE__"
NO_SUBTYPE_LABEL = "No HPA diagnostic subtype provided"


def row_diagnostic_subtype_key(row):
    """Return the validated HPA diagnostic subtype key used for subtype-level selection."""
    subtype = valid_hpa_diagnostic_subtype(
        row.get("DiagnosticCategory") or row.get("CancerType") or row.get("Diagnosis") or "",
        section_category=row.get("SectionCategory") or row.get("Category") or row.get("CategoryFolder") or "",
    )
    return subtype or NO_SUBTYPE_KEY


def selected_rows_from_bulk_inventory(bulk_inv, selection):
    """Return selected rows from gene -> antibody -> category inventory.

    Selection supports two levels:
      1) category-level selection, preserving the original behavior;
      2) optional diagnostic subtype-level selection for cancer categories.

    If no subtype filter is present for a selected category, all rows in that
    category are included. This keeps CLI and older GUI behavior compatible.
    """
    categories_by_gene_ab = selection.get("categories_by_gene_antibody", {})
    subtypes_by_gene_ab_category = selection.get("subtypes_by_gene_antibody_category", {})
    rows = []
    seen = set()
    for gene, inv in bulk_inv.items():
        for ab_id, categories in inv.items():
            selected_categories = categories_by_gene_ab.get((gene, ab_id), set())
            for cat, payload in categories.items():
                if cat not in selected_categories:
                    continue

                subtype_filter = subtypes_by_gene_ab_category.get((gene, ab_id, cat))

                for item in payload.get("items", []):
                    if subtype_filter is not None:
                        subtype_key = row_diagnostic_subtype_key(item)
                        if subtype_key not in subtype_filter:
                            continue

                    key = (
                        item.get("Gene"),
                        item.get("HPASection"),
                        item.get("AntibodyID"),
                        item.get("CategoryFolder") or item.get("CancerFolder"),
                        item.get("ImageLink"),
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    rows.append(item)
    return rows


def write_summary_csv(path, rows):
    """Write rows to a CSV file with stable fieldnames and Excel-safe text."""
    if not rows:
        return

    preferred = [
        "Gene", "HPASection", "AntibodyID",
        "CategoryFolder", "CategoryType", "Category",
        "SectionCategory", "SectionCategoryFolder",
        "DiagnosticCategory", "DiagnosticCategoryFolder",
        "CancerFolder", "StoredRelativePath", "DownloadStatus", "FailureReason",
        "FileSizeBytes", "SHA256", "DownloadDate", "SoftwareVersion",
        "ImageName", "ImageLink", "ImageKey",
        "PatientID", "Gender", "Age",
        "Tissue", "TissueCode",
        "Diagnosis", "DiagnosisCode",
        "CancerType", "CancerCode", "LocationCodes",
        "AnnotationType", "AntibodyStaining", "Intensity", "Quantity",
        "Location", "AnnotationSummary",
    ]
    all_keys = []
    for r in rows:
        for k in r.keys():
            if k not in all_keys:
                all_keys.append(k)
    fieldnames = [k for k in preferred if k in all_keys] + [k for k in all_keys if k not in preferred]

    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(sanitize_csv_row(r) for r in rows)




def compute_file_integrity(path):
    """Return file size and SHA256 checksum for reproducible dataset audits."""
    path = Path(path)
    size = path.stat().st_size
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return size, h.hexdigest()


def flatten_bulk_inventory_rows(bulk_inv):
    """Flatten gene -> antibody -> category inventory into preview rows."""
    rows = []
    seen = set()
    for gene, inv in bulk_inv.items():
        for ab_id, categories in inv.items():
            for cat, payload in categories.items():
                for item in payload.get("items", []):
                    row = dict(item)
                    key = (
                        row.get("Gene"), row.get("HPASection"), row.get("AntibodyID"),
                        row.get("CategoryFolder") or row.get("CancerFolder"),
                        row.get("ImageLink"), row.get("ImageKey"),
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    rows.append(row)
    return rows


def summarize_rows_by_gene_antibody_category(rows):
    """Return compact counts for Markdown reports."""
    summary = {}
    for row in rows or []:
        key = (
            row.get("Gene") or "UnknownGene",
            row.get("AntibodyID") or "UnknownAntibody",
            row.get("HPASection") or "ihc",
            row.get("Category") or row.get("CategoryFolder") or row.get("CancerFolder") or "Unknown",
            row.get("DiagnosticCategory") or "",
        )
        summary[key] = summary.get(key, 0) + 1
    return summary


def write_json(path, data):
    """Write UTF-8 JSON with stable indentation."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def make_reproducibility_manifest(input_urls, img_ext, preview_total, completed_rows, failed_rows, operation):
    """Build a machine-readable manifest describing the run."""
    completed_rows = completed_rows or []
    failed_rows = failed_rows or []
    genes = sorted({r.get("Gene") for r in completed_rows + failed_rows if r.get("Gene")})
    antibodies = sorted({r.get("AntibodyID") for r in completed_rows + failed_rows if r.get("AntibodyID")})
    sections = sorted({r.get("HPASection") for r in completed_rows + failed_rows if r.get("HPASection")})
    return {
        "software": "HPA IHC Image Downloader",
        "version": __version__,
        "author": __author__,
        "github_url": __github_url__,
        "doi": __doi__,
        "run_datetime": datetime.now().isoformat(timespec="seconds"),
        "operation": operation,
        "input_urls": list(input_urls or []),
        "image_format": img_ext,
        "output_structure": "gene/antibody/HPA category/HPA diagnostic subtype when available",
        "preview_total_items": int(preview_total or 0),
        "downloaded_or_available_items": len(completed_rows),
        "failed_items": len(failed_rows),
        "genes": genes,
        "antibodies": antibodies,
        "hpa_sections": sections,
        "python_version": sys.version.replace("\n", " "),
        "platform": platform.platform(),
    }


def write_citation_helper(path, manifest):
    """Write a human-readable citation and methods helper file."""
    path = Path(path)
    urls = manifest.get("input_urls") or []
    urls_text = "\n".join(f"- {u}" for u in urls) if urls else "- Not recorded"
    text = f"""# Citation and methods helper

## Suggested software citation

Rodríguez-Rojas J. HPA IHC Image Downloader. Version {manifest.get('version', __version__)}. {__github_url__}. DOI: {__doi__}.

## Suggested methods text

Immunohistochemistry image records were retrieved from the Human Protein Atlas using HPA IHC Image Downloader v{manifest.get('version', __version__)}. The software was used to preview, select, download, organize, and document HPA image records by gene, antibody, HPA section, HPA category, and diagnostic subtype when available. The download run exported structured metadata tables, a reproducibility manifest, file-size information, and SHA256 checksums for integrity verification.

## Run summary

- Run date/time: {manifest.get('run_datetime')}
- Operation: {manifest.get('operation')}
- Image format: {manifest.get('image_format')}
- Previewed items: {manifest.get('preview_total_items')}
- Downloaded or already available items: {manifest.get('downloaded_or_available_items')}
- Failed items: {manifest.get('failed_items')}
- Software version: {manifest.get('version')}
- Python: {manifest.get('python_version')}
- Platform: {manifest.get('platform')}

## Source URLs

{urls_text}

## BibTeX template

```bibtex
@software{{rodriguez_rojas_hpa_ihc_image_downloader,
  title = {{HPA IHC Image Downloader}},
  author = {{Rodríguez-Rojas, José}},
  version = {{{manifest.get('version', __version__)}}},
  doi = {{{__doi__.replace('https://doi.org/', '')}}},
  url = {{{__github_url__}}},
  year = {{{datetime.now().year}}}
}}
```
"""
    path.write_text(text, encoding="utf-8")


def write_markdown_report(path, manifest, completed_rows=None, failed_rows=None, preview_rows=None):
    """Write a compact audit report in Markdown."""
    completed_rows = completed_rows or []
    failed_rows = failed_rows or []
    preview_rows = preview_rows or []
    path = Path(path)

    lines = []
    lines.append("# HPA IHC Image Downloader report")
    lines.append("")
    lines.append("## Run summary")
    lines.append("")
    lines.append(f"- Software version: {manifest.get('version')}")
    lines.append(f"- Run date/time: {manifest.get('run_datetime')}")
    lines.append(f"- Operation: {manifest.get('operation')}")
    lines.append(f"- Image format: {manifest.get('image_format')}")
    lines.append(f"- Previewed items: {manifest.get('preview_total_items')}")
    lines.append(f"- Downloaded or already available items: {manifest.get('downloaded_or_available_items')}")
    lines.append(f"- Failed items: {manifest.get('failed_items')}")
    lines.append(f"- Platform: {manifest.get('platform')}")
    lines.append("")

    rows_for_summary = completed_rows or preview_rows
    summary = summarize_rows_by_gene_antibody_category(rows_for_summary)
    if summary:
        lines.append("## Counts by gene, antibody, section, category, and diagnostic subtype")
        lines.append("")
        lines.append("| Gene | Antibody | Section | Category | Diagnostic subtype | Count |")
        lines.append("|---|---|---|---|---|---:|")
        for (gene, ab, sec, cat, diag), count in sorted(summary.items()):
            lines.append(f"| {gene} | {ab} | {sec} | {cat} | {diag or '-'} | {count} |")
        lines.append("")

    if failed_rows:
        lines.append("## Failed downloads")
        lines.append("")
        lines.append("| Gene | Antibody | Category | Image name | Reason | URL |")
        lines.append("|---|---|---|---|---|---|")
        for row in failed_rows:
            lines.append(
                f"| {row.get('Gene','')} | {row.get('AntibodyID','')} | "
                f"{row.get('Category') or row.get('CategoryFolder') or ''} | "
                f"{row.get('ImageName','')} | {row.get('FailureReason','')} | {row.get('ImageLink','')} |"
            )
        lines.append("")

    lines.append("## Source URLs")
    lines.append("")
    for u in manifest.get("input_urls") or []:
        lines.append(f"- {u}")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_preview_outputs(root_dir, bulk_inv, input_urls=None, img_ext=".tif", log=None):
    """Export preview inventory before downloading, useful for reproducible discovery."""
    root_dir = Path(root_dir)
    rows = flatten_bulk_inventory_rows(bulk_inv)
    if not rows:
        return None
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    preview_dir = root_dir / f"HPA_preview_inventory_{stamp}"
    preview_dir.mkdir(parents=True, exist_ok=True)
    write_summary_csv(preview_dir / "preview_inventory.csv", rows)
    manifest = make_reproducibility_manifest(input_urls, img_ext, len(rows), [], [], operation="preview")
    write_json(preview_dir / "preview_manifest.json", manifest)
    write_markdown_report(preview_dir / "preview_report.md", manifest, preview_rows=rows)
    write_citation_helper(preview_dir / "citation_and_methods_helper.md", manifest)
    if log:
        log(f"Preview inventory exported: {preview_dir}")
    return preview_dir

def write_summary_outputs(root_dir, rows, log=None):
    """Write per-gene, per-section, and bulk CSV summaries."""
    if not rows:
        return

    rows_by_gene = {}
    rows_by_section = {}
    for row in rows:
        gene = clean_text_value(row.get("Gene") or "UnknownGene", max_len=80) or "UnknownGene"
        sec = clean_text_value(row.get("HPASection") or "ihc", max_len=40).lower() or "ihc"
        rows_by_gene.setdefault(gene, []).append(row)
        rows_by_section.setdefault(sec, []).append(row)

    for gene, gene_rows in rows_by_gene.items():
        gene_dir = root_dir / safe_folder_name(gene)
        gene_dir.mkdir(parents=True, exist_ok=True)
        write_summary_csv(gene_dir / f"{gene}_all_antibodies_summary.csv", gene_rows)

        gene_rows_by_section = {}
        for row in gene_rows:
            sec = clean_text_value(row.get("HPASection") or "ihc", max_len=40).lower() or "ihc"
            gene_rows_by_section.setdefault(sec, []).append(row)
        for sec, sec_rows in gene_rows_by_section.items():
            write_summary_csv(gene_dir / f"{gene}_{sec}_summary.csv", sec_rows)
            if log:
                log(f"Summary CSV saved: {gene}/{gene}_{sec}_summary.csv")

    # Bulk CSVs at the selected output root. These are useful when several URLs/genes are processed together.
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    bulk_dir = root_dir / f"HPA_bulk_summary_{stamp}"
    bulk_dir.mkdir(parents=True, exist_ok=True)
    write_summary_csv(bulk_dir / "bulk_all_genes_summary.csv", rows)
    for sec, sec_rows in rows_by_section.items():
        write_summary_csv(bulk_dir / f"bulk_{sec}_summary.csv", sec_rows)
    if log:
        log(f"Bulk summary CSVs saved in: {bulk_dir}")



def write_download_report_outputs(root_dir, completed_rows, failed_rows=None, input_urls=None, img_ext=None, preview_total=None, log=None):
    """
    Write publication-oriented download outputs:
    - downloaded_successfully.csv
    - failed_downloads.csv only if failures exist
    - download_manifest.json
    - download_report.md
    - citation_and_methods_helper.md
    """
    completed_rows = completed_rows or []
    failed_rows = failed_rows or []
    if not completed_rows and not failed_rows:
        return None

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_dir = Path(root_dir) / f"HPA_download_report_{stamp}"
    report_dir.mkdir(parents=True, exist_ok=True)

    if completed_rows:
        write_summary_csv(report_dir / "downloaded_successfully.csv", completed_rows)

    if failed_rows:
        write_summary_csv(report_dir / "failed_downloads.csv", failed_rows)

    manifest = make_reproducibility_manifest(
        input_urls=input_urls,
        img_ext=img_ext,
        preview_total=preview_total if preview_total is not None else len(completed_rows) + len(failed_rows),
        completed_rows=completed_rows,
        failed_rows=failed_rows,
        operation="download",
    )
    write_json(report_dir / "download_manifest.json", manifest)
    write_markdown_report(report_dir / "download_report.md", manifest, completed_rows=completed_rows, failed_rows=failed_rows)
    write_citation_helper(report_dir / "citation_and_methods_helper.md", manifest)

    if log:
        log(f"Download report saved in: {report_dir}")
        log(f"Manifest saved: {report_dir / 'download_manifest.json'}")
        if failed_rows:
            log(f"Failed download CSV saved: {report_dir / 'failed_downloads.csv'}")
    return report_dir


def download_from_bulk_inventory(root_dir, bulk_inv, img_ext, selection, progress_cb=None, log=None, input_urls=None):
    """Download selected images from a bulk inventory."""
    root_dir.mkdir(parents=True, exist_ok=True)
    download_list = selected_rows_from_bulk_inventory(bulk_inv, selection)
    if not download_list:
        raise RuntimeError("There are no selected items to download.")

    total = len(download_list)
    if log:
        log(f"Total unique items to download: {total}")

    completed_rows = []
    failed_rows = []
    for i, row in enumerate(download_list, start=1):
        gene = safe_folder_name(row.get("Gene") or "UnknownGene")
        ab_id = row.get("AntibodyID") or "UnknownAntibody"
        category_folder = row.get("CategoryFolder") or row.get("CancerFolder") or "Unknown"

        image_link = image_link_for_format(row["ImageLink"], img_ext)
        row["ImageLink"] = image_link

        category_dir = root_dir / gene / ab_id / safe_folder_name(category_folder)

        # Optional second-level folder for cancer diagnosis/histology inside the
        # broader cancer category. Example:
        # FGFR2/CAB010886/colorectal cancer/Adenocarcinoma/ID_1506_1.tif
        subfolder = valid_hpa_diagnostic_subtype(
            row.get("DownloadSubfolder") or row.get("DiagnosticCategoryFolder") or row.get("DiagnosticCategory") or "",
            section_category=row.get("SectionCategory") or row.get("Category") or category_folder,
        )
        if (row.get("HPASection") or "").lower() == "cancer" and subfolder:
            category_dir = category_dir / safe_folder_name(subfolder)

        category_dir.mkdir(parents=True, exist_ok=True)

        image_path = category_dir / row["ImageName"]
        row["StoredRelativePath"] = str(image_path.relative_to(root_dir)).replace("\\", "/")

        # A zero-byte file is likely a partial/failed previous download, so it
        # should not be treated as a valid existing image.
        if image_path.exists():
            try:
                if image_path.stat().st_size == 0:
                    image_path.unlink()
                    if log:
                        log(f"⚠ Removed empty existing file before retry: {image_path.name}")
            except Exception:
                pass

        already_exists = image_path.exists()
        download_ok = True
        failure_reason = ""

        if not already_exists:
            download_ok = download_image_with_retry(image_link, image_path, log=log)
            if not download_ok:
                failure_reason = "Download failed after retry attempts"
        else:
            if log:
                log(f"📁 Already exists: {row.get('StoredRelativePath') or image_path.name}")

        if download_ok and image_path.exists():
            completed_row = dict(row)
            completed_row["DownloadStatus"] = "already_exists" if already_exists else "downloaded"
            completed_row["FailureReason"] = ""
            completed_row["DownloadDate"] = datetime.now().isoformat(timespec="seconds")
            completed_row["SoftwareVersion"] = __version__
            try:
                fsize, sha256 = compute_file_integrity(image_path)
                completed_row["FileSizeBytes"] = fsize
                completed_row["SHA256"] = sha256
            except Exception as e:
                completed_row["FileSizeBytes"] = ""
                completed_row["SHA256"] = ""
                completed_row["FailureReason"] = f"Could not compute file integrity metadata: {e}"
            completed_rows.append(completed_row)
        else:
            failed_row = dict(row)
            failed_row["DownloadStatus"] = "failed"
            failed_row["FailureReason"] = failure_reason or "File was not present after download attempt"
            failed_row["DownloadDate"] = datetime.now().isoformat(timespec="seconds")
            failed_row["SoftwareVersion"] = __version__
            failed_row["FileSizeBytes"] = ""
            failed_row["SHA256"] = ""
            failed_rows.append(failed_row)
            if log:
                log(f"⚠ Download failed and was written to failure report: {image_link}")

        if progress_cb:
            progress_cb(i, total)

    if completed_rows:
        write_summary_outputs(root_dir, completed_rows, log=log)
    write_download_report_outputs(root_dir, completed_rows, failed_rows, input_urls=input_urls, img_ext=img_ext, preview_total=total, log=log)

    if log:
        log(f"Done. Output folder: {root_dir}")
        log(f"Successful/available files: {len(completed_rows)}")
        if failed_rows:
            log(f"Failed files: {len(failed_rows)}")

def download_from_inventory(root_dir, gene_name, inv, img_ext, selection, progress_cb=None, log=None):
    """Backward-compatible wrapper for single-gene downloads."""
    old = selection.get("cancers_by_antibody", {})
    converted = {"categories_by_gene_antibody": {(gene_name, ab): cats for ab, cats in old.items()}}
    download_from_bulk_inventory(root_dir, {gene_name: inv}, img_ext, converted, progress_cb=progress_cb, log=log)


# ---------------------------------------------------------------------
# Graphical user interface
# ---------------------------------------------------------------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("HPA IHC Image Downloader")
        icon_path = resource_path("assets/icons/hpa_jjrr_icon.ico")
        try:
            self.iconbitmap(str(icon_path))
        except Exception:
            pass
        self.geometry("1280x820")
        self.minsize(1120, 720)

        self.inv = None                  # Backward compatibility: current single-gene inventory
        self.gene_name = None             # Backward compatibility: current single gene
        self.bulk_inv = None              # gene -> antibody -> category inventory
        self.total_items = 0
        self.preview_img_ext = None
        self.preview_urls = []

        self.msg_q = queue.Queue()
        self.output_dir = ROOT_DIR
        self.worker_threads = []
        self._closing = False
        self._poll_after_id = None
        self._build_ui()
        self.protocol("WM_DELETE_WINDOW", self.on_close)

        self.log(f"HPA IHC Image Downloader v{__version__}")
        self.log(f"DOI: {__doi__}")
        self.log(f"Output directory: {self.output_dir}")
        self.log("Tip: paste one HPA URL and click Preview, or use Bulk URLs for several links/genes.")
        self._poll_queue()

    def _build_ui(self):
        menubar = tk.Menu(self)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Select Output Directory...", command=self.select_output_directory)
        file_menu.add_command(label="Open Output Directory", command=self.open_output_directory)
        file_menu.add_command(label="Bulk URLs...", command=self.open_bulk_urls_dialog)
        file_menu.add_command(label="Clear preview", command=self.clear_preview)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.on_close)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="User Guide", command=self.show_help)
        help_menu.add_separator()
        help_menu.add_command(label="About", command=self.show_about)

        menubar.add_cascade(label="File", menu=file_menu)
        menubar.add_cascade(label="Help", menu=help_menu)
        self.config(menu=menubar)

        top = ttk.Frame(self, padding=10)
        top.pack(fill="x")

        url_row = ttk.Frame(top)
        url_row.pack(fill="x", pady=(0, 6))
        ttk.Label(url_row, text="URL:").pack(side="left")
        self.url_var = tk.StringVar()
        url_entry = ttk.Entry(url_row, textvariable=self.url_var)
        url_entry.pack(side="left", padx=8, fill="x", expand=True)

        ctrl_row = ttk.Frame(top)
        ctrl_row.pack(fill="x")
        ttk.Label(ctrl_row, text="Format:").pack(side="left")
        self.ext_var = tk.StringVar(value=".tif")
        ttk.Combobox(ctrl_row, textvariable=self.ext_var, values=[".tif", ".jpg"], width=6, state="readonly").pack(side="left", padx=(4, 12))

        self.preview_btn = ttk.Button(ctrl_row, text="Preview", command=self.on_preview)
        self.preview_btn.pack(side="left", padx=4)
        self.bulk_btn = ttk.Button(ctrl_row, text="Bulk URLs", command=self.open_bulk_urls_dialog)
        self.bulk_btn.pack(side="left", padx=4)
        self.download_btn = ttk.Button(ctrl_row, text="Download", command=self.on_download, state="disabled")
        self.download_btn.pack(side="left", padx=4)

        out_row = ttk.Frame(top)
        out_row.pack(fill="x", pady=(6, 0))
        ttk.Label(out_row, text="Output:").pack(side="left")
        self.output_var = tk.StringVar(value=str(self.output_dir))
        out_entry = ttk.Entry(out_row, textvariable=self.output_var, state="readonly")
        out_entry.pack(side="left", fill="x", expand=True, padx=(8, 0))

        main = ttk.PanedWindow(self, orient="horizontal")
        main.pack(fill="both", expand=True, padx=10, pady=(0, 10))

        left = ttk.Frame(main, padding=(0, 0, 6, 0))
        main.add(left, weight=3)
        ttk.Label(left, text="Preview (gene → antibody → cancer/tissue category)").pack(anchor="w")

        tree_box = ttk.Frame(left)
        tree_box.pack(fill="both", expand=True, pady=6)
        self.tree = ttk.Treeview(tree_box, columns=("count",), show="tree headings")
        self.tree.heading("#0", text="Item")
        self.tree.heading("count", text="Count")
        self.tree.column("#0", width=420, minwidth=260, stretch=True)
        self.tree.column("count", width=90, minwidth=70, anchor="e", stretch=False)
        tree_y = ttk.Scrollbar(tree_box, orient="vertical", command=self.tree.yview)
        tree_x = ttk.Scrollbar(tree_box, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=tree_y.set, xscrollcommand=tree_x.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        tree_y.grid(row=0, column=1, sticky="ns")
        tree_x.grid(row=1, column=0, sticky="ew")
        tree_box.rowconfigure(0, weight=1)
        tree_box.columnconfigure(0, weight=1)

        right = ttk.Frame(main, padding=(6, 0, 0, 0))
        main.add(right, weight=2)
        ttk.Label(right, text="Selection").pack(anchor="w")

        self.sel_canvas = tk.Canvas(right, highlightthickness=0)
        self.sel_scroll = ttk.Scrollbar(right, orient="vertical", command=self.sel_canvas.yview)
        self.sel_canvas.configure(yscrollcommand=self.sel_scroll.set)
        self.sel_scroll.pack(side="right", fill="y")
        self.sel_canvas.pack(side="left", fill="both", expand=True, pady=6)

        self.sel_frame = ttk.Frame(self.sel_canvas)
        self.sel_canvas.create_window((0, 0), window=self.sel_frame, anchor="nw")
        self.sel_frame.bind("<Configure>", lambda e: self.sel_canvas.configure(scrollregion=self.sel_canvas.bbox("all")))
        self.sel_canvas.bind_all("<MouseWheel>", self._on_mousewheel)

        self.gene_vars = {}
        self.ab_vars = {}
        self.category_vars = {}
        self.subtype_vars = {}
        self.cancer_vars = self.category_vars  # compatibility with older method names

        bottom = ttk.Frame(self, padding=10)
        bottom.pack(fill="x")
        self.progress = ttk.Progressbar(bottom, mode="determinate")
        self.progress.pack(fill="x")
        self.status_var = tk.StringVar(value="Ready.")
        ttk.Label(bottom, textvariable=self.status_var).pack(anchor="w", pady=(6, 0))

        log_frame = ttk.Frame(self, padding=(10, 0, 10, 10))
        log_frame.pack(fill="both", expand=False)
        ttk.Label(log_frame, text="Log").pack(anchor="w")
        self.log_text = tk.Text(log_frame, height=9, wrap="word")
        self.log_text.pack(fill="both", expand=True, pady=6)
        self.log_text.configure(state="disabled")

    def open_bulk_urls_dialog(self):
        """Open a larger, friendlier dialog for managing several HPA URLs."""
        win = tk.Toplevel(self)
        win.title("Bulk URL manager")
        win.geometry("920x680")
        win.minsize(820, 580)
        win.transient(self)
        win.grab_set()

        frame = ttk.Frame(win, padding=12)
        frame.pack(fill="both", expand=True)
        frame.rowconfigure(2, weight=1)
        frame.columnconfigure(0, weight=1)

        intro = (
            "Add one or more Human Protein Atlas URLs. General /cancer URLs will be expanded "
            "to available cancer types during Preview. You can mix genes and sections."
        )
        ttk.Label(frame, text=intro, wraplength=860).grid(row=0, column=0, sticky="ew", pady=(0, 8))

        add_box = ttk.LabelFrame(frame, text="Add URLs")
        add_box.grid(row=1, column=0, sticky="ew", pady=(0, 10))
        add_box.columnconfigure(0, weight=1)

        self.bulk_entry_var = tk.StringVar()
        entry = ttk.Entry(add_box, textvariable=self.bulk_entry_var)
        entry.grid(row=0, column=0, sticky="ew", padx=8, pady=8)

        url_list_var = tk.Variable(value=[])

        def get_list():
            return list(url_list_var.get())

        def set_list(values):
            # Keep order but remove duplicates and empty strings.
            seen = set()
            clean = []
            for u in values:
                u = (u or "").strip()
                if not u or u in seen:
                    continue
                seen.add(u)
                clean.append(u)
            url_list_var.set(clean)
            count_var.set(f"{len(clean)} URL(s) ready")

        def add_urls_from_text(raw):
            found = parse_bulk_urls(raw)
            if not found:
                messagebox.showwarning("No valid URLs", "No Human Protein Atlas URLs were detected.", parent=win)
                return
            set_list(get_list() + found)

        def add_entry_url():
            raw = self.bulk_entry_var.get().strip()
            if not raw:
                return
            add_urls_from_text(raw)
            self.bulk_entry_var.set("")
            entry.focus_set()

        ttk.Button(add_box, text="Add", command=add_entry_url).grid(row=0, column=1, padx=(0, 8), pady=8)

        paste_box = ttk.Frame(add_box)
        paste_box.grid(row=1, column=0, columnspan=2, sticky="ew", padx=8, pady=(0, 8))
        paste_box.columnconfigure(0, weight=1)
        paste_txt = tk.Text(paste_box, height=5, wrap="word")
        paste_txt.grid(row=0, column=0, sticky="ew")

        paste_buttons = ttk.Frame(paste_box)
        paste_buttons.grid(row=0, column=1, sticky="ns", padx=(8, 0))
        ttk.Button(paste_buttons, text="Add pasted URLs", command=lambda: (add_urls_from_text(paste_txt.get("1.0", "end")), paste_txt.delete("1.0", "end"))).pack(fill="x", pady=(0, 4))

        def paste_from_clipboard():
            try:
                clip = win.clipboard_get()
            except Exception:
                clip = ""
            if clip:
                paste_txt.insert("end", clip)
        ttk.Button(paste_buttons, text="Paste clipboard", command=paste_from_clipboard).pack(fill="x")

        list_box = ttk.LabelFrame(frame, text="URLs to preview")
        list_box.grid(row=2, column=0, sticky="nsew")
        list_box.rowconfigure(0, weight=1)
        list_box.columnconfigure(0, weight=1)

        lb = tk.Listbox(list_box, listvariable=url_list_var, selectmode="extended", activestyle="dotbox")
        lb.grid(row=0, column=0, sticky="nsew", padx=(8, 0), pady=8)
        yscroll = ttk.Scrollbar(list_box, orient="vertical", command=lb.yview)
        yscroll.grid(row=0, column=1, sticky="ns", pady=8)
        xscroll = ttk.Scrollbar(list_box, orient="horizontal", command=lb.xview)
        xscroll.grid(row=1, column=0, sticky="ew", padx=(8, 0), pady=(0, 8))
        lb.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)

        list_actions = ttk.Frame(list_box)
        list_actions.grid(row=0, column=2, sticky="ns", padx=8, pady=8)

        def remove_selected():
            current = get_list()
            selected = set(lb.curselection())
            set_list([u for i, u in enumerate(current) if i not in selected])

        def clear_all_urls():
            set_list([])

        def load_urls_from_file():
            path = filedialog.askopenfilename(
                title="Load URL list",
                filetypes=[("Text or CSV files", "*.txt *.csv"), ("All files", "*.*")],
                parent=win,
            )
            if not path:
                return
            try:
                raw = Path(path).read_text(encoding="utf-8", errors="ignore")
                add_urls_from_text(raw)
            except Exception as e:
                messagebox.showerror("Could not read file", str(e), parent=win)

        ttk.Button(list_actions, text="Remove selected", command=remove_selected).pack(fill="x", pady=(0, 4))
        ttk.Button(list_actions, text="Clear list", command=clear_all_urls).pack(fill="x", pady=(0, 4))
        ttk.Button(list_actions, text="Load .txt/.csv", command=load_urls_from_file).pack(fill="x")

        bottom = ttk.Frame(frame)
        bottom.grid(row=3, column=0, sticky="ew", pady=(10, 0))
        bottom.columnconfigure(0, weight=1)
        count_var = tk.StringVar(value="0 URL(s) ready")
        ttk.Label(bottom, textvariable=count_var).grid(row=0, column=0, sticky="w")

        def do_preview():
            urls = get_list()
            if not urls:
                # Allow direct use of the single URL box or paste box without pressing Add first.
                candidate_text = "\n".join([self.bulk_entry_var.get(), paste_txt.get("1.0", "end")])
                urls = parse_bulk_urls(candidate_text)
            if not urls:
                messagebox.showwarning("No valid URLs", "Add at least one valid Human Protein Atlas URL.", parent=win)
                return
            bad = [u for u in urls if not is_valid_hpa_url(u)]
            if bad:
                messagebox.showwarning("Invalid URLs", "Some URLs do not appear to be Human Protein Atlas URLs.", parent=win)
                return
            win.destroy()
            self.on_bulk_preview(urls)

        ttk.Button(bottom, text="Cancel", command=win.destroy).grid(row=0, column=1, sticky="e", padx=(8, 0))
        ttk.Button(bottom, text="Preview URLs", command=do_preview).grid(row=0, column=2, sticky="e", padx=(8, 0))

        # Preload current URL(s), if available.
        current = self.preview_urls or ([self.url_var.get().strip()] if self.url_var.get().strip() else [])
        set_list(current)
        entry.bind("<Return>", lambda e: add_entry_url())
        entry.focus_set()

    def select_output_directory(self):
        folder = filedialog.askdirectory(title="Select output directory", initialdir=str(self.output_dir))
        if folder:
            self.output_dir = Path(folder)
            if hasattr(self, "output_var"):
                self.output_var.set(str(self.output_dir))
            self.log(f"Output directory set to: {self.output_dir}")
            self.set_status("Output directory updated.")

    def open_output_directory(self):
        """Open the current output directory in the operating system file browser."""
        try:
            self.output_dir.mkdir(parents=True, exist_ok=True)
            if sys.platform.startswith("win"):
                import os
                os.startfile(str(self.output_dir))
            elif sys.platform == "darwin":
                import subprocess
                subprocess.Popen(["open", str(self.output_dir)])
            else:
                import subprocess
                subprocess.Popen(["xdg-open", str(self.output_dir)])
        except Exception as e:
            messagebox.showerror("Could not open output folder", str(e), parent=self)

    def _on_mousewheel(self, event):
        """Scroll the selection panel with the mouse wheel when the pointer is over the app."""
        try:
            widget = self.winfo_containing(event.x_root, event.y_root)
            if widget is not None and str(widget).startswith(str(self.sel_canvas)):
                self.sel_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        except Exception:
            pass

    def log(self, msg):
        self.log_text.configure(state="normal")
        self.log_text.insert("end", msg + "\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def set_status(self, s):
        self.status_var.set(s)

    def show_help(self):
        """Show simple usage instructions and output-structure information."""
        help_win = tk.Toplevel(self)
        help_win.title("Help - User Guide")
        help_win.geometry("760x620")
        help_win.minsize(680, 520)
        help_win.transient(self)
        help_win.grab_set()

        frame = ttk.Frame(help_win, padding=16)
        frame.pack(fill="both", expand=True)

        ttk.Label(
            frame,
            text="HPA IHC Image Downloader - User Guide",
            font=("Segoe UI", 12, "bold")
        ).pack(anchor="w", pady=(0, 10))

        txt = tk.Text(frame, wrap="word", height=26)
        txt.pack(fill="both", expand=True)
        txt.insert("end", """Purpose
This tool previews, selects, downloads, organizes, and documents immunohistochemistry (IHC) images from Human Protein Atlas (HPA) cancer and normal tissue pages.

Basic use
1. Paste one HPA gene URL in the URL box.
2. Choose the image format: .tif for maximum quality or .jpg for smaller files.
3. Click Preview.
4. Review the available marker, antibody, HPA category, and diagnostic subtype information.
5. Select the antibodies, categories, and cancer diagnostic subtypes you want to download.
6. Click Download selected.
7. Use File > Open Output Directory to inspect the downloaded files and reports.

Bulk use
1. Open File > Bulk URLs.
2. Add or paste several HPA URLs, or load a .txt/.csv file.
3. Click Preview URLs.
4. Select the desired items and download.

Output folder
Use File > Select Output Directory to change where files are saved.
Use File > Open Output Directory to open the current output folder.

Output organization
Downloaded images are stored using a fixed biomarker-centered structure:

Marker / Antibody / HPA category / HPA diagnostic subtype when available

Example:
FGFR2 / CAB010886 / lung cancer / Adenocarcinoma / ID_1506_1.tif

For normal tissue pages, HPA diagnostic subtype selection is not shown because diagnostic subtypes are specific to cancer pages.

Reports and reproducibility files
The tool exports structured CSV files and reports, including:
- preview_inventory.csv
- downloaded_successfully.csv
- failed_downloads.csv, only when failures occur
- download_manifest.json
- download_report.md
- citation_and_methods_helper.md

SHA256
SHA256 is a file integrity checksum. If the image content changes, the SHA256 value changes. This helps verify that downloaded files are unchanged and reproducible.
""")
        txt.configure(state="disabled")

        ttk.Button(frame, text="Close", command=help_win.destroy).pack(anchor="e", pady=(10, 0))


    def show_about(self):
        about_win = tk.Toplevel(self)
        about_win.title("About")
        about_win.resizable(False, False)
        about_win.transient(self)
        about_win.grab_set()
        frame = ttk.Frame(about_win, padding=16)
        frame.pack(fill="both", expand=True)
        ttk.Label(frame, text="HPA IHC Image Downloader", font=("Segoe UI", 12, "bold")).pack(anchor="w", pady=(0, 10))
        ttk.Label(frame, text=f"Version: {__version__}").pack(anchor="w")
        ttk.Label(frame, text=f"Author: {__author__}").pack(anchor="w", pady=(0, 10))
        ttk.Label(frame, text="GitHub:", font=("Segoe UI", 9, "bold")).pack(anchor="w")
        github_link = ttk.Label(frame, text=__github_url__, foreground="blue", cursor="hand2")
        github_link.pack(anchor="w", pady=(0, 8))
        github_link.bind("<Button-1>", lambda e: webbrowser.open(__github_url__))
        ttk.Label(frame, text="DOI:", font=("Segoe UI", 9, "bold")).pack(anchor="w")
        doi_link = ttk.Label(frame, text=__doi__, foreground="blue", cursor="hand2")
        doi_link.pack(anchor="w", pady=(0, 12))
        doi_link.bind("<Button-1>", lambda e: webbrowser.open(__doi__))
        ttk.Button(frame, text="Close", command=about_win.destroy).pack(anchor="e")

    def clear_preview(self):
        self.inv = None
        self.gene_name = None
        self.bulk_inv = None
        self.total_items = 0
        self.preview_img_ext = None
        self.preview_urls = []
        self._clear_tree()
        self._clear_selection_panel()
        self.progress.configure(value=0, maximum=1)
        self.download_btn.configure(state="disabled", text="Download")
        self.set_status("Preview cleared.")
        self.log("Preview cleared.")

    def _clear_tree(self):
        for item in self.tree.get_children():
            self.tree.delete(item)

    def _clear_selection_panel(self):
        for child in self.sel_frame.winfo_children():
            child.destroy()
        self.gene_vars.clear()
        self.ab_vars.clear()
        self.category_vars.clear()
        self.subtype_vars.clear()

    def _diagnostic_counts_for_payload(self, payload):
        """Return subtype counts for a category using only HPA-provided labels."""
        counts = {}
        direct = 0
        for row in payload.get("items", []):
            subtype = valid_hpa_diagnostic_subtype(
                row.get("DiagnosticCategory") or row.get("CancerType") or row.get("Diagnosis") or "",
                section_category=row.get("SectionCategory") or row.get("Category") or "",
            )
            if subtype:
                counts[subtype] = counts.get(subtype, 0) + 1
            else:
                direct += 1
        return counts, direct

    def _populate_tree_and_selection(self):
        self._clear_tree()
        self._clear_selection_panel()
        if not self.bulk_inv:
            return

        n_genes = len(self.bulk_inv)
        ttk.Label(self.sel_frame, text=f"Preview: {n_genes} gene(s), {self.total_items} item(s)").pack(anchor="w", pady=(0, 8))
        util = ttk.Frame(self.sel_frame)
        util.pack(fill="x", pady=(0, 8))
        ttk.Button(util, text="Select all", command=self._select_all).pack(side="left")
        ttk.Button(util, text="Clear all", command=self._select_none).pack(side="left", padx=6)

        root_id = self.tree.insert("", "end", text=f"Preview: {n_genes} gene(s)", values=(self.total_items,))
        self.tree.item(root_id, open=True)

        for gene, inv in sorted(self.bulk_inv.items(), key=lambda kv: kv[0].lower()):
            gene_count = sum(payload["count"] for cats in inv.values() for payload in cats.values())
            gene_node = self.tree.insert(root_id, "end", text=f"Gene: {gene}", values=(gene_count,))
            self.tree.item(gene_node, open=(n_genes == 1))

            gene_var = tk.BooleanVar(value=True)
            self.gene_vars[gene] = gene_var
            gene_row = ttk.Frame(self.sel_frame)
            gene_row.pack(fill="x", pady=(8, 2))
            ttk.Checkbutton(gene_row, text=f"{gene} ({gene_count})", variable=gene_var, command=lambda g=gene: self._toggle_gene(g)).pack(anchor="w")

            gene_box = ttk.Frame(self.sel_frame, padding=(14, 0, 0, 0))
            gene_box.pack(fill="x")

            for ab_id, categories in sorted(inv.items(), key=lambda kv: kv[0].lower()):
                ab_count = sum(v["count"] for v in categories.values())
                ab_node = self.tree.insert(gene_node, "end", text=f"Antibody: {ab_id}", values=(ab_count,))
                self.tree.item(ab_node, open=False)

                ab_var = tk.BooleanVar(value=True)
                self.ab_vars[(gene, ab_id)] = ab_var
                ttk.Checkbutton(
                    gene_box,
                    text=f"{ab_id} ({ab_count})",
                    variable=ab_var,
                    command=lambda g=gene, a=ab_id: self._toggle_antibody(g, a)
                ).pack(anchor="w", pady=(4, 1))

                cat_box = ttk.Frame(gene_box, padding=(18, 0, 0, 0))
                cat_box.pack(fill="x")
                for cat, payload in sorted(categories.items(), key=lambda kv: kv[0].lower()):
                    c_count = payload["count"]
                    cat_node = self.tree.insert(ab_node, "end", text=cat, values=(c_count,))

                    # Show cancer diagnostic/histology subtypes under the HPA cancer
                    # category for transparency. These subtype labels are not invented;
                    # they come from HPA metadata/table labels. If missing, we show a
                    # small direct-to-category count instead of creating an Unknown folder.
                    subtype_counts, direct_count = self._diagnostic_counts_for_payload(payload)
                    if subtype_counts:
                        self.tree.item(cat_node, open=False)
                        for subtype, s_count in sorted(subtype_counts.items(), key=lambda kv: kv[0].lower()):
                            self.tree.insert(cat_node, "end", text=f"Subtype: {subtype}", values=(s_count,))
                    if direct_count and subtype_counts:
                        self.tree.insert(cat_node, "end", text="No HPA diagnostic subtype provided", values=(direct_count,))

                    c_var = tk.BooleanVar(value=True)
                    self.category_vars[(gene, ab_id, cat)] = c_var
                    if subtype_counts:
                        subtype_label = f"{cat} ({c_count}; {len(subtype_counts)} subtype(s))"
                    else:
                        subtype_label = f"{cat} ({c_count})"
                    ttk.Checkbutton(
                        cat_box,
                        text=subtype_label,
                        variable=c_var,
                        command=lambda g=gene, a=ab_id, c=cat: self._toggle_category(g, a, c),
                    ).pack(anchor="w")

                    # Let users select only specific HPA diagnostic/histology
                    # subtypes while keeping the broader cancer category visible.
                    # This is intentionally only based on labels parsed from HPA;
                    # no synthetic subtype folders are created.
                    is_cancer_payload = any(
                        (row.get("HPASection") or "").lower() == "cancer"
                        for row in payload.get("items", [])
                    )

                    # Show selectable subtype checkboxes only for cancer pages.
                    # Normal tissue categories do not have HPA diagnostic subtypes,
                    # so showing "No HPA diagnostic subtype provided" there would
                    # be redundant and confusing.
                    show_subtype_controls = is_cancer_payload and bool(subtype_counts)
                    if show_subtype_controls:
                        subtype_box = ttk.Frame(cat_box, padding=(18, 0, 0, 0))
                        subtype_box.pack(fill="x")
                        for subtype, s_count in sorted(subtype_counts.items(), key=lambda kv: kv[0].lower()):
                            st_var = tk.BooleanVar(value=True)
                            self.subtype_vars[(gene, ab_id, cat, subtype)] = st_var
                            ttk.Checkbutton(
                                subtype_box,
                                text=f"↳ {subtype} ({s_count})",
                                variable=st_var,
                            ).pack(anchor="w")
                        if direct_count:
                            st_var = tk.BooleanVar(value=True)
                            self.subtype_vars[(gene, ab_id, cat, NO_SUBTYPE_KEY)] = st_var
                            ttk.Checkbutton(
                                subtype_box,
                                text=f"↳ {NO_SUBTYPE_LABEL} ({direct_count})",
                                variable=st_var,
                            ).pack(anchor="w")

        self.tree.see(root_id)

    def _toggle_gene(self, gene):
        state = self.gene_vars[gene].get()
        for (g, ab), var in self.ab_vars.items():
            if g == gene:
                var.set(state)
        for (g, ab, cat), var in self.category_vars.items():
            if g == gene:
                var.set(state)
        for (g, ab, cat, subtype), var in self.subtype_vars.items():
            if g == gene:
                var.set(state)

    def _toggle_antibody(self, gene, antibody_id):
        state = self.ab_vars[(gene, antibody_id)].get()
        for (g, ab, cat), var in self.category_vars.items():
            if g == gene and ab == antibody_id:
                var.set(state)
        for (g, ab, cat, subtype), var in self.subtype_vars.items():
            if g == gene and ab == antibody_id:
                var.set(state)

    def _toggle_category(self, gene, antibody_id, category):
        """Toggle all diagnostic subtype checkboxes under one category."""
        state = self.category_vars[(gene, antibody_id, category)].get()
        for (g, ab, cat, subtype), var in self.subtype_vars.items():
            if g == gene and ab == antibody_id and cat == category:
                var.set(state)

    def _select_all(self):
        for v in self.gene_vars.values():
            v.set(True)
        for v in self.ab_vars.values():
            v.set(True)
        for v in self.category_vars.values():
            v.set(True)
        for v in self.subtype_vars.values():
            v.set(True)

    def _select_none(self):
        for v in self.gene_vars.values():
            v.set(False)
        for v in self.ab_vars.values():
            v.set(False)
        for v in self.category_vars.values():
            v.set(False)
        for v in self.subtype_vars.values():
            v.set(False)

    def _get_selection(self):
        categories_by_gene_ab = {}
        subtypes_by_gene_ab_category = {}

        for (gene, ab_id), ab_var in self.ab_vars.items():
            selected_categories = set()

            for (g, ab, cat), var in self.category_vars.items():
                if g != gene or ab != ab_id:
                    continue

                selected_subtypes = {
                    subtype
                    for (sg, sab, scat, subtype), st_var in self.subtype_vars.items()
                    if sg == gene and sab == ab_id and scat == cat and st_var.get()
                }

                # Include a category when either the category checkbox is on or
                # at least one subtype under it is selected. The subtype filter
                # is applied only when subtype checkboxes exist for that category.
                has_subtype_controls = any(
                    sg == gene and sab == ab_id and scat == cat
                    for (sg, sab, scat, subtype) in self.subtype_vars.keys()
                )

                if var.get() or selected_subtypes:
                    selected_categories.add(cat)
                    if has_subtype_controls:
                        subtypes_by_gene_ab_category[(gene, ab_id, cat)] = selected_subtypes

            if selected_categories:
                categories_by_gene_ab[(gene, ab_id)] = selected_categories

        return {
            "categories_by_gene_antibody": categories_by_gene_ab,
            "subtypes_by_gene_antibody_category": subtypes_by_gene_ab_category,
        }

    def _workers_alive(self):
        """Return True if a background preview/download thread is still running."""
        self.worker_threads = [t for t in self.worker_threads if t.is_alive()]
        return bool(self.worker_threads)

    def _start_thread(self, target):
        """Start a non-daemon worker so Tcl/Tk is not destroyed while Python threads still exist."""
        t = threading.Thread(target=target, daemon=False)
        self.worker_threads.append(t)
        t.start()
        return t

    def on_close(self):
        """Close safely without deleting Tcl handlers from the wrong thread."""
        if self._workers_alive():
            messagebox.showinfo(
                "Process still running",
                "A preview or download is still running. Please wait until it finishes before closing the app.",
                parent=self,
            )
            return
        self._closing = True
        if self._poll_after_id is not None:
            try:
                self.after_cancel(self._poll_after_id)
            except Exception:
                pass
            self._poll_after_id = None
        try:
            self.grab_release()
        except Exception:
            pass
        self.destroy()

    def _start_preview_worker(self, urls):
        self.preview_btn.configure(state="disabled")
        self.bulk_btn.configure(state="disabled")
        self.download_btn.configure(state="disabled")
        self.set_status("Building preview...")
        self.progress.configure(value=0, maximum=1)
        self.log("=== PREVIEW ===")
        self.log(f"URLs: {len(urls)}")
        for u in urls:
            self.log(f"- {u}")
        self.log(f"Format: {self.ext_var.get()}")

        def log_cb(m):
            self.msg_q.put(("log", m))

        def worker():
            try:
                bulk_inv, total = build_bulk_preview_inventory(urls, img_ext=self.ext_var.get(), log=log_cb)
                self.msg_q.put(("preview_ok", bulk_inv, total, urls))
            except Exception as e:
                self.msg_q.put(("error", f"Preview failed: {e}"))

        self._start_thread(worker)

    def on_preview(self):
        url = self.url_var.get().strip()
        if not url:
            messagebox.showwarning("Missing URL", "Paste a Human Protein Atlas URL.")
            return
        if not is_valid_hpa_url(url):
            messagebox.showwarning("Invalid URL", "This does not appear to be a valid Human Protein Atlas URL.")
            return
        self._start_preview_worker([url])

    def on_bulk_preview(self, urls):
        bad = [u for u in urls if not is_valid_hpa_url(u)]
        if bad:
            messagebox.showwarning("Invalid URLs", "Some URLs do not appear to be Human Protein Atlas URLs.")
            return
        self._start_preview_worker(urls)

    def on_download(self):
        if not self.bulk_inv:
            messagebox.showinfo("No preview available", "Run Preview first to build the selection tree.")
            return
        if self.preview_img_ext != self.ext_var.get():
            messagebox.showwarning("Format changed", "The output format was changed after preview. Please run Preview again.")
            return

        selection = self._get_selection()
        sel_rows = selected_rows_from_bulk_inventory(self.bulk_inv, selection)
        if not sel_rows:
            messagebox.showwarning("Nothing selected", "Select at least one antibody/category to download.")
            return

        self.preview_btn.configure(state="disabled")
        self.bulk_btn.configure(state="disabled")
        self.download_btn.configure(state="disabled")
        self.set_status("Downloading...")
        self.log("=== DOWNLOAD ===")
        self.progress.configure(value=0, maximum=max(len(sel_rows), 1))

        def progress_cb(i, total):
            self.msg_q.put(("progress", i, total))

        def log_cb(m):
            self.msg_q.put(("log", m))

        def worker():
            try:
                download_from_bulk_inventory(
                    self.output_dir,
                    self.bulk_inv,
                    self.ext_var.get(),
                    selection,
                    progress_cb=progress_cb,
                    log=log_cb,
                    input_urls=self.preview_urls,
                )
                self.msg_q.put(("download_ok",))
            except Exception as e:
                self.msg_q.put(("error", f"Download failed: {e}"))

        self._start_thread(worker)

    def _poll_queue(self):
        if self._closing or not self.winfo_exists():
            return
        try:
            while True:
                msg = self.msg_q.get_nowait()
                kind = msg[0]
                if kind == "log":
                    self.log(msg[1])
                elif kind == "progress":
                    i, total = msg[1], msg[2]
                    self.progress.configure(value=i, maximum=total)
                    self.set_status(f"Downloading... {i}/{total}")
                elif kind == "preview_ok":
                    _, bulk_inv, total, urls = msg
                    self.bulk_inv = bulk_inv
                    self.total_items = total
                    self.preview_urls = urls
                    self.preview_img_ext = self.ext_var.get()
                    # Backward compatibility for old single-gene state.
                    if len(bulk_inv) == 1:
                        self.gene_name = next(iter(bulk_inv.keys()))
                        self.inv = bulk_inv[self.gene_name]
                        self.url_var.set(urls[0])
                    else:
                        self.gene_name = None
                        self.inv = None
                    self._populate_tree_and_selection()
                    try:
                        write_preview_outputs(self.output_dir, self.bulk_inv, input_urls=self.preview_urls, img_ext=self.preview_img_ext, log=self.log)
                    except Exception as e:
                        self.log(f"⚠ Could not export preview inventory: {e}")
                    self.set_status(f"Preview ready: {len(bulk_inv)} gene(s), {total} items")
                    self.log(f"Preview OK: {len(bulk_inv)} gene(s) | {total} items")
                    self.preview_btn.configure(state="normal")
                    self.bulk_btn.configure(state="normal")
                    self.download_btn.configure(state="normal", text="Download selected")
                elif kind == "download_ok":
                    self.set_status("Download completed.")
                    self.log("✅ Download completed.")
                    self.preview_btn.configure(state="normal")
                    self.bulk_btn.configure(state="normal")
                    self.download_btn.configure(state="normal", text="Download selected")
                elif kind == "error":
                    self.set_status("Error.")
                    self.log("❌ " + msg[1])
                    messagebox.showerror("Error", msg[1])
                    self.preview_btn.configure(state="normal")
                    self.bulk_btn.configure(state="normal")
                    self.download_btn.configure(state="normal" if self.bulk_inv else "disabled", text="Download selected" if self.bulk_inv else "Download")
        except queue.Empty:
            pass
        if not self._closing and self.winfo_exists():
            self._poll_after_id = self.after(100, self._poll_queue)



# ---------------------------------------------------------------------
# Command-line interface for reproducible workflows
# ---------------------------------------------------------------------

def read_url_file(path):
    """Read URLs from a text/CSV file using the same parser as the GUI."""
    raw = Path(path).read_text(encoding="utf-8", errors="ignore")
    return parse_bulk_urls(raw)


def run_cli(argv=None):
    """Run the downloader from command line without opening the GUI."""
    parser = argparse.ArgumentParser(
        description="Preview and download HPA IHC images with reproducible metadata outputs."
    )
    parser.add_argument("--urls", nargs="*", default=[], help="One or more HPA URLs.")
    parser.add_argument("--url-file", help="Text/CSV file containing HPA URLs, one per line or mixed with labels.")
    parser.add_argument("--output", default=str(ROOT_DIR), help="Output directory.")
    parser.add_argument("--format", choices=[".tif", ".jpg", "tif", "jpg"], default=".tif", help="Image format to request.")
    parser.add_argument("--preview-only", action="store_true", help="Only export preview inventory; do not download images.")
    parser.add_argument("--download", action="store_true", help="Download all previewed items. If omitted, --preview-only controls behavior.")
    args = parser.parse_args(argv)

    urls = []
    urls.extend(parse_bulk_urls("\n".join(args.urls or [])))
    if args.url_file:
        urls.extend(read_url_file(args.url_file))
    # De-duplicate while preserving order.
    seen = set()
    urls = [u for u in urls if not (u in seen or seen.add(u))]

    if not urls:
        raise SystemExit("No valid Human Protein Atlas URLs were provided. Use --urls or --url-file.")

    img_ext = args.format
    if not img_ext.startswith("."):
        img_ext = "." + img_ext

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    def cli_log(msg):
        print(msg, flush=True)

    cli_log(f"HPA IHC Image Downloader v{__version__}")
    cli_log(f"Output directory: {output_dir}")
    cli_log(f"URLs: {len(urls)}")
    bulk_inv, total = build_bulk_preview_inventory(urls, img_ext=img_ext, log=cli_log)
    write_preview_outputs(output_dir, bulk_inv, input_urls=urls, img_ext=img_ext, log=cli_log)

    if args.preview_only and not args.download:
        cli_log("Preview-only mode completed.")
        return 0

    # CLI default when --preview-only is not used: download everything found.
    selection = {"categories_by_gene_antibody": {}}
    for gene, inv in bulk_inv.items():
        for ab_id, categories in inv.items():
            selection["categories_by_gene_antibody"][(gene, ab_id)] = set(categories.keys())

    def progress(i, n):
        cli_log(f"Progress: {i}/{n}")

    download_from_bulk_inventory(
        output_dir,
        bulk_inv,
        img_ext,
        selection,
        progress_cb=progress,
        log=cli_log,
        input_urls=urls,
    )
    cli_log("CLI run completed.")
    return 0

if __name__ == "__main__":
    # If command-line arguments are provided, run in reproducible CLI mode.
    # Double-clicking or running without arguments still opens the GUI.
    if len(sys.argv) > 1:
        raise SystemExit(run_cli(sys.argv[1:]))

    try:
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass
    app = App()
    app.mainloop()