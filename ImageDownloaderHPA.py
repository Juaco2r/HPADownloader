# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 23:25:37 2025

@author: juaco
"""


import os
import re
import csv
import time
import threading
import queue
import requests
from bs4 import BeautifulSoup
from html import unescape
from pathlib import Path
import tkinter as tk
from tkinter import ttk, messagebox

BASE_IMG_HOST = "https://images.proteinatlas.org"
USER_DOWNLOAD_DIR = Path.home() / "Downloads"
ROOT_DIR = USER_DOWNLOAD_DIR / "HPA Images"


# -------------------------
# Tu lógica (con mínimos cambios)
# -------------------------

def ensure_root():
    ROOT_DIR.mkdir(parents=True, exist_ok=True)

def download_html(url):
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.text

def extract_gene_name(soup):
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
    text = soup.get_text(" ", strip=True)
    ids = sorted(set(re.findall(r"\b(HPA\d{6}|CAB\d{6})\b", text)))
    return ids

def clean_html_title(html_title):
    txt = unescape(html_title)
    txt = re.sub(r"</?b>", "", txt)
    txt = re.sub(r"<br\s*/?>", "\n", txt)
    return txt.strip()

def parse_metadata_from_title(title_text):
    meta = {
        "Gender": None,
        "Age": None,
        "PatientID": None,
        "LocationsRaw": [],
        "CancerType": None,
        "CancerCode": None,
        "AntibodyStaining": None,
        "Intensity": None,
        "Quantity": None,
        "Location": None,
    }

    lines = [l.strip() for l in title_text.split("\n") if l.strip()]

    if lines:
        m = re.match(r"(Male|Female)\s*,\s*age\s*(\d+)", lines[0])
        if m:
            meta["Gender"] = m.group(1)
            meta["Age"] = m.group(2)

    for line in lines[1:]:
        m = re.match(r"Patient id:\s*(\d+)", line, re.IGNORECASE)
        if m:
            meta["PatientID"] = m.group(1)
            continue

        handled = False
        for key, field in [
            ("Antibody staining:", "AntibodyStaining"),
            ("Intensity:", "Intensity"),
            ("Quantity:", "Quantity"),
            ("Location:", "Location"),
        ]:
            if line.startswith(key):
                meta[field] = line.split(":", 1)[1].strip()
                handled = True
                break
        if handled:
            continue

        code_match = re.search(r"\((M-\d+|T-\w+)\)", line)
        if code_match:
            code = code_match.group(1)
            if code.startswith("M-"):
                meta["CancerType"] = line.rsplit("(", 1)[0].strip()
                meta["CancerCode"] = code
            else:
                meta["LocationsRaw"].append(line.strip())

    return meta

def normalize_cancer_folder(cancer_type):
    if cancer_type is None:
        return "Unknown"
    ct = cancer_type.replace("NOS", "").strip()
    ct = ct.replace(",", "")
    return ct or "Unknown"

def download_image_with_retry(image_link, image_path, max_retries=3, log=None):
    for attempt in range(max_retries):
        try:
            if log: log(f"Descargando: {image_path.name}")
            r = requests.get(image_link, timeout=20)
            r.raise_for_status()
            with open(image_path, "wb") as out:
                out.write(r.content)
            if log: log(f"✓ OK: {image_path.name}")
            return True
        except Exception as e:
            if log: log(f"✗ Error intento {attempt+1}: {e}")
            if attempt < max_retries - 1:
                time.sleep(1)
            else:
                if log: log(f"✗ Falló definitivamente: {image_link}")
                return False


# -------------------------
# Preview: construir inventario (sin descargar)
# -------------------------

def build_preview_inventory(url, img_ext=".jpg"):
    """
    Retorna:
      gene_name: str
      inv: dict[antibody_id][cancer_folder] = {
           'count': int,
           'items': list of dicts (metadata + image_link + image_name sugerido)
      }
      total_items: int
    """
    html = download_html(url)
    soup = BeautifulSoup(html, "html.parser")

    gene_name = extract_gene_name(soup)
    antibody_ids = extract_antibody_ids(soup)
    if not antibody_ids:
        raise RuntimeError("No se encontraron IDs de anticuerpo en la página.")

    rows = soup.find_all("tr")

    inv = {}
    total = 0

    for antibody_id in antibody_ids:
        antibody_number = re.sub(r'^[HCAPB]*0*', '', antibody_id)

        inv.setdefault(antibody_id, {})
        image_counters = {}  # para sugerir nombres consistentes

        for tr in rows:
            ths = tr.find_all("th")
            tds = tr.find_all("td")
            if len(ths) < 1 or len(tds) < 1:
                continue

            for col_idx, th in enumerate(ths):
                cancer_label = th.get_text(" ", strip=True)
                if not cancer_label:
                    continue
                if col_idx >= len(tds):
                    continue

                td = tds[col_idx]
                cancer_divs = td.find_all("div", class_="cancerAnnoations")
                if not cancer_divs:
                    continue

                for div in cancer_divs:
                    a_tags = div.find_all("a", href=True)
                    for a in a_tags:
                        img = a.find("img")
                        if not img or "title" not in img.attrs:
                            continue

                        href = a["href"]
                        if antibody_number not in href:
                            continue

                        if href.startswith("//"):
                            href = "https:" + href
                        elif href.startswith("/"):
                            href = BASE_IMG_HOST + href.replace("/images_static", "")

                        image_link = href
                        if img_ext == ".tif":
                            image_link = image_link.replace(".jpg", ".tif")

                        title_text = clean_html_title(img["title"])
                        meta = parse_metadata_from_title(title_text)

                        cancer_type = meta["CancerType"] if meta["CancerType"] else cancer_label
                        cancer_folder = normalize_cancer_folder(cancer_type)

                        patient_id = meta["PatientID"] if meta["PatientID"] else "unknown"
                        counter_key = (cancer_folder, patient_id)
                        idx = image_counters.get(counter_key, 0) + 1
                        image_counters[counter_key] = idx

                        image_name = f"ID_{patient_id}_{idx}{img_ext}"

                        inv[antibody_id].setdefault(cancer_folder, {"count": 0, "items": []})
                        inv[antibody_id][cancer_folder]["count"] += 1
                        inv[antibody_id][cancer_folder]["items"].append({
                            "Gene": gene_name,
                            "AntibodyID": antibody_id,
                            "CancerFolder": cancer_folder,
                            "ImageName": image_name,
                            "ImageLink": image_link,
                            "PatientID": patient_id,
                            "Gender": meta["Gender"],
                            "Age": meta["Age"],
                            "CancerType": cancer_type,
                            "CancerCode": meta["CancerCode"],
                            "LocationCodes": "; ".join(meta["LocationsRaw"]) if meta["LocationsRaw"] else "",
                            "AntibodyStaining": meta["AntibodyStaining"],
                            "Intensity": meta["Intensity"],
                            "Quantity": meta["Quantity"],
                            "Location": meta["Location"],
                        })
                        total += 1

    return gene_name, inv, total


# -------------------------
# Download: usa inventario y selección
# -------------------------

def download_from_inventory(gene_name, inv, img_ext, selection, progress_cb=None, log=None):
    """
    selection:
      {
        'antibodies': set([...]) or empty -> all
        'cancers_by_antibody': dict[antibody_id] = set([...]) or empty -> all per antibody
      }
    """
    ensure_root()
    gene_dir = ROOT_DIR / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)

    antibodies_sel = selection.get("antibodies", set())
    cancers_by_ab = selection.get("cancers_by_antibody", {})

    # construir lista final de items a descargar
    download_list = []
    for ab_id, cancers in inv.items():
        if antibodies_sel and ab_id not in antibodies_sel:
            continue

        cancers_sel = cancers_by_ab.get(ab_id, set())
        for cancer_folder, payload in cancers.items():
            if cancers_sel and cancer_folder not in cancers_sel:
                continue
            download_list.extend(payload["items"])

    total = len(download_list)
    if total == 0:
        raise RuntimeError("No hay nada seleccionado para descargar (o la página no tiene items).")

    if log:
        log(f"Total a descargar: {total}")

    # para CSV global por gen
    all_rows_for_gene = []

    for i, row in enumerate(download_list, start=1):
        ab_id = row["AntibodyID"]
        cancer_folder = row["CancerFolder"]

        antibody_dir = gene_dir / ab_id
        antibody_dir.mkdir(parents=True, exist_ok=True)

        cancer_dir = antibody_dir / cancer_folder
        cancer_dir.mkdir(parents=True, exist_ok=True)

        image_path = cancer_dir / row["ImageName"]
        if not image_path.exists():
            download_image_with_retry(row["ImageLink"], image_path, log=log)
        else:
            if log: log(f"📁 Ya existe: {row['ImageName']}")

        all_rows_for_gene.append(row)

        # progreso
        if progress_cb:
            progress_cb(i, total)

    # CSV resumen por gen
    if all_rows_for_gene:
        gene_summary_path = gene_dir / f"{gene_name}_all_antibodies_summary.csv"
        fieldnames = list(all_rows_for_gene[0].keys())
        with open(gene_summary_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows_for_gene)

    if log:
        log(f"🎉 Listo. Carpeta: {gene_dir}")


# -------------------------
# GUI
# -------------------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("HPA Cancer/IHC Downloader (ligero)")
        self.geometry("1100x650")
        self.minsize(980, 600)

        self.inv = None
        self.gene_name = None
        self.total_items = 0

        self.msg_q = queue.Queue()

        self._build_ui()
        self._poll_queue()

    def _build_ui(self):
        # ===== TOP AREA =====
        top = ttk.Frame(self, padding=10)
        top.pack(fill="x")
    
        # --- Row 1: URL ---
        url_row = ttk.Frame(top)
        url_row.pack(fill="x", pady=(0, 6))
    
        ttk.Label(url_row, text="URL:").pack(side="left")
    
        self.url_var = tk.StringVar()
        url_entry = ttk.Entry(url_row, textvariable=self.url_var)
        url_entry.pack(side="left", padx=8, fill="x", expand=True)
    
        # --- Row 2: Controls ---
        ctrl_row = ttk.Frame(top)
        ctrl_row.pack(fill="x")
    
        ttk.Label(ctrl_row, text="Formato:").pack(side="left")
    
        self.ext_var = tk.StringVar(value=".jpg")
        ttk.Combobox(
            ctrl_row,
            textvariable=self.ext_var,
            values=[".jpg", ".tif"],
            width=6,
            state="readonly"
        ).pack(side="left", padx=(4, 12))
    
        self.preview_btn = ttk.Button(
            ctrl_row, text="Preview", command=self.on_preview
        )
        self.preview_btn.pack(side="left", padx=4)
    
        self.download_btn = ttk.Button(
            ctrl_row,
            text="Descargar TODO / Selección",
            command=self.on_download,
            state="disabled"
        )
        self.download_btn.pack(side="left", padx=4)

        # Main split
        main = ttk.PanedWindow(self, orient="horizontal")
        main.pack(fill="both", expand=True, padx=10, pady=(0,10))

        # Left: Tree preview
        left = ttk.Frame(main, padding=(0,0,6,0))
        main.add(left, weight=3)

        ttk.Label(left, text="Preview (árbol)").pack(anchor="w")

        self.tree = ttk.Treeview(left, columns=("count",), show="tree headings")
        self.tree.heading("#0", text="Elemento")
        self.tree.heading("count", text="#")
        self.tree.column("count", width=70, anchor="e")
        self.tree.pack(fill="both", expand=True, pady=6)

        # Right: Selection panel
        right = ttk.Frame(main, padding=(6,0,0,0))
        main.add(right, weight=2)

        ttk.Label(right, text="Selección").pack(anchor="w")

        self.sel_canvas = tk.Canvas(right, highlightthickness=0)
        self.sel_scroll = ttk.Scrollbar(right, orient="vertical", command=self.sel_canvas.yview)
        self.sel_canvas.configure(yscrollcommand=self.sel_scroll.set)
        self.sel_scroll.pack(side="right", fill="y")
        self.sel_canvas.pack(side="left", fill="both", expand=True, pady=6)

        self.sel_frame = ttk.Frame(self.sel_canvas)
        self.sel_canvas.create_window((0,0), window=self.sel_frame, anchor="nw")
        self.sel_frame.bind("<Configure>", lambda e: self.sel_canvas.configure(scrollregion=self.sel_canvas.bbox("all")))

        self.ab_vars = {}         # antibody -> BooleanVar
        self.cancer_vars = {}     # (antibody, cancer) -> BooleanVar

        # Bottom: progress + log
        bottom = ttk.Frame(self, padding=10)
        bottom.pack(fill="x")

        self.progress = ttk.Progressbar(bottom, mode="determinate")
        self.progress.pack(fill="x")

        self.status_var = tk.StringVar(value="Listo.")
        ttk.Label(bottom, textvariable=self.status_var).pack(anchor="w", pady=(6,0))

        log_frame = ttk.Frame(self, padding=(10,0,10,10))
        log_frame.pack(fill="both", expand=False)
        ttk.Label(log_frame, text="Log").pack(anchor="w")
        self.log_text = tk.Text(log_frame, height=8, wrap="word")
        self.log_text.pack(fill="both", expand=True, pady=6)
        self.log_text.configure(state="disabled")

    def log(self, msg):
        self.log_text.configure(state="normal")
        self.log_text.insert("end", msg + "\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def set_status(self, s):
        self.status_var.set(s)

    def _clear_tree(self):
        for item in self.tree.get_children():
            self.tree.delete(item)

    def _clear_selection_panel(self):
        for child in self.sel_frame.winfo_children():
            child.destroy()
        self.ab_vars.clear()
        self.cancer_vars.clear()

    def _populate_tree_and_selection(self):
        self._clear_tree()
        self._clear_selection_panel()

        if not self.inv or not self.gene_name:
            return

        root_id = self.tree.insert("", "end", text=f"Gene: {self.gene_name}", values=(self.total_items,))
        self.tree.item(root_id, open=True)

        ttk.Label(self.sel_frame, text=f"Gene: {self.gene_name} ({self.total_items} items)").pack(anchor="w", pady=(0,8))

        # Botones utilitarios selección
        util = ttk.Frame(self.sel_frame)
        util.pack(fill="x", pady=(0,8))
        ttk.Button(util, text="Marcar todo", command=self._select_all).pack(side="left")
        ttk.Button(util, text="Desmarcar todo", command=self._select_none).pack(side="left", padx=6)

        for ab_id, cancers in self.inv.items():
            ab_count = sum(v["count"] for v in cancers.values())
            ab_node = self.tree.insert(root_id, "end", text=f"Antibody: {ab_id}", values=(ab_count,))
            self.tree.item(ab_node, open=False)

            # checkbox anticuerpo
            ab_var = tk.BooleanVar(value=True)
            self.ab_vars[ab_id] = ab_var

            ab_row = ttk.Frame(self.sel_frame)
            ab_row.pack(fill="x", pady=(6,2))
            ttk.Checkbutton(
                ab_row, text=f"{ab_id} ({ab_count})", variable=ab_var,
                command=lambda a=ab_id: self._toggle_antibody(a)
            ).pack(anchor="w")

            # cancers dentro
            cancers_box = ttk.Frame(self.sel_frame, padding=(18,0,0,0))
            cancers_box.pack(fill="x")

            for cancer_folder, payload in sorted(cancers.items(), key=lambda kv: kv[0].lower()):
                c_count = payload["count"]
                self.tree.insert(ab_node, "end", text=cancer_folder, values=(c_count,))

                c_var = tk.BooleanVar(value=True)
                self.cancer_vars[(ab_id, cancer_folder)] = c_var
                ttk.Checkbutton(cancers_box, text=f"{cancer_folder} ({c_count})", variable=c_var).pack(anchor="w")

        self.tree.see(root_id)

    def _toggle_antibody(self, antibody_id):
        # si desmarcas anticuerpo, desmarca sus cancers; si marcas, marca todo
        state = self.ab_vars[antibody_id].get()
        for (ab, cancer), var in self.cancer_vars.items():
            if ab == antibody_id:
                var.set(state)

    def _select_all(self):
        for v in self.ab_vars.values():
            v.set(True)
        for v in self.cancer_vars.values():
            v.set(True)

    def _select_none(self):
        for v in self.ab_vars.values():
            v.set(False)
        for v in self.cancer_vars.values():
            v.set(False)

    def _get_selection(self):
        antibodies = {ab for ab, v in self.ab_vars.items() if v.get()}
        cancers_by_ab = {}
        for (ab, cancer), v in self.cancer_vars.items():
            if v.get():
                cancers_by_ab.setdefault(ab, set()).add(cancer)
        return {"antibodies": antibodies, "cancers_by_antibody": cancers_by_ab}

    def on_preview(self):
        url = self.url_var.get().strip()
        if not url:
            messagebox.showwarning("Falta URL", "Pega una URL de Protein Atlas.")
            return

        self.preview_btn.configure(state="disabled")
        self.download_btn.configure(state="disabled")
        self.set_status("Haciendo preview...")
        self.progress.configure(value=0, maximum=1)
        self.log("=== PREVIEW ===")
        self.log(f"URL: {url}")
        self.log(f"Formato: {self.ext_var.get()}")

        def worker():
            try:
                gene, inv, total = build_preview_inventory(url, img_ext=self.ext_var.get())
                self.msg_q.put(("preview_ok", gene, inv, total))
            except Exception as e:
                self.msg_q.put(("error", f"Preview falló: {e}"))

        threading.Thread(target=worker, daemon=True).start()

    def on_download(self):
        if not self.inv or not self.gene_name:
            messagebox.showinfo("Sin preview", "Haz Preview primero para construir la selección.")
            return

        url = self.url_var.get().strip()
        if not url:
            messagebox.showwarning("Falta URL", "Pega una URL de Protein Atlas.")
            return

        selection = self._get_selection()

        self.preview_btn.configure(state="disabled")
        self.download_btn.configure(state="disabled")
        self.set_status("Descargando...")
        self.log("=== DOWNLOAD ===")

        # calcular total seleccionado para progreso
        sel_total = 0
        antibodies_sel = selection["antibodies"]
        cancers_by_ab = selection["cancers_by_antibody"]

        for ab_id, cancers in self.inv.items():
            if antibodies_sel and ab_id not in antibodies_sel:
                continue
            cancers_sel = cancers_by_ab.get(ab_id, set())
            for cancer_folder, payload in cancers.items():
                if cancers_sel and cancer_folder not in cancers_sel:
                    continue
                sel_total += payload["count"]

        if sel_total == 0:
            # si quedó todo desmarcado, asumimos "descarga todo"
            # (alternativa: avisar; yo lo hago amigable)
            self.log("Nada seleccionado → descargando TODO.")
            selection = {"antibodies": set(), "cancers_by_antibody": {}}
            sel_total = self.total_items

        self.progress.configure(value=0, maximum=max(sel_total, 1))

        def progress_cb(i, total):
            self.msg_q.put(("progress", i, total))

        def log_cb(m):
            self.msg_q.put(("log", m))

        def worker():
            try:
                download_from_inventory(
                    self.gene_name, self.inv, self.ext_var.get(), selection,
                    progress_cb=progress_cb, log=log_cb
                )
                self.msg_q.put(("download_ok",))
            except Exception as e:
                self.msg_q.put(("error", f"Descarga falló: {e}"))

        threading.Thread(target=worker, daemon=True).start()

    def _poll_queue(self):
        try:
            while True:
                msg = self.msg_q.get_nowait()
                kind = msg[0]

                if kind == "log":
                    self.log(msg[1])

                elif kind == "progress":
                    i, total = msg[1], msg[2]
                    self.progress.configure(value=i, maximum=total)
                    self.set_status(f"Descargando... {i}/{total}")

                elif kind == "preview_ok":
                    _, gene, inv, total = msg
                    self.gene_name = gene
                    self.inv = inv
                    self.total_items = total
                    self._populate_tree_and_selection()
                    self.set_status(f"Preview listo: {gene} ({total} items)")
                    self.log(f"Preview OK: {gene} | {total} items")
                    self.preview_btn.configure(state="normal")
                    self.download_btn.configure(state="normal")

                elif kind == "download_ok":
                    self.set_status("Descarga terminada.")
                    self.log("✅ Descarga terminada.")
                    self.preview_btn.configure(state="normal")
                    self.download_btn.configure(state="normal")

                elif kind == "error":
                    self.set_status("Error.")
                    self.log("❌ " + msg[1])
                    messagebox.showerror("Error", msg[1])
                    self.preview_btn.configure(state="normal")
                    # habilita download solo si había preview previa
                    self.download_btn.configure(state="normal" if self.inv else "disabled")

        except queue.Empty:
            pass

        self.after(100, self._poll_queue)


if __name__ == "__main__":
    try:
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass

    app = App()
    app.mainloop()
