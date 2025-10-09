#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, time, shutil, subprocess
from urllib.parse import quote
import requests

BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"
OUT_DIR = "biom_data"
TARGET_PER_CLASS = 200
REQ_TIMEOUT = 60
SLEEP = 0.05
MAX_SAMPLE_PAGES     = 120
MAX_RUNS_PER_SAMPLE  = 12
MAX_ANALYSES_PER_RUN = 12
CLASS_TIME_LIMIT     = None  # None чтобы снять лимит

BIOMES = {
    "forest":       "root:Environmental:Terrestrial:Soil:Forest soil",
}

session = requests.Session()
session.headers.update({"Accept": "application/json"})
os.makedirs(OUT_DIR, exist_ok=True)

def get_json(url, params=None, retries=5):
    for a in range(retries):
        try:
            r = session.get(url, params=params, timeout=REQ_TIMEOUT)
            if r.status_code == 429:
                time.sleep(3 + a); continue
            r.raise_for_status()
            return r.json()
        except Exception as e:
            if a == retries - 1: raise
            time.sleep(1 + a)

def iter_pages(url, params=None):
    while url:
        data = get_json(url, params=params)
        if not data: break
        yield data
        url = (data.get("links") or {}).get("next")

def ensure_unique_path(dirpath, base_name):
    path = os.path.join(dirpath, base_name)
    if not os.path.exists(path): return path
    stem, ext = os.path.splitext(base_name); i = 2
    while True:
        cand = os.path.join(dirpath, f"{stem}__{i}{ext}")
        if not os.path.exists(cand): return cand
        i += 1

def download_file(file_url, out_path, retries=5):
    for a in range(retries):
        try:
            with session.get(file_url, stream=True, timeout=REQ_TIMEOUT) as r:
                if r.status_code == 429:
                    time.sleep(3 + a); continue
                r.raise_for_status()
                with open(out_path, "wb") as f:
                    for chunk in r.iter_content(1048576):
                        if chunk: f.write(chunk)
            return True
        except Exception as e:
            if a == retries - 1:
                print(f"[!] FAIL {file_url} -> {out_path}: {e}", file=sys.stderr)
                return False
            time.sleep(1 + a)

def find_biom_or_tsv(downloads_json):
    """
    FIX: Брать только реальные ссылки из downloads.
    Возвращает (biom_url, biom_alias, tsv_url, tsv_alias).
    """
    biom_url = biom_alias = tsv_url = tsv_alias = None
    for item in (downloads_json.get("data") or []):
        attrs = item.get("attributes") or {}
        alias = (attrs.get("alias") or "").strip()
        link  = (item.get("links") or {}).get("self")
        if not alias or not link:
            continue
        low = alias.lower()
        if low.endswith(".biom"):
            biom_url, biom_alias = link, alias
            break  # приоритет .biom
        # запасной вариант: любые OTU-таблицы в TSV
        if "otu" in low and low.endswith(".tsv") and tsv_url is None:
            tsv_url, tsv_alias = link, alias
    return biom_url, biom_alias, tsv_url, tsv_alias

def locate_biom_cli():
    cli = shutil.which("biom")
    if cli: return cli
    guess = os.path.join(sys.exec_prefix, "Scripts", "biom.exe")
    return guess if os.path.exists(guess) else None

def convert_tsv_to_biom(tsv_path, biom_out_path):
    biom_cli = locate_biom_cli()
    if biom_cli:
        cmd = [biom_cli, "convert", "-i", tsv_path, "-o", biom_out_path,
               "--table-type=OTU table", "--to-hdf5"]
        subprocess.check_call(cmd)
        return True
    # python fallback
    try:
        import pandas as pd
        from biom.table import Table
        df = pd.read_csv(tsv_path, sep="\t", index_col=0)
        table = Table(df.values,
                      observation_ids=list(df.index.astype(str)),
                      sample_ids=list(df.columns.astype(str)))
        table.to_hdf5(biom_out_path, generated_by="mgnify-collector")
        return True
    except Exception as e:
        print(f"[!] TSV→BIOM fallback не удался: {e}", file=sys.stderr)
        return False

def harvest_class(class_name, lineage, n_target=TARGET_PER_CLASS):
    print(f"\n=== {class_name.upper()} ===")
    start_ts = time.time(); saved = 0
    class_dir = os.path.join(OUT_DIR, class_name); os.makedirs(class_dir, exist_ok=True)
    page_url = f"{BASE}/biomes/{quote(lineage, safe='')}/samples"

    page_i = 0
    for page in iter_pages(page_url):
        page_i += 1
        print(f"[{class_name}] page={page_i} saved={saved}")
        if saved >= n_target: break
        if CLASS_TIME_LIMIT and (time.time() - start_ts > CLASS_TIME_LIMIT):
            print(f"[{class_name}] Достигнут лимит времени."); break
        if page_i > MAX_SAMPLE_PAGES:
            print(f"[{class_name}] Достигнут лимит страниц samples."); break

        for samp in (page.get("data") or []):
            if saved >= n_target: break
            runs_link = (((samp.get("relationships") or {}).get("runs") or {}).get("links") or {}).get("related")
            if not runs_link: continue
            runs = get_json(runs_link) or {}
            run_cnt = 0

            for run in (runs.get("data") or []):
                run_cnt += 1
                if run_cnt > MAX_RUNS_PER_SAMPLE: break
                if saved >= n_target: break

                analyses_link = (((run.get("relationships") or {}).get("analyses") or {}).get("links") or {}).get("related")
                if not analyses_link: continue
                analyses = get_json(analyses_link) or {}
                an_cnt = 0

                for an in (analyses.get("data") or []):
                    an_cnt += 1
                    if an_cnt > MAX_ANALYSES_PER_RUN: break
                    if saved >= n_target: break

                    dl_link = (((an.get("relationships") or {}).get("downloads") or {}).get("links") or {}).get("related")
                    if not dl_link: continue

                    downloads = get_json(dl_link) or {}
                    biom_url, biom_alias, tsv_url, tsv_alias = find_biom_or_tsv(downloads)

                    # DEBUG/лог: что реально нашли
                    if not biom_url and not tsv_url:
                        print(f"[{class_name}] skip analysis {an.get('id')} — нет .biom и OTU .tsv")
                        continue

                    # 1) BIOM напрямую
                    if biom_url:
                        out_name = f"{class_name}_{saved+1}.biom"
                        out_path = ensure_unique_path(class_dir, out_name)
                        print(f"[{class_name}] download BIOM: {biom_alias} ← {biom_url}")
                        if download_file(biom_url, out_path):
                            saved += 1
                            print(f"[{class_name}] {saved}/{n_target}: {os.path.basename(out_path)} (BIOM)")
                        time.sleep(SLEEP)
                        if saved >= n_target: break
                        continue

                    # 2) OTU TSV → BIOM
                    if tsv_url:
                        tmp_tsv = os.path.join(class_dir, f"__tmp_{class_name}_{saved+1}.tsv")
                        out_name = f"{class_name}_{saved+1}.biom"
                        out_path = ensure_unique_path(class_dir, out_name)
                        print(f"[{class_name}] download TSV:  {tsv_alias} ← {tsv_url}")
                        if download_file(tsv_url, tmp_tsv):
                            try:
                                if convert_tsv_to_biom(tmp_tsv, out_path):
                                    saved += 1
                                    print(f"[{class_name}] {saved}/{n_target}: {os.path.basename(out_path)} (from TSV)")
                            finally:
                                try: os.remove(tmp_tsv)
                                except OSError: pass
                            time.sleep(SLEEP)
                            if saved >= n_target: break
                        continue
        # конец страницы
    print(f"Итого для {class_name}: {saved} файлов.")
    return saved

def main():
    total = 0
    for cname, lineage in BIOMES.items():
        try:
            total += harvest_class(cname, lineage, TARGET_PER_CLASS)
        except subprocess.CalledProcessError as e:
            print(f"[!] Ошибка конвертации для {cname}: {e}", file=sys.stderr)
        except KeyboardInterrupt:
            print("\n[!] Прервано пользователем."); break
        except Exception as e:
            print(f"[!] Ошибка на классе {cname}: {e}", file=sys.stderr)
    print(f"\nГотово: всего скачано {total} BIOM-файлов. Папка: {OUT_DIR}")

if __name__ == "__main__":
    main()



