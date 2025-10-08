#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, time, sys, re
import requests
from urllib.parse import quote

BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"
OUT_DIR = "biom_data"
os.makedirs(OUT_DIR, exist_ok=True)

# 5 целевых почвенных классов → lineage в MGnify
BIOMES = {
    "forest":       "root:Environmental:Terrestrial:Soil:Forest soil",
    "wetland":      "root:Environmental:Terrestrial:Soil:Wetlands",
    "grassland":    "root:Environmental:Terrestrial:Soil:Grasslands",
    "desert":       "root:Environmental:Terrestrial:Soil:Desert",
    "agricultural": "root:Environmental:Terrestrial:Soil:Agricultural",
}

TARGET_PER_CLASS = 100
REQ_TIMEOUT = 30
SLEEP = 0.2  # небольшой бэкофф, чтобы не душить API

session = requests.Session()
session.headers.update({"Accept": "application/json"})

def get_json(url, params=None):
    for attempt in range(5):
        try:
            r = session.get(url, params=params, timeout=REQ_TIMEOUT)
            if r.status_code == 429:
                # rate limit – подождать подольше
                time.sleep(3 + attempt)
                continue
            r.raise_for_status()
            return r.json()
        except Exception as e:
            if attempt == 4:
                raise
            time.sleep(1 + attempt)
    return None

def iter_pages(url, params=None):
    """Итерация по страницам JSON:API (links.next)."""
    while url:
        data = get_json(url, params=params)
        if data is None: break
        yield data
        links = data.get("links", {})
        url = links.get("next")

def pick_biom_download(downloads):
    """Вернуть (url, alias) первого доступного .biom из списка downloads."""
    for item in downloads.get("data", []):
        attrs = item.get("attributes", {})
        alias = attrs.get("alias", "")
        # MGnify помечает BIOM так: *_OTU_TABLE_JSON.biom или *_HDF5.biom
        if alias.lower().endswith(".biom"):
            link = item.get("links", {}).get("self")
            if link:
                return link, alias
    return None, None

def ensure_unique_path(dirpath, base_name):
    """forest_1.biom ... forest_100.biom, без перезаписи."""
    path = os.path.join(dirpath, base_name)
    if not os.path.exists(path):
        return path
    # если вдруг имя занято, добавим суффикс
    stem, ext = os.path.splitext(base_name)
    i = 2
    while True:
        candidate = os.path.join(dirpath, f"{stem}__{i}{ext}")
        if not os.path.exists(candidate):
            return candidate
        i += 1

def download_file(file_url, out_path):
    for attempt in range(5):
        try:
            with session.get(file_url, stream=True, timeout=REQ_TIMEOUT) as r:
                if r.status_code == 429:
                    time.sleep(3 + attempt)
                    continue
                r.raise_for_status()
                with open(out_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1048576):
                        if chunk:
                            f.write(chunk)
            return True
        except Exception as e:
            if attempt == 4:
                print(f"[!] FAIL {file_url} -> {out_path}: {e}", file=sys.stderr)
                return False
            time.sleep(1 + attempt)
    return False

def harvest_for_biome(class_name, lineage, n_target=TARGET_PER_CLASS):
    print(f"\n=== {class_name.upper()} ===")
    saved = 0
    page_url = f"{BASE}/biomes/{quote(lineage, safe='')}/samples"

    # Папка для класса
    class_dir = os.path.join(OUT_DIR, class_name)
    os.makedirs(class_dir, exist_ok=True)

    # Перебираем все sample-ы данного биома
    for page in iter_pages(page_url):
        for samp in page.get("data", []):
            if saved >= n_target:
                break

            sid = samp.get("id")  # sample accession (e.g. SRS..., ERS...)
            runs_link = samp.get("relationships", {}).get("runs", {}).get("links", {}).get("related")
            if not runs_link:
                continue

            # Просматриваем все runs этого sample (лучше брать amplicon – там BIOM точно есть)
            runs = get_json(runs_link)
            for run in runs.get("data", []):
                if saved >= n_target:
                    break

                exp_type = run.get("attributes", {}).get("experiment-type", "")
                # Приоритет ампликонных данных (SSU/ITS/LSU) – там всегда есть OTU BIOM
                if exp_type not in ("amplicon", "SSU", "ITS", "LSU"):
                    continue

                analyses_link = run.get("relationships", {}).get("analyses", {}).get("links", {}).get("related")
                if not analyses_link:
                    continue

                analyses = get_json(analyses_link)
                for an in analyses.get("data", []):
                    if saved >= n_target:
                        break
                    an_id = an.get("id")
                    dl_link = an.get("relationships", {}).get("downloads", {}).get("links", {}).get("related")
                    if not dl_link:
                        continue

                    downloads = get_json(dl_link)
                    file_url, alias = pick_biom_download(downloads)
                    if not file_url:
                        continue

                    # Имя по ТЗ: forest_1.biom, ...
                    out_name = f"{class_name}_{saved+1}.biom"
                    out_path = ensure_unique_path(class_dir, out_name)
                    ok = download_file(file_url, out_path)
                    if ok:
                        saved += 1
                        print(f"[{class_name}] {saved}/{n_target}: {out_name}")
                        time.sleep(SLEEP)
                    if saved >= n_target:
                        break

            if saved >= n_target:
                break

        # маленькая пауза между страницами
        time.sleep(SLEEP)

    print(f"Итого для {class_name}: {saved} файлов.\n")
    return saved

def main():
    total = 0
    for cname, lineage in BIOMES.items():
        try:
            total += harvest_for_biome(cname, lineage, TARGET_PER_CLASS)
        except Exception as e:
            print(f"[!] Ошибка на классе {cname}: {e}", file=sys.stderr)
    print(f"Готово. Всего скачано: {total} BIOM-файлов. См. папку: {OUT_DIR}")

if __name__ == "__main__":
    main()
