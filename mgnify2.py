#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, time, shutil, subprocess, re, json, hashlib
from urllib.parse import quote
import requests

# ================== НАСТРОЙКИ ==================
BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"
OUT_DIR = "biom_data"
TARGET_PER_CLASS = 100
REQ_TIMEOUT = 60
SLEEP = 0.05

# Лимиты обхода
MAX_SAMPLE_PAGES     = 120
MAX_RUNS_PER_SAMPLE  = 12
MAX_ANALYSES_PER_RUN = 12
CLASS_TIME_LIMIT     = 20*60*10  # None → снять лимит

# ДОБАВЬ СВОИ КЛАССЫ СЮДА
BIOMES = {
    "tropical_rainforest": "root:Environmental:Terrestrial:Soil:Tropical rainforest",
}

# ===============================================
session = requests.Session()
session.headers.update({"Accept": "application/json"})
os.makedirs(OUT_DIR, exist_ok=True)

IDX_RE = re.compile(r"^(?P<cls>[a-z_]+)_(?P<idx>\d+)\.biom$", re.IGNORECASE)

# -------- state --------
def _state_path(class_dir): return os.path.join(class_dir, ".state.json")

def load_state(class_dir):
    if os.path.exists(_state_path(class_dir)):
        try:
            data = json.load(open(_state_path(class_dir), "r", encoding="utf-8"))
            return {
                "seen_links": set(data.get("seen_links") or []),
                "hash_to_name": dict(data.get("hash_to_name") or {}),   # старый формат, пусть будет
                "biom_sig_to_name": dict(data.get("biom_sig_to_name") or {}),
            }
        except Exception:
            pass
    return {"seen_links": set(), "hash_to_name": {}, "biom_sig_to_name": {}}

def save_state(class_dir, state):
    data = {
        "seen_links": sorted(state["seen_links"]),
        "hash_to_name": state["hash_to_name"],
        "biom_sig_to_name": state["biom_sig_to_name"],
    }
    with open(_state_path(class_dir), "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)

# -------- helpers --------
def next_index_for_class(class_dir: str, class_name: str) -> int:
    max_idx = 0
    if not os.path.isdir(class_dir): return 1
    for fname in os.listdir(class_dir):
        m = IDX_RE.match(fname)
        if not m: continue
        if m.group("cls").lower() != class_name.lower(): continue
        try:
            max_idx = max(max_idx, int(m.group("idx")))
        except ValueError:
            pass
    return max_idx + 1

def get_json(url, params=None, retries=5):
    for a in range(retries):
        try:
            r = session.get(url, params=params, timeout=REQ_TIMEOUT)
            if r.status_code == 429:
                time.sleep(3 + a); continue
            r.raise_for_status()
            return r.json()
        except Exception:
            if a == retries - 1: raise
            time.sleep(1 + a)

def iter_pages(url, params=None):
    while url:
        data = get_json(url, params=params)
        if not data: break
        yield data
        url = (data.get("links") or {}).get("next")

def _try_url_variants(url: str):
    variants = [url]
    if url.endswith(".bio"):
        variants.append(url + "m")
    if "/file/" in url:
        variants.append(url.replace("/file/", "/download/"))
    if "/download/" in url:
        variants.append(url.replace("/download/", "/file/"))
    seen, out = set(), []
    for v in variants:
        if v not in seen:
            out.append(v); seen.add(v)
    return out

def download_file_atomic(url, out_path, retries=5):
    tmp_path = out_path + ".part"
    for candidate in _try_url_variants(url):
        for a in range(retries):
            try:
                with session.get(candidate, stream=True, timeout=REQ_TIMEOUT) as r:
                    if r.status_code == 429:
                        time.sleep(3 + a); continue
                    if r.status_code == 404:
                        break
                    r.raise_for_status()
                    with open(tmp_path, "wb") as f:
                        for chunk in r.iter_content(1048576):
                            if chunk: f.write(chunk)
                os.replace(tmp_path, out_path)
                return True, candidate
            except Exception as e:
                try:
                    if os.path.exists(tmp_path): os.remove(tmp_path)
                except Exception:
                    pass
                if a == retries - 1:
                    print(f"[!] FAIL {candidate} -> {out_path}: {e}", file=sys.stderr)
                time.sleep(1 + a)
    return False, None

# ---------- BIOM-сигнатура (контентная) ----------
def compute_biom_signature(path):
    """
    Контентная сигнатура BIOM:
    - загружаем через biom.load_table()
    - получаем dense DF (obs × samples)
    - сортируем индексы (строки/столбцы)
    - округляем до целых (OTU/ASV счётчики)
    - sha256( row_ids || col_ids || matrix_bytes )
    """
    try:
        from biom import load_table
        import numpy as np
    except Exception as e:
        # biom не установлен — вернём None, чтобы не падать
        print("[warn] biom не установлен, контентная дедупликация отключена.", file=sys.stderr)
        return None

    tbl = load_table(path)  # авто JSON/HDF5
    df = tbl.to_dataframe(dense=True)

    # нормализуем имена
    def dec(x): 
        if isinstance(x, bytes):
            try: return x.decode()
            except Exception: return str(x)
        return str(x)
    df.index = df.index.map(dec)
    df.columns = df.columns.map(dec)

    # сортируем
    df = df.sort_index(axis=0).sort_index(axis=1)

    # к целым
    arr = df.values
    if issubclass(arr.dtype.type, (float,)):
        arr = (arr + 0.0)  # make sure float
        arr = np.rint(arr).astype("int64", copy=False)
    else:
        arr = arr.astype("int64", copy=False)

    h = hashlib.sha256()
    for rid in df.index:
        h.update(rid.encode("utf-8")); h.update(b"\0")
    h.update(b"|")
    for cid in df.columns:
        h.update(cid.encode("utf-8")); h.update(b"\0")
    h.update(b"|")
    h.update(arr.tobytes(order="C"))
    return h.hexdigest()

# ---------- поиск .biom / OTU .tsv ----------
def find_biom_or_tsv(downloads_json):
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
        if "otu" in low and low.endswith(".tsv") and tsv_url is None:
            tsv_url, tsv_alias = link, alias
    return biom_url, biom_alias, tsv_url, tsv_alias

# ---------- BIOM CLI / fallback ----------
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

# ---------- сбор одного класса (resume + антидубликаты) ----------
def harvest_class(class_name, lineage, n_target=TARGET_PER_CLASS):
    print(f"\n=== {class_name.upper()} ===")
    start_ts = time.time()
    class_dir = os.path.join(OUT_DIR, class_name)
    os.makedirs(class_dir, exist_ok=True)

    state = load_state(class_dir)

    # предварительно посчитаем сигнатуры для уже существующих файлов (если state пуст)
    if not state["biom_sig_to_name"]:
        for fname in sorted(os.listdir(class_dir)):
            if not fname.lower().endswith(".biom"): continue
            if not fname.startswith(class_name + "_"): continue
            fpath = os.path.join(class_dir, fname)
            sig = compute_biom_signature(fpath)
            if sig:
                state["biom_sig_to_name"].setdefault(sig, fname)
        save_state(class_dir, state)

    next_idx = next_index_for_class(class_dir, class_name)
    saved = next_idx - 1
    remain = max(0, n_target - saved)
    print(f"[resume] уже есть: {saved}, цель: {n_target}, докачать: {remain}")
    if remain == 0:
        print(f"Итого для {class_name}: {saved} файлов (цель достигнута).")
        return saved

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

                    # анти-дуб по ссылкам (чтобы не ходить повторно)
                    candidate_link = biom_url or tsv_url
                    if candidate_link and candidate_link in state["seen_links"]:
                        continue

                    # ------- 1) прямой BIOM -------
                    if biom_url:
                        out_name = f"{class_name}_{next_idx}.biom"
                        out_path = os.path.join(class_dir, out_name)

                        # качаем во временный .part, но с расширением .biom для совместимости с loader
                        tmp_path = out_path + ".tmp"
                        ok, used = download_file_atomic(biom_url, tmp_path)
                        if not ok:
                            state["seen_links"].add(candidate_link)  # больше не пытаться
                            save_state(class_dir, state)
                            continue

                        # контентная сигнатура
                        sig = compute_biom_signature(tmp_path)
                        if sig and sig in state["biom_sig_to_name"]:
                            # дубликат — удаляем и НЕ увеличиваем индекс
                            try: os.remove(tmp_path)
                            except OSError: pass
                            state["seen_links"].add(candidate_link)
                            save_state(class_dir, state)
                            # print(f"[{class_name}] duplicate content → drop {out_name}")
                            continue

                        # уникальный — переименовываем
                        os.replace(tmp_path, out_path)
                        if sig:
                            state["biom_sig_to_name"][sig] = out_name
                        state["seen_links"].add(candidate_link)
                        save_state(class_dir, state)

                        saved += 1
                        print(f"[{class_name}] {saved}/{n_target}: {out_name} (BIOM)")
                        next_idx += 1
                        time.sleep(SLEEP)
                        if saved >= n_target: break
                        continue

                    # ------- 2) TSV → BIOM -------
                    if tsv_url:
                        tmp_tsv  = os.path.join(class_dir, f"__tmp_{class_name}_{next_idx}.tsv")
                        ok, used = download_file_atomic(tsv_url, tmp_tsv)
                        if not ok:
                            state["seen_links"].add(candidate_link)
                            save_state(class_dir, state)
                            continue

                        out_name = f"{class_name}_{next_idx}.biom"
                        out_path = os.path.join(class_dir, out_name)

                        try:
                            if convert_tsv_to_biom(tmp_tsv, out_path):
                                # сигнатура уже по BIOM
                                sig = compute_biom_signature(out_path)
                                if sig and sig in state["biom_sig_to_name"]:
                                    # дубликат — удалить
                                    try: os.remove(out_path)
                                    except OSError: pass
                                    state["seen_links"].add(candidate_link)
                                    save_state(class_dir, state)
                                    # print(f"[{class_name}] duplicate content (from TSV) → drop {out_name}")
                                else:
                                    if sig:
                                        state["biom_sig_to_name"][sig] = out_name
                                    state["seen_links"].add(candidate_link)
                                    save_state(class_dir, state)
                                    saved += 1
                                    print(f"[{class_name}] {saved}/{n_target}: {out_name} (from TSV)")
                                    next_idx += 1
                                    time.sleep(SLEEP)
                        finally:
                            try: os.remove(tmp_tsv)
                            except OSError: pass
        # конец страницы

    print(f"Итого для {class_name}: {saved} файлов.")
    return saved

def main():
    total = 0
    try:
        for cname, lineage in BIOMES.items():
            try:
                total += harvest_class(cname, lineage, TARGET_PER_CLASS)
            except subprocess.CalledProcessError as e:
                print(f"[!] Ошибка конвертации для {cname}: {e}", file=sys.stderr)
            except KeyboardInterrupt:
                print("\n[!] Прервано пользователем.")
                break
            except Exception as e:
                print(f"[!] Ошибка на классе {cname}: {e}", file=sys.stderr)
    finally:
        print(f"\nГотово: всего скачано {total} BIOM-файлов. Папка: {OUT_DIR}")

if __name__ == "__main__":
    main()
