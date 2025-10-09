# list_soil_biomes_turbo.py
# pip install requests

import requests, time, sys, json, collections
from urllib.parse import quote

BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"
ROOT = "root:Environmental:Terrestrial:Soil"

# --- настройки скорости/ограничений ---
TIMEOUT_S         = 30          # HTTP-таймаут на запрос
GLOBAL_TIME_LIMIT = 120         # общий лимит на обход (сек)
MAX_DEPTH         = 2           # глубина: 0=Soil, 1=дети Soil, 2=внуки
MAX_PAGES_PER_NODE= 5           # страниц children на один узел
PAGE_SIZE         = 100         # просим больше объектов на страницу (сервер может игнорить)
# --------------------------------------

session = requests.Session()
session.headers.update({"Accept":"application/json"})

def get_json(url, params=None):
    r = session.get(url, params=params, timeout=TIMEOUT_S)
    r.raise_for_status()
    return r.json()

def iter_children_pages(lineage):
    url = f"{BASE}/biomes/{quote(lineage, safe='')}/children"
    params = {"page_size": PAGE_SIZE}
    pages = 0
    while url and pages < MAX_PAGES_PER_NODE:
        data = get_json(url, params=params)
        yield data
        url = (data.get("links") or {}).get("next")
        params = None
        pages += 1

def crawl_soil(max_depth=MAX_DEPTH, global_time_limit=GLOBAL_TIME_LIMIT):
    start = time.time()
    seen = set([ROOT])
    # сохраняем: name -> {'lineage':..., 'samples': int}
    biomes = {}

    # добавим сам ROOT (одним запросом, чтобы знать его имя/счётчик)
    try:
        root_obj = get_json(f"{BASE}/biomes/{quote(ROOT, safe='')}")
        ra = (root_obj.get("data") or {}).get("attributes") or {}
        biomes[ra.get("biome-name","Soil")] = {"lineage": ROOT, "samples": ra.get("samples-count")}
    except Exception:
        pass

    Q = collections.deque([(ROOT, 0)])
    requests_done = 0

    while Q:
        if global_time_limit and time.time() - start > global_time_limit:
            print(f"[!] Прервано по лимиту времени {global_time_limit}s", file=sys.stderr)
            break

        lineage, depth = Q.popleft()
        if depth >= max_depth:
            continue

        try:
            for page in iter_children_pages(lineage):
                requests_done += 1
                for item in (page.get("data") or []):
                    attrs = item.get("attributes") or {}
                    child_lin = attrs.get("lineage")
                    name = attrs.get("biome-name") or (child_lin.split(":")[-1] if child_lin else None)
                    samples = attrs.get("samples-count")
                    if not child_lin or child_lin in seen:
                        continue
                    seen.add(child_lin)
                    # запишем сразу без доп. запросов
                    biomes[name] = {"lineage": child_lin, "samples": samples}
                    # пойдём глубже
                    Q.append((child_lin, depth+1))
        except requests.HTTPError as e:
            # не валимся на 404/500: просто пропустим узел
            print(f"[skip] children of {lineage} -> {e}", file=sys.stderr)
            continue
        except Exception as e:
            print(f"[skip] {lineage} -> {e}", file=sys.stderr)
            continue

    # упорядочим по числу сэмплов
    order = sorted(
        [(n, v["samples"], v["lineage"]) for n,v in biomes.items() if isinstance(v.get("samples"), int)],
        key=lambda x: (-x[1], x[0])
    )
    return order, biomes, {"requests": requests_done, "elapsed_s": round(time.time()-start,2)}

def emit_dict(biomes, min_samples=50):
    # делаем красивый python-словарь для вставки в код
    lines = ["BIOMES = {"]
    used = set()
    def slug(s):
        return (s.lower().replace(" ","_").replace("-","_").replace("/","_")
                .replace("(","").replace(")","").replace("__","_"))
    for name, info in sorted(biomes.items(), key=lambda kv: (-(kv[1]["samples"] or 0), kv[0])):
        s = info["samples"]
        if s is None or s < min_samples: 
            continue
        key = slug(name)
        i = 2
        while key in used:
            key = f"{key}_{i}"; i += 1
        used.add(key)
        lines.append(f'    "{key}": "{info["lineage"]}",  # {name} — samples={s}')
    lines.append("}")
    return "\n".join(lines)

if __name__ == "__main__":
    top, biomes, stats = crawl_soil()
    print(f"Done in {stats['elapsed_s']}s, requests={stats['requests']}. Узлов: {len(biomes)}\n")

    # покажем топ-20
    print("Топ-20 Soil-биомов по числу сэмплов:")
    for name, s, lin in top[:20]:
        print(f"- {name:35s}  samples={s:5d}  | {lin}")
    print("\nСловарь для вставки (>=50 samples):\n")
    print(emit_dict({n: {"lineage": l, "samples": s} for n,s,l in top}, min_samples=50))

