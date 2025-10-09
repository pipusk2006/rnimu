import requests
from urllib.parse import quote

BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"
BIOMES = {

    # === малорастительные и экстремальные среды ===
    "desert":       "root:Environmental:Terrestrial:Soil:Desert",
    "shrubland":    "root:Environmental:Terrestrial:Soil:Shrubland",
    "permafrost":   "root:Environmental:Terrestrial:Soil:Permafrost",
    "tundra":       "root:Environmental:Terrestrial:Soil:Tundra",

    # === влажные и заболоченные ===
    "wetland":      "root:Environmental:Terrestrial:Soil:Wetlands",
    "peatland":     "root:Environmental:Terrestrial:Soil:Peatland",

    # === антропогенные и загрязнённые ===
    "contaminated": "root:Environmental:Terrestrial:Soil:Contaminated",
    "pasture":      "root:Environmental:Terrestrial:Soil:Pasture",

    # === менее распространённые, но встречаются ===
    "steppe":       "root:Environmental:Terrestrial:Soil:Steppe",
    "savanna":      "root:Environmental:Terrestrial:Soil:Savanna",
}

for name, lineage in BIOMES.items():
    url = f"{BASE}/biomes/{quote(lineage, safe='')}"
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        obj = r.json().get("data", {}).get("attributes", {})
        print(f"{name:12s} → samples-count = {obj.get('samples-count')}  (lineage: {lineage})")
    except Exception as e:
        print(f"{name:12s} → ERROR: {e}")
