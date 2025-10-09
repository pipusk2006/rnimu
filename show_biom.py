# Если пакета нет, раскомментируй следующую строку:
# !pip install biom-format pandas matplotlib

from biom import load_table
import pandas as pd
import numpy as np
import os

PATH = "biom_data/tropical_rainforest/tropical_rainforest_2.biom"   # поправь путь при необходимости

def load_biom(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Не найден файл: {path}")
    table = load_table(path)  # авто-детект JSON/HDF5
    return table

def table_info(table):
    n_obs, n_samp = table.shape
    sample_ids = list(table.ids(axis='sample'))
    print(f"Форма таблицы: {n_obs} наблюдений (OTU/ASV) × {n_samp} сэмплов")
    print("Первые сэмплы:", sample_ids[:5])

def biom_to_df(table, dense=True):
    """DataFrame: строки — наблюдения (OTU/ASV), столбцы — сэмплы."""
    df = table.to_dataframe(dense=dense)
    # Приведём имена наблюдений к строке (иногда бывают байты)
    df.index = df.index.map(lambda x: x.decode() if isinstance(x, bytes) else x)
    df.columns = df.columns.map(lambda x: x.decode() if isinstance(x, bytes) else x)
    return df

def extract_taxonomy(table):
    """
    Возвращает Series с таксономией для каждого observation_id.
    Пытается вытащить 'taxonomy' из метаданных (список рангов или строка).
    Если таксономии нет — подставляет сам observation_id.
    """
    obs_ids = list(table.ids(axis='observation'))
    meta = table.metadata(axis='observation')
    taxa = []

    for oid, m in zip(obs_ids, meta if meta is not None else [None]*len(obs_ids)):
        name = None
        if isinstance(m, dict):
            tx = m.get('taxonomy') or m.get('Taxonomy') or m.get('lineage')
            if tx is not None:
                # В MGnify taxonomy бывает списком рангов или строкой с ';'
                if isinstance(tx, (list, tuple)):
                    # берём самый “глубокий” ненулевой ранг
                    tx_clean = [str(t).strip() for t in tx if t and str(t).strip() and str(t).strip() != 'NA']
                    name = tx_clean[-1] if tx_clean else None
                else:
                    # строка: 'k__Bacteria; p__...; g__...; s__...'
                    parts = [p.strip() for p in str(tx).split(';') if p.strip()]
                    name = parts[-1] if parts else None
        # если не получилось — fallback на id наблюдения
        if not name:
            name = oid.decode() if isinstance(oid, bytes) else str(oid)
        taxa.append(name)

    # чистим префиксы рангов (k__, p__, g__ и т.п.)
    def strip_rank_prefix(s):
        s = str(s)
        if '__' in s[:4]:
            return s.split('__', 1)[1] or s
        return s
    taxa = pd.Series([strip_rank_prefix(t) for t in taxa], index=[oid.decode() if isinstance(oid, bytes) else oid for oid in obs_ids], name='taxonomy')
    return taxa

def top_taxa_overall(table, topn=20):
    df = biom_to_df(table, dense=True)
    taxa = extract_taxonomy(table)
    # сумма по всем сэмплам для каждого наблюдения
    totals = df.sum(axis=1)
    # агрегируем по таксон-имени
    agg = totals.groupby(taxa).sum().sort_values(ascending=False)
    rel = agg / agg.sum()
    out = pd.DataFrame({"abundance": agg, "relative": rel})
    return out.head(topn)

def top_taxa_for_sample(table, sample_id, topn=20):
    df = biom_to_df(table, dense=True)
    if sample_id not in df.columns:
        raise KeyError(f"Сэмпл {sample_id} не найден. Доступные: {list(df.columns)[:5]} ...")
    taxa = extract_taxonomy(table)
    counts = df[sample_id]
    agg = counts.groupby(taxa).sum().sort_values(ascending=False)
    rel = agg / agg.sum()
    out = pd.DataFrame({"abundance": agg, "relative": rel})
    return out.head(topn)

# --------- запуск ----------
table = load_biom(PATH)
table_info(table)

# Топ таксонов по всей таблице (агрегировано по всем сэмплам)
overall_top = top_taxa_overall(table, topn=20)
print("\nТоп-20 таксонов (общая абунд.):")
print(overall_top)

# Если хочешь — топ по конкретному сэмплу:
samples = list(table.ids(axis='sample'))
if samples:
    sample_top = top_taxa_for_sample(table, samples[0], topn=15)
    print(f"\nТоп-15 таксонов для сэмпла {samples[0]}:")
    print(sample_top)

# (опционально) простой барчарт — раскомментируй в ноутбуке
# import matplotlib.pyplot as plt
# overall_top['relative'].head(15).plot(kind='bar')
# plt.ylabel('Relative abundance')
# plt.title('Top-15 taxa (overall)')
# plt.tight_layout()
# plt.show()
