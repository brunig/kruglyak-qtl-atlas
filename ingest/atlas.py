"""Shared pieces of the atlas build: genome tables, number formatting, and the merge of
per-dataset fragments into docs/data.json.

A fragment is a dict of columnar tables whose cross-table references are row indices
local to that dataset (gene indices are global: every dataset uses the same gene table).
`assemble()` rebases those indices and writes one docs/data.json for all datasets.
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCES = os.path.join(ROOT, "sources")
DOCS = os.path.join(ROOT, "docs")
FRAG_DIR = os.path.join(ROOT, "build", "fragments")

SCHEMA = "atlas-1"
BUILD = "S288C R64 (UCSC sacCer3)"
HEAT_BIN = 50000

# Table order in data.json, and which columns of each table point into which table.
TABLES = ["datasets", "sources", "parents", "crosses", "traits", "qtl", "jeff", "causal", "cand", "vexp"]
REFS = {
    "sources": {"ds": "datasets"},
    "parents": {"ds": "datasets"},
    "crosses": {"ds": "datasets", "p1": "parents", "p2": "parents", "g1": "parents"},
    "traits": {"ds": "datasets"},
    "qtl": {"ds": "datasets", "t": "traits", "x": "crosses", "up": "parents"},
    "jeff": {"q": "qtl", "x": "crosses"},
    "causal": {"ds": "datasets", "t": "traits", "sp": "parents", "q1": "qtl", "q2": "qtl"},
    "cand": {"ds": "datasets"},
    "vexp": {"t": "traits", "x": "crosses"},
}


def chrom_lengths():
    out = {}
    with open(os.path.join(SOURCES, "sgd", "chrom_lengths_R64.tsv")) as f:
        rows = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
    for c, n in rows[1:]:
        out[c] = int(n)
    return out


def genes():
    cols = {k: [] for k in ["c", "a", "e", "s", "id", "n", "k"]}
    with open(os.path.join(SOURCES, "sgd", "genes_R64-5-1.tsv")) as f:
        rows = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
    for r in rows[1:]:
        r = r + [""] * (7 - len(r))
        cols["c"].append(r[0]); cols["a"].append(int(r[1])); cols["e"].append(int(r[2]))
        cols["s"].append(r[3]); cols["id"].append(r[4]); cols["n"].append(r[5]); cols["k"].append(r[6])
    return cols


def r6(x):
    """Round a float to 6 significant digits for the normalized tables; None stays None."""
    if x is None:
        return None
    x = float(x)
    if x != x:
        return None
    return float(f"{x:.6g}")


def nlen(table):
    return len(next(iter(table.values()))) if table else 0


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def save_fragment(dataset_id, frag):
    write_json(os.path.join(FRAG_DIR, dataset_id + ".json"), frag)


def assemble():
    """Merge every fragment in build/fragments/ (sorted by file name) into docs/data.json."""
    names = sorted(f for f in os.listdir(FRAG_DIR) if f.endswith(".json"))
    out = {t: {} for t in TABLES}
    offset = {t: 0 for t in TABLES}
    for name in names:
        with open(os.path.join(FRAG_DIR, name), encoding="utf-8") as f:
            frag = json.load(f)
        for t in TABLES:
            tab = frag.get(t)
            if not tab:
                continue
            base = nlen(out[t])
            for k, col in tab.items():
                ref = REFS.get(t, {}).get(k)
                if ref:
                    col = [v + offset[ref] if isinstance(v, int) and v >= 0 else v for v in col]
                if k not in out[t]:
                    out[t][k] = [None] * base
                out[t][k].extend(col)
            n = nlen(tab)
            for k in out[t]:
                if k not in tab:
                    out[t][k].extend([None] * n)
        for t in TABLES:
            offset[t] = nlen(out[t])
    data = {"schema": SCHEMA, "build": BUILD, "chroms": chrom_lengths(), "hb": HEAT_BIN}
    data.update(out)
    data["genes"] = genes()
    write_json(os.path.join(DOCS, "data.json"), data)
    return data
