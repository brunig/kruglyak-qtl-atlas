"""Ingest Foss et al. 2007 (Nat Genet 39:1369): genetic basis of proteome variation in BY x RM,
protein levels by label-free LC-MS.

    python3 -I ingest/foss2007.py          # from the repository root

Reads only from sources/foss2007/ (and the shared SGD tables). The main text and Online Methods
are paywalled; the Supplementary Information PDF is the only data file. It tabulates just the six
local linkages of Supplementary Table 1 (peptide, protein, linked markers, LOD, marker distance from
the middle of the gene). Genome-wide counts are in the Supplementary Note; the four protein
hotspots are only plotted (Figure 4). The six linkages go into a lazily loaded per-gene table
("mspqtl", shared with Foss 2011's mass-spec pQTL track). Markers are named, not positioned
(ORF-named markers like YJR103W_1873 and chromosome-arm markers like gGR12_1350), so each linkage is
placed by its gene's R64-5-1 position: interval = gene midpoint +/- the table's marker distance (the
marker lies at one end; the side is not given). Deterministic: same inputs, same bytes out.
"""
import collections
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "foss2007", "f07"
SRC = os.path.join(atlas.SOURCES, DS_ID)
RAW = os.path.join(SRC, "raw")
PDF = "41588_2007_BFng200722_MOESM17_ESM.pdf"
ST1 = "../extracted/supp_table1.tsv"
TF = "../extracted/text_facts.tsv"
HB = atlas.HEAT_BIN
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]


def tsv(rel):
    """Rows of an extracted TSV with their 1-based line number in the file (line 1 = comment, line 2 = header)."""
    with open(os.path.join(RAW, rel)) as f:
        lines = f.read().splitlines()
    head = lines[1].split("\t")
    return head, [(i + 1, l.split("\t")) for i, l in enumerate(lines) if i >= 2]


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def marker_chrom(m, G, gid):
    """Chromosome of a marker name: ORF-named markers by their SGD R64-5-1 gene; gXLnn / gXRnn markers by the letter
    (A = chrI ... P = chrXVI), the naming the Brem et al. marker sets use. None if neither applies."""
    name = m.rsplit("_", 1)[0]
    if name in gid:
        return G["c"][gid[name]]
    if len(name) == 5 and name[0] == "g" and name[2] in "LR" and "A" <= name[1] <= "P":
        return "chr" + ROMAN[ord(name[1]) - ord("A")]
    return None


def main():
    sha = manifest()
    for rel in (ST1, TF):
        sha[rel] = hashlib.sha256(open(os.path.join(RAW, rel), "rb").read()).hexdigest()
    C = collections.OrderedDict
    T = {k: C() for k in atlas.TABLES}

    def add(table, **kw):
        tab = T[table]
        n = atlas.nlen(tab)
        for k, v in kw.items():
            if k not in tab:
                tab[k] = [None] * n
            tab[k].append(v)
        for k in tab:
            if k not in kw:
                tab[k].append(None)
        return n

    add("sources", code="f07:ST1", ds=0, file=PDF, sheet="page 4, Supplementary Table 1 (rows = text lines on the page; extracted copy in ../extracted/supp_table1.tsv, column pdf_line)", hrow=1, sha256=sha[PDF])
    add("sources", code="f07:TF", ds=0, file=TF, sheet="text_facts", hrow=2, sha256=sha[TF])
    st1h, st1 = tsv(ST1)
    tfh, tf = tsv(TF)
    tfrow = {r[0]: ln for ln, r in tf}

    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1=None, s2=None, p18=None, org=None)
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1=None, s2=None, p18=None, org=None)
    add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=None, dip=None, mm=None,
        pv=f"f07:TF#{tfrow['design']}")
    add("traits", id=f"{SC}-protein-abundance-mass-spec", ds=0, raw="protein abundance (label-free LC-MS)",
        name="Protein abundance (mass spec, 2007)", cond=None, bt=None, ck="protein abundance", cat="expression",
        catb="atlas assignment: peptide ion-current levels by label-free LC-MS in each segregant (Foss et al. 2007)")

    # ---- Supplementary Table 1 -> lazy "mspqtl" -------------------------------------------------
    G = atlas.genes()
    gid = {}
    for i, s in enumerate(G["id"]):
        gid.setdefault(s, i)
    h = {c: i for i, c in enumerate(st1h)}
    L = {c: [] for c in ["id", "g", "c", "a", "e", "p", "lod", "ef", "loc", "pv", "mk", "dist", "pep"]}
    dens = collections.Counter()
    orig = collections.defaultdict(dict)
    for ln, r in st1:
        orf = r[h["protein"]]
        g = gid[orf]
        c = G["c"][g]
        mid = (G["a"][g] + G["e"][g]) // 2
        dist = float(r[h["marker distance from middle of gene"]])
        markers = [m for m in r[h["linkage marker"]].split(";") if m]
        same = [m for m in markers if marker_chrom(m, G, gid) == c]
        unknown = [m for m in markers if marker_chrom(m, G, gid) is None]
        local = same[0] if same else (unknown[0] if len(unknown) == 1 else None)
        a, e = max(1, int(round(mid - dist))), min(atlas.chrom_lengths()[c], int(round(mid + dist)))
        lid = f"{SC}-p-{int(r[h['pdf_line']])}"
        for k, v in (("id", lid), ("g", g), ("c", c), ("a", a), ("e", e), ("p", mid), ("lod", atlas.r6(float(r[h["lod score"]]))),
                     ("ef", None), ("loc", 1), ("pv", f"f07:ST1#{int(r[h['pdf_line']])}"), ("mk", local), ("dist", dist), ("pep", r[h["peptide"]])):
            L[k].append(v)
        for b in range(a // HB, e // HB + 1):
            dens[("protein, local", c, b)] += 1
        orig[c][lid] = r
    chroms = list(atlas.chrom_lengths())
    for (k, c, b), n in sorted(dens.items(), key=lambda kv: (kv[0][0], chroms.index(kv[0][1]), kv[0][2])):
        add("dens", ds=0, k=k, c=c, b=b, n=n)
    n = len(L["id"])
    add("lazy", ds=0, name="mspqtl", file=f"v/{SC}_mspqtl.json", n=n,
        desc=f"The {n} local protein linkages of Supplementary Table 1 (marker < 20 kb from the gene). Placed by the gene's R64-5-1 position: "
             "interval = gene midpoint ± the table's marker distance; the marker lies at one end, side not given.")

    add("datasets", id=DS_ID, sc=SC, label="Foss 2007",
        cite="Foss EJ, Radulovic D, Shaffer SA, Ruderfer DM, Bedalov A, Goodlett DR, Kruglyak L. 2007. Genetic basis of proteome variation in yeast. Nat Genet 39:1369–1375",
        doi="10.1038/ng.2007.22", year=2007,
        design="segregant panel BY × RM; peptide levels by label-free LC-MS of total unfractionated protein in each segregant; linkage of peptide and transcript levels",
        build="Not stated, and no marker positions are tabulated (markers are named: ORF-based like YJR103W_1873, or chromosome-arm segments like gGR12_1350). "
              "Each local linkage is placed by its gene's SGD R64-5-1 coordinates; the table's marker-to-gene distances are bp-scale and consistent with R64 neighbours (e.g. YLR350W–YLR354C).",
        method="Linkage of peptide abundance (and of transcript abundance) to markers; the detection procedure is in the paywalled Online Methods. "
               "The Supplementary Note counts linkages at LOD > 3 for 163 peptides with a unique protein match and transcript data.",
        fdr=None,
        interval="Not reported. Atlas placement: gene midpoint ± the marker distance from Supplementary Table 1 (the marker lies at one end).",
        effect="Not reported in the supplement; Supplementary Table 1 gives no direction or size.",
        nq=0,
        notes=f"1 cross (segregant count in the paywalled main text, not in the bundle). Only the {n} local linkages are tabulated; genome-wide, 71 of 163 peptides had at least one "
              "linkage at LOD > 3 (93 distinct linkages) and 110 of 163 transcripts (140); 3 of 4 protein hotspots are also transcript hotspots (Figure 4 only, not tabulated).")

    atlas.save_fragment(DS_ID, T)

    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_mspqtl.json"), L)
    for c in chroms:
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{c}.json"), {"cols": {"f07:ST1": st1h}, "rows": {"f07:ST1": dict(orig.get(c, {}))}})
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"),
                     {"f07:ST1": {"cols": st1h, "rows": [r for _, r in st1]}, "f07:TF": {"cols": tfh, "rows": [r for _, r in tf]}})
    data = atlas.assemble()
    print(f"{DS_ID}: {n} local protein linkages; {atlas.nlen(T['dens'])} density bins; data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
