"""Ingest Ehrenreich et al. 2010 (Nature 464:1039), X-QTL in extremely large BY x RM segregant pools.

    python3 -I ingest/ehrenreich2010.py          # from the repository root

Reads only from sources/: the Supplementary Information PDF in sources/ehrenreich2010/raw/
(Tables S1 and S2 are parsed straight from its page 1, and every text line is kept),
sources/ehrenreich2010/extracted/text_facts.tsv (statements quoted verbatim from the abstract
and the SI figure legends, each with its location), and the SGD tables in sources/sgd/.
Writes build/fragments/ehrenreich2010.json, docs/v/e10_*.json and, through ingest/atlas.py,
docs/data.json. Deterministic.

What the files hold: Table S1 lists the 17 chemicals and their selective doses; Table S2 gives
the only tabulated loci, the two known 4-NQO loci (RAD5, MKT1) with the causal gene's position,
the estimated X-QTL peak and their distance. The other detected loci are drawn in Figures
S5A-S5Q but not tabulated, and the main text is behind the paywall, so they are not ingested
(OPEN_QUESTIONS.md, Ehrenreich 2010 section, Q1).

Genome build: not stated. Table S2's gene spans lie inside the R64-5-1 gene features of RAD5 and
MKT1 and its peaks fall within R64 chromosome lengths, so positions are used as given, with no
transform (Q3).
"""
import collections
import hashlib
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "ehrenreich2010", "e10"
SRC = os.path.join(atlas.SOURCES, "ehrenreich2010")
RAW = os.path.join(SRC, "raw")
SI = "41586_2010_BFnature08923_MOESM285_ESM.pdf"
TF = "../extracted/text_facts.tsv"   # path relative to raw/, as sources.file is
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
S1_COLS = ["Chemical", "Selective dose", "Source of chemical"]
S2_COLS = ["Gene", "Chromosome", "Position (bp)", "Estimated Peak Location (bp)", "Distance from Causal Gene (bp)"]
S1_ROW = re.compile(r"^(?P<chem>.+?) (?P<dose>[\d.]+ ?(?:ug/ml|mg/ml|ng/ml|uM|mM|M|%)) (?P<src>.+)$")
S2_GENE = re.compile(r"^(?P<gene>[A-Z0-9]+) (?P<chr>[IVX]+) (?P<peak>\d+) (?P<dist>\d+)$")
# display names and condition keys (Bloom 2019 style: lower-case compound, dose without spaces)
NAME = {"4-nitroquinoline": "4-NQO", "SDS": "SDS"}   # Figure S5A calls 4-nitroquinoline "4-NQO"
CK = {"4-nitroquinoline": "4-nqo"}


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def pdf_pages():
    import pdfplumber
    with pdfplumber.open(os.path.join(RAW, SI)) as pdf:
        return [(pg.extract_text() or "").split("\n") for pg in pdf.pages]


def unit(d):
    return d.replace("ug/ml", "ug/mL").replace("mg/ml", "mg/mL").replace("ng/ml", "ng/mL")


def show_dose(d):
    d = unit(d)
    return re.sub(r"^([\d.]+)\s*(?=[A-Za-z])", r"\1 ", d).replace(" ug/", " µg/").replace(" uM", " µM")


def slug(s):
    return re.sub(r"-+", "-", re.sub(r"[^A-Za-z0-9]+", "-", s)).strip("-")


def main():
    sha = manifest()
    sha[TF] = hashlib.sha256(open(os.path.join(SRC, "extracted", "text_facts.tsv"), "rb").read()).hexdigest()
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

    srcs = {"e10:S1": (SI, "page 1, Table S1 (text lines)", 1), "e10:S2": (SI, "page 1, Table S2 (text lines)", 1),
            "e10:TF": (TF, "text_facts", 1), "e10:SI": (SI, "all pages (text lines)", None)}
    for code, (f, sheet, hrow) in srcs.items():
        add("sources", code=code, ds=0, file=f, sheet=sheet, hrow=hrow, sha256=sha[f])

    # ---- text facts ---------------------------------------------------------------------
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        tfr = [l.rstrip("\n").split("\t") for l in f]
    tf_head = tfr[0]
    tf = [(i + 2, dict(zip(tf_head, r))) for i, r in enumerate(tfr[1:])]
    row_of = lambda kind, contains: next(n for n, d in tf if d["kind"] == kind and contains in d["quote"])

    # ---- parents and the one cross --------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY4716", s1="BY4716 (named in the Figure S6 legend)",
             s2=None, p18=None, org=None)
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1=None, s2=None, p18=None, org=None)
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=None,
            dip=None, mm=None, pv=f"e10:TF#{row_of('genotyping', 'BY and RM allele-specific')}")

    # ---- Table S1 and S2 from page 1 ----------------------------------------------------------
    pages = pdf_pages()
    p1 = pages[0]
    i0 = p1.index("Chemical Selective dose Source of chemical")
    s1rows = []
    for k in range(i0 + 1, len(p1)):
        if p1[k].startswith("Table S1."):
            break
        m = S1_ROW.match(p1[k])
        if not m:
            raise SystemExit(f"Table S1 line did not parse: {p1[k]!r}")
        s1rows.append((len(s1rows) + 2, [m.group("chem"), m.group("dose"), m.group("src")]))
    j0 = next(k for k, l in enumerate(p1) if l.startswith("(bp) Location (bp) Gene (bp)"))
    s2rows, k = [], j0 + 1
    while k < len(p1) and not p1[k].startswith("Table S2."):
        lo = re.fullmatch(r"(\d+)-", p1[k]); m = S2_GENE.match(p1[k + 1]); hi = re.fullmatch(r"(\d+)", p1[k + 2])
        if not (lo and m and hi):
            raise SystemExit(f"Table S2 lines did not parse: {p1[k:k + 3]!r}")
        s2rows.append((len(s2rows) + 2, [m.group("gene"), m.group("chr"), f"{lo.group(1)}-{hi.group(1)}",
                                         int(m.group("peak")), int(m.group("dist"))]))
        k += 3

    # ---- traits: the 17 chemicals, plus mitochondrial function (abstract) ----------------------
    tidx = {}
    for rn, (chem, dose, _src) in s1rows:
        name = (NAME.get(chem) or chem[0].upper() + chem[1:]) + " " + show_dose(dose)
        ck = (CK.get(chem) or chem.lower()) + " " + unit(dose).replace(" ", "")
        cat, catb = "drug/chemical", "atlas assignment: chemical resistance selection (Table S1)"
        tidx[chem] = add("traits", id=f"{SC}-{slug(chem)}", ds=0, raw=f"{chem};{dose}", name=name, cond=dose, bt=None,
                         ck=ck, cat=cat, catb=catb)
    mito = add("traits", id=f"{SC}-mitochondrial-function", ds=0, raw="mitochondrial function", name="Mitochondrial function",
               cond=None, bt=None, ck="mitochondrial function", cat="cellular phenotype",
               catb="atlas assignment: the abstract's 'mitochondrial function' (Mitotracker Red staining, Figure S6); no existing category fits")

    G = atlas.genes()
    gname = {}
    for gi, n in enumerate(G["n"]):
        if n:
            gname.setdefault(n, gi)

    # ---- QTL: the two loci of Table S2 (4-NQO) -------------------------------------------------
    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    nqo = tidx["4-nitroquinoline"]
    for rn, r in s2rows:
        gene, chrr, span, peak, dist = r
        chrom = "chr" + chrr
        g = gname[gene]
        qid = f"{SC}-x-{rn}"
        q = add("qtl", id=qid, ds=0, t=nqo, x=x, m="x", c=chrom, a=peak, e=peak, p=peak, cm=None, lod=None, q=None,
                ef=None, ve=None, lv=None, af=None, cc=None, pa=None, anc=None, daf=None, up=-1, ups=None, dns=None,
                xs=None, k="known 4-NQO locus (Table S2)", g=g, gf=None, gt="v", pv=f"e10:S2#{rn}")
        lo, hi = (int(v) for v in span.split("-"))
        cid = f"{SC}-c-{rn}"
        add("causal", id=cid, ds=0, t=nqo, g=g, orf=G["id"][g], name=gene, c=chrom, a=lo, e=hi, pcs=None, sp=None,
            eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier="v", pv=f"e10:S2#{rn}")
        orig[chrom]["e10:S2"][qid] = r
        orig[chrom]["e10:S2"][cid] = r

    # HAP1 and mitochondrial function: Figure S6 compares BY4716 with BY4716 carrying a functional HAP1
    hr = row_of("candidate", "functional HAP1")
    g = gname["HAP1"]
    cid = f"{SC}-c-tf{hr}"
    add("causal", id=cid, ds=0, t=mito, g=g, orf=G["id"][g], name="HAP1", c=G["c"][g], a=G["a"][g], e=G["e"][g],
        pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=-1, q2=-1, tier="x", pv=f"e10:TF#{hr}")
    orig[G["c"][g]]["e10:TF"][cid] = [dict(tf)[hr][k] for k in tf_head]

    cz = T["causal"]
    rank = {"x": 0, "v": 1, "5": 2, "c": 3, "20": 4}
    per_gene = C()
    for i in range(atlas.nlen(cz)):
        per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier=min((cz["tier"][i] for i in ii), key=rank.get), fdr=None,
            nt=len({cz["t"][i] for i in ii}), n=len(ii))

    add("datasets", id=DS_ID, sc=SC, label="Ehrenreich 2010",
        cite="Ehrenreich IM, Torabi N, Jia Y, Kent J, Martis S, Shapiro JA, Gresham D, Caudy AA, Kruglyak L. 2010. "
             "Dissection of genetically complex traits with extremely large pools of yeast segregants. Nature 464:1039-1042",
        doi="10.1038/nature08923", year=2010,
        design="X-QTL (bulk segregant): very large BY × RM segregant pools (more than 10^5 haploids) selected under 17 chemicals, "
               "pools genotyped on BY/RM allele-specific microarrays; also mitochondrial function (Mitotracker Red)",
        build="not stated. Table S2's gene spans for RAD5 and MKT1 lie inside the R64-5-1 gene features and its peaks are inside "
              "R64 chromosome lengths; positions are used as given, no transform.",
        method="Allele frequencies in selected versus control pools from allele-specific microarray probes; loci called at a global "
               "FDR of 0.05 (Figure S5 legends). The main text describing the statistic is behind the paywall.",
        fdr=0.05,
        interval="None reported. Only Table S2 tabulates loci, each as an estimated peak position; drawn as points.",
        effect="Not reported for the X-QTL. (Figure S4 reports effects of the 14 4-NQO loci in 96 separately genotyped segregants as a plot; RAD5 ~40% of variance.)",
        nq=atlas.nlen(T["qtl"]),
        notes=f"1 cross (BY × RM), pools not counted, {len(s1rows)} chemicals plus mitochondrial function. Only the two known 4-NQO loci "
              "(RAD5, MKT1; Table S2) have coordinates in the files; the loci for the other traits are drawn in Figures S5A-S5Q "
              "but not tabulated, and the main text is paywalled. HAP1 for mitochondrial function from Figure S6.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy files: verbatim rows ------------------------------------------------------------
    for chrom in atlas.chrom_lengths():
        out = {"cols": {"e10:S2": S2_COLS, "e10:TF": tf_head},
               "rows": {"e10:S2": dict(orig[chrom].get("e10:S2", {})), "e10:TF": dict(orig[chrom].get("e10:TF", {}))}}
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"), out)
    extra = {"e10:S1": {"cols": S1_COLS, "rows": [r for _, r in s1rows]},
             "e10:S2": {"cols": S2_COLS, "rows": [r for _, r in s2rows]},
             "e10:TF": {"cols": tf_head, "rows": [[d[k] for k in tf_head] for _, d in tf]},
             "e10:SI": {"cols": ["page", "line", "text"],
                        "rows": [[pi + 1, li + 1, l] for pi, pg in enumerate(pages) for li, l in enumerate(pg)]}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: 1 cross, {len(tidx) + 1} traits ({len(s1rows)} chemicals + mitochondrial function), "
          f"{atlas.nlen(T['qtl'])} tabulated loci, {atlas.nlen(T['causal'])} candidate rows ({atlas.nlen(T['cand'])} genes); "
          f"data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
