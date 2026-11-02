"""Ingest Ronald et al. 2005 (PLoS Genet 1:e25): local regulatory variation ("self-linkage") for
5,727 transcripts in 112 BY x RM segregants, with allele-specific expression in the hybrid for 93 genes.

    python3 -I ingest/ronald2005.py          # from the repository root

Reads only from sources/ronald2005/ (and the shared SGD tables). Table S1 gives, for every tested
ORF, the test at the single marker closest to the gene's start site: a self-linkage flag (p < 0.012,
FDR < 5%), log2 fold change BY/RM with a bootstrap CI, -ln p, and 1- and 2-LOD support intervals in
base pairs of the 2005 SGD assembly. The 1,428 self-linkages are local eQTL by construction (the
tested marker is the one nearest the gene), so each is placed at the start site of its gene on
SGD R64-5-1 (atlas placement; the paper's marker positions and support intervals are kept verbatim,
not converted). They go into a lazily loaded per-gene table "eqtl" (sharing the mRNA eQTL track),
with per-50-kb density rows. Deterministic: same inputs, same bytes out.
"""
import collections
import hashlib
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "ronald2005", "r05"
SRC = os.path.join(atlas.SOURCES, DS_ID)
RAW = os.path.join(SRC, "raw")
P = "journal.pgen.0010025."
HB = atlas.HEAT_BIN
TF = "../extracted/text_facts.tsv"


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def val(v):
    """Verbatim cell: 'NA' and '' become null, numbers stay numbers, everything else stays text."""
    if v in ("", "NA"):
        return None
    try:
        x = float(v)
        return int(x) if x.is_integer() and "." not in v and "e" not in v.lower() else x
    except ValueError:
        return v


def table(name):
    """(raw header cells, header without empty cells, [(row number, cells)]). Row 1 is the header."""
    with open(os.path.join(RAW, P + name), encoding="utf-8") as f:
        lines = [l.rstrip("\r\n") for l in f]
    raw = lines[0].split("\t")
    head = [h for h in raw if h != ""]
    rows = [(i + 2, l.split("\t")) for i, l in enumerate(lines[1:]) if l.strip()]
    return raw, head, rows


def facts():
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        lines = [l.rstrip("\n") for l in f if not l.startswith("#")]
    return lines[0].split("\t"), [(i + 2, l.split("\t")) for i, l in enumerate(lines[1:])]


def main():
    sha = manifest()
    sha[TF] = hashlib.sha256(open(os.path.join(RAW, TF), "rb").read()).hexdigest()
    C = collections.OrderedDict
    T = {k: C() for k in atlas.TABLES}

    def add(table_, **kw):
        tab = T[table_]
        n = atlas.nlen(tab)
        for k, v in kw.items():
            if k not in tab:
                tab[k] = [None] * n
            tab[k].append(v)
        for k in tab:
            if k not in kw:
                tab[k].append(None)
        return n

    s1raw, s1cols, s1rows = table("st001.txt")
    s2raw, s2cols, s2rows = table("st002.txt")
    add("sources", code="r05:S1", ds=0, file=P + "st001.txt", sheet="Table S1 (tab-separated text)", hrow=1, sha256=sha[P + "st001.txt"])
    add("sources", code="r05:S2", ds=0, file=P + "st002.txt", sheet="Table S2 (tab-separated text)", hrow=1, sha256=sha[P + "st002.txt"])
    add("sources", code="r05:TF", ds=0, file=TF, sheet="text_facts", hrow=1, sha256=sha[TF])
    tfh, tfr = facts()
    fact_row = {r[0]: rn for rn, r in tfr}

    # ---- parents, cross, trait ---------------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1="BY4716 (isogenic to S288C)", s2=None, p18=None, org="laboratory")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1=None, s2="RM11-1a (wild vineyard strain)", p18=None, org="vineyard")
    add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=112, dip=None, mm=None,
        pv=f"r05:TF#{fact_row['design']}")
    add("traits", id=f"{SC}-transcript-levels-local", ds=0, raw="transcript levels (self-linkage)", name="Transcript levels (local eQTL)",
        cond=None, bt=None, ck="transcript abundance", cat="expression",
        catb="atlas assignment: microarray transcript levels of Brem & Kruglyak 2005, tested for linkage at the marker nearest each gene")

    # ---- self-linkages (Table S1) -> lazy "eqtl" + density -----------------------------------
    G = atlas.genes()
    gid = {}
    for i, s in enumerate(G["id"]):
        gid.setdefault(s, i)
    ase = {r[0]: r for _, r in s2rows}
    cols = ["id", "g", "c", "a", "e", "p", "lod", "ef", "loc", "pv", "up", "nlp", "ci_l", "ci_r", "ase_ef", "ase_p"]
    L = {c: [] for c in cols}
    dens = collections.Counter()
    orig = collections.defaultdict(dict)
    unplaced = []
    for rn, r in s1rows:
        if r[2] != "1":
            continue
        orf = r[0]
        g = gid.get(orf, -1)
        if g < 0:
            unplaced.append(orf)
            continue
        c = G["c"][g]
        start = G["a"][g] if G["s"][g] != "-" else G["e"][g]   # start site: the marker tested is the one nearest it
        ef = float(r[3])
        a2 = ase.get(orf)
        lid = f"{SC}-e-{rn}"
        row = {"id": lid, "g": g, "c": c, "a": start, "e": start, "p": start, "lod": None, "ef": atlas.r6(ef), "loc": 1,
               "pv": f"r05:S1#{rn}", "up": "BY" if ef > 0 else "RM" if ef < 0 else None,
               "nlp": atlas.r6(float(r[6]) / math.log(10)), "ci_l": atlas.r6(float(r[4])), "ci_r": atlas.r6(float(r[5])),
               "ase_ef": atlas.r6(float(a2[6])) if a2 else None, "ase_p": atlas.r6(float(a2[9])) if a2 else None}
        for k in cols:
            L[k].append(row[k])
        dens[(c, start // HB)] += 1
        orig[c][lid] = [val(v) for v in r]
    chroms = list(atlas.chrom_lengths())
    for (c, b), n in sorted(dens.items(), key=lambda kv: (chroms.index(kv[0][0]), kv[0][1])):
        add("dens", ds=0, k="local", c=c, b=b, n=n)
    nl = len(L["id"])
    add("lazy", ds=0, name="eqtl", file=f"v/{SC}_eqtl.json", n=nl,
        desc=f"The {nl} self-linkages of Table S1 (p < 0.012, FDR < 5%), each a local eQTL tested at the marker nearest the gene and "
             "drawn at the gene's start site on R64-5-1; effect = log2 fold change BY/RM; no LOD (-log10 p in column nlp).")

    # ---- dataset row ---------------------------------------------------------------------------
    nby = sum(1 for u in L["up"] if u == "BY")
    add("datasets", id=DS_ID, sc=SC, label="Ronald 2005",
        cite="Ronald J, Brem RB, Whittle J, Kruglyak L. 2005. Local regulatory variation in Saccharomyces cerevisiae. PLoS Genet 1:e25",
        doi="10.1371/journal.pgen.0010025", year=2005,
        design="segregant panel BY4716 × RM11-1a, 112 segregants (expression and genotypes of Brem & Kruglyak 2005); each of 5,727 transcripts tested for linkage at the single marker closest to its gene's start site; allele-specific expression (TaqMan) in the BY/RM hybrid for 77 self-linking genes and 16 controls",
        build="Table S1's support intervals are base pairs on the 2005 SGD assembly and are kept verbatim, not converted. The atlas places each self-linkage at its gene's start site on SGD R64-5-1, because the test is by definition at the marker nearest that gene (atlas placement, not a published coordinate).",
        method="Linkage test at one marker per transcript (as in Brem et al. 2002), permutation p-values; p < 0.012 corresponds to FDR < 0.05 over 5,727 transcripts (1,428 self-linkages, about 1,357 expected true). Effect = log2 of mean expression with the BY allele over mean with the RM allele, 95% CI by bootstrap over segregants.",
        fdr=0.05,
        interval="Point at the gene's start site on R64-5-1 (atlas placement). The paper's 1- and 2-LOD support intervals (2005 coordinates) are in the source rows.",
        effect="Log2 fold change BY/RM at the marker nearest the gene; positive = the BY allele raises expression.",
        nq=atlas.nlen(T["qtl"]),
        notes=f"1 cross, 112 segregants; 5,727 transcripts tested, {nl} self-linkages loaded per gene on demand ({nby} with the BY allele higher, {nl - nby} with RM higher). "
              f"ASE in the hybrid for 93 genes: 44 of 77 self-linking genes nominally significant, 2 of 16 controls (Table S2; columns ase_ef, ase_p). "
              + (f"{len(unplaced)} self-linking ORFs are not in the R64-5-1 gene table and are not placed: {', '.join(unplaced)}." if unplaced else "Every self-linking ORF is in the R64-5-1 gene table."))

    atlas.save_fragment(DS_ID, T)

    # ---- lazy, verbatim and extra files ---------------------------------------------------------
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_eqtl.json"), L)
    for c in chroms:
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{c}.json"), {"cols": {"r05:S1": s1cols}, "rows": {"r05:S1": dict(orig.get(c, {}))}})
    extra = {"r05:S1": {"cols": s1cols, "header_raw": s1raw, "rows": [[val(v) for v in r] for _, r in s1rows]},
             "r05:S2": {"cols": s2cols, "rows": [[val(v) for v in r] for _, r in s2rows]},
             "r05:TF": {"cols": tfh, "rows": [r for _, r in tfr]}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: {nl} self-linkages placed ({len(unplaced)} unplaced), {atlas.nlen(T['dens'])} density bins; "
          f"data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
