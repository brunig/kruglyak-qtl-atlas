"""Ingest Storey, Akey & Kruglyak 2005 (PLoS Biol 3:e267): two-locus joint linkage of genome-wide
expression in 112 BY x RM segregants.

    python3 -I ingest/storey2005.py          # from the repository root

Reads only from sources/storey2005/ (and the shared SGD tables). The supporting tables are PDFs;
sources/storey2005/extracted/extract.py turned them into TSVs that keep every printed column plus
the PDF page and line, and recorded each PDF's SHA-256, which this script checks against MANIFEST.tsv.

Table S1 lists the 170 traits with joint linkage to a pair of loci at FDR <= 10%. Each locus is a
single marker position, so the 340 loci go into a lazily loaded per-gene table (lazy name "eqtl",
docs/v/s05_eqtl.json, shared with the other eQTL datasets' mRNA track), two rows per trait
(lk = 1 primary, lk = 2 secondary), each with the pair's joint q-value; `loc` follows the paper's
cis flag. There are no qtl rows: a locus pair per transcript is not a trait-level QTL, and the
paper tabulates no single-locus results. Table S2 (10 pair-wise bins with >= 3 linkages and their
putative functions) is kept verbatim. Positions are used as given (R64-compatible; OPEN_QUESTIONS).
Deterministic: same inputs, same bytes out.
"""
import collections
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "storey2005", "s05"
SRC = os.path.join(atlas.SOURCES, DS_ID)
RAW = os.path.join(SRC, "raw")
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
HB = atlas.HEAT_BIN
EX = {"s05:S1": ("../extracted/s1_table.tsv", "journal.pbio.0030267.st001.pdf", "Table S1 (PDF), extracted"),
      "s05:S2": ("../extracted/s2_table.tsv", "journal.pbio.0030267.st002.pdf", "Table S2 (PDF), extracted"),
      "s05:TF": ("../extracted/text_facts.tsv", None, "text_facts")}


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def tsv(rel, pdf, sha):
    """Read an extracted TSV; check the source PDF checksum it records. Rows are numbered as in a
    spreadsheet whose row 1 is the header (the leading comment line is not counted)."""
    lines = open(os.path.join(RAW, rel)).read().splitlines()
    assert lines[0].startswith("#"), rel
    if pdf:
        assert f"sha256={sha[pdf]}" in lines[0], f"{rel} was extracted from a different {pdf}"
    head = lines[1].split("\t")
    return head, [(i + 2, l.split("\t")) for i, l in enumerate(lines[2:])]


def chrom(n):
    return "chr" + ROMAN[int(n) - 1]


def main():
    sha = manifest()
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

    tabs = {}
    for code, (rel, pdf, sheet) in EX.items():
        tabs[code] = tsv(rel, pdf, sha)
        add("sources", code=code, ds=0, file=rel, sheet=sheet, hrow=1,
            sha256=hashlib.sha256(open(os.path.join(RAW, rel), "rb").read()).hexdigest())
    s1h, s1 = tabs["s05:S1"]
    s2h, s2 = tabs["s05:S2"]
    tfh, tf = tabs["s05:TF"]
    design_row = [rn for rn, r in tf if r[0] == "design"][0]

    # ---- parents, cross, trait --------------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY4716", s1="BY4716 (isogenic to S288C)", s2=None, p18=None, org="laboratory")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM11-1a", s1=None, s2="RM11-1a (wild isolate)", p18=None, org=None)
    add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY4716 × RM11-1a", p1=by, p2=rm, g1=None, n=112, dip=None, mm=None,
        pv=f"s05:TF#{design_row}")
    add("traits", id=f"{SC}-transcript-levels-two-locus", ds=0, raw="gene expression traits, two-locus joint linkage",
        name="Transcript levels (two-locus joint linkage)", cond=None, bt=None, ck="transcript abundance", cat="expression",
        catb="atlas assignment: microarray expression of 6,216 ORFs in 112 segregants (Storey et al. 2005)")

    # ---- Table S1 -> lazy eQTL table (two loci per trait) + density ---------------------------
    G = atlas.genes()
    gid = {}
    for i, s in enumerate(G["id"]):
        gid.setdefault(s, i)
    h = {c: i for i, c in enumerate(s1h)}
    L = {c: [] for c in ["id", "g", "c", "a", "e", "p", "lod", "ef", "loc", "lk", "q", "pv"]}
    dens = collections.Counter()
    orig = collections.defaultdict(dict)
    for rn, r in s1:
        cis = int(r[h["Cis"]])
        g = gid.get(r[h["Trait"]], -1)
        for lk in (1, 2):
            c = chrom(r[h[f"Locus {lk} Chr"]])
            p = int(r[h[f"Locus {lk} Pos"]])
            loc = 1 if cis == lk else 0
            lid = f"{SC}-e-{rn}-{lk}"
            for k, v in (("id", lid), ("g", g), ("c", c), ("a", p), ("e", p), ("p", p), ("lod", None), ("ef", None), ("loc", loc),
                         ("lk", lk), ("q", atlas.r6(float(r[h["q-value Joint"]]))), ("pv", f"s05:S1#{rn}")):
                L[k].append(v)
            dens[("local" if loc else "distant", c, p // HB)] += 1
            orig[c][lid] = list(r)
    chroms = list(atlas.chrom_lengths())
    for (k, c, b), n in sorted(dens.items(), key=lambda kv: (kv[0][0], chroms.index(kv[0][1]), kv[0][2])):
        add("dens", ds=0, k=k, c=c, b=b, n=n)
    n = len(L["id"])
    add("lazy", ds=0, name="eqtl", file=f"v/{SC}_eqtl.json", n=n,
        desc=f"The {n} loci of Table S1: for each of 170 transcripts with joint linkage to a pair of loci at FDR <= 10%, "
             "the primary (lk = 1) and secondary (lk = 2) locus as single marker positions, with the pair's joint q-value. "
             "No LOD, effect or interval is tabulated.")

    # ---- dataset row --------------------------------------------------------------------------
    ncis = sum(1 for _, r in s1 if r[h["Cis"]] != "0")
    nmiss = len({r[h["Trait"]] for _, r in s1 if r[h["Trait"]] not in gid})
    add("datasets", id=DS_ID, sc=SC, label="Storey 2005",
        cite="Storey JD, Akey JM, Kruglyak L. 2005. Multiple locus linkage analysis of genomewide expression in yeast. PLoS Biol 3:e267",
        doi="10.1371/journal.pbio.0030267", year=2005,
        design="segregant panel BY4716 × RM11-1a, 112 segregants; microarray expression of 6,216 ORFs; sequential two-locus (multiple locus) linkage scan with epistasis",
        build="S288C as printed in 2005. Every marker lies inside its R64 chromosome, and the trait genes' printed coordinates match their R64-5-1 positions to within the extent of the gene model (no per-chromosome shift), so positions are used as R64 with no transform; agreement is at gene/marker resolution.",
        method="For each trait, the primary locus is the best single marker and the secondary locus the best marker given the primary (allowing a locus1 × locus2 interaction); significance by permutation; Bayesian joint-linkage probabilities, with FDR among traits called jointly linked estimated from them. Table S1: joint linkage at FDR <= 10%.",
        fdr=0.10,
        interval="None: each locus is a single marker position (a = e = p).",
        effect="Not reported: Table S1 gives positions, joint q-value and a cis flag only, no effect size or direction.",
        nq=atlas.nlen(T["qtl"]),
        notes=f"1 cross, 112 segregants, 3,312 markers, 6,216 expression traits. 170 traits with joint two-locus linkage at FDR <= 10% "
              f"({n} loci, on demand per gene; {ncis} traits with a cis locus by the paper's flag). {nmiss} trait ORFs are not in the R64-5-1 gene table and keep g = -1. "
              "Table S2's 10 pair-wise 50-kb bins with three or more linkages (mating, amino acid and mitochondria metabolism, mitochondrial TCA cycle, response to stress, ...) are kept verbatim, not drawn.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy, verbatim and extra files --------------------------------------------------------
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_eqtl.json"), L)
    for c in chroms:
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{c}.json"), {"cols": {"s05:S1": s1h}, "rows": {"s05:S1": dict(orig.get(c, {}))}})
    extra = {"s05:S1": {"cols": s1h, "rows": [r for _, r in s1]},
             "s05:S2": {"cols": s2h, "rows": [r for _, r in s2]},
             "s05:TF": {"cols": tfh, "rows": [r for _, r in tf]}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: 170 traits, {n} loci ({sum(L['loc'])} local), {atlas.nlen(T['dens'])} density bins; "
          f"data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
