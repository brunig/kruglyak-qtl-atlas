"""Ingest Torabi & Kruglyak 2012 (PLoS Genet 8:e1002546): QTL revealed by increased translational readthrough.

    python3 -I ingest/torabi2012.py          # from the repository root

Reads only from sources/: the PLOS supporting information in sources/torabi2012/raw/, the article page
sources/torabi2012/article.html, and sources/torabi2012/extracted/ (Tables S1-S3 written cell for cell
from the .doc files, and statements quoted verbatim from the article, both by extract_sources.py).
Writes build/fragments/torabi2012.json, docs/v/t12_*.json and, through ingest/atlas.py, docs/data.json.
Deterministic.

Design: X-QTL in BY × RM. For each of nine conditions, a wildtype segregant pool and a pool from parents
that both carry sup35 C653R (raised stop-codon readthrough) were selected and compared; loci whose
allele-frequency skew differs between the two pools are "readthrough-dependent". The paper reports 18
such locus-condition pairs (ten distinct loci), drawn in Figure 2 and Figure S2 only. The text names one
locus, on the right arm of chromosome XIII, for diamide and hydrogen peroxide, with SKY1 confirmed as
causal by allele replacement. Those two rows are ingested, with a chromosome and no coordinates; the rest
is recorded as not tabulated.

Genome build: none applies (no coordinates). SKY1's position comes from the SGD R64-5-1 table.
"""
import collections
import hashlib
import html
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "torabi2012", "t12"
SRC = os.path.join(atlas.SOURCES, "torabi2012")
EXT = os.path.join(SRC, "extracted")
# sources.file is relative to raw/, as in the other datasets
TF, S1, S2, S3, ART = "../extracted/text_facts.tsv", "../extracted/table_s1.tsv", "../extracted/table_s2.tsv", "../extracted/table_s3.tsv", "../article.html"
DISPLAY = {"CoCl2": "Cobalt chloride", "H2O2": "Hydrogen peroxide", "E6-berbamine": "E6-berbamine"}


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def tsv(name):
    with open(os.path.join(SRC, "raw", name)) as f:
        lines = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
    first = 1 + sum(1 for l in open(os.path.join(SRC, "raw", name)) if l.startswith("#"))
    return lines[0], [(first + 1 + i, r) for i, r in enumerate(lines[1:])]   # file line numbers


def article_lines():
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace").read()
    m = re.search(r'<div class="article-text"[^>]*>(.*)', t, re.S); t = m.group(1) if m else t
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    t = re.sub(r"<h[1-6][^>]*>", "\n## ", t); t = re.sub(r"</(p|h[1-6]|div|li|tr)>", "\n", t)
    t = html.unescape(re.sub(r"<[^>]+>", " ", t))
    return [[i + 1, re.sub(r"\s+", " ", l).strip()] for i, l in enumerate(x for x in t.split("\n") if x.strip())]


def main():
    sha = manifest()
    for f in (TF, S1, S2, S3, ART):
        sha[f] = hashlib.sha256(open(os.path.join(SRC, "raw", f), "rb").read()).hexdigest()
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

    srcs = {"t12:TF": (TF, "text_facts", 1), "t12:S1": (S1, "Table S1 (from s007.doc)", 3),
            "t12:S2": (S2, "Table S2 (from s008.doc)", 3), "t12:S3": (S3, "Table S3 (from s009.doc)", 3),
            "t12:ART": (ART, "article page (text lines)", None)}
    for code, (f, sheet, hrow) in srcs.items():
        add("sources", code=code, ds=0, file=f, sheet=sheet, hrow=hrow, sha256=sha[f])

    tf_head, tf = tsv(TF)
    TFD = {n: dict(zip(tf_head, r)) for n, r in tf}

    def fact(topic):
        return next(n for n, r in tf if dict(zip(tf_head, r))["topic"] == topic)

    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY4724", s1="BY4724 (a laboratory strain)", s2=None, p18=None, org="laboratory strain")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM11_1a", s1=None, s2="RM11_1a (a wine strain)", p18=None, org="wine strain")
    pools = f"t12:TF#{fact('pools')}"
    xwt = add("crosses", id=f"{SC}-WT", ds=0, code="WT", name="BY × RM, wildtype SUP35", p1=by, p2=rm, g1=None, n=None, dip=None, mm=None, pv=pools)
    xsup = add("crosses", id=f"{SC}-sup35", ds=0, code="sup35", name="BY × RM, both parents sup35 C653R", p1=by, p2=rm, g1=None, n=None, dip=None, mm=None, pv=pools)

    s1_head, s1 = tsv(S1)
    tidx = {}
    for ln, r in s1:
        cond, rng, fig1, xconc, unit, source = r
        if xconc in ("", "NA"):
            continue
        dose = xconc + ("" if unit in ("V/V", "W/V", "NA") else " " + unit)
        if unit == "V/V" and not xconc.endswith("%"):
            dose = xconc + " %"
        base = DISPLAY.get(cond, cond[0].upper() + cond[1:])
        ckname = {"CoCl2": "cobalt chloride", "H2O2": "hydrogen peroxide"}.get(cond, cond.lower().replace("-", " "))
        tidx[cond] = add("traits", id=f"{SC}-{re.sub(r'[^A-Za-z0-9]+', '-', cond).strip('-')}", ds=0, raw=f"{cond};{xconc} {unit}",
                         name=f"{base} {dose}", cond=f"{xconc} {unit}", bt=None, ck=f"{ckname} {xconc.replace(' ', '')}{'' if unit in ('V/V', 'W/V') else unit}",
                         cat="drug/chemical", catb="atlas assignment: agent added to YPD selection plates (Methods, X-QTL); dose and unit as printed in Table S1")
    assert len(tidx) == 9, tidx

    G = atlas.genes()
    sky1 = next(i for i, n in enumerate(G["n"]) if n == "SKY1")

    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    rows = [("Diamide", fact("chrXIII")), ("H2O2", fact("hydrogen peroxide"))]
    for k, (cond, frow) in enumerate(rows):
        qid = f"{SC}-x-{k + 1}"
        q = add("qtl", id=qid, ds=0, t=tidx[cond], x=xsup, m="x", c="chrXIII", a=None, e=None, p=None, cm=None, lod=None, q=None, ef=None,
                ve=None, lv=None, af=None, cc=None, pa=None, anc=None, daf=None, up=rm, ups=None, dns=None, xs=None,
                k="readthrough-dependent locus (sup35 vs wildtype pool), right arm of chrXIII; position drawn in Figures 2 and S3 only",
                g=sky1, gf=None, gt="x", pv=f"t12:TF#{frow}")
        orig["chrXIII"]["t12:TF"][qid] = [TFD[frow][h] for h in tf_head]
        cid = f"{SC}-c-{k + 1}"
        add("causal", id=cid, ds=0, t=tidx[cond], g=sky1, orf=G["id"][sky1], name="SKY1", c=G["c"][sky1], a=G["a"][sky1], e=G["e"][sky1],
            pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier="x", pv=f"t12:TF#{fact('SKY1')}")
        orig["chrXIII"]["t12:TF"][cid] = [TFD[fact("SKY1")][h] for h in tf_head]
    add("cand", ds=0, g=sky1, tier="x", fdr=None, nt=2, n=2)

    add("datasets", id=DS_ID, sc=SC, label="Torabi 2012",
        cite="Torabi N, Kruglyak L. 2012. Genetic Basis of Hidden Phenotypic Variation Revealed by Increased Translational Readthrough in Yeast. PLoS Genet 8:e1002546",
        doi="10.1371/journal.pgen.1002546", year=2012,
        design="X-QTL (bulk segregant) in BY × RM: for each of nine conditions, a wildtype pool and a pool from parents both carrying "
               "sup35 C653R (increased stop-codon readthrough), two biological replicates each, genotyped on an allele-specific array (~18,000 SNPs)",
        build="none applies: no locus coordinates are tabulated (chromosome only). SKY1's position is from SGD R64-5-1.",
        method="X-QTL allele-frequency skews in selected vs control pools per cross (FDR 5%); readthrough-dependent loci called from the "
               "difference between wildtype and sup35 skews with a Savitzky-Golay smoothing approach (FDR 5%) and by a t-test "
               "(Bonferroni-corrected p < 0.05).",
        fdr=0.05,
        interval="Not tabulated: loci are drawn in Figure 2 and Figure S2 only. The two rows here carry a chromosome and no coordinates.",
        effect="Not signed numerically. Direction from the text: at the chromosome XIII locus the RM allele is favoured in the sup35 pool "
               "(qtl.up = RM); SKY1's RM allele increases growth on diamide and hydrogen peroxide under sup35 C653R.",
        nq=atlas.nlen(T["qtl"]),
        notes="2 crosses (BY × RM wildtype; BY × RM with sup35 C653R in both parents), 9 conditions. The paper reports 18 readthrough-"
              "dependent locus-condition pairs (ten distinct loci), drawn only; the atlas lists the one locus the text names "
              "(chromosome XIII, diamide and hydrogen peroxide) with SKY1 confirmed by allele replacement (tier x). Raw arrays: GEO GSE33817.")

    atlas.save_fragment(DS_ID, T)

    for chrom in atlas.chrom_lengths():
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"),
                         {"cols": {"t12:TF": tf_head}, "rows": {"t12:TF": dict(orig[chrom].get("t12:TF", {}))}})
    extra = {"t12:TF": {"cols": tf_head, "rows": [r for _, r in tf]},
             "t12:ART": {"cols": ["line", "text"], "rows": article_lines()}}
    for code, f in (("t12:S1", S1), ("t12:S2", S2), ("t12:S3", S3)):
        h, rs = tsv(f)
        extra[code] = {"cols": h, "rows": [r for _, r in rs]}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: 2 crosses, {len(tidx)} traits, {atlas.nlen(T['qtl'])} QTL (chromosome only), SKY1 tier x; "
          f"data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
