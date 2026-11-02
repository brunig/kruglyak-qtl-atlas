"""Ingest Dimitrov et al. 2009 (Genetics 183:365), mitochondrial genome instability (petite frequency) in BY x RM.

    python3 -I ingest/dimitrov2009.py          # from the repository root

Reads only from sources/: the PubMed Central full-text page sources/dimitrov2009/raw/PMC2746160.html
(Tables 1-4 are parsed from it), sources/dimitrov2009/extracted/text_facts.tsv (statements quoted
verbatim from that page, with section), and the SGD tables in sources/sgd/. Writes
build/fragments/dimitrov2009.json, docs/v/d09_*.json and, through ingest/atlas.py, docs/data.json.
Deterministic.

What the source holds: 122 BY4716 x RM11-1a segregants scored for median petite frequency, linkage
tested at 2957 markers with a Wilcoxon-Mann-Whitney test (statistic -ln p, "linkage likelihood"; 5%
genome-wide threshold 9.3 by permutation). Three QTL: locus 1 on chrXIV, loci 2 and 3 on chrXV. The
paper gives, per locus, the window between the two markers with the highest linkage likelihood
(Figure 3 legend); those windows are the QTL intervals here, and the peak is drawn at the window
midpoint because the single peak marker is not stated (labelled in qtl.k). Linkage statistics are
only plotted. At all three loci the RM allele is associated with lower petite frequency (qtl.up = BY).
Allele replacement confirms SAL1 and MKT1 (locus 1), CAT5 (locus 2) and MIP1 (locus 3): tier x.

Genome build: marker positions are as printed (2009). MKT1, SAL1 and MIP1 lie inside their loci's
windows on SGD R64-5-1 coordinates, so positions are used as R64 with no transform; CAT5 lies 12.5 kb
right of the locus-2 window but inside the 97-kbp region the paper analysed around that peak.
"""
import collections
import hashlib
import html as htmlmod
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "dimitrov2009", "d09"
SRC = os.path.join(atlas.SOURCES, "dimitrov2009")
RAW = os.path.join(SRC, "raw")
PAGE = "PMC2746160.html"
TF = "../extracted/text_facts.tsv"   # relative to raw/, as sources.file is


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def cell(s):
    return re.sub(r"\s+", " ", htmlmod.unescape(re.sub(r"<[^>]+>", "", s))).strip()


def page_tables():
    """Tables 1-4 from the PMC page: [(caption, [[cells], ...])] in page order, every row verbatim."""
    t = open(os.path.join(RAW, PAGE), encoding="utf-8", errors="replace").read()
    out = []
    for m in re.finditer(r"<table.*?</table>", t, re.S):
        before = cell(t[max(0, m.start() - 600):m.start()])
        cap = re.findall(r"(TABLE \d+\.\s*[^.]*)", before)
        rows = [[cell(c) for c in re.findall(r"<t[hd][^>]*>(.*?)</t[hd]>", r, re.S)] for r in re.findall(r"<tr.*?</tr>", m.group(0), re.S)]
        out.append((cap[-1].strip() if cap else "", rows))
    return out


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

    tables = page_tables()
    assert len(tables) == 4, len(tables)
    srcs = {"d09:TF": (TF, "text_facts", 1)}
    for k, (cap, _) in enumerate(tables):
        srcs[f"d09:T{k + 1}"] = (PAGE, cap or f"table {k + 1}", 1)
    for code, (f, sheet, hrow) in srcs.items():
        add("sources", code=code, ds=0, file=f, sheet=sheet, hrow=hrow, sha256=sha[f])

    # ---- text facts -------------------------------------------------------------------------
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        tfr = [l.rstrip("\n").split("\t") for l in f]
    tf_head = tfr[0]
    tf = [(i + 2, dict(zip(tf_head, r))) for i, r in enumerate(tfr[1:])]
    TFD = dict(tf)

    def fact(kind, contains):
        return next(n for n, d in tf if d["kind"] == kind and contains in d["quote"])

    # ---- parents, cross ---------------------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY4716", s1="BY4716 (MATα lys2Δ0)", s2=None, p18=None,
             org="laboratory strain (S288C derivative)")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM11-1a", s1=None,
             s2="RM11-1a (MATa leu2Δ0 ura3Δ0 ho∷loxP-kanMX4-loxP)", p18=None, org="vineyard isolate")
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=122,
            dip=None, mm=None, pv=f"d09:TF#{fact('design', '122 segregants')}")

    # ---- trait ----------------------------------------------------------------------------------
    t = add("traits", id=f"{SC}-petite-frequency", ds=0, raw="median petite frequency",
            name="Petite frequency (spontaneous mtDNA loss)", cond=None, bt=None, ck="petite frequency",
            cat="cellular phenotype",
            catb="atlas assignment: mitochondrial genome instability, as for Ehrenreich 2010's mitochondrial function")

    G = atlas.genes()
    gname = {}
    for gi, n in enumerate(G["n"]):
        if n:
            gname.setdefault(n, gi)

    # ---- QTL ----------------------------------------------------------------------------------
    win = [fact("window", "449,639"), fact("window", "519,776"), fact("window", "937,036")]
    LOCI = [("chrXIV", 449639, 486861, win[0]), ("chrXV", 519776, 546197, win[1]), ("chrXV", 937036, 945781, win[2])]
    dir_row, loci_row = fact("direction", "inheritance of the RM alleles"), fact("loci", "Locus 1 is on chromosome XIV")
    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    qix = []
    for k, (chrom, a, e, row) in enumerate(LOCI):
        qid = f"{SC}-w-{k + 1}"
        q = add("qtl", id=qid, ds=0, t=t, x=x, m="w", c=chrom, a=a, e=e, p=(a + e) // 2, cm=None, lod=None, q=None, ef=None,
                ve=None, lv=None, af=None, cc=None, pa=None, anc=None, daf=None, up=by, ups=None, dns=None, xs=None,
                k=f"locus {k + 1}: window between the two markers with the highest linkage likelihood; peak drawn at its midpoint",
                g=-1, gf=None, gt=None, pv=f"d09:TF#{row}")
        orig[chrom]["d09:TF"][qid] = [TFD[row][h] for h in tf_head]
        qix.append(q)

    # candidate genes confirmed by allele replacement (tier x); MKT1's BY allele acts in the opposite direction
    GENES = [("SAL1", 0, fact("gene", "SAL1 is a quantitative trait gene"), "increase"),
             ("MKT1", 0, fact("gene", "BY MKT1-30G strain was 82%"), "decrease"),
             ("CAT5", 1, fact("gene", "CAT5 is a QTG"), "increase"),
             ("MIP1", 2, fact("gene", "MIP1 is a QTG"), "increase")]
    for gene, li, row, eff in GENES:
        g = gname[gene]
        q = qix[li]
        cid = f"{SC}-c-{G['id'][g]}"
        add("causal", id=cid, ds=0, t=t, g=g, orf=G["id"][g], name=gene, c=G["c"][g], a=G["a"][g], e=G["e"][g],
            pcs=None, sp=by, eff=eff, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier="x", pv=f"d09:TF#{row}")
        orig[G["c"][g]]["d09:TF"][cid] = [TFD[row][h] for h in tf_head]
        if T["qtl"]["g"][q] < 0:   # the locus's card shows its first confirmed gene; all are in causal
            T["qtl"]["g"][q], T["qtl"]["gt"][q] = g, "x"
    for i in range(atlas.nlen(T["causal"])):
        add("cand", ds=0, g=T["causal"]["g"][i], tier="x", fdr=None, nt=1, n=1)

    add("datasets", id=DS_ID, sc=SC, label="Dimitrov 2009",
        cite="Dimitrov LN, Brem RB, Kruglyak L, Gottschling DE. 2009. Polymorphisms in Multiple Genes Contribute to the "
             "Spontaneous Mitochondrial Genome Instability of Saccharomyces cerevisiae S288C Strains. Genetics 183:365-383",
        doi="10.1534/genetics.109.104497", year=2009,
        design="segregant panel, BY4716 × RM11-1a (122 segregants, 2957 markers from Brem and Kruglyak 2005); median petite "
               "frequency (spontaneous loss of functional mtDNA) per segregant",
        build="as printed (2009); MKT1, SAL1 and MIP1 lie inside their loci's windows on SGD R64-5-1 coordinates, so positions "
              "are used as R64 with no transform",
        method="Single-marker test: segregants split by BY or RM allele at each marker, median petite frequencies compared by "
               "Wilcoxon-Mann-Whitney; statistic -ln p ('linkage likelihood'); 5% genome-wide threshold 9.3 from permutations. "
               "Candidate genes tested by allele replacement in BY.",
        fdr=None,
        interval="Window between the two markers with the highest linkage likelihood at each locus (Figure 3 legend); "
                 "the peak is drawn at the window midpoint because the single peak marker is not stated. Linkage statistics "
                 "are only plotted (Figure 3A).",
        effect="No effect sizes in a table. Direction from the text: at all three QTL the RM allele is associated with lower "
               "petite frequency, so the BY allele raises it (qtl.up). By allele replacement, the BY alleles of SAL1, CAT5 and "
               "MIP1 raise petite frequency and the BY allele of MKT1 (MKT1-30D) lowers it (causal.eff, for the shared BY parent).",
        nq=atlas.nlen(T["qtl"]),
        notes="1 cross (BY × RM), 122 segregants, 1 trait (petite frequency). Three QTL (chrXIV; chrXV twice). SAL1, MKT1, "
              "CAT5 and MIP1 confirmed by allele replacement; together the BY alleles of SAL1, CAT5 and MIP1 explain ~94% of "
              "the BY-RM difference when both carry MKT1-30G. The paper's supporting information (genetics.org) was not "
              "reachable.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy files ------------------------------------------------------------------------------
    for chrom in atlas.chrom_lengths():
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"),
                         {"cols": {"d09:TF": tf_head}, "rows": {"d09:TF": dict(orig[chrom].get("d09:TF", {}))}})
    extra = {"d09:TF": {"cols": tf_head, "rows": [[d[h] for h in tf_head] for _, d in tf]}}
    for k, (cap, rows) in enumerate(tables):
        extra[f"d09:T{k + 1}"] = {"cols": rows[0], "rows": rows[1:], "caption": cap}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: 1 cross (122 segregants), 1 trait, {atlas.nlen(T['qtl'])} QTL, {atlas.nlen(T['causal'])} confirmed genes, "
          f"{len(tables)} tables kept; data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
