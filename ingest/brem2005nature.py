"""Ingest Brem, Storey, Whittle & Kruglyak 2005 (Nature 436:701), genetic interactions affecting transcript levels.

    python3 -I ingest/brem2005nature.py          # from the repository root

Reads only from sources/: sources/brem2005nature/extracted/text_facts.tsv (statements quoted verbatim from
the PubMed Central author manuscript raw/PMC1409747.html and the Nature page article.html, written and checked
by extracted/extract_facts.py), the MANIFEST, and the SGD tables in sources/sgd/. Writes
build/fragments/brem2005nature.json, docs/v/b05n_*.json and, through ingest/atlas.py, docs/data.json.
Deterministic.

What the source holds: a two-stage, two-locus linkage search over the transcript levels of 112 BY x RM
segregants (data of Brem & Kruglyak 2005 / ref. 7). Significant locus pairs for 225 transcripts at FDR 5%,
estimated involvement of locus pairs in 57% of transcripts, and one named pair: a primary locus near MAT
and a secondary locus containing GPA1 (S469I in BY), linked by 17 transcripts and confirmed with
engineered strains carrying the four MAT x GPA1 allele combinations. No table of locus pairs, positions or
linkage statistics exists in the article or its supplement (Figure 1 is a histogram in 20 kb bins).

Atlas rows: the two loci of the MAT-GPA1 pair, as chromosome-only QTL (a/e/p null: the paper gives no
coordinates, and placing them at the genes would be an atlas reading), qtl.k naming each as primary or
secondary locus of the pair and qtl.cc = 17 linking transcripts. GPA1 (YHR005C) and the MAT locus
(represented by MATALPHA1, YCR040W, the S288C reference gene at MAT) are tier x: the interaction was
confirmed with engineered strains. All counts are kept in datasets.notes and tested.
"""
import collections
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "brem2005nature", "b05n"
SRC = os.path.join(atlas.SOURCES, "brem2005nature")
TF = "../extracted/text_facts.tsv"   # relative to raw/, as sources.file is


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
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

    add("sources", code="b05n:TF", ds=0, file=TF, sheet="text_facts", hrow=1, sha256=sha[TF])

    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        tfr = [l.rstrip("\n").split("\t") for l in f]
    head = tfr[0]
    tf = [(i + 2, dict(zip(head, r))) for i, r in enumerate(tfr[1:])]   # spreadsheet-style row numbers: header is row 1
    TFD = dict(tf)

    def fact(kind, contains):
        return next(n for n, d in tf if d["kind"] == kind and contains in d["quote"])

    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1=None, s2=None, p18=None, org="laboratory strain (S288C derivative)")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1=None, s2=None, p18=None, org="wild strain")
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=112,
            dip=None, mm=None, pv=f"b05n:TF#{fact('design', '112 segregants')}")
    t = add("traits", id=f"{SC}-transcript-levels", ds=0, raw="transcript levels (two-stage two-locus search)",
            name="Transcript levels (two-locus linkage)", cond=None, bt=None, ck="transcript levels", cat="expression",
            catb="atlas assignment: transcript levels of all genes, summarised as one trait")

    G = atlas.genes()
    gix = {g: i for i, g in enumerate(G["id"])}
    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    pair17 = fact("pair", "leaving 17 transcripts")
    LOCI = [("chrIII", "MAT-GPA1 locus pair: primary locus, near MAT (position not given)", fact("pair", "near MAT"), "YCR040W"),
            ("chrVIII", "MAT-GPA1 locus pair: secondary locus, contains GPA1 (position not given)", fact("pair", "secondary locus contained"), "YHR005C")]
    qix = []
    for k, (chrom, kind, row, orf) in enumerate(LOCI):
        qid = f"{SC}-w-{k + 1}"
        q = add("qtl", id=qid, ds=0, t=t, x=x, m="w", c=chrom, a=None, e=None, p=None, cm=None, lod=None, q=None, ef=None,
                ve=None, lv=None, af=None, cc=17, pa=None, anc=None, daf=None, up=-1, ups=None, dns=None, xs=None, k=kind,
                g=gix[orf], gf=None, gt="x", pv=f"b05n:TF#{row}")
        orig[chrom]["b05n:TF"][qid] = [TFD[row][h] for h in head]
        qix.append(q)
    eng = fact("engineered", "four possible combinations")
    for q, orf, name in ((qix[0], "YCR040W", "MAT (MATALPHA1)"), (qix[1], "YHR005C", "GPA1")):
        g = gix[orf]
        cid = f"{SC}-c-{orf}"
        add("causal", id=cid, ds=0, t=t, g=g, orf=orf, name=G["n"][g] or name, c=G["c"][g], a=G["a"][g], e=G["e"][g],
            pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=qix[1] if q == qix[0] else qix[0],
            tier="x", pv=f"b05n:TF#{eng}")
        orig[G["c"][g]]["b05n:TF"][cid] = [TFD[eng][h] for h in head]
        add("cand", ds=0, g=g, tier="x", fdr=None, nt=1, n=1)

    add("datasets", id=DS_ID, sc=SC, label="Brem 2005 (Nature)",
        cite="Brem RB, Storey JD, Whittle J, Kruglyak L. 2005. Genetic interactions between polymorphisms that affect gene "
             "expression in yeast. Nature 436:701-703",
        doi="10.1038/nature03865", year=2005,
        design="segregant panel, BY × RM (112 segregants; genotypes and expression data of ref. 7); two-stage search for "
               "interacting locus pairs underlying each transcript level; MAT-GPA1 interaction tested in engineered strains",
        build="none: no positions are given; the two loci are named by gene (MAT, GPA1), whose R64-5-1 coordinates the "
              "candidate-gene rows carry",
        method="For each transcript, primary QTL = marker with the strongest Wilcoxon rank-sum statistic; segregants split by "
               "inheritance at it and each subgroup scanned for a secondary locus. Joint significance by empirical-Bayes "
               "posterior probabilities (Storey et al.), FDR 5%; interaction tested by an F-test on t = ax + by + cxy + d.",
        fdr=0.05,
        interval="None reported: loci are not given positions; the MAT-GPA1 pair is named by gene.",
        effect="No effect sizes in a table; Figure 2 and Supplementary Figure S1 plot expression by MAT-GPA1 genotype "
               "(RM is MATa, BY is MATα; GPA1 S469I is the BY allele).",
        nq=atlas.nlen(T["qtl"]),
        notes="1 cross (BY × RM), 112 segregants. Locus pairs estimated to be involved in 57% of transcripts; significant "
              "pairs for 225 transcripts at FDR 5%; interaction term p < 0.05 for 65% of them (13% of 547 independently "
              "mapped pairs); 74 (33%) of the 225 secondary loci detectable by a genome scan, 67% missed. The largest pair, "
              "MAT-GPA1, links 17 transcripts and was confirmed with engineered strains. No table of locus pairs is published; "
              "the atlas holds the MAT-GPA1 pair only.")

    atlas.save_fragment(DS_ID, T)
    for chrom in atlas.chrom_lengths():
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"),
                         {"cols": {"b05n:TF": head}, "rows": {"b05n:TF": dict(orig[chrom].get("b05n:TF", {}))}})
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"),
                     {"b05n:TF": {"cols": head, "rows": [[d[h] for h in head] for _, d in tf]}})
    data = atlas.assemble()
    print(f"{DS_ID}: 1 cross (112 segregants), 1 trait, {atlas.nlen(T['qtl'])} QTL (MAT-GPA1 pair, chromosome only), "
          f"{atlas.nlen(T['causal'])} confirmed genes; data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
