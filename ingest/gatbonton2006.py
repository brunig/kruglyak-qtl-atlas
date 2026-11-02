"""Ingest Gatbonton et al. 2006 (PLoS Genet 2:e35), telomere length as a quantitative trait.

    python3 -I ingest/gatbonton2006.py          # from the repository root

Reads only from sources/: sources/gatbonton2006/raw/ (Table S1 as XLS; Figure S1, Tables 1-2 and the
correction's Table 1 as TIFF images, not parsed), the article and correction pages, and the extractions
in sources/gatbonton2006/extracted/ (text_facts.tsv: quotes checked verbatim against the pages;
table_s1.tsv: Table S1 cell for cell), plus the SGD tables in sources/sgd/. Writes
build/fragments/gatbonton2006.json, docs/v/g06_*.json and, through ingest/atlas.py, docs/data.json.
Deterministic.

What the files hold: the paper has two parts. (1) A genome-wide survey of deletion strains (72 short,
80 long telomere genes; Table 1, corrected in PLoS Genet 2:e104): context only, an image table, not
drawn. (2) Linkage mapping of telomere length in 122 BY x RM segregants: two QTL, on chromosome 12
(LOD 7.4, 25% of variance, RM allele longer) and chromosome 13 (LOD 4.3, 13%). Positions are given only
in Table 2 (an image) and Figure 6C, so the QTL rows have a chromosome and no coordinates. Candidate
genes are the five the text names in the mapped regions (tier c); Table S1 lists all genes it
considered, kept verbatim.

Genome build: no coordinates are tabulated, so none applies. Gene positions come from SGD R64-5-1.
"""
import collections
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "gatbonton2006", "g06"
SRC = os.path.join(atlas.SOURCES, "gatbonton2006")
RAW = os.path.join(SRC, "raw")
S1 = "journal.pgen.0020035.st001.xls"
TF = "../extracted/text_facts.tsv"     # relative to raw/, as sources.file is
ART, COR = "../article.html", "../correction.html"
CAND_TEXT = ["EST1", "VPS34", "RPS28B", "ARV1", "VPS9"]   # the genes the Results name in the mapped regions


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def sha_of(rel):
    return hashlib.sha256(open(os.path.normpath(os.path.join(RAW, rel)), "rb").read()).hexdigest()


def s1_rows():
    import xlrd
    sh = xlrd.open_workbook(os.path.join(RAW, S1)).sheet_by_name("Sheet1")
    rows = []
    for i in range(sh.nrows):
        v = sh.row_values(i)
        if any(str(x).strip() for x in v):
            rows.append((i + 1, [None if x == "" else x for x in v]))
    return rows


def main():
    sha = manifest()
    for rel in (TF, ART, COR):
        sha[rel] = sha_of(rel)
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

    srcs = {"g06:TF": (TF, "text_facts", 1), "g06:S1": (S1, "Sheet1", 9),
            "g06:ART": (ART, "article page (text)", None), "g06:COR": (COR, "correction page (text)", None)}
    for code, (f, sheet, hrow) in srcs.items():
        add("sources", code=code, ds=0, file=f, sheet=sheet, hrow=hrow, sha256=sha[f])

    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        tfr = [l.rstrip("\n").split("\t") for l in f]
    tf_head = tfr[0]
    tf = [(i + 2, dict(zip(tf_head, r))) for i, r in enumerate(tfr[1:])]
    TFD = dict(tf)

    def fact(kind, contains):
        return next(n for n, d in tf if d["kind"] == kind and contains in d["quote"])

    def frow(n):
        return [TFD[n][h] for h in tf_head]

    # ---- parents, cross, trait -------------------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY4716", s1="BY4716 MATa lys2 (S288C derivative)", s2=None,
             p18=None, org="laboratory strain")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM11-1a", s1=None,
             s2="RM11-1a MATa leu2 ura3 HO::KanMx (from vineyard isolate Bb32(2))", p18=None, org="vineyard strain")
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=122, dip=None, mm=None,
            pv=f"g06:TF#{fact('design', 'The 122 segregants')}")
    t = add("traits", id=f"{SC}-telomere-length", ds=0, raw="telomere length", name="Telomere length", cond=None, bt=None,
            ck="telomere length", cat="cellular phenotype",
            catb="atlas assignment: telomere length (Southern blot of terminal restriction fragments) is a cellular phenotype")

    G = atlas.genes()
    gname = {}
    for gi, n in enumerate(G["n"]):
        if n:
            gname.setdefault(n, gi)

    # ---- QTL: chromosome only ------------------------------------------------------------------
    lod_row, ve_row, dir_row = fact("lod", "7.4 and 4.3"), fact("ve", "25% and 13%"), fact("direction", "RM11 allele conferred longer")
    QTL = [("chrXII", 7.4, 0.25, rm, "chromosome 12 telomere-length QTL (position in Table 2 image and Figure 6C only)"),
           ("chrXIII", 4.3, 0.13, -1, "chromosome 13 telomere-length QTL (position in Table 2 image and Figure 6C only)")]
    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    qidx = {}
    for k, (chrom, lod, ve, up, kind) in enumerate(QTL):
        qid = f"{SC}-w-{k + 1}"
        qidx[chrom] = add("qtl", id=qid, ds=0, t=t, x=x, m="w", c=chrom, a=None, e=None, p=None, cm=None, lod=lod, q=None, ef=None,
                          ve=ve, lv=None, af=None, cc=None, pa=None, anc=None, daf=None, up=up, ups=None, dns=None, xs=None,
                          k=kind, g=-1, gf=None, gt=None, pv=f"g06:TF#{lod_row}")
        orig[chrom]["g06:TF"][qid] = frow(lod_row)

    # ---- candidate genes: the five named in the text, on the QTL chromosome --------------------
    cand_row = fact("candidate", "EST1, VPS34, RPS28B, ARV1, and VPS9")
    for name in CAND_TEXT:
        g = gname[name]
        chrom = G["c"][g]
        q = qidx[chrom]
        cid = f"{SC}-c-{G['id'][g]}"
        add("causal", id=cid, ds=0, t=t, g=g, orf=G["id"][g], name=name, c=chrom, a=G["a"][g], e=G["e"][g],
            pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier="c", pv=f"g06:TF#{cand_row}")
        orig[chrom]["g06:TF"][cid] = frow(cand_row)
    cz = T["causal"]
    for i in range(atlas.nlen(cz)):
        add("cand", ds=0, g=cz["g"][i], tier="c", fdr=None, nt=1, n=1)
    # qtl.g: the text names several genes per locus; leave qtl.g empty rather than pick one (causal rows carry them)

    add("datasets", id=DS_ID, sc=SC, label="Gatbonton 2006",
        cite="Gatbonton T, Imbesi M, Nelson M, Akey JM, Ruderfer DM, Kruglyak L, Simon JA, Bedalov A. 2006. Telomere length as a "
             "quantitative trait: genome-wide survey and genetic mapping of telomere length-control genes in yeast. PLoS Genet 2:e35 "
             "(correction: PLoS Genet 2:e104)",
        doi="10.1371/journal.pgen.0020035", year=2006,
        design="segregant panel, BY4716 × RM11-1a, 122 segregants genotyped at 3,312 microarray markers; telomere length by Southern "
               "blot. Also a genome-wide survey of deletion strains (context only, not drawn).",
        build="none applies: QTL positions are given only in Table 2 (an image) and Figure 6C, so no coordinates are ingested. "
              "Gene positions come from SGD R64-5-1.",
        method="Single-locus linkage with R/qtl over 3,312 markers; significance from 100 permutations of segregant identities "
               "(no permuted LOD above 4; genome-wide corrected p < 0.01 for both loci).",
        fdr=None,
        interval="No interval is tabulated: the loci are named by chromosome in the text, with positions only in an image table "
                 "and a figure. QTL rows have a chromosome and no coordinates.",
        effect="Not signed numerically. Direction from the text: at the chromosome 12 locus the RM11 allele confers longer "
               "telomeres (qtl.up); the chromosome 13 direction is in Table 2 (image) only and is not recorded.",
        nq=atlas.nlen(T["qtl"]),
        notes="1 cross (BY × RM), 122 segregants, 1 trait (telomere length). Two QTL: chromosome 12 (LOD 7.4, 25% of variance) and "
              "chromosome 13 (LOD 4.3, 13%); the abstract says they account for 30%–35% of the variation between the strains. "
              "Candidate genes named in the text: EST1, VPS34, RPS28B, ARV1 (chrXII) and VPS9 (chrXIII), tier c; Table S1's full "
              "gene list is kept verbatim. Deletion survey: 72 short and 80 long telomere genes (Table 1 is an image; the "
              "correction prints the complete Table 1). Listed, not drawn.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy files ------------------------------------------------------------------------------
    for chrom in atlas.chrom_lengths():
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"),
                         {"cols": {"g06:TF": tf_head}, "rows": {"g06:TF": dict(orig[chrom].get("g06:TF", {}))}})
    rows = s1_rows()
    extra = {"g06:TF": {"cols": tf_head, "rows": [frow(n) for n, _ in tf]},
             "g06:S1": {"cols": ["row", "col_A", "col_B", "col_C", "col_D", "col_E"],
                        "rows": [[rn] + v[:5] + [None] * (5 - len(v[:5])) for rn, v in rows]}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: 1 cross (122 segregants), 1 trait, {atlas.nlen(T['qtl'])} QTL (chromosome only), "
          f"{atlas.nlen(T['causal'])} candidate rows; data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
