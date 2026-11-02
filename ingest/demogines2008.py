"""Ingest Demogines et al. 2008 (PLoS Genet 4:e1000123), DNA-damage sensitivity QTL in BY x RM.

    python3 -I ingest/demogines2008.py          # from the repository root

Reads only from sources/: the article page sources/demogines2008/article.html, the statements
quoted verbatim from it in sources/demogines2008/extracted/text_facts.tsv, and the SGD tables in
sources/sgd/. The four supporting-information files (figures) are inventoried, not parsed.
Writes build/fragments/demogines2008.json, docs/v/d08_*.json and, through ingest/atlas.py,
docs/data.json. Deterministic.

What the files hold. The article has no QTL table; its Results text gives:
- linkage in 121 genotyped BY x RM segregants (R/qtl, non-parametric; 1.5-LOD support intervals
  whose bounds are not printed): MMS and 4-NQO at chrXII marker YLR034C_1989 (LOD 16.6, 10.1),
  bleomycin at YBR138C_275 (LOD 5.2), caffeine at YBR161W_293 (LOD 5.0);
- bulk-segregant analysis of 40 resistant vs 40 sensitive third-backcross segregants (tiling
  arrays): chrVII 433,546-551,683 bp (peak 463,019) and chrXIV 437,935-473,640 bp (peak 467,209).
RAD5 (MMS, 4-NQO) and MKT1 (4-NQO) are confirmed experimentally (tier x).

Placement. The bulk-segregant coordinates are used as printed (BY reference genome of 2008): the
peaks fall inside KAP122 and MKT1 on R64-5-1, the genes the text places at them, so no transform.
Linkage markers are named ORF_offset. A linkage QTL is drawn over the ORF that contains its peak
marker (a, e = SGD R64-5-1 ORF bounds; this is not the support interval, which is not printed),
with the peak p read as `offset` bp from the start codon (most defensible reading: it puts the
chrXII peak ~1.6 kb from RAD5 and ~12 kb from AAT2, as the text's "~2 kb" and "~12.5 kb").
"""
import collections
import hashlib
import html
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "demogines2008", "d08"
SRC = os.path.join(atlas.SOURCES, "demogines2008")
TF = "../extracted/text_facts.tsv"   # relative to raw/, as sources.file is
ART = "../article.html"


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def article_lines():
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    m = re.search(r'<div class="article-text"', t)
    t = t[m.start():] if m else t
    t = re.sub(r"<h[1-6][^>]*>", "\n## ", t)
    t = re.sub(r"</?(p|div|h[1-6]|li|tr|table|br|section)[^>]*>", "\n", t)
    t = html.unescape(re.sub(r"<[^>]+>", "", t))
    t = re.sub(r"[ \t ]+", " ", t)
    return [[i + 1, l.strip()] for i, l in enumerate(x for x in t.split("\n") if x.strip())]


def main():
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

    add("sources", code="d08:TF", ds=0, file=TF, sheet="text_facts", hrow=1,
        sha256=sha(os.path.join(SRC, "extracted", "text_facts.tsv")))
    add("sources", code="d08:ART", ds=0, file=ART, sheet="article page (text lines)", hrow=None,
        sha256=sha(os.path.join(SRC, "article.html")))

    with open(os.path.join(SRC, "extracted", "text_facts.tsv"), encoding="utf-8") as f:
        tfr = [l.rstrip("\n").split("\t") for l in f]
    tf_head = tfr[0]
    tf = [(i + 2, dict(zip(tf_head, r))) for i, r in enumerate(tfr[1:])]
    TFD = dict(tf)

    def fact(kind, contains):
        return next(n for n, d in tf if d["kind"] == kind and contains in d["quote"])

    def quote(row):
        return TFD[row]["quote"]

    # ---- parents, crosses ------------------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1="BY (S288c background)", s2=None, p18=None, org=None)
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1=None, s2="RM", p18=None, org=None)
    assert "121 segregants" in quote(fact("method", "121 segregants")) and "123 haploid" in quote(fact("design", "123"))
    f1 = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM (F1 segregants)", p1=by, p2=rm, g1=None, n=121,
             dip=None, mm=None, pv=f"d08:TF#{fact('method', '121 segregants')}")
    bk = add("crosses", id=f"{SC}-BK3", ds=0, code="BK3", name="third backcross of BY::RAD5-RM × RM (pools)", p1=by, p2=rm,
             g1=None, n=80, dip=None, mm=None, pv=f"d08:TF#{fact('bsa_design', '40 from each set')}")

    # ---- traits ------------------------------------------------------------------------------
    dose = fact("doses", "0.03% MMS")
    catb = "atlas assignment: DNA-damaging agent on YPD agar (Materials and Methods, Plating Assays); growth scored by spot dilution"
    TR = {}
    for key, raw, name, cond, ck in (("mms", "MMS;0.03%", "MMS 0.03%", "0.03% (vol/vol)", "mms 0.03%"),
                                     ("nqo", "4-NQO;0.25 ug/ml", "4-NQO 0.25 µg/mL", "0.25 µg/ml (0.30 µg/ml also used)", "4nqo 0.25ug/mL"),
                                     ("bleo", "bleomycin;10 ug/ml", "Bleomycin 10 µg/mL", "10 µg/ml", "bleomycin 10ug/mL"),
                                     ("caff", "caffeine;10 mM", "Caffeine 10 mM", "10 mM", "caffeine 10mM")):
        TR[key] = add("traits", id=f"{SC}-{key}", ds=0, raw=raw, name=name, cond=cond, bt=None, ck=ck,
                      cat="drug/chemical", catb=catb)

    G = atlas.genes()
    gi_by_id = {g: i for i, g in enumerate(G["id"])}
    gi_by_name = {}
    for i, n in enumerate(G["n"]):
        if n:
            gi_by_name.setdefault(n, i)

    def marker(m):
        orf, off = m.split("_")
        g = gi_by_id[orf]
        a, e, s = G["a"][g], G["e"][g], G["s"][g]
        off = int(off)
        assert off < e - a + 1, m     # the offset lies inside the ORF
        p = a + off if s == "+" else e - off
        return G["c"][g], a, e, p

    # ---- QTL -----------------------------------------------------------------------------------
    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    qrows = []
    LINK = (("mms", "YLR034C_1989", 16.6, by, fact("qtl_mms", "16.6"), "RAD5"),
            ("nqo", "YLR034C_1989", 10.1, by, fact("qtl_4nqo", "10.1"), "RAD5"),
            ("bleo", "YBR138C_275", 5.2, -1, fact("qtl_bleo", "5.2"), None),
            ("caff", "YBR161W_293", 5.0, -1, fact("qtl_caff", "5.0"), None))
    for k, (tk, mk, lod, up, row, gene) in enumerate(LINK):
        assert mk in quote(row) and str(lod) in quote(row)
        c, a, e, p = marker(mk)
        qid = f"{SC}-w-{k + 1}"
        q = add("qtl", id=qid, ds=0, t=TR[tk], x=f1, m="w", c=c, a=a, e=e, p=p, cm=None, lod=lod, q=None, ef=None, ve=None,
                lv=mk, af=None, cc=None, pa=None, anc=None, daf=None, up=up, ups=None, dns=None, xs=None,
                k=f"linkage peak marker {mk}; drawn over its ORF (1.5-LOD support interval not printed)",
                g=gi_by_name[gene] if gene else -1, gf=None, gt="x" if gene else None, pv=f"d08:TF#{row}")
        qrows.append((q, gene, row))
        orig[c]["d08:TF"][qid] = [TFD[row][h] for h in tf_head]
    bsa = fact("bsa_regions", "433,546")
    for k, (c, a, e, p, up, thr) in enumerate((("chrVII", 433546, 551683, 463019, -1, ">1"),
                                               ("chrXIV", 437935, 473640, 467209, rm, ">1.5"))):
        for v in (a, e, p):
            assert f"{v:,}" in quote(bsa)
        qid = f"{SC}-x-{k + 1}"
        gene = "MKT1" if c == "chrXIV" else None
        q = add("qtl", id=qid, ds=0, t=TR["nqo"], x=bk, m="x", c=c, a=a, e=e, p=p, cm=None, lod=None, q=None, ef=None, ve=None,
                lv=None, af=None, cc=None, pa=None, anc=None, daf=None, up=up, ups=None, dns=None, xs=None,
                k=f"bulk-segregant region (smoothed log10 intensity difference {thr}), third-backcross pools",
                g=gi_by_name[gene] if gene else -1, gf=None, gt="x" if gene else None, pv=f"d08:TF#{bsa}")
        qrows.append((q, gene, bsa))
        orig[c]["d08:TF"][qid] = [TFD[bsa][h] for h in tf_head]

    # ---- candidate genes ------------------------------------------------------------------------
    ev = {"RAD5": {TR["mms"]: fact("rad5_mms", "replaced with the BY copy"), TR["nqo"]: fact("rad5_4nqo", "RM copy of RAD5")},
          "MKT1": {TR["nqo"]: fact("mkt1", "MKT1-BY allele contributed")}}
    for q, gene, _ in qrows:
        if not gene:
            continue
        t = T["qtl"]["t"][q]
        g = gi_by_name[gene]
        row = ev[gene][t]
        cid = f"{SC}-c-{T['qtl']['id'][q]}-{G['id'][g]}"
        add("causal", id=cid, ds=0, t=t, g=g, orf=G["id"][g], name=G["n"][g], c=G["c"][g], a=G["a"][g], e=G["e"][g],
            pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier="x", pv=f"d08:TF#{row}")
        orig[G["c"][g]]["d08:TF"][cid] = [TFD[row][h] for h in tf_head]
    cz = T["causal"]
    per_gene = C()
    for i in range(atlas.nlen(cz)):
        per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier="x", fdr=None, nt=len({cz["t"][i] for i in ii}), n=len(ii))

    add("datasets", id=DS_ID, sc=SC, label="Demogines 2008",
        cite="Demogines A, Smith E, Kruglyak L, Alani E. 2008. Identification and Dissection of a Complex DNA Repair "
             "Sensitivity Phenotype in Baker's Yeast. PLoS Genet 4:e1000123",
        doi="10.1371/journal.pgen.1000123", year=2008,
        design="segregant panel, BY × RM (123 segregants phenotyped, 121 genotyped at 2,956 markers) on MMS, 4-NQO, bleomycin "
               "and caffeine; plus bulk-segregant analysis (tiling arrays) of 40 4-NQO-resistant vs 40 sensitive third-backcross "
               "segregants of BY::RAD5-RM × RM",
        build="Bulk-segregant coordinates as printed (BY reference genome, 2008); their peaks fall inside KAP122 and MKT1 on "
              "R64-5-1, the genes the text places there, so no transform. Linkage peaks are named by marker (ORF_offset) and "
              "drawn over the ORF on R64-5-1.",
        method="Linkage: R/qtl, non-parametric model, 121 segregants; significance by 4,000 permutations at p = 0.0125 "
               "(0.05 / 4 traits): LOD 3.9 (MMS), 3.8 (4-NQO), 3.8 (bleomycin), 3.7 (caffeine). Bulk-segregant analysis: "
               "log10 intensity difference between resistant and sensitive pools per variant probe, smoothed; regions where "
               "the smoothed value exceeds 1 (chrVII) or 1.5 (chrXIV).",
        fdr=None,
        interval="Linkage QTL: 1.5-LOD support intervals in the paper, bounds not printed; drawn over the ORF containing the peak "
                 "marker, with the peak read as the marker's offset from the start codon. Bulk-segregant QTL: the printed region.",
        effect="No signed effects. Traits are growth (resistance) on each agent; qtl.up names the parent whose allele raises "
               "resistance where the text says: BY at the RAD5 locus (MMS, 4-NQO), RM at the MKT1 region (4-NQO).",
        nq=atlas.nlen(T["qtl"]),
        notes="2 crosses (F1 segregants, third-backcross pools), 4 traits, 6 QTL: 4 linkage peaks (MMS and 4-NQO at the chrXII "
              "RAD5 locus, bleomycin and caffeine on chrII) and 2 bulk-segregant regions for 4-NQO (chrVII, chrXIV/MKT1). RAD5 "
              "confirmed by allele replacement; MKT1 by deletion and allele-specific plasmid complementation (homologous "
              "replacement failed). KAP122 and PDR1 tested at the chrVII region and did not suppress.")

    atlas.save_fragment(DS_ID, T)

    for chrom in atlas.chrom_lengths():
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"),
                         {"cols": {"d08:TF": tf_head}, "rows": {"d08:TF": dict(orig[chrom].get("d08:TF", {}))}})
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"),
                     {"d08:TF": {"cols": tf_head, "rows": [[d[h] for h in tf_head] for _, d in tf]},
                      "d08:ART": {"cols": ["line", "text"], "rows": article_lines()}})

    data = atlas.assemble()
    print(f"{DS_ID}: 2 crosses, {atlas.nlen(T['traits'])} traits, {atlas.nlen(T['qtl'])} QTL, "
          f"{atlas.nlen(T['causal'])} candidate rows ({atlas.nlen(T['cand'])} genes); "
          f"data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
