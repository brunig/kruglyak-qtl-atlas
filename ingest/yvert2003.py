"""Ingest Yvert et al. 2003 (Nat Genet 35:57): trans-acting regulatory variation in 86 BY x RM
segregants.

    python3 -I ingest/yvert2003.py          # from the repository root

Reads only from sources/yvert2003/ (and the shared SGD tables). The main text is paywalled;
the two supplementary PDFs are extracted to TSVs by sources/yvert2003/extracted/extract.py:
  - Supplementary Table 2 summary: 13 groups of messages linking to loci with widespread
    transcriptional effects (trans-eQTL hotspots), each with its linkage bin as printed
    (e.g. II:560000), message count and putative regulator  -> 13 qtl rows
  - Supplementary Table 2 member lists: 1,265 messages, one per (group, gene)
    -> lazy per-gene table "eqtl" (docs/v/y03_eqtl.json), each a distant linkage to its group's bin,
       with per-50-kb density rows for the heatmap
  - Supplementary Table 1: 593 co-expression clusters -> kept verbatim (not drawn)
Positions are the printed linkage-bin labels (20-kb bins; a group may merge up to four bins), used
as given: every named regulator lies on the bin's chromosome within 48 kb on R64-5-1. Points
(a = e = p), no transform. Deterministic: same inputs, same bytes out.
"""
import collections
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "yvert2003", "y03"
SRC = os.path.join(atlas.SOURCES, DS_ID)
RAW = os.path.join(SRC, "raw")
EXT = "../extracted/"
SI1 = "41588_2003_BFng1222_MOESM1_ESM.pdf"
SI2 = "41588_2003_BFng1222_MOESM2_ESM.pdf"
HB = atlas.HEAT_BIN
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]

# Regulators the supplement names per group, with the evidence it gives (quotes in text_facts.tsv).
# tier x = "we confirmed" (positional cloning and functional assays, per the abstract); c = named as putative/candidate.
CANDIDATES = [
    (2, "AMN1", "x", "Group 2 contains the cluster of 12 genes shown in Figure 1 of the main text, whose regulator we confirmed to be Amn1."),
    (2, "MAK5", "c", "A good candidate is the RNA helicase MAK5"),
    (4, "LEU2", "c", "Group 4 contains 100 messages linking to a region near the gene LEU2 (an auxotrophic marker deleted in the RM strain)."),
    (6, "URA3", "c", "Group 6 contains 19 messages linking to the URA3 locus (another auxotrophic marker deleted in the RM strain)."),
    (7, "GPA1", "x", "Group 7 contains the cluster of 7 pheromone response genes shown in Figure 1 of the main text, whose regulator we confirmed to be GPA1."),
    (8, "HAP1", "c", "We showed previously that our BY strain carries a Ty insertion in HAP1"),
    (9, "SIR3", "c", "The messages in Group 9 link to a region on chromosome XII containing SIR3."),
    (13, "CAT5", "c", "A good candidate in the linkage region is CAT5"),
]
PRIMARY = {2: "AMN1", 4: "LEU2", 6: "URA3", 7: "GPA1", 8: "HAP1", 9: "SIR3", 13: "CAT5"}


def tsv(name):
    """Return (header, [(file line number, row)]) for an extracted TSV (line 1 is a source comment)."""
    with open(os.path.join(SRC, "extracted", name)) as f:
        lines = f.read().split("\n")
    head = lines[1].split("\t")
    rows = [(k + 1, l.split("\t")) for k, l in enumerate(lines) if k >= 2 and l != ""]
    return head, rows


def sha_of(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


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

    gs_h, gs_rows = tsv("groups_summary.tsv")
    gm_h, gm_rows = tsv("group_members.tsv")
    cl_h, cl_rows = tsv("clusters.tsv")
    tf_h, tf_rows = tsv("text_facts.tsv")
    for code, f, sheet in (("y03:SI1", SI1, "Supplementary Table 1 (PDF)"), ("y03:SI2", SI2, "Supplementary Table 2 (PDF)")):
        add("sources", code=code, ds=0, file=f, sheet=sheet, hrow=0, sha256=sha[f])
    for code, f, sheet in (("y03:GS", "groups_summary.tsv", "Supplementary Table 2 summary, extracted"),
                           ("y03:GM", "group_members.tsv", "Supplementary Table 2 group members, extracted"),
                           ("y03:CL", "clusters.tsv", "Supplementary Table 1 clusters, extracted"),
                           ("y03:TF", "text_facts.tsv", "quoted sentences")):
        add("sources", code=code, ds=0, file=EXT + f, sheet=sheet, hrow=2, sha256=sha_of(os.path.join(SRC, "extracted", f)))

    def fact_row(quote):
        return next(rn for rn, r in tf_rows if r[1] == quote)

    # ---- parents, cross, trait ---------------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1="BY (laboratory strain)", s2=None, p18=None, org="laboratory")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1="RM (wild isolate)", s2=None, p18=None, org=None)
    xq = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=86, dip=None, mm=None,
             pv=f"y03:TF#{fact_row('86 segregants: 1 culture each')}")
    t = add("traits", id=f"{SC}-trans-eQTL-hotspots", ds=0, raw="groups of messages linking to loci with widespread transcriptional effects",
            name="Expression hotspots (trans-eQTL)", cond=None, bt=None, ck="expression hotspots", cat="expression",
            catb="atlas assignment: Supplementary Table 2 groups of transcripts linking to one locus (trans-acting hotspots)")

    # ---- 13 hotspot groups -> qtl rows ---------------------------------------------------------
    G = atlas.genes()
    gid, gname = {}, {}
    for i, s in enumerate(G["id"]):
        gid.setdefault(s, i)
        if G["n"][i]:
            gname.setdefault(G["n"][i], i)
    h = {c: i for i, c in enumerate(gs_h)}
    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    qidx, gpos, extra_lines = {}, {}, []
    for rn, r in gs_rows:
        if not r[h["group"]]:
            extra_lines.append((rn, r))   # continuation line of the printed table (second function / regulator)
            continue
        grp = int(r[h["group"]])
        c = "chr" + r[h["chrom_roman"]]
        p = int(r[h["position"]])
        reg = PRIMARY.get(grp)
        g = gname.get(reg, -1) if reg else -1
        tier = next((tr for gg, nm, tr, _ in CANDIDATES if gg == grp and nm == reg), None)
        qid = f"{SC}-w-{rn}"
        qidx[grp] = add("qtl", id=qid, ds=0, t=t, x=xq, m="w", c=c, a=p, e=p, p=p, lod=None, q=None, ef=None, ve=None,
                        lv=None, af=None, cc=int(r[h["n_msgs"]]), pa=None, anc=None, daf=None, up=-1,
                        g=g, gf=None, gt=tier if g >= 0 else None,
                        k=f"trans-eQTL hotspot (group {grp}; linkage bin as printed)", pv=f"y03:GS#{rn}")
        gpos[grp] = (c, p)
        orig[c]["y03:GS"][qid] = r
    # second regulator line under group 2 ("rRNA processing MAK5") stays verbatim with the group's row
    for rn, r in extra_lines:
        orig["chrII"]["y03:GS"][f"{SC}-w-{rn}"] = r

    # ---- candidate genes ------------------------------------------------------------------------
    for grp, name, tier, quote in CANDIDATES:
        g = gname[name]
        add("causal", id=f"{SC}-c-{grp}-{name}", ds=0, t=t, g=g, orf=G["id"][g], name=name, c=G["c"][g], a=G["a"][g], e=G["e"][g],
            pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=qidx[grp], q2=-1, tier=tier,
            pv=f"y03:TF#{fact_row(quote)}")
    per_gene = C()
    for i, g in enumerate(T["causal"]["g"]):
        per_gene.setdefault(g, []).append(i)
    rank = {"x": 0, "v": 1, "5": 2, "c": 3, "20": 4}
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier=min((T["causal"]["tier"][i] for i in ii), key=rank.get), fdr=None,
            nt=len({T["causal"]["t"][i] for i in ii}), n=len(ii))

    # ---- 1,265 member messages -> lazy eqtl table + density ------------------------------------
    hm = {c: i for i, c in enumerate(gm_h)}
    L = {k: [] for k in ["id", "g", "c", "a", "e", "p", "lod", "ef", "loc", "grp", "pv"]}
    dens = collections.Counter()
    missing = []
    for rn, r in gm_rows:
        grp = int(r[hm["group"]])
        c, p = gpos[grp]
        orf = r[hm["orf"]]
        g = gid.get(orf, -1)
        if g < 0:
            missing.append(orf)
        eid = f"{SC}-e-{rn}"
        for k, v in (("id", eid), ("g", g), ("c", c), ("a", p), ("e", p), ("p", p), ("lod", None), ("ef", None),
                     ("loc", 0), ("grp", grp), ("pv", f"y03:GM#{rn}")):
            L[k].append(v)
        dens[("distant", c, p // HB)] += 1
        orig[c]["y03:GM"][eid] = r
    chroms = list(atlas.chrom_lengths())
    for (k, c, b), n in sorted(dens.items(), key=lambda kv: (kv[0][0], chroms.index(kv[0][1]), kv[0][2])):
        add("dens", ds=0, k=k, c=c, b=b, n=n)
    add("lazy", ds=0, name="eqtl", file=f"v/{SC}_eqtl.json", n=len(L["id"]),
        desc=f"The {len(L['id'])} messages of Supplementary Table 2, each linking to one of the 13 trans-acting groups "
             "(cis-linkages removed by the authors): regulated gene and the group's linkage bin as printed. No LOD, interval or direction is tabulated.")

    # ---- dataset row --------------------------------------------------------------------------
    add("datasets", id=DS_ID, sc=SC, label="Yvert 2003",
        cite="Yvert G, Brem RB, Whittle J, Akey JM, Foss E, Smith EN, Mackelprang R, Kruglyak L. 2003. Trans-acting regulatory variation in Saccharomyces cerevisiae and the role of transcription factors. Nat Genet 35:57–64",
        doi="10.1038/ng1222", year=2003,
        design="segregant panel BY × RM, 86 segregants plus parents; microarray expression (two dye-swapped hybridizations per sample); linkage of transcript levels; trans-acting loci found as 20-kb bins with excess linkages",
        build="S288C; not stated (2003). Positions are the printed linkage-bin labels; every named regulator (AMN1, GPA1, LEU2, URA3, HAP1, SIR3, CAT5) lies on the bin's chromosome within 48 kb on R64-5-1, consistent with R64 at bin resolution. Used as given, no transform.",
        method="Main-text methods are paywalled. Supplementary Table 2: genome divided into 20-kb bins; messages linking to each bin counted after removing cis-linkages; bins with more than 12 messages significant at p < 0.01 by random reassignment (22 bins), adjacent bins merged into 13 groups.",
        fdr=None,
        interval="Linkage bin label as printed (20-kb bins; Group 2 merges four bins, Group 11 three, Groups 3, 7, 8, 12 two). Drawn as points at the label: the table does not say which merged bin the label marks.",
        effect="Not reported: no LOD, effect size or direction is tabulated for the groups or their messages.",
        nq=atlas.nlen(T["qtl"]),
        notes=f"1 cross, 86 segregants; 13 trans-acting groups (hotspots) with {len(L['id'])} member messages, loaded per gene on demand; "
              f"593 co-expression clusters kept verbatim. Regulators named by the paper: AMN1 and GPA1 confirmed (tier x); MAK5, LEU2, URA3, HAP1, SIR3, CAT5 putative (tier c); "
              f"Group 5 maps to the MAT locus (no single gene). {len(missing)} member ORFs are not in the R64-5-1 gene table (kept with g = -1).")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy, verbatim and extra files -------------------------------------------------------
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_eqtl.json"), L)
    for c in chroms:
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{c}.json"),
                         {"cols": {"y03:GS": gs_h, "y03:GM": gm_h},
                          "rows": {"y03:GS": dict(orig[c]["y03:GS"]), "y03:GM": dict(orig[c]["y03:GM"])}})
    extra = {"y03:CL": {"cols": cl_h, "rows": [r for _, r in cl_rows]},
             "y03:TF": {"cols": tf_h, "rows": [r for _, r in tf_rows]}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: {atlas.nlen(T['qtl'])} hotspot groups, {len(L['id'])} member messages ({len(missing)} ORFs not in R64-5-1: "
          f"{', '.join(sorted(set(missing))[:8])}{'…' if len(set(missing)) > 8 else ''}), {len(cl_rows)} clusters, "
          f"{atlas.nlen(T['dens'])} density bins; data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
