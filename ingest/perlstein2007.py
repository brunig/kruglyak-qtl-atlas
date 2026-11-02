"""Ingest Perlstein EO et al. 2007 (Nat Genet 39:496), drug-response QTL in BY x RM segregants.

    python3 -I ingest/perlstein2007.py          # from the repository root

Reads only from sources/: the supplementary files in sources/perlstein2007/raw/ (Nature
supplementary information; the main text and Methods are paywalled), their line-for-line
flattenings in sources/perlstein2007/extracted/table_s*.tsv (made by extract_tables.py from the
PDFs), the quoted sentences in extracted/text_facts.tsv, and the SGD tables in sources/sgd/.
Writes build/fragments/perlstein2007.json, docs/v/p07_*.json and, through ingest/atlas.py,
docs/data.json. Deterministic.

Design (abstract): 104 genotyped BY x RM segregants grown with ~100 small molecules ("SMPs"); growth
read as OD600 at time points; linkage analysis. Supplementary Table 3 lists all 219 compound ×
dose × time-point linkages (marker, chromosome, position in what is most defensibly cM, LOD);
Table 4 the 123 collapsed linkages; Table 5 the same collapsed linkages (124 rows; the abstract's
"124 distinct linkages") with base-pair peak and confidence interval. Table 5 is therefore the
drawn table: one qtl row per Table 5 row. Its trait codes ("29_D16HighAT.48") are labelled with
compound, dose and time point from the Table 4 row they align to (same order, same marker, LOD
within rounding); the two Table 5 rows with no Table 4 partner are labelled from Table 3 by marker,
LOD and time point. Traits are compound + dose; the time point goes in qtl.k.

Table 6 gives the eight QTL hotspots (chromosome, base-pair position as tabulated, number of linked
SMPs, and the SMP names): one point qtl row each under a "Drug-response hotspots" trait.

Candidate genes: LEU2 (deletion) for niguldipine and PHO84 (coding polymorphism) for "two
polychlorinated phenols" are experimentally verified (abstract) -> tier x. The only polychlorinated
phenols in Table 1 are dichlorophene and pentachlorophenol, both linked on chrXIII near PHO84, so
PHO84 is attached to those two QTL (reading, labelled). SWH1 (chrI hotspot) and MKS1 (chrXIV
hotspots) are favoured in the Supplementary Note as literature associations "which must ultimately
be verified" -> tier c.

Genome build: not stated. Table 5's base-pair positions sit a median ~260 bp (max ~4.5 kb) from the
R64-5-1 midpoint of the gene that names each marker, every peak lies inside its interval and every
interval inside R64; used as R64 with no transform. Direction of effect is not tabulated.
"""
import collections
import csv
import hashlib
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "perlstein2007", "p07"
SRC = os.path.join(atlas.SOURCES, "perlstein2007")
EXT = os.path.join(SRC, "extracted")
PDF = "41588_2007_BFng1991_MOESM2{}_ESM.pdf"
TABLES = {  # code -> (raw file, flattened TSV in extracted/, caption)
    "p07:S1": (PDF.format("36"), "table_s1.tsv", "Supplementary Table 1: complete list of SMPs"),
    "p07:S3": (PDF.format("38"), "table_s3.tsv", "Supplementary Table 3: complete list of SMP/linkages"),
    "p07:S4": (PDF.format("39"), "table_s4.tsv", "Supplementary Table 4: collapsed list of SMP/linkages"),
    "p07:S5": (PDF.format("40"), "table_s5.tsv", "Supplementary Table 5: collapsed list of confidence intervals"),
    "p07:S6": (PDF.format("41"), "table_s6_hotspots.tsv", "Supplementary Table 6: QTL hotspots (chromosome, position, number of SMPs)"),
    "p07:S6S": (PDF.format("41"), "table_s6_hotspot_smps.tsv", "Supplementary Table 6: SMPs listed under each hotspot (columns by word position)"),
}
TF = "../extracted/text_facts.tsv"   # relative to raw/, as sources.file is
ROM = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def tsv(name):
    """Header and [(line number in the file, row)] for an extracted table (line 1 is a comment, line 2 the header)."""
    with open(os.path.join(EXT, name), encoding="utf-8") as f:
        lines = f.read().split("\n")
    rows = [l.split("\t") for l in lines]
    start = 1 if lines[0].startswith("#") else 0
    head = rows[start]
    return head, [(i + 1, r) for i, r in enumerate(rows) if i > start and r != [""]]


def chrom(n):
    return "chr" + ROM[int(n) - 1]


def unit_ascii(u):
    return "uM" if u == "µM" else u


def base_marker(m):
    return m[:-2] if m.endswith(".1") else m


def same_linkage(a, b):
    """Table 4/3 row a and Table 5 row b describe the same linkage: marker (ignoring '.1'), chromosome, LOD within rounding."""
    return base_marker(a["marker"]) == base_marker(b["marker"]) and a["chrom"] == b["chrom"] and abs(float(a["lod"]) - float(b["lod"])) <= 0.55


def code_time(code):
    m = re.search(r"AT\.(\d+)$", code)
    return m.group(1) if m else None


def align(s5, s4, s3):
    """Label each Table 5 row with a Table 4 row (in order), else a Table 3 row (marker, LOD, time)."""
    out, j = [], 0
    for line5, r in s5:
        t5 = code_time(r["trait"])
        ok = lambda row: same_linkage(row, r) and (t5 is None or row["time_h"] == t5)
        k = next((k for k in range(j, min(j + 3, len(s4))) if ok(s4[k][1])), None)
        if k is not None:
            out.append((line5, r, "p07:S4", s4[k][0], s4[k][1]))
            j = k + 1
            continue
        hits = [(l3, x) for l3, x in s3 if ok(x)]
        if len(hits) != 1:
            raise SystemExit(f"cannot label Table 5 row {line5}: {r}")
        out.append((line5, r, "p07:S3", hits[0][0], hits[0][1]))
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

    for code, (raw, flat, cap) in TABLES.items():
        add("sources", code=code, ds=0, file=raw, sheet=f"{cap} (flattened line for line in extracted/{flat})", hrow=2, sha256=sha[raw])
    add("sources", code="p07:TF", ds=0, file=TF, sheet="sentences quoted from article.html and the supplementary PDFs", hrow=1,
        sha256=hashlib.sha256(open(os.path.join(EXT, "text_facts.tsv"), "rb").read()).hexdigest())

    with open(os.path.join(EXT, "text_facts.tsv"), encoding="utf-8") as f:
        tfl = f.read().split("\n")
    tf_head = tfl[0].split("\t")
    tf = [(i + 1, l.split("\t")) for i, l in enumerate(tfl) if i > 0 and l]

    def fact(contains):
        return next(n for n, r in tf if contains in r[3])

    # ---- parents, cross ---------------------------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1=None, s2=None, p18=None, org=None)
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1=None, s2=None, p18=None, org=None)
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=104,
            dip=None, mm=None, pv=f"p07:TF#{fact('104 genotyped segregants')}")

    # ---- tables --------------------------------------------------------------------------------------
    heads, rows = {}, {}
    for code, (_, flat, _) in TABLES.items():
        h, rs = tsv(flat)
        heads[code], rows[code] = h, rs
    as_dict = lambda code: [(l, dict(zip(heads[code], r))) for l, r in rows[code]]
    s3, s4, s5 = as_dict("p07:S3"), as_dict("p07:S4"), as_dict("p07:S5")
    labelled = align(s5, s4, s3)

    # ---- traits: compound + dose -----------------------------------------------------------------------
    catb = "atlas assignment: small molecule from the SMP panel (Supplementary Table 1), growth read as OD600"
    tidx = {}

    def trait(cmp, dose, unit):
        key = (cmp, dose, unit)
        if key not in tidx:
            tidx[key] = add("traits", id=f"{SC}-" + re.sub(r"[^A-Za-z0-9]+", "-", f"{cmp}-{dose}{unit_ascii(unit)}").strip("-"), ds=0,
                            raw=f"{cmp};{dose} {unit}", name=f"{cmp} {dose} {unit}", cond=f"{dose} {unit}", bt=None,
                            ck=f"{cmp.lower()} {dose}{unit_ascii(unit)}", cat="drug/chemical", catb=catb)
        return tidx[key]

    G = atlas.genes()
    gname = {}
    for gi, n in enumerate(G["n"]):
        if n:
            gname.setdefault(n, gi)

    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    qrows = []
    for line5, r, lcode, lline, lab in labelled:
        c = chrom(r["chrom"])
        t = trait(lab["compound"], lab["dose"], lab["unit"])
        qid = f"{SC}-w-{line5}"
        q = add("qtl", id=qid, ds=0, t=t, x=x, m="w", c=c, a=round(float(r["ci_start"])), e=round(float(r["ci_end"])),
                p=round(float(r["pos_bp"])), cm=None, lod=atlas.r6(float(r["lod"])), q=None, ef=None, ve=None, lv=r["marker"],
                af=None, cc=None, pa=None, anc=None, daf=None, up=-1, ups=None, dns=None, xs=None,
                k=f"OD600 at T={lab['time_h']} h (Table 5 code {r['trait']}; labelled from {lcode} line {lline})",
                g=-1, gf=None, gt=None, pv=f"p07:S5#{line5}")
        qrows.append((q, lab["compound"], c))
        orig[c]["p07:S5"][qid] = rows["p07:S5"][[l for l, _ in rows["p07:S5"]].index(line5)][1]

    # ---- hotspots (Table 6) ------------------------------------------------------------------------------
    th = add("traits", id=f"{SC}-hotspots", ds=0, raw="QTL hotspots (Supplementary Table 6)", name="Drug-response hotspots",
             cond=None, bt=None, ck=None, cat="drug/chemical", catb="atlas assignment: QTL hotspots of the SMP panel")
    smps = collections.defaultdict(list)
    for l, r in as_dict("p07:S6S"):
        smps[r["column"]].append(r["smp"])
    hot = {}
    for l, r in as_dict("p07:S6"):
        c = chrom(r["chrom"])
        pos = int(r["pos_bp"])
        qid = f"{SC}-h-{r['column']}"
        hot[r["column"]] = add("qtl", id=qid, ds=0, t=th, x=x, m="w", c=c, a=pos, e=pos, p=pos, cm=None, lod=None, q=None, ef=None,
                               ve=None, lv=None, af=None, cc=int(r["n_smp"]), pa=None, anc=None, daf=None, up=-1, ups=None, dns=None,
                               xs=None, k=f"drug-response hotspot (Table 6 column {r['column']}): {r['n_smp']} linked SMPs; position as tabulated",
                               g=-1, gf=None, gt=None, pv=f"p07:S6#{l}")
        orig[c]["p07:S6"][qid] = list(rows["p07:S6"][[ll for ll, _ in rows["p07:S6"]].index(l)][1]) + [", ".join(smps[r["column"]])]

    # ---- candidate genes ----------------------------------------------------------------------------------
    def cand(q, gene, tier, fct, t):
        g = gname[gene]
        cid = f"{SC}-c-{T['qtl']['id'][q]}-{gene}"
        add("causal", id=cid, ds=0, t=t, g=g, orf=G["id"][g], name=gene, c=G["c"][g], a=G["a"][g], e=G["e"][g], pcs=None, sp=None,
            eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier=tier, pv=f"p07:TF#{fct}")
        T["qtl"]["g"][q], T["qtl"]["gt"][q] = g, tier
        orig[G["c"][g]]["p07:TF"][cid] = tf[fct - 2][1]

    f_leu2, f_pho84 = fact("deletion of LEU2"), fact("polymorphism in the inorganic phosphate transporter PHO84")
    for q, cmp, c in qrows:
        if cmp == "niguldipine" and c == "chrIII":
            cand(q, "LEU2", "x", f_leu2, T["qtl"]["t"][q])
        if cmp in ("dichlorophene", "pentachlorophenol") and c == "chrXIII":
            cand(q, "PHO84", "x", f_pho84, T["qtl"]["t"][q])
    cand(hot["1"], "SWH1", "c", fact("We favor SWH1"), th)
    for col in ("5", "6"):
        cand(hot[col], "MKS1", "c", fact("we favor MKS1"), th)
    cz = T["causal"]
    per_gene = C()
    for i in range(atlas.nlen(cz)):
        per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier=min((cz["tier"][i] for i in ii), key={"x": 0, "c": 1}.get), fdr=None,
            nt=len({cz["t"][i] for i in ii}), n=len(ii))

    nq = atlas.nlen(T["qtl"])
    add("datasets", id=DS_ID, sc=SC, label="Perlstein 2007",
        cite="Perlstein EO, Ruderfer DM, Roberts DC, Schreiber SL, Kruglyak L. 2007. Genetic basis of individual differences in the response to small-molecule drugs in yeast. Nat Genet 39:496",
        doi="10.1038/ng1991", year=2007,
        design="segregant panel, BY × RM, 104 genotyped segregants grown with a panel of small molecules (SMPs); growth read as OD600 "
               "at time points; QTL per compound, dose and time point",
        build="not stated. Supplementary Table 5's base-pair positions sit within ~4.5 kb (median ~260 bp) of the R64-5-1 midpoint of "
              "the gene that names each marker, and every interval lies inside R64: used as R64 with no transform. Tables 3 and 4 "
              "give a second position column at ~3 kb per unit, most defensibly cM (not used).",
        method="Linkage analysis in 104 segregants (abstract); the Methods are paywalled, so the test statistic, LOD threshold and FDR "
               "are not in the bundle. Tables 3-5 list the linked markers with peak LOD; Table 5 gives base-pair confidence intervals "
               "(definition not stated in the supplement).",
        fdr=None,
        interval="Supplementary Table 5 confidence interval (ci_start-ci_end, base pairs); definition not stated in the bundle. "
                 "Hotspots: a single tabulated position.",
        effect="Not tabulated: no effect size or direction is given in the supplement.",
        nq=nq,
        notes="1 cross (BY × RM), 104 segregants, 92 SMPs listed in Table 1 (the abstract says 100); 124 collapsed linkages to "
              "84 compound names (83 compounds, counting rapamycin in two media once) drawn from Table 5 with their intervals; "
              "all 219 compound × dose × time-point linkages of Table 3 kept verbatim. 8 drug-response hotspots (Table 6). "
              "LEU2 (niguldipine) and PHO84 (two polychlorinated phenols, read as dichlorophene and pentachlorophenol) verified "
              "experimentally; SWH1 and MKS1 favoured at hotspots as unverified literature associations.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy files ------------------------------------------------------------------------------------------
    cols = {"p07:S5": heads["p07:S5"], "p07:S6": heads["p07:S6"] + ["smps (Table 6 column)"], "p07:TF": tf_head}
    for c in atlas.chrom_lengths():
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{c}.json"),
                         {"cols": cols, "rows": {code: dict(orig[c].get(code, {})) for code in cols}})
    extra = {"p07:TF": {"cols": tf_head, "rows": [r for _, r in tf]},
             "p07:ALIGN": {"cols": ["table5_line", "table5_trait", "label_source", "label_line", "compound", "dose", "unit", "time_h"],
                           "rows": [[l5, r["trait"], lc, ll, lab["compound"], lab["dose"], lab["unit"], lab["time_h"]]
                                    for l5, r, lc, ll, lab in labelled]}}
    for code in ("p07:S1", "p07:S3", "p07:S4", "p07:S6S"):
        extra[code] = {"cols": heads[code], "rows": [r for _, r in rows[code]]}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: 1 cross (104 segregants), {len(tidx)} compound-dose traits + hotspots, {len(qrows)} QTL (Table 5) + "
          f"{len(hot)} hotspots, {atlas.nlen(T['causal'])} candidate rows; data.json has {len(data['datasets']['id'])} dataset(s)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
