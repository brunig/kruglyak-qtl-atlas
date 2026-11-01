"""Ingest Albert et al. 2018 (eLife 7:e35471), eQTL in 1,012 BY x RM segregants.

    python3 -I ingest/albert2018.py          # from the repository root

Reads only from sources/: the unpacked eLife bundle in sources/albert2018/raw/ (zip checksum
in sources/albert2018/SOURCE_ZIP.tsv), sources/albert2018/extracted/text_facts.tsv (statements
transcribed from the article, each with its quote) and the SGD tables in sources/sgd/.
Writes build/fragments/albert2018.json, docs/v/a18_*.json and, through ingest/atlas.py,
docs/data.json. Deterministic.

Display design (chosen by the user): the 102 trans-eQTL hotspots are ordinary QTL rows; the
36,498 individual eQTL go to the lazily loaded docs/v/a18_eqtl.json (table `lazy`), with
per-50-kb counts for the heatmap in table `dens`. No trait per gene.

Semantics (joshsbloom/eQTL_BYxRM at commit ad83786, code/eQTL_BYxRM1000_fx_stranded.R, and the
article): genotypes are coded -1 = BY, +1 = RM (README; Source data 3 caption), so `r`, the
Pearson correlation between residual expression and genotype at the peak marker, is positive
when the RM allele raises expression. eQTL intervals are 1.5-LOD drops extended over markers in
perfect LD. `cis` is the paper's local eQTL: interval overlaps the gene padded by 1 kb upstream
and 200 bp downstream (and, in the code, peak within 40 kb). Hotspot intervals are the central
95% of 1,000 bootstrap peaks; `bootstrapPeak` is the median bootstrap position. Marker positions
are sacCer3 (README: "postitions are based on the S.Cerevisiae SacCer3 genome build"), i.e. R64.
"""
import collections
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "albert2018", "a18"
SRC = os.path.join(atlas.SOURCES, "albert2018")
RAW = os.path.join(SRC, "raw")
P = "elife_poa_e35471_"
TF = "../extracted/text_facts.tsv"   # relative to raw/, as sources.file is
HB = atlas.HEAT_BIN

# provenance code -> (file, sheet)
SHEETS = {
    "a18:S2": ("Source_data_2.xlsx", "SI_Data_02_covariates"),
    "a18:S4": ("Source_data_4.xlsx", "SI_Data_04_eQTL"),
    "a18:S5A": ("Source_data_5.xlsx", "heritability_A"),
    "a18:S5AA": ("Source_data_5.xlsx", "heritability_A_AA"),
    "a18:S6": ("Source_data_6.xlsx", "SI_Data_06_geneFeatures"),
    "a18:S7": ("Source_data_7.xlsx", "SI_Data_07_ASE_localeQTL"),
    "a18:S8": ("Source_data_8.xlsx", "hotspotOverview"),
    "a18:S8L": ("Source_data_8.xlsx", "legend"),
    "a18:S9": ("Source_data_9.xlsx", "SI_Data_09_hotspotEffects"),
    "a18:S11": ("Source_data_11.xlsx", "TFBS_allResults4paper_171009"),
    "a18:S12BP": ("Source_data_12.xlsx", "BiologicalProcess"),
    "a18:S12MF": ("Source_data_12.xlsx", "MolecularFunction"),
    "a18:S12CC": ("Source_data_12.xlsx", "CellularProcess"),
    "a18:S13E": ("Source_data_13.xlsx", "pQTL_effects_at_eQTL"),
    "a18:S13P": ("Source_data_13.xlsx", "eQTL_effects_at_pQTL"),
    "a18:S14A": ("Source_data_14.xlsx", "between_additive_loci.scan"),
    "a18:S14M": ("Source_data_14.xlsx", "marginal2d.scan"),
    "a18:S14F": ("Source_data_14.xlsx", "full2d.scan"),
    "a18:S14L": ("Source_data_14.xlsx", "legend"),
}
EXTRA = ["a18:S2", "a18:S5A", "a18:S5AA", "a18:S6", "a18:S7", "a18:S8L", "a18:S12BP", "a18:S12MF", "a18:S12CC",
         "a18:S13E", "a18:S13P", "a18:S14A", "a18:S14M", "a18:S14F", "a18:S14L"]


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


_books = {}


def read(code):
    """(column names, [(spreadsheet row, values)]) for a provenance code; header on row 1."""
    f, sheet = SHEETS[code]
    if f not in _books:
        import openpyxl
        wb = openpyxl.load_workbook(os.path.join(RAW, P + f), read_only=True, data_only=True)
        _books[f] = {ws.title: [list(r) for r in ws.iter_rows(values_only=True)] for ws in wb.worksheets}
    rows = _books[f][sheet]
    width = max(len(r) for r in rows)
    rows = [r + [None] * (width - len(r)) for r in rows]
    cols = ["" if v is None else str(v) for v in rows[0]]
    return cols, [(i + 1, r) for i, r in enumerate(rows) if i >= 1]


def verbatim(r):
    return [v.item() if hasattr(v, "item") else v for v in r]


def mpos(m):
    """'chrI:48890_C/T' -> 48890."""
    return int(str(m).split(":")[1].split("_")[0])


def mchr(m):
    return str(m).split(":")[0]


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

    for code, (f, sheet) in SHEETS.items():
        add("sources", code=code, ds=0, file=P + f, sheet=sheet, hrow=1, sha256=sha[P + f])
    add("sources", code="a18:TF", ds=0, file=TF, sheet="text_facts", hrow=1, sha256=sha[TF])

    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        tfr = [l.rstrip("\n").split("\t") for l in f]
    tf_head = tfr[0]
    tf = [(i + 2, dict(zip(tf_head, r))) for i, r in enumerate(tfr[1:])]
    design_row = [rn for rn, d in tf if d["kind"] == "design"][0]

    # ---- parents, cross, trait ---------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1="BY (laboratory strain; S288C reference background)", s2=None, p18=None, org=None)
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1="RM (wine strain)", s2=None, p18=None, org=None)
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=by, n=1012,
            dip=None, mm=None, pv=f"a18:TF#{design_row}")
    t = add("traits", id=f"{SC}-trans-eQTL-hotspots", ds=0, raw="trans-eQTL hotspots", name="Expression hotspots (trans-eQTL)",
            cond=None, bt=None, ck="expression hotspots", cat="expression",
            catb="atlas assignment: loci affecting the mRNA levels of many genes (Albert et al. 2018, Source data 8)")

    G = atlas.genes()
    gid, gname = {}, {}
    for i, (s, n) in enumerate(zip(G["id"], G["n"])):
        gid.setdefault(s, i)
        if n:
            gname.setdefault(n, i)

    orig = collections.defaultdict(lambda: collections.defaultdict(dict))

    # ---- hotspots (Source data 8) ------------------------------------------------------
    hcols, hrows = read("a18:S8")
    hi = {c: i for i, c in enumerate(hcols)}
    hq = []
    for rn, r in hrows:
        chrom = r[hi["chromosome"]]
        a, e, p = int(r[hi["bootstrapIntervalLeft"]]), int(r[hi["bootstrapIntervalRight"]]), int(r[hi["bootstrapPeak"]])
        qid = f"{SC}-w-{rn}"
        q = add("qtl", id=qid, ds=0, t=t, x=x, m="w", c=chrom, a=a, e=e, p=p, lod=None, q=None, ef=None, ve=None,
                lv=r[hi["hotspotMarker"]], af=None, cc=int(r[hi["numberEQTLInHotspot"]]), pa=None, anc=None, daf=None,
                up=-1, ups=None, dns=None, xs=None, k="trans-eQTL hotspot", g=-1, gf=None, gt=None, pv=f"a18:S8#{rn}")
        hq.append((q, chrom, a, e, r))
        orig[chrom]["a18:S8"][qid] = verbatim(r)

    # ---- candidate genes named in the text, placed on the hotspot whose interval contains them
    rank = {"x": 0, "c": 1}
    unplaced = []
    for rn, d in tf:
        if d["kind"] != "candidate":
            continue
        g = gname[d["gene"]]
        hits = [(q, r) for q, chrom, a, e, r in hq if chrom == G["c"][g] and a <= G["e"][g] and e >= G["a"][g]]
        if not hits:
            unplaced.append(d["gene"]); continue
        for q, r in hits:
            assert d["gene"] in str(r[hi["allGenesInInterval"]]).split(";"), d["gene"]
            cid = f"{SC}-c-{rn}-{T['qtl']['id'][q].rsplit('-', 1)[1]}"
            add("causal", id=cid, ds=0, t=t, g=g, orf=G["id"][g], name=d["gene"], c=G["c"][g], a=G["a"][g], e=G["e"][g],
                pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier=d["tier"], pv=f"a18:TF#{rn}")
            cur = T["qtl"]["gt"][q]
            if cur is None or rank[d["tier"]] < rank[cur]:
                T["qtl"]["g"][q] = g; T["qtl"]["gt"][q] = d["tier"]
            orig[G["c"][g]]["a18:TF"][cid] = [d[k] for k in tf_head]
    cz = T["causal"]
    per_gene = C()
    for i in range(atlas.nlen(cz)):
        per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier=min((cz["tier"][i] for i in ii), key=rank.get), fdr=None,
            nt=len({cz["t"][i] for i in ii}), n=len(ii))

    # ---- the 36,498 eQTL (Source data 4): lazy file + density ----------------------------
    ecols, erows = read("a18:S4")
    ei = {c: i for i, c in enumerate(ecols)}
    L = {k: [] for k in ["id", "g", "c", "a", "e", "p", "lod", "ef", "loc", "pv"]}
    dens = collections.Counter()
    nogene = set()
    for rn, r in erows:
        gene, pm = r[ei["gene"]], r[ei["pmarker"]]
        chrom = r[ei["chr"]]
        a, e, p = mpos(r[ei["CI.l"]]), mpos(r[ei["CI.r"]]), mpos(pm)
        loc = 1 if r[ei["cis"]] is True else 0
        eid = f"{SC}-e-{rn}"
        g = gid.get(gene, -1)
        if g < 0:
            nogene.add(gene)
        for k, v in (("id", eid), ("g", g), ("c", chrom), ("a", a), ("e", e), ("p", p), ("lod", atlas.r6(r[ei["LOD"]])),
                     ("ef", atlas.r6(r[ei["r"]])), ("loc", loc), ("pv", f"a18:S4#{rn}")):
            L[k].append(v)
        for b in range(a // HB, e // HB + 1):
            dens[("local" if loc else "distant", chrom, b)] += 1
        orig[chrom]["a18:S4"][eid] = verbatim(r)
    chroms = list(atlas.chrom_lengths())
    for (k, chrom, b), n in sorted(dens.items(), key=lambda kv: (kv[0][0], chroms.index(kv[0][1]), kv[0][2])):
        add("dens", ds=0, k=k, c=chrom, b=b, n=n)
    add("lazy", ds=0, name="eqtl", file=f"v/{SC}_eqtl.json", n=len(L["id"]),
        desc="All 36,498 eQTL at FDR 5% (Source data 4): the gene whose expression is affected, peak marker, 1.5-LOD interval, LOD, r (positive = RM allele raises expression) and whether the eQTL is local.")

    # ---- dataset row --------------------------------------------------------------------
    nloc = sum(L["loc"])
    add("datasets", id=DS_ID, sc=SC, label="Albert 2018",
        cite="Albert FW, Bloom JS, Siegel J, Day L, Kruglyak L. 2018. Genetics of trans-regulatory variation in gene expression. eLife 7:e35471",
        doi="10.7554/eLife.35471", year=2018,
        design="segregant panel; 1,012 BY × RM segregants, RNA-seq (5,720 expressed genes) and low-coverage genotyping (42,052 markers, 11,530 after collapsing perfect LD); eQTL per gene, plus multivariate mapping of trans-eQTL hotspots",
        build="S288C R64 (sacCer3): the analysis code's README states marker positions are on the SacCer3 build; all positions fall inside R64 chromosomes. No transform.",
        method="Per gene and chromosome, forward stepwise linkage on residual log2(TPM+0.5) (batch, growth, large background eQTL on other chromosomes and a polygenic term regressed out); LOD = -n ln(1-r^2)/(2 ln 10); permutation FDR 5% (1,000 permutations). "
               "Hotspots: multivariate (SVD) linkage per chromosome with a 99% permutation threshold, refined, then 1,000 bootstraps.",
        fdr=0.05,
        interval="eQTL (lazy file): 1.5-LOD drop, extended to markers in perfect LD. Hotspots: central 95% of 1,000 bootstrap peak positions; the peak is the median bootstrap position.",
        effect="eQTL: r, the Pearson correlation between residual expression and genotype coded -1 (BY) / +1 (RM) at the peak marker; positive = the RM allele raises expression. Hotspots: no effect size in the overview table (per-gene hotspot effects are in Source data 9).",
        nq=atlas.nlen(T["qtl"]),
        notes=f"1 cross, 1,012 segregants, {len(L['id']):,} eQTL for {len(set(r[ei['gene']] for _, r in erows)):,} genes ({nloc:,} local, {len(L['id']) - nloc:,} distant), "
              f"{len(hq)} trans-eQTL hotspots. Hotspot rows carry the number of eQTL whose interval overlaps the hotspot in qtl.cc. "
              "Candidate genes are those the paper names for a hotspot: tier x for the causal genes it calls experimentally confirmed in earlier work (MKT1, HAP1, IRA2, GPA1), tier c for the others; not claimed causal beyond the paper.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy files ---------------------------------------------------------------------------
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_eqtl.json"), L)
    for chrom in chroms:
        out = {"cols": {"a18:S8": hcols, "a18:S4": ecols, "a18:TF": tf_head},
               "rows": {code: dict(orig[chrom].get(code, {})) for code in ("a18:S8", "a18:S4", "a18:TF")}}
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"), out)
    extra = {}
    for code in EXTRA:
        cols, rows = read(code)
        extra[code] = {"cols": cols, "rows": [verbatim(r) for _, r in rows]}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)
    for code, name in (("a18:S9", "hotspot_effects"), ("a18:S11", "tfbs")):
        cols, rows = read(code)
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{name}.json"), {code: {"cols": cols, "rows": [verbatim(r) for _, r in rows]}})

    data = atlas.assemble()
    size = os.path.getsize(os.path.join(atlas.DOCS, "v", f"{SC}_eqtl.json"))
    print(f"{DS_ID}: {len(hq)} hotspots, {len(L['id'])} eQTL ({nloc} local), {atlas.nlen(T['dens'])} density bins, "
          f"{atlas.nlen(T['causal'])} candidate rows ({atlas.nlen(T['cand'])} genes; not inside any hotspot: {unplaced}); "
          f"eQTL genes not in the R64-5-1 table: {len(nogene)}; a18_eqtl.json {size:,} bytes; "
          f"data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
