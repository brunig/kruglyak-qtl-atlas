"""Ingest Boocock et al. 2025 (eLife 13:RP95566), single-cell ("one-pot") eQTL in three yeast crosses.

    python3 -I ingest/boocock2025.py          # from the repository root

Reads only from sources/: the eLife supplementary files in sources/boocock2025/raw/ (fetched
from cdn.elifesciences.org; the article page refused plain requests, so the article text is the
PMC full text, sources/boocock2025/article_pmc.xml, from NCBI E-utilities), the verbatim quotes in
sources/boocock2025/extracted/text_facts.tsv, and the SGD tables in sources/sgd/. Writes
build/fragments/boocock2025.json, docs/v/b25_*.json and, through ingest/atlas.py, docs/data.json.
Deterministic.

Design (as Albert 2018): the 25 trans-eQTL hotspots and the 20 cell-cycle occupancy QTL are
ordinary QTL rows; the individual eQTL (significant local eQTL, FDR < 5%, and all distant eQTL)
go to the lazily loaded docs/v/b25_sceqtl.json (lazy table name "sceqtl"), with per-50-kb counts
in table `dens`. Every tested local row is kept verbatim in docs/v/b25_<chrom>.json.

Semantics (article methods and the supplementary captions):
- Crosses (Supplementary file 1, Table S4): A = BY x RM, B = YJM145 x YPS163, C = CBS2888 x YJM981.
  Its column `effect_directions` gives "-BY,+RM", "-YJM145,+YPS163", "-CBS2888,+YJM981": a positive
  estimate means the second parent's allele raises expression. Each lazy row carries `up`, the
  name of the parent whose allele raises expression for that row.
- Local eQTL (Supplementary file 2): the closest marker to each transcript, with estimate, LOD and
  FDR for every tested transcript; significant = FDR < 0.05 (text: "We mapped 1031 local eQTLs at
  an FDR of 5%"). Distant eQTL (Supplementary file 3): peak marker on another chromosome, FDR < 5%.
  No eQTL interval is tabulated, so lazy rows are points at the marker (a = e = p).
- Hotspots: 50-kb bins (Supplementary file 3, column `hotspot bin`) merged into the regions named
  by the fine-mapping workbooks of Supplementary file 4 (`<bin>_<n eQTL>_finemapping.xlsx` inside
  FileS4_VOR.tar.gz). a/e = the merged bin range; p = the peak marker of the strongest distant
  eQTL in the hotspot (most defensible reading: the paper's hotspot CI, from the top-20 eQTL CIs,
  is not tabulated); qtl.cc = number of distant eQTL in the hotspot.
- Cell-cycle occupancy QTL (Supplementary file 7): 1.5-LOD intervals, estimate, LOD, stage. One
  trait per stage.
- Build: Cellranger with "the S2888C reference genome (SGD, R64-2-1)", i.e. R64. No transform.
"""
import collections
import hashlib
import io
import os
import re
import sys
import tarfile
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "boocock2025", "b25"
SRC = os.path.join(atlas.SOURCES, "boocock2025")
RAW = os.path.join(SRC, "raw")
P = "elife-95566-"
TF = "../extracted/text_facts.tsv"   # relative to raw/, as sources.file is
S4M = "../extracted/supp4_members.tsv"   # member list of FileS4_VOR.tar.gz (line 1 a comment, line 2 the header)
HB = atlas.HEAT_BIN
CROSSES = ["A", "B", "C"]

# provenance code -> (file, sheet)
SHEETS = {
    "b25:S1T1": ("supp1-v1.xlsx", "TableS1. Yeast Strains"),
    "b25:S1T2": ("supp1-v1.xlsx", "Table S2. Plasmids"),
    "b25:S1T3": ("supp1-v1.xlsx", "Table S3. Primers"),
    "b25:S1T4": ("supp1-v1.xlsx", "Table S4. Dataset Summary"),
    "b25:S1T5": ("supp1-v1.xlsx", "Table S5. Variance Decomp"),
    "b25:S1T6": ("supp1-v1.xlsx", "Table S6. Local eQTLs"),
    "b25:S2A": ("supp2-v1.xlsx", "S1 (Cross A)"), "b25:S2B": ("supp2-v1.xlsx", "S2 (Cross B)"), "b25:S2C": ("supp2-v1.xlsx", "S3 (Cross C)"),
    "b25:S3A": ("supp3-v1.xlsx", "S1 (Cross A)"), "b25:S3B": ("supp3-v1.xlsx", "S2 (Cross B)"), "b25:S3C": ("supp3-v1.xlsx", "S3 (Cross C)"),
    "b25:S5A": ("supp5-v1.xlsx", "S1 (Cross A)"), "b25:S5B": ("supp5-v1.xlsx", "S2 (Cross B)"), "b25:S5C": ("supp5-v1.xlsx", "S3 (Cross C)"),
    "b25:S6A": ("supp6-v1.xlsx", "S1 (Cross A)"), "b25:S6B": ("supp6-v1.xlsx", "S2 (Cross B)"), "b25:S6C": ("supp6-v1.xlsx", "S3 (Cross C)"),
    "b25:S7A": ("supp7-v1.xlsx", "S1 (Cross A)"), "b25:S7B": ("supp7-v1.xlsx", "S2 (Cross B)"), "b25:S7C": ("supp7-v1.xlsx", "S3 (Cross C)"),
    "b25:S9": ("supp9-v1.xlsx", "S1"),
}
EXTRA = ["b25:S1T1", "b25:S1T2", "b25:S1T3", "b25:S1T4", "b25:S1T5", "b25:S1T6", "b25:S5A", "b25:S5B", "b25:S5C",
         "b25:S6A", "b25:S6B", "b25:S6C", "b25:S9"]
STAGES = ["G1", "G1/S", "S", "G2/M", "M/G1"]


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


_books = {}


def cell(v):
    """Spreadsheet booleans were exported as the text '=TRUE()' / '=FALSE()'; read them as booleans."""
    if isinstance(v, str) and v.upper() in ("=TRUE()", "=FALSE()"):
        return v.upper() == "=TRUE()"
    return v


def read(code):
    """(column names, header row number, [(spreadsheet row, values)]) for a provenance code.
    Some sheets open with a title/caption row; the header is the first row with 2+ filled cells."""
    f, sheet = SHEETS[code]
    if f not in _books:
        import openpyxl
        wb = openpyxl.load_workbook(os.path.join(RAW, P + f), read_only=True, data_only=True)
        _books[f] = {ws.title: [list(r) for r in ws.iter_rows(values_only=True)] for ws in wb.worksheets}
    rows = _books[f][sheet]
    used = max((j + 1 for r in rows for j, v in enumerate(r) if v is not None), default=0)
    rows = [(r + [None] * used)[:used] for r in rows]
    hi = next(i for i, r in enumerate(rows) if sum(v is not None for v in r) >= 2)
    cols = ["" if v is None else str(v) for v in rows[hi]]
    data = [(i + 1, [cell(v) for v in r]) for i, r in enumerate(rows) if i > hi and any(v is not None for v in r)]
    return cols, hi + 1, data


def mpos(m):
    """'chrI_147785_C_G_839' -> 147785."""
    return int(str(m).split("_")[1])


def mchr(m):
    return str(m).split("_")[0]


def hotspot_files():
    """Fine-mapping workbook names inside Supplementary file 4 (a zip holding FileS4_VOR.tar.gz):
    [(cross, chrom, start, end, n eQTL, member path)] for crosses A, B, C."""
    with zipfile.ZipFile(os.path.join(RAW, P + "supp4-v1.zip")) as z:
        tgz = z.read("FileS4_VOR.tar.gz")
    out = []
    with tarfile.open(fileobj=io.BytesIO(tgz), mode="r:gz") as t:
        for m in t.getnames():
            mm = re.fullmatch(r"fine_mapping/([ABC])/(chr[IVX]+):(\d+)-(\d+)_(\d+)_finemapping\.xlsx", m)
            if mm:
                out.append((mm.group(1), mm.group(2), int(mm.group(3)), int(mm.group(4)), int(mm.group(5)), m))
    return sorted(out, key=lambda r: (r[0], list(atlas.chrom_lengths()).index(r[1]), r[2]))


def main():
    sha = manifest()
    sha[TF] = hashlib.sha256(open(os.path.join(SRC, "extracted", "text_facts.tsv"), "rb").read()).hexdigest()
    sha[S4M] = hashlib.sha256(open(os.path.join(SRC, "extracted", "supp4_members.tsv"), "rb").read()).hexdigest()
    with open(os.path.join(SRC, "extracted", "supp4_members.tsv")) as f:
        s4line = {l.split("\t")[0]: i + 1 for i, l in enumerate(f.read().splitlines())}
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

    hrow = {}
    for code, (f, sheet) in SHEETS.items():
        cols, hr, _ = read(code)
        hrow[code] = hr
        add("sources", code=code, ds=0, file=P + f, sheet=sheet, hrow=hr, sha256=sha[P + f])
    add("sources", code="b25:S4", ds=0, file=S4M, sheet="members of FileS4_VOR.tar.gz in elife-95566-supp4-v1.zip", hrow=2, sha256=sha[S4M])
    add("sources", code="b25:TF", ds=0, file=TF, sheet="text_facts", hrow=1, sha256=sha[TF])

    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        tfr = [l.rstrip("\n").split("\t") for l in f]
    tf_head = tfr[0]
    tf = [(i + 2, dict(zip(tf_head, r))) for i, r in enumerate(tfr[1:])]
    tfrow = {(d["kind"], d["key"]): rn for rn, d in tf}

    # ---- parents and crosses (Supplementary file 1, Table S4: parent1, parent2, effect_directions)
    _, _, t4 = read("b25:S1T4")
    t4c, _, _ = read("b25:S1T4")
    ti = {c: i for i, c in enumerate(t4c)}
    pidx, xidx, up_of = {}, {}, {}
    ORG = {}
    for rn, r in t4:
        if r[ti["type"]] != "eQTL" or r[ti["cross"]] not in CROSSES:
            continue
        cx, p1, p2, dirs = r[ti["cross"]], r[ti["parent1"]], r[ti["parent2"]], r[ti["effect_directions"]]
        assert dirs == f"-{p1},+{p2}", (cx, dirs)
        for p in (p1, p2):
            if p not in pidx:
                pidx[p] = add("parents", id=f"{SC}-{p}", ds=0, name=p, code=p, s1=None, s2=None, p18=None, org=ORG.get(p))
        xidx[cx] = add("crosses", id=f"{SC}-{cx}", ds=0, code=cx, name=f"Cross {cx} ({p1} × {p2})", p1=pidx[p1], p2=pidx[p2], g1=pidx[p1],
                       n=None, dip=None, mm=None, pv=f"b25:S1T4#{rn}")
        up_of[cx] = (p1, p2)

    t_hot = add("traits", id=f"{SC}-sc-trans-eQTL-hotspots", ds=0, raw="trans-eQTL hotspots (single-cell)", name="Expression hotspots (single-cell trans-eQTL)",
                cond=None, bt=None, ck="expression hotspots", cat="expression",
                catb="atlas assignment: loci affecting the mRNA levels of many genes (Boocock et al. 2025, Supplementary files 3 and 4)")
    t_occ = {}
    for st in STAGES:
        t_occ[st] = add("traits", id=f"{SC}-occupancy-{st.replace('/', '-')}", ds=0, raw=f"cell-cycle occupancy {st}", name=f"Cell-cycle occupancy, {st}",
                        cond=st, bt=None, ck=f"cell-cycle occupancy {st.lower()}", cat="cellular phenotype",
                        catb="atlas assignment: fraction of cells assigned to a cell-cycle stage (Supplementary file 7)")

    G = atlas.genes()
    gid, gname = {}, {}
    for i, (s, n) in enumerate(zip(G["id"], G["n"])):
        gid.setdefault(s, i)
        if n:
            gname.setdefault(n, i)
    chroms = list(atlas.chrom_lengths())
    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    cols_of = {}

    # ---- eQTL: lazy rows (significant local + all distant), every tested local row verbatim
    L = {k: [] for k in ["id", "g", "c", "a", "e", "p", "lod", "ef", "loc", "q", "xid", "up", "cci", "hs", "pv"]}
    dens = collections.Counter()
    nogene = set()
    nlocal = {}
    distant = {cx: [] for cx in CROSSES}   # (lod, pos, chrom, bin, eid) for hotspot peaks
    for cx in CROSSES:
        p1, p2 = up_of[cx]
        for code, loc in ((f"b25:S2{cx}", 1), (f"b25:S3{cx}", 0)):
            cols, _, rows = read(code)
            cols_of[code] = cols
            ci = {c: i for i, c in enumerate(cols)}
            mk = "closest marker" if loc else "peak marker"
            nsig = 0
            for rn, r in rows:
                m = r[ci[mk]]
                chrom, pos = mchr(m), mpos(m)
                eid = f"{SC}-e{cx}{loc}-{rn}"
                orig[chrom][code][eid] = r
                fdr = r[ci["FDR (one-pot)"]]
                if fdr is None or fdr >= 0.05:
                    continue
                nsig += 1
                tr = r[ci["transcript"]]
                g = gid.get(tr, -1)
                if g < 0:
                    nogene.add(tr)
                ef = r[ci["estimate (one-pot)"]]
                hs = r[ci["hotspot bin"]] if not loc and r[ci["eQTL in hotspot"]] is True else None
                for k, v in (("id", eid), ("g", g), ("c", chrom), ("a", pos), ("e", pos), ("p", pos), ("lod", atlas.r6(r[ci["LOD (one-pot)"]])),
                             ("ef", atlas.r6(ef)), ("loc", loc), ("q", atlas.r6(fdr)), ("xid", f"{SC}-{cx}"),
                             ("up", p2 if ef > 0 else p1 if ef < 0 else None), ("cci", 1 if r[ci["has cell-cycle interaction"]] is True else 0),
                             ("hs", hs), ("pv", f"{code}#{rn}")):
                    L[k].append(v)
                dens[("local" if loc else "distant", chrom, pos // HB)] += 1
                if not loc:
                    distant[cx].append((r[ci["LOD (one-pot)"]], pos, chrom, hs, eid))
            nlocal[(cx, loc)] = nsig
    for (k, chrom, b), n in sorted(dens.items(), key=lambda kv: (kv[0][0], chroms.index(kv[0][1]), kv[0][2])):
        add("dens", ds=0, k=k, c=chrom, b=b, n=n)
    add("lazy", ds=0, name="sceqtl", file=f"v/{SC}_sceqtl.json", n=len(L["id"]),
        desc="Single-cell eQTL at FDR 5% in three crosses (Supplementary files 2 and 3): local eQTL at the transcript's closest marker and distant eQTL at their peak marker; positions only (no interval tabulated); `up` names the parent whose allele raises expression; `xid` is the cross.")

    # ---- hotspots (Supplementary file 4 workbook names; bins and counts from Supplementary file 3)
    hfiles = hotspot_files()
    hq = {}
    for cx, chrom, a, e, n, member in hfiles:
        inh = [d for d in distant[cx] if d[3] and d[2] == chrom and a <= int(d[3].split(":")[1].split("-")[0]) and int(d[3].split("-")[1]) <= e]
        assert len(inh) == n, (cx, chrom, a, e, n, len(inh))
        best = max(inh)  # strongest distant eQTL: (lod, pos, ...)
        qid = f"{SC}-h-{cx}-{chrom}-{a}"
        q = add("qtl", id=qid, ds=0, t=t_hot, x=xidx[cx], m="w", c=chrom, a=a, e=e, p=best[1], lod=None, q=None, ef=None, ve=None,
                lv=None, af=None, cc=n, pa=None, anc=None, daf=None, up=-1, ups=None, dns=None, xs=None,
                k="trans-eQTL hotspot (single-cell)", g=-1, gf=None, gt=None, pv=f"b25:S4#{s4line[member]}")
        hq[(cx, chrom, a, e)] = q
        orig[chrom]["b25:S4"][qid] = [member, cx, chrom, a, e, n]

    # ---- cell-cycle occupancy QTL (Supplementary file 7)
    oq = []
    for cx in CROSSES:
        code = f"b25:S7{cx}"
        cols, _, rows = read(code)
        cols_of[code] = cols
        ci = {c: i for i, c in enumerate(cols)}
        p1, p2 = up_of[cx]
        for rn, r in rows:
            st = r[ci["cell-cycle"]]
            m = r[ci["marker"]]
            chrom, pos = r[ci["chromosome"]], mpos(m)
            ef = r[ci["estimate (cell-cycle effect)"]]
            qid = f"{SC}-o-{cx}-{rn}"
            q = add("qtl", id=qid, ds=0, t=t_occ[st], x=xidx[cx], m="w", c=chrom, a=int(r[ci["C.I. left"]]), e=int(r[ci["C.I. right"]]), p=pos,
                    lod=atlas.r6(r[ci["LOD (cell-cycle effect)"]]), q=None, ef=atlas.r6(ef), ve=None, lv=m, af=None, cc=None, pa=None, anc=None, daf=None,
                    up=pidx[p2] if ef > 0 else pidx[p1], ups=None, dns=None, xs=None, k="cell-cycle occupancy QTL", g=-1, gf=None, gt=None, pv=f"{code}#{rn}")
            oq.append((q, cx, chrom, r[ci["QTL bin"]]))
            orig[chrom][code][qid] = r

    # ---- candidate genes named in the text
    rank = {"x": 0, "c": 1}
    CAND = [  # (gene, text-fact key, tier, cross, target: "hotspot" or "occupancy")
        ("GPA1", ("candidate", "GPA1"), "x", "B", "hotspot"),
        ("CYR1", ("candidate", "CYR1"), "c", "C", "hotspot"),
        ("ZAP1", ("candidate", "ZAP1"), "c", "A", "hotspot"),
        ("MKT1", ("candidate", "MKT1"), "c", "A", "occupancy"),
        ("GPA1", ("candidate", "GPA1"), "x", "B", "occupancy"),
        ("CYR1", ("candidate", "CYR1"), "c", "C", "occupancy"),
    ]
    unplaced = []
    for gene, key, tier, cx, target in CAND:
        g = gname[gene]
        if target == "hotspot":
            hits = [q for (hx, chrom, a, e), q in hq.items() if hx == cx and chrom == G["c"][g] and a <= G["e"][g] and e >= G["a"][g]]
        else:
            hits = [q for q, ox, chrom, _ in oq if ox == cx and chrom == G["c"][g] and T["qtl"]["a"][q] <= G["e"][g] and T["qtl"]["e"][q] >= G["a"][g]]
        if not hits:
            unplaced.append((gene, cx, target)); continue
        rn = [r for r, d in tf if (d["kind"], d["key"]) == key and (target == "hotspot") == ("occupancy" not in d["value"])][0]
        for q in hits:
            cid = f"{SC}-c-{gene}-{T['qtl']['id'][q].split('-', 1)[1]}"
            add("causal", id=cid, ds=0, t=T["qtl"]["t"][q], g=g, orf=G["id"][g], name=gene, c=G["c"][g], a=G["a"][g], e=G["e"][g],
                pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier=tier, pv=f"b25:TF#{rn}")
            cur = T["qtl"]["gt"][q]
            if cur is None or rank[tier] < rank[cur]:
                T["qtl"]["g"][q] = g; T["qtl"]["gt"][q] = tier
            orig[G["c"][g]]["b25:TF"][cid] = [dict(tf)[rn][k] for k in tf_head]
    cz = T["causal"]
    per_gene = C()
    for i in range(atlas.nlen(cz)):
        per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier=min((cz["tier"][i] for i in ii), key=rank.get), fdr=None, nt=len({cz["t"][i] for i in ii}), n=len(ii))

    # ---- dataset row --------------------------------------------------------------------
    nl = sum(L["loc"])
    add("datasets", id=DS_ID, sc=SC, label="Boocock 2025",
        cite="Boocock J, Alexander N, Alamo Tapia L, Walter-McNeill L, Patel SP, Munugala C, Bloom JS, Kruglyak L. 2025. Single-cell eQTL mapping in yeast reveals a tradeoff between growth and reproduction. eLife 13:RP95566",
        doi="10.7554/eLife.95566", year=2025,
        design="one-pot single-cell eQTL: pools of MATa segregants from three crosses (A = BY × RM, 27,744 cells; B = YJM145 × YPS163, 44,784 cells; C = CBS2888 × YJM981, 6,595 cells), scRNA-seq with genotypes inferred per cell by an HMM; cell-cycle stage assigned per cell",
        build="S288C R64: reads were processed with 'the S2888C reference genome (SGD, R64-2-1)'; all marker positions fall inside R64 chromosomes. No transform.",
        method="Count-based (negative binomial) model per transcript with cell-cycle stage as covariate; LOD from the likelihood ratio statistic / (2 ln 10); FDR 5% from permutations. Local eQTL: closest marker to each transcript. Distant eQTL: peak marker per transcript and chromosome, genes on a different chromosome. Hotspots: 50-kb bins enriched for distant eQTL. Cell-cycle occupancy QTL: allele over-representation per stage, 1.5-LOD intervals.",
        fdr=0.05,
        interval="eQTL (lazy file): not tabulated; drawn at the marker. Hotspots: the merged 50-kb bins named in Supplementary file 4. Cell-cycle occupancy QTL: 1.5-LOD drop (Supplementary file 7).",
        effect="Per cross, positive = the second-named parent's allele raises the trait (Supplementary file 1, Table S4, effect_directions: A '-BY,+RM', B '-YJM145,+YPS163', C '-CBS2888,+YJM981'); each lazy eQTL row names that parent in `up`. Hotspots: no effect size tabulated.",
        nq=atlas.nlen(T["qtl"]),
        notes=f"3 crosses; {len(L['id']):,} eQTL at FDR 5% ({nl:,} local, {len(L['id']) - nl:,} distant) loaded on demand, {len(hq)} trans-eQTL hotspots and {len(oq)} cell-cycle occupancy QTL rows. "
              "Segregant counts are not given (each cell is a segregant; cell counts are in the design). Candidate genes are those the text names: GPA1 (tier x, allele replacement in Supplementary file 9), CYR1, ZAP1, MKT1 (tier c); not claimed causal beyond the paper.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy files ---------------------------------------------------------------------------
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_sceqtl.json"), L)
    codes = [f"b25:S2{c}" for c in CROSSES] + [f"b25:S3{c}" for c in CROSSES] + [f"b25:S7{c}" for c in CROSSES] + ["b25:S4", "b25:TF"]
    cols_of["b25:S4"] = ["member", "cross", "chromosome", "start", "end", "n distant eQTL"]
    cols_of["b25:TF"] = tf_head
    for chrom in chroms:
        out = {"cols": {c: cols_of[c] for c in codes}, "rows": {c: dict(orig[chrom].get(c, {})) for c in codes}}
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"), out)
    extra = {}
    for code in EXTRA:
        cols, _, rows = read(code)
        extra[code] = {"cols": cols, "rows": [r for _, r in rows]}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    size = os.path.getsize(os.path.join(atlas.DOCS, "v", f"{SC}_sceqtl.json"))
    print(f"{DS_ID}: eQTL local/distant per cross {dict(sorted(nlocal.items()))}; {len(hq)} hotspots; {len(oq)} occupancy QTL rows; "
          f"{atlas.nlen(T['causal'])} candidate rows (unplaced: {unplaced}); genes not in R64-5-1 table: {len(nogene)}; "
          f"b25_sceqtl.json {size:,} bytes; data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
