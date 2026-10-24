"""Ingest Bloom et al. 2019 (eLife 8:e49212), the 16-parent round-robin segregant panel.

    python3 -I ingest/bloom2019.py          # from the repository root

Reads only from sources/ (the unzipped eLife bundle in sources/bloom2019/raw/, the
code-derived tables in sources/bloom2019/code_tables/, and the SGD tables in
sources/sgd/). Writes build/fragments/bloom2019.json, docs/v/b19_*.json, and, through
ingest/atlas.py, docs/data.json. Deterministic: same inputs, same bytes out.

Genome build: the source coordinates are S288C R64 / UCSC sacCer3 (the analysis code
loads BSgenome.Scerevisiae.UCSC.sacCer3 and its chromosome lengths equal the R64 lengths
in sources/sgd/chrom_lengths_R64.tsv). No coordinate transform is applied.

Column semantics come from https://github.com/joshsbloom/yeast-16-parents (commit
c913c9a); see SCHEMA.md and OPEN_QUESTIONS.md.
"""
import collections
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "bloom2019", "b19"
SRC = os.path.join(atlas.SOURCES, "bloom2019")
RAW = os.path.join(SRC, "raw")
PFX = "elife_poa_e49212_"

FILES = {
    "F1D1": PFX + "Figure_1_source_data_1.xls",
    "F2D1": PFX + "Figure_2_source_data_1.xls",
    "F3D1": PFX + "Figure_3_source_data_1.xls",
    "F3D2": PFX + "Figure_3_source_data_2.xls",
    "F4D1": PFX + "Figure_4_source_data_1.xls",
    "SF1": PFX + "Supplementary_file_1.xlsx",
}
# provenance code -> (file code, sheet, 1-based header row)
SHEETS = {
    "F1D1:x": ("F1D1", "Crosses and Strains", 1),
    "F1D1:ph": ("F1D1", "Phenotypes", 1),
    "F2D1:v": ("F2D1", "Additive QTL Model", 1),
    "F2D1:vc3": ("F2D1", "Multiple VC model (3 epistasis components)", 1),
    "F2D1:vc1": ("F2D1", "Multiple VC Model (1 epistasis component)", 1),
    "F3D1:pc": ("F3D1", "joint model (effect per cross)", 1),
    "F3D1:j": ("F3D1", "joint model (average effect)", 1),
    "F3D1:w": ("F3D1", "within-cross model", 1),
    "F3D2:af2": ("F3D2", "2 components", 1),
    "F3D2:af7": ("F3D2", "7 components", 1),
    "F3D2:af7p": ("F3D2", "7 components private-only", 1),
    "F4D1:cg": ("F4D1", "Causal_Genes", 2),
    "F4D1:bp": ("F4D1", "GO_Biological_Process", 1),
    "F4D1:mf": ("F4D1", "GO_Molecular_Function", 1),
    "F4D1:cc": ("F4D1", "GO_Cellular_Compartment", 1),
    "SF1:kr": ("SF1", "Sheet1", 1),
}

# Trait categories. The Phenotypes sheet of Figure 1-source data 1 groups the media into a
# compound table, a condition table, a "Carbon Sources" section and a "pH" section; the
# category below names that section. Where the sheet does not settle it, the basis says so.
PH = "Figure 1-source data 1, Phenotypes sheet: "
CATEGORY = {
    "6-azauracil": ("drug/chemical", PH + "compound table"),
    "Cadmium_Chloride": ("drug/chemical", PH + "compound table"),
    "Caffeine": ("drug/chemical", PH + "compound table"),
    "Cobalt_Chloride": ("drug/chemical", PH + "compound table"),
    "Congo_Red": ("drug/chemical", PH + "compound table"),
    "Diamide": ("drug/chemical", PH + "compound table"),
    "EGTA": ("drug/chemical", PH + "compound table"),
    "Fluconazole": ("drug/chemical", PH + "compound table"),
    "Formamide": ("drug/chemical", PH + "compound table"),
    "Lithium_Chloride": ("drug/chemical", PH + "compound table"),
    "Methotrexate": ("drug/chemical", PH + "compound table"),
    "Neomycin": ("drug/chemical", PH + "compound table"),
    "Paraquat": ("drug/chemical", PH + "compound table"),
    "SDS": ("drug/chemical", PH + "compound table"),
    "Tunicamycin": ("drug/chemical", PH + "compound table"),
    "Zeocin": ("drug/chemical", PH + "compound table"),
    "Copper_Sulfate": ("drug/chemical", PH + "condition table, added to YP 2% glucose"),
    "Magnesium_Chloride": ("drug/chemical", PH + "condition table, added to YP 2% glucose"),
    "Manganese_Sulfate": ("drug/chemical", PH + "condition table, added to YP 2% glucose"),
    "EtOH_Glucose": ("drug/chemical", PH + "condition table ('Ethanol with Glucose', YP 2% glucose)"),
    "EtOH": ("carbon source", PH + "condition table ('Ethanol NO glucose')"),
    "Glycerol": ("carbon source", PH + "condition table ('Glycerol 3%', YP, no carbon source listed)"),
    "Fructose": ("carbon source", PH + "Carbon Sources section"),
    "Galactose": ("carbon source", PH + "Carbon Sources section"),
    "Lactate": ("carbon source", PH + "Carbon Sources section"),
    "Maltose": ("carbon source", PH + "Carbon Sources section"),
    "Mannose": ("carbon source", PH + "Carbon Sources section"),
    "Raffinose": ("carbon source", PH + "Carbon Sources section"),
    "Sucrose": ("carbon source", PH + "Carbon Sources section"),
    "Trehalose": ("carbon source", PH + "Carbon Sources section"),
    "Xylose": ("carbon source", PH + "Carbon Sources section"),
    "Sorbitol": ("osmotic stress", "atlas assignment (Sorbitol is not in the Phenotypes sheet)"),
}
DISPLAY = {"6-azauracil": "6-Azauracil", "EtOH": "Ethanol", "EtOH_Glucose": "Ethanol + glucose",
           "Congo_Red": "Congo red", "SDS": "SDS", "EGTA": "EGTA", "YPD": "YPD", "YNB": "YNB"}


def num(v):
    """xlrd stores every number as float; give integral values back as int."""
    if isinstance(v, float) and v.is_integer() and abs(v) < 1e15:
        return int(v)
    return v


def verbatim(row):
    return [None if v == "" else num(v) for v in row]


def fnum(v):
    if v in ("", None, "NA"):
        return None
    return float(v)


def sstr(v):
    v = num(v)
    return "" if v is None else str(v)


def flag(v):
    if isinstance(v, bool):
        return v
    return str(v).upper() == "TRUE"


def marker_pos(m):
    return int(str(m).split("_")[1])


def marker_chr(m):
    return str(m).split("_")[0]


_books = {}


def read(code):
    """Return (column names, [(spreadsheet row, values)]) for a provenance code."""
    fcode, sheet, hrow = SHEETS[code]
    path = os.path.join(RAW, FILES[fcode])
    if path not in _books:
        if path.endswith(".xls"):
            import xlrd
            b = xlrd.open_workbook(path)
            _books[path] = {s.name: [s.row_values(i) for i in range(s.nrows)] for s in b.sheets()}
        else:
            import openpyxl
            b = openpyxl.load_workbook(path, read_only=True, data_only=True)
            _books[path] = {ws.title: [["" if v is None else v for v in r] for r in ws.iter_rows(values_only=True)]
                            for ws in b.worksheets}
    rows = _books[path][sheet]
    width = max(len(r) for r in rows)
    rows = [list(r) + [""] * (width - len(r)) for r in rows]
    cols = [sstr(v) for v in rows[hrow - 1]]
    return cols, [(i + 1, r) for i, r in enumerate(rows) if i >= hrow]


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def tsv(path):
    with open(path) as f:
        rows = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def slug(s):
    return re.sub(r"-+", "-", re.sub(r"[^A-Za-z0-9]+", "-", s)).strip("-")


def trait_name(raw):
    comp, cond, _ = (raw.split(";") + ["", ""])[:3]
    base = DISPLAY.get(comp) or comp.replace("_", " ")
    base = base[0].upper() + base[1:].lower() if comp not in DISPLAY else base
    c = cond
    if comp == "YPD" and re.fullmatch(r"\d+", c):
        c = c + " °C"
    elif c.lower().startswith("ph"):
        c = "pH " + c[2:]
    else:
        c = re.sub(r"^([\d.]+)(?=[A-Za-z%])", r"\1 ", c).replace(" u", " µ").replace(" %", "%")
        c = c.replace("ug/mL", "µg/mL")
    return (base + " " + c).strip()


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

    # ---- sources -------------------------------------------------------------------
    for code, (fcode, sheet, hrow) in SHEETS.items():
        add("sources", code=code, ds=0, file=FILES[fcode], sheet=sheet, hrow=hrow, sha256=sha[FILES[fcode]])

    # ---- parents and crosses (Figure 1-source data 1 + the code's parent labels) ---
    xcols, xrows = read("F1D1:x")
    code_rows = {r["cross"]: (r["parent1_code"], r["parent2_code"])
                 for r in tsv(os.path.join(SRC, "code_tables", "crosses_to_parents.tsv"))}
    pidx, s1, s2, p18 = {}, {}, {}, {}
    for rn, r in xrows:
        cross = sstr(r[4])
        c1, c2 = code_rows[cross]
        s1[c1] = r[2]; s2[c2] = r[3]
        p18[c1] = sstr(r[1]) or None
        if c1 not in pidx:
            pidx[c1] = len(pidx)
    for c in pidx:
        name = re.sub(r"(xa|a|x)$", "", c)
        add("parents", id=f"{SC}-{name}", ds=0, name=name, code=c, s1=s1.get(c), s2=s2.get(c), p18=p18.get(c), org=None)

    # Which parent is genotype 1 in each cross's within-cross coding. Not recorded in the
    # code or bundle; recovered from the Causal_Genes sheet, where the code labelled each
    # shared-parent pair 'increase'/'decrease' from sign(beta) and that parent's genotype
    # column (QTL_causality.R lines 377-383). See OPEN_QUESTIONS Q2.
    ccols, crows = read("F4D1:cg")
    votes = collections.defaultdict(collections.Counter)
    seen = set()
    for rn, r in crows:
        key = (r[0], r[23], r[48])
        if key in seen:
            continue
        seen.add(key)
        d = -1 if r[10] == "increase" else 1
        for cross, beta in ((sstr(r[31]), r[26]), (sstr(r[56]), r[51])):
            sgn = 1 if beta > 0 else -1
            other = [p for p in code_rows[cross] if p != r[9]][0]
            first = r[9] if sgn == d else other
            votes[cross][first] += 1
    g1 = {}
    for cross, (c1, c2) in code_rows.items():
        g1[cross] = votes[cross].most_common(1)[0][0]
        assert g1[cross] == min(c1, c2), (cross, votes[cross])  # every cross: the alphabetically first label
    xidx = {}
    for rn, r in xrows:
        cross = sstr(r[4])
        c1, c2 = code_rows[cross]
        xidx[cross] = add("crosses", id=f"{SC}-{cross}", ds=0, code=cross, p1=pidx[c1], p2=pidx[c2], g1=pidx[g1[cross]],
                          n=num(r[6]), dip=r[0], mm=r[5], pv=f"F1D1:x#{rn}")

    # ---- QTL tables ----------------------------------------------------------------
    jcols, jrows = read("F3D1:j")
    wcols, wrows = read("F3D1:w")
    pccols, pcrows = read("F3D1:pc")
    traits_raw = sorted({r[0] for _, r in jrows} | {r[0] for _, r in wrows} | {r[0] for _, r in pcrows}, key=str.lower)
    tidx = {}
    for raw in traits_raw:
        comp, cond, bt = (raw.split(";") + ["", ""])[:3]
        if comp in CATEGORY:
            cat, catb = CATEGORY[comp]
        elif comp == "YPD" and cond:
            cat, catb = "temperature", PH + "condition table ('YPD " + cond + "')"
        elif comp == "YNB" and cond.lower().startswith("ph"):
            cat, catb = "pH", PH + "pH section"
        elif comp in ("YPD", "YNB"):
            cat, catb = "reference medium", PH + "condition table ('" + comp + "')"
        else:
            raise SystemExit("no category for trait " + raw)
        tidx[raw] = add("traits", id=f"{SC}-{slug(raw)}", ds=0, raw=raw, name=trait_name(raw), cond=cond, bt=bt, cat=cat, catb=catb)

    prev = {r["name"] for r in tsv(os.path.join(SRC, "code_tables", "previously_identified.tsv"))}
    G = atlas.genes()
    gidx = {}
    for i, g in enumerate(G["id"]):
        gidx.setdefault(g, i)

    def tier(name, fdr):
        if name in prev:
            return "v"
        return "5" if fdr < 0.05 else "20"

    orig = collections.defaultdict(lambda: collections.defaultdict(dict))  # chrom -> code -> id -> row
    qkey = {}
    # joint model (average effect): one row per joint QTL
    for rn, r in jrows:
        trait, peak = r[0], r[1]
        ef = fnum(r[6])
        resolved = flag(r[7])
        g = gidx[r[8]] if resolved else -1
        gf = fnum(r[10]) if resolved else None
        cc = num(r[3])
        qid = f"{SC}-j-{rn}"
        q = add("qtl", id=qid, ds=0, t=tidx[trait], x=-1, m="j", c=marker_chr(peak), a=marker_pos(peak), e=marker_pos(peak),
                p=marker_pos(peak), lod=None, q=None, ef=atlas.r6(ef), ve=None, lv=peak, af=atlas.r6(fnum(r[2])),
                cc=cc, pa=None if cc is None else cc // 2, anc=None if r[4] == "" else num(r[4]), daf=atlas.r6(fnum(r[5])),
                up=-1, g=g, gf=atlas.r6(gf), gt=tier(r[9], gf) if resolved else None, pv=f"F3D1:j#{rn}")
        qkey[("j", trait, peak)] = q
        orig[marker_chr(peak)]["F3D1:j"][qid] = verbatim(r)
    # within-cross model
    other = {}
    for cross, (c1, c2) in code_rows.items():
        other[cross] = pidx[c2] if g1[cross] == c1 else pidx[c1]
    for rn, r in wrows:
        trait, chrom, cross = r[0], r[1], sstr(r[3])
        beta = fnum(r[13])
        qid = f"{SC}-w-{rn}"
        q = add("qtl", id=qid, ds=0, t=tidx[trait], x=xidx[cross], m="w", c=chrom, a=marker_pos(r[10]), e=marker_pos(r[11]),
                p=num(r[2]), lod=atlas.r6(fnum(r[8])), q=atlas.r6(fnum(r[7])), ef=atlas.r6(beta), ve=atlas.r6(fnum(r[12])),
                lv=str(r[9]).rsplit("_", 1)[0], af=None, cc=None, pa=None, anc=None, daf=None,
                up=other[cross] if beta > 0 else pidx[g1[cross]], g=-1, gf=None, gt=None, pv=f"F3D1:w#{rn}")
        qkey[("w", trait, cross, r[9])] = q
        orig[chrom]["F3D1:w"][qid] = verbatim(r)
    # joint model effects per cross
    for rn, r in pcrows:
        trait, cross, peak = r[0], sstr(r[1]), r[2]
        q = qkey[("j", trait, peak)]
        add("jeff", q=q, x=xidx[cross], b=atlas.r6(fnum(r[3])), ve=atlas.r6(fnum(r[4])), cn=num(r[6]), pv=f"F3D1:pc#{rn}")
        orig[marker_chr(peak)]["F3D1:pc"].setdefault(T["qtl"]["id"][q], []).append(verbatim(r))

    # ---- candidate genes -------------------------------------------------------------
    best = {}  # qtl index -> (fdr, -pcs, causal index)
    unmatched = []
    for rn, r in crows:
        trait, orf, name = r[0], r[1], sstr(r[2]) or None
        fdr = float(r[11])
        q1 = qkey.get(("w", trait, sstr(r[31]), r[20]), -1)
        q2 = qkey.get(("w", trait, sstr(r[56]), r[45]), -1)
        cid = f"{SC}-c-{rn}"
        k = add("causal", id=cid, ds=0, t=tidx[trait], g=gidx[orf], orf=orf, name=name, c=r[3], a=num(r[4]), e=num(r[5]),
                pcs=atlas.r6(r[8]), sp=pidx[r[9]], eff=r[10], fdr=atlas.r6(fdr), jppc=atlas.r6(r[12]), jv=r[13],
                q1=q1, q2=q2, tier=tier(name, fdr), pv=f"F4D1:cg#{rn}")
        for q, pm, cross, qv in ((q1, r[20], sstr(r[31]), r[18]), (q2, r[45], sstr(r[56]), r[43])):
            if q < 0:
                unmatched.append((rn, trait, cross, pm, qv))
                continue
            cand = (fdr, -float(r[8]), k)
            if q not in best or cand < best[q]:
                best[q] = cand
        orig[r[3]]["F4D1:cg"][cid] = verbatim(r)
    for q, (fdr, _, k) in best.items():
        T["qtl"]["g"][q] = T["causal"]["g"][k]
        T["qtl"]["gf"][q] = atlas.r6(fdr)
        T["qtl"]["gt"][q] = T["causal"]["tier"][k]
    rank = {"v": 0, "5": 1, "20": 2}
    per_gene = C()
    cz = T["causal"]
    for i in range(atlas.nlen(cz)):
        per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier=min((cz["tier"][i] for i in ii), key=rank.get), fdr=min(cz["fdr"][i] for i in ii),
            nt=len({cz["t"][i] for i in ii}), n=len(ii))

    # ---- variance explained ------------------------------------------------------------
    vcols, vrows = read("F2D1:v")
    for rn, r in vrows:
        add("vexp", t=tidx[r[0]], x=xidx[sstr(r[1])], wve=atlas.r6(fnum(r[2])), jve=atlas.r6(fnum(r[3])),
            h2=atlas.r6(fnum(r[4])), h2se=atlas.r6(fnum(r[5])), re=atlas.r6(fnum(r[6])), rese=atlas.r6(fnum(r[7])), pv=f"F2D1:v#{rn}")

    # ---- dataset row ---------------------------------------------------------------
    nseg = sum(T["crosses"]["n"])
    add("datasets", id=DS_ID, sc=SC, label="Bloom 2019",
        cite="Bloom JS et al. 2019. eLife 8:e49212",
        doi="10.7554/eLife.49212", year=2019,
        design="segregant panel; 16 biparental crosses in a round robin (each parent in two crosses); within-cross and joint (all crosses) models",
        build="S288C R64 (UCSC sacCer3). The analysis code loads BSgenome.Scerevisiae.UCSC.sacCer3 and its chromosome lengths match R64; variant names are chr_pos_ref_alt on that reference. No coordinate transform.",
        method="Within-cross: forward stepwise regression on scaled genotypes with G'Sell FDR control (10,000 permutations per step), then 1.5-LOD-drop intervals around each QTL refitted with the others in the model. "
               "Joint: per chromosome, forward stepwise search over variants biallelic in the 1,011-isolate panel, phenotypes residualized for QTL on other chromosomes and a polygenic background, structured permutations within cross (1,000), same FDR procedure.",
        fdr=0.05,
        interval="Within-cross QTL: 1.5-LOD-drop interval (left and right markers). Joint QTL: the lead variant only; the source gives no interval.",
        effect="Standardized: SD of the cross-scaled colony-radius phenotype per SD of genotype. Within-cross: positive = genotype-2 parent's allele (crosses.g1 is genotype 1) increases colony radius. Joint: positive = alternate (non-S288C) allele increases colony radius; mean over the crosses where the variant segregates.",
        nq=atlas.nlen(T["qtl"]),
        notes=f"{len(pidx)} parents, {len(xidx)} crosses, {nseg:,} segregants, {len(tidx)} traits. "
              "Candidate genes: assigned to pairs of within-cross QTL that share a parent, by fine-mapping posterior (PPC) within 25 kb of each peak; tiers are the paper's local FDR (5%, 20%) and the analysis code's list of previously identified genes.")

    frag = {k: v for k, v in T.items()}
    atlas.save_fragment(DS_ID, frag)

    # ---- lazy files --------------------------------------------------------------------
    for chrom in atlas.chrom_lengths():
        out = {"cols": {}, "rows": {}}
        for code in ("F3D1:j", "F3D1:w", "F3D1:pc", "F4D1:cg"):
            out["cols"][code] = read(code)[0]
            out["rows"][code] = dict(orig[chrom].get(code, {}))
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"), out)
    extra = {}
    for code in ("F1D1:ph", "F2D1:vc3", "F2D1:vc1", "F3D2:af2", "F3D2:af7", "F3D2:af7p", "F4D1:bp", "F4D1:mf", "F4D1:cc", "SF1:kr"):
        cols, rows = read(code)
        extra[code] = {"cols": cols, "rows": [verbatim(r) for _, r in rows]}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()

    # ---- report ------------------------------------------------------------------------
    q = T["qtl"]
    print(f"{DS_ID}: {len(pidx)} parents, {len(xidx)} crosses, {nseg} segregants, {len(tidx)} traits", file=sys.stderr)
    print(f"  qtl: {q['m'].count('j')} joint, {q['m'].count('w')} within-cross; jeff {atlas.nlen(T['jeff'])}; "
          f"causal {atlas.nlen(T['causal'])} rows, {atlas.nlen(T['cand'])} genes; vexp {atlas.nlen(T['vexp'])}", file=sys.stderr)
    print(f"  genotype-1 parent per cross: " + ", ".join(f"{c}={g1[c]}" for c in code_rows), file=sys.stderr)
    print(f"  causal-sheet QTL not in the within-cross table: {len(unmatched)}", file=sys.stderr)
    for u in unmatched:
        print("    row %d %s cross %s %s q=%s" % u, file=sys.stderr)
    print(f"  docs/data.json: {os.path.getsize(os.path.join(atlas.DOCS, 'data.json')):,} bytes, {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
