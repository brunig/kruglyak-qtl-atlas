"""Ingest Treusch et al. 2015 (PLOS Genetics 11:e1004913): X-QTL in a 12-strain round-robin cross.

    python3 -I ingest/treusch2015.py          # from the repository root

Reads only from sources/ (the supporting-information files in sources/treusch2015/raw/ and
the SGD tables in sources/sgd/). Writes build/fragments/treusch2015.json, docs/v/t15_*.json
and, through ingest/atlas.py, docs/data.json. Deterministic.

Genome build: S288C R64 (sacCer3). Methods: "We created a set of genetic markers for each
cross by aligning parent strain sequences to the sacCer3 reference genome". Chromosome names
in the tables are chr01..chr16 and are renamed chrI..chrXVI; positions are not transformed.

Semantics, from the article (Methods, "Allele frequency measurements and QTL mapping", and
the S4/S8 captions): LOD scores come from MULTIPOOL comparing allele frequencies of a selected
pool with its YPD control pool (S3, S4), or of the MATa pool with the MATalpha pool (S8);
QTL are called at LOD > 5 with 2-LOD drop intervals ("Lower/upper bound"). A round-robin QTL
had to pass LOD 5 in both the MATa and MATalpha selections of a cross; S4 groups those per-cross
QTL across crosses into regions and gives one combined LOD, max-LOD position and interval per
region. X-QTL give no effect sizes; S7 lists strains carrying beneficial or deleterious alleles.
See OPEN_QUESTIONS.md, Treusch 2015 section.
"""
import collections
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "treusch2015", "t15"
SRC = os.path.join(atlas.SOURCES, DS_ID)
RAW = os.path.join(SRC, "raw")
P = "journal.pgen.1004913."
# provenance code -> (file, sheet, 1-based header row)
SHEETS = {
    "t15:S1": (P + "s005.xlsx", "Sheet1", 3),
    "t15:S2": (P + "s006.xlsx", "Sheet1", 3),
    "t15:S3": (P + "s007.xlsx", "Sheet1", 3),
    "t15:S4": (P + "s008.xlsx", "Sheet1", 3),
    "t15:S5": (P + "s009.xlsx", "GO Analyses", 5),
    "t15:S6": (P + "s010.xlsx", "ENA CN", 3),
    "t15:S7": (P + "s011.xlsx", "RR QTL Alleles", 3),
    "t15:S8": (P + "s012.xlsx", "MAT RR QTL", 3),
    "t15:S9": (P + "s013.xlsx", "Plasmids", 3),
    "t15:S10": (P + "s014.xlsx", "Primers", 3),
}
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]

# Display names shared with Bloom 2019 where the strain background is the same, so strains link
# across datasets by name. Keys are the source spellings upper-cased.
DISPLAY = {"273614X": "273614", "CLIB219": "CLIB219", "NC-02": "NC-02"}

# Trait per source condition string. Category: caffeine is a drug; high NaCl is labelled osmotic
# stress (atlas assignment; the paper: "High salt conditions lead to activation of the HOG
# signaling cascade"). S8 traits are MATa-vs-MATalpha comparisons, a different quantity.
TRAITS = [
    # raw, display name, cond, ck, category, category basis
    ("Caffeine", "Caffeine 15 + 20 mM (combined)", "15 mM and 20 mM combined", "caffeine 15-20mM", "drug/chemical", "atlas assignment"),
    ("NaCl", "NaCl 0.5 + 1 M (combined)", "0.5 M and 1 M combined", "nacl 0.5-1M", "osmotic stress", "atlas assignment (HOG pathway; article Results)"),
    ("15mM Caffeine", "Caffeine 15 mM", "15mM", "caffeine 15mM", "drug/chemical", "atlas assignment"),
    ("1M NaCl", "NaCl 1 M", "1M", "nacl 1M", "osmotic stress", "atlas assignment (HOG pathway; article Results)"),
    ("0.5M & 1M NaCl", "NaCl 0.5 + 1 M", "0.5M & 1M", "nacl 0.5-1M", "osmotic stress", "atlas assignment (HOG pathway; article Results)"),
    ("MAT:YPD", "MATa vs MATα, YPD", "YPD", None, "reference medium", "atlas assignment (S8: MATa vs MATalpha pools)"),
    ("MAT:0.5M NaCl", "MATa vs MATα, NaCl 0.5 M", "0.5M", None, "osmotic stress", "atlas assignment (S8: MATa vs MATalpha pools)"),
    ("MAT:1M NaCl", "MATa vs MATα, NaCl 1 M", "1M", None, "osmotic stress", "atlas assignment (S8: MATa vs MATalpha pools)"),
    ("MAT:15mM Caffeine", "MATa vs MATα, caffeine 15 mM", "15mM", None, "drug/chemical", "atlas assignment (S8: MATa vs MATalpha pools)"),
    ("MAT:20mM Caffeine", "MATa vs MATα, caffeine 20 mM", "20mM", None, "drug/chemical", "atlas assignment (S8: MATa vs MATalpha pools)"),
]


def num(v):
    if isinstance(v, float) and v.is_integer() and abs(v) < 1e15:
        return int(v)
    return v


def verbatim(row):
    return [None if v in ("", None) else num(v) for v in row]


def chrom(s):
    m = re.fullmatch(r"chr0?(\d+)", str(s).strip())
    return "chr" + ROMAN[int(m.group(1)) - 1]


_books = {}


def read(code):
    """(column names, [(spreadsheet row, values, bold)]) for every row after the header."""
    import openpyxl
    f, sheet, hrow = SHEETS[code]
    path = os.path.join(RAW, f)
    if path not in _books:
        _books[path] = openpyxl.load_workbook(path)
    ws = _books[path][sheet]
    width = ws.max_column
    rows = []
    for r in ws.iter_rows(min_row=1, max_row=ws.max_row, max_col=width):
        rows.append(([c.value for c in r], bool(r[0].font and r[0].font.bold)))
    cols = ["" if v is None else str(v) for v in rows[hrow - 1][0]]
    while cols and cols[-1] == "":
        cols.pop()
    out = [(i + 1, vals[:max(len(cols), 1)] if code not in ("t15:S7",) else vals, b)
           for i, (vals, b) in enumerate(rows) if i >= hrow]
    return cols, out


def full_sheet(code):
    """Every cell of the sheet, row by row, for the verbatim extra file."""
    import openpyxl
    f, sheet, _ = SHEETS[code]
    path = os.path.join(RAW, f)
    if path not in _books:
        _books[path] = openpyxl.load_workbook(path)
    ws = _books[path][sheet]
    return [verbatim([c.value for c in r]) for r in ws.iter_rows()]


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def nonempty(v):
    return v not in (None, "")


def main():
    sha = manifest()
    T = {k: collections.OrderedDict() for k in atlas.TABLES}

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

    for code, (f, sheet, hrow) in SHEETS.items():
        add("sources", code=code, ds=0, file=f, sheet=sheet, hrow=hrow, sha256=sha[f])

    G = atlas.genes()
    gname = {}
    for i, (sid, n) in enumerate(zip(G["id"], G["n"])):
        gname.setdefault(sid.upper(), i)
        if n:
            gname.setdefault(n.upper(), i)

    # ---- parents: the 12 round-robin strains (S2) plus BY and RM (S3, BY/RM cross) ----------
    s1cols, s1rows = read("t15:S1")
    s1_label = {}   # strain (upper, without _suffix) -> S1 label and whether it is bold (RR parent)
    for rn, vals, bold in s1rows:
        if nonempty(vals[0]):
            lab = str(vals[0])
            s1_label.setdefault(re.sub(r"_.*$", "", lab).upper(), (lab, bold))
    s7cols, s7rows = read("t15:S7")
    origin = {}
    for rn, vals, bold in s7rows:
        if len(vals) > 8 and nonempty(vals[6]) and vals[6] != "Strain":
            origin[str(vals[6]).upper()] = f"{vals[7]}: {vals[8]}"
    s2cols, s2rows = read("t15:S2")
    pidx, pcode = {}, []

    def parent(src_name):
        key = str(src_name).strip().upper()
        if key not in pidx:
            name = DISPLAY.get(key, str(src_name).strip())
            lab = s1_label.get(key, (None, False))[0]
            pidx[key] = add("parents", id=f"{SC}-{name}", ds=0, name=name, code=str(src_name).strip(),
                            s1=lab, s2=None, p18=None, org=origin.get(key))
        return pidx[key]

    xidx = {}
    for rn, vals, bold in s2rows:
        if not nonempty(vals[0]):
            continue
        n = num(vals[0])
        px, pa = parent(vals[1]), parent(vals[2])
        xidx[str(n)] = add("crosses", id=f"{SC}-{n}", ds=0, code=str(n), name=f"Cross {n}", p1=px, p2=pa, g1=None,
                           n=None, dip=None, mm=None, pv=f"t15:S2#{rn}")
    # BY x RM: the replicate experiments (S3). BY and RM are named in the Methods; the cross is
    # evidenced by S3's first data row.
    s3cols, s3rows = read("t15:S3")
    first_s3 = next(rn for rn, vals, b in s3rows if nonempty(vals[0]))
    xby = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=parent("BY"), p2=parent("RM"), g1=None,
              n=None, dip="yLK1993", mm=None, pv=f"t15:S3#{first_s3}")
    xidx["BYxRM"] = xby

    # ---- traits ------------------------------------------------------------------------
    tidx = {}
    for raw, name, cond, ck, cat, catb in TRAITS:
        slug = re.sub(r"-+", "-", re.sub(r"[^A-Za-z0-9]+", "-", raw)).strip("-")
        tidx[raw] = add("traits", id=f"{SC}-{slug}", ds=0, raw=raw, name=name, cond=cond, bt=None, ck=ck, cat=cat, catb=catb)

    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    T_c = []  # candidate rows to add after QTL

    def candidates(text):
        """Split a 'Candidate QTG' cell into (token, gene index or -1, confirmed?)."""
        out = []
        if not nonempty(text):
            return out
        for tok in re.split(r",|&", str(text)):
            tok = tok.strip()
            if not tok:
                continue
            conf = "(confirmed)" in tok
            name = tok.replace("(confirmed)", "").strip()
            out.append((name, gname.get(name.upper(), -1), conf))
        return out

    # ---- S7: beneficial / deleterious strains per grouped QTL --------------------------
    s7 = {}
    for rn, vals, bold in s7rows:
        if nonempty(vals[0]) and str(vals[0]) != "QTL":
            ups = [parent(s) for s in str(vals[1]).split(",") if s.strip()] if nonempty(vals[1]) else []
            dns = [parent(s) for s in str(vals[3]).split(",") if s.strip()] if nonempty(vals[3]) else []
            s7[str(vals[0])] = (ups, dns, rn, verbatim(vals[:5]))

    # ---- S4: grouped round-robin QTL ---------------------------------------------------
    s4cols, s4rows = read("t15:S4")
    cond = None
    for rn, vals, bold in s4rows:
        if not nonempty(vals[1]):
            continue
        if nonempty(vals[0]):
            cond = str(vals[0])
        c = chrom(vals[2])
        xs = [xidx[str(int(float(x)))] for x in str(vals[7]).split(",") if x.strip()]
        ups, dns, s7rn, s7row = s7.get(str(vals[1]), (None, None, None, None))
        cands = candidates(vals[13])
        g = next((gi for _, gi, _ in cands if gi >= 0), -1)
        tier = None
        if g >= 0:
            tier = "x" if any(conf for _, gi, conf in cands if gi == g) else "c"
        qid = f"{SC}-x-S4-{rn}"
        q = add("qtl", id=qid, ds=0, t=tidx[cond], x=-1, xs=xs, m="x", k="grouped round-robin QTL", c=c,
                a=num(vals[5]), e=num(vals[6]), p=num(vals[4]), lod=atlas.r6(vals[3]), q=None, ef=None, ve=None,
                lv=None, af=None, cc=len(xs), pa=None, anc=None, daf=None, up=-1, ups=ups, dns=dns,
                g=g, gf=None, gt=tier, pv=f"t15:S4#{rn}")
        row = verbatim(vals)
        orig[c]["t15:S4"][qid] = row
        if s7row is not None:
            orig[c]["t15:S7"][qid] = s7row
        for name, gi, conf in cands:
            T_c.append((tidx[cond], gi, name, conf, q, f"t15:S4#{rn}", c))

    # ---- S3: replicate-confirmed BY/RM QTL ---------------------------------------------
    for rn, vals, bold in s3rows:
        if not nonempty(vals[0]):
            continue
        c = chrom(vals[2])
        qid = f"{SC}-x-S3-{rn}"
        add("qtl", id=qid, ds=0, t=tidx[str(vals[1])], x=xby, xs=None, m="x", k="replicate-confirmed BY/RM QTL", c=c,
            a=num(vals[5]), e=num(vals[6]), p=num(vals[4]), lod=atlas.r6(vals[3]), q=None, ef=None, ve=None,
            lv=None, af=None, cc=None, pa=None, anc=None, daf=None, up=-1, ups=None, dns=None,
            g=-1, gf=None, gt=None, pv=f"t15:S3#{rn}")
        orig[c]["t15:S3"][qid] = verbatim(vals)

    # ---- S8: MAT-dependent QTL ---------------------------------------------------------
    s8cols, s8rows = read("t15:S8")
    for rn, vals, bold in s8rows:
        if not (nonempty(vals[0]) and str(vals[0]).startswith("Cross ")):
            continue   # the footnote row ("QTL identified in multiple conditions.") is kept in t15_extra.json
        n = str(vals[0]).split()[1]
        c = chrom(vals[2])
        cands = candidates(vals[7])
        g = next((gi for _, gi, _ in cands if gi >= 0), -1)
        qid = f"{SC}-x-S8-{rn}"
        q = add("qtl", id=qid, ds=0, t=tidx["MAT:" + str(vals[1])], x=xidx[n], xs=None, m="x",
                k="MAT-linked" + (" (bold in S8: found in multiple conditions)" if bold else ""), c=c,
                a=num(vals[5]), e=num(vals[6]), p=num(vals[4]), lod=atlas.r6(vals[3]), q=None, ef=None, ve=None,
                lv=None, af=None, cc=None, pa=None, anc=None, daf=None, up=-1, ups=None, dns=None,
                g=g, gf=None, gt="c" if g >= 0 else None, pv=f"t15:S8#{rn}")
        orig[c]["t15:S8"][qid] = verbatim(vals) + [bold]
        for name, gi, conf in cands:
            T_c.append((tidx["MAT:" + str(vals[1])], gi, name, conf, q, f"t15:S8#{rn}", c))

    # ---- candidate genes named by the paper ----------------------------------------------
    # tier "x": the paper confirmed the gene by allele replacement ("(confirmed)" in S4);
    # tier "c": a candidate gene the paper names for the QTL (S4, S8 "Candidate QTG").
    # Tokens that are not SGD gene names (e.g. "ENA locus") get g = -1 and keep their text.
    k = 0
    for t, gi, name, conf, q, pv, c in T_c:
        k += 1
        tier = "x" if conf else "c"
        add("causal", id=f"{SC}-c-{k}", ds=0, t=t, g=gi, orf=G["id"][gi] if gi >= 0 else None, name=name,
            c=G["c"][gi] if gi >= 0 else c, a=G["a"][gi] if gi >= 0 else None, e=G["e"][gi] if gi >= 0 else None,
            pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier=tier, pv=pv)
    cz = T["causal"]
    per_gene = collections.OrderedDict()
    for i in range(atlas.nlen(cz)):
        if cz["g"][i] >= 0:
            per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier="x" if any(cz["tier"][i] == "x" for i in ii) else "c", fdr=None,
            nt=len({cz["t"][i] for i in ii}), n=len(ii))

    # ---- dataset row ---------------------------------------------------------------------
    q = T["qtl"]
    add("datasets", id=DS_ID, sc=SC, label="Treusch 2015",
        cite="Treusch S, Albert FW, Bloom JS, Kotenko IE, Kruglyak L. 2015. PLoS Genet 11(1):e1004913",
        doi="10.1371/journal.pgen.1004913", year=2015,
        design="X-QTL (bulk segregant) in a 12-strain round-robin cross (each strain in two crosses); MATa and MATα pools of each cross as replicate selections; replicate BY × RM experiments",
        build="S288C R64 (sacCer3): 'aligning parent strain sequences to the sacCer3 reference genome' (Methods). Tables use chr01..chr16, renamed chrI..chrXVI; no coordinate transform.",
        method="Pools of 500,000 FACS-sorted MATa or MATα segregants grown on YPD, 0.5 or 1 M NaCl, 15 or 20 mM caffeine and bulk-sequenced; MULTIPOOL LOD for the selected pool against its YPD control (S3, S4) or MATa against MATα pools (S8); QTL at LOD > 5. "
               "Round-robin QTL had to pass LOD 5 in both mating-type selections of a cross (155 per-cross QTL) and were grouped across crosses into regions (S4), combining doses of the same condition.",
        fdr=None,
        interval="2-LOD drop intervals (Methods; S4 caption). Grouped QTL (S4): one interval per region, from the mean of replicated per-cross intervals; per-cross intervals are not reported.",
        effect="No effect sizes. Direction is reported only for grouped QTL, as the strains carrying beneficial and deleterious alleles (S7), in qtl.ups and qtl.dns.",
        nq=atlas.nlen(q),
        notes=f"{len(pidx)} parent strains (12 round-robin parents plus BY and RM), {len(xidx)} crosses (12 round-robin and BY × RM), {len(tidx)} traits. "
              f"{q['k'].count('grouped round-robin QTL')} grouped round-robin QTL regions, {q['k'].count('replicate-confirmed BY/RM QTL')} replicate-confirmed BY/RM QTL, "
              f"{sum(1 for v in q['k'] if v.startswith('MAT-linked'))} MAT-dependent QTL. QTL called at LOD 5, not at an FDR. Candidate genes are the paper's 'Candidate QTG' column; WHI2 and TOR1 were confirmed by allele replacement.")

    atlas.save_fragment(DS_ID, dict(T))

    # ---- lazy files ------------------------------------------------------------------------
    cols = {"t15:S3": s3cols, "t15:S4": s4cols, "t15:S7": s7cols[:5], "t15:S8": s8cols + ["bold (source formatting)"]}
    for ch in atlas.chrom_lengths():
        out = {"cols": cols, "rows": {code: dict(orig[ch].get(code, {})) for code in cols}}
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{ch}.json"), out)
    extra = {}
    for code in SHEETS:
        extra[code] = {"cols": read(code)[0], "rows": full_sheet(code)}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: {len(pidx)} parents, {len(xidx)} crosses, {len(tidx)} traits, {atlas.nlen(q)} QTL "
          f"(S4 {q['k'].count('grouped round-robin QTL')}, S3 {q['k'].count('replicate-confirmed BY/RM QTL')}, "
          f"S8 {sum(1 for v in q['k'] if v.startswith('MAT-linked'))}); causal {atlas.nlen(cz)} rows, {atlas.nlen(T['cand'])} genes; "
          f"unmapped candidate tokens: {sorted({cz['name'][i] for i in range(atlas.nlen(cz)) if cz['g'][i] < 0})}", file=sys.stderr)
    print(f"  docs/data.json: {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
