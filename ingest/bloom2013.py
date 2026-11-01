"""Ingest Bloom et al. 2013 (Nature 494:234, doi 10.1038/nature11867), QTL for 46 traits in BY x RM segregants.

    python3 -I ingest/bloom2013.py          # from the repository root

Reads only from sources/: the supplementary tables in sources/bloom2013/raw/ (Supplementary
Tables 1-3) and the SGD tables in sources/sgd/. Writes build/fragments/bloom2013.json,
docs/v/b13_*.json and, through ingest/atlas.py, docs/data.json. Deterministic.

What the bundle settles and what it does not (OPEN_QUESTIONS.md, Bloom 2013 section):
- Supplementary Table 3 is the QTL table: trait, fraction of phenotypic variance explained,
  chromosome (1-16), peak position (bp), LOD, 1.5-LOD confidence interval (bp), and the genes
  under the interval. It gives no effect direction, so qtl.ef and qtl.up stay empty.
- The Online Methods are paywalled on the article page, so how QTL were detected and how the
  table's LOD was computed are not in the bundle. Some LODs are far below any detection
  threshold (minimum 0.0001); the most defensible reading is that the column is each QTL's
  LOD in the final multiple-QTL model, conditional on the trait's other QTL.
- Genome build: the table's coordinates agree with S288C R64. Every SGD systematic name the
  table lists under an interval lies within 5 kb of that interval on R64-5-1 coordinates
  (tests/test_bloom2013.py checks this). No transform.
"""
import collections
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "bloom2013", "b13"
SRC = os.path.join(atlas.SOURCES, "bloom2013")
RAW = os.path.join(SRC, "raw")
S1 = "41586_2013_BFnature11867_MOESM86_ESM.xls"
S2 = "41586_2013_BFnature11867_MOESM87_ESM.xls"
S3 = "41586_2013_BFnature11867_MOESM88_ESM.xls"
SHEETS = {"b13:S1": (S1, "SuppTable1"), "b13:S2": (S2, "SuppTable2"), "b13:S3": (S3, "SuppTable3")}
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]

# Categories. This bundle has no media table, so they follow the Bloom 2019 assignment for the
# same compounds; temperature, pH and the reference media are read from the trait name.
ANALOGY = "atlas assignment by analogy with Bloom 2019 (this bundle has no media table)"
CARBON = {"Ethanol", "Galactose", "Lactate", "Lactose", "Maltose", "Mannose", "Raffinose", "Trehalose", "Xylose"}
OSMOTIC = {"Sorbitol"}
KEEP_CASE = {"4NQO", "SDS", "YPD", "YNB", "E6_Berbamine"}


def num(v):
    if isinstance(v, float) and v.is_integer() and abs(v) < 1e15:
        return int(v)
    return v


def verbatim(row):
    return [None if v == "" else num(v) for v in row]


_books = {}


def book(f):
    if f not in _books:
        import xlrd
        _books[f] = xlrd.open_workbook(os.path.join(RAW, f), formatting_info=True)
    return _books[f]


def read(code):
    """(column names, [(spreadsheet row, values)]) for a provenance code; header on row 1."""
    f, sheet = SHEETS[code]
    sh = book(f).sheet_by_name(sheet)
    rows = [sh.row_values(i) for i in range(sh.nrows)]
    return [str(num(v)) for v in rows[0]], [(i + 1, r) for i, r in enumerate(rows) if i >= 1]


def shown(f, sheet, r, c):
    """A cell as Excel displays it: percent-formatted numbers become 'n%' (0.02 with format '0%' -> '2%')."""
    b = book(f); sh = b.sheet_by_name(sheet); cell = sh.cell(r, c)
    if cell.ctype == 2:
        fmt = b.format_map[b.xf_list[cell.xf_index].format_key].format_str
        if fmt.endswith("%"):
            dec = len(fmt.split(".")[1]) - 1 if "." in fmt else 0
            return f"{cell.value * 100:.{dec}f}%"
        return str(num(cell.value))
    return str(cell.value).strip()


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as fh:
        for line in fh.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def slug(s):
    return re.sub(r"-+", "-", re.sub(r"[^A-Za-z0-9]+", "-", s)).strip("-")


def base_name(comp):
    if comp in KEEP_CASE:
        return comp.replace("_", " ")
    words = comp.split("_")
    return " ".join([words[0]] + [w.lower() for w in words[1:]])


def trait_fields(raw, dose):
    """Display name, condition, condition key and category for a Table S1/S3 trait string."""
    comp, _, mod = raw.partition(":")
    if comp == "YPD" and mod:
        t = mod.rstrip("Cc")
        return f"YPD {t} °C", mod, f"ypd {t}", "temperature", "trait name ('" + raw + "')"
    if comp == "YNB" and mod:
        return f"YNB pH {mod[2:]}", mod, f"ynb {mod.lower()}", "pH", "trait name ('" + raw + "')"
    if comp in ("YPD", "YNB"):
        return comp, "", comp.lower(), "reference medium", "trait name ('" + raw + "')"
    disp_dose = re.sub(r"(?<=[\d ])u(?=[gmM])", "µ", dose)
    name = (base_name(comp) + (" " + disp_dose if dose else "")).strip()
    # condition key in Bloom 2019's style (compound lower case, dose without spaces, unit case kept; "/ml" written "/mL")
    ck = (comp.replace("_", " ").lower() + (" " + dose.replace(" ", "").replace("/ml", "/mL") if dose else "")).strip()
    cat = "carbon source" if comp in CARBON else "osmotic stress" if comp in OSMOTIC else "drug/chemical"
    return name, dose, ck, cat, ANALOGY


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

    for code, (f, sheet) in SHEETS.items():
        add("sources", code=code, ds=0, file=f, sheet=sheet, hrow=1, sha256=sha[f])

    s1cols, s1rows = read("b13:S1")
    s3cols, s3rows = read("b13:S3")
    s2cols, s2rows = read("b13:S2")

    # one cross; segregant count is not stated in the bundle (Methods paywalled): use the largest
    # number of segregants with at least one measurement for any trait (Table S1, column D)
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1=None, s2=None, p18=None, org=None)
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1=None, s2=None, p18=None, org=None)
    nmax_row, nmax = max(((rn, num(r[3])) for rn, r in s1rows), key=lambda t: (t[1], -t[0]))
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=nmax,
            dip=None, mm=None, pv=f"b13:S1#{nmax_row}")

    # traits, in Table S1 order
    tidx = {}
    for rn, r in s1rows:
        raw = r[0]
        dose = shown(S1, "SuppTable1", rn - 1, 1)
        name, cond, ck, cat, catb = trait_fields(raw, dose)
        tidx[raw] = add("traits", id=f"{SC}-{slug(raw)}", ds=0, raw=raw, name=name, cond=cond, bt=str(num(r[2])),
                        ck=ck, cat=cat, catb=catb)
        add("vexp", t=tidx[raw], x=x, wve=atlas.r6(r[12]), jve=None, h2=atlas.r6(r[7]), h2se=atlas.r6(r[8]),
            re=None, rese=None, pv=f"b13:S1#{rn}")

    # QTL
    orig = collections.defaultdict(dict)
    for rn, r in s3rows:
        chrom = "chr" + ROMAN[int(r[2]) - 1]
        qid = f"{SC}-w-{rn}"
        add("qtl", id=qid, ds=0, t=tidx[r[0]], x=x, xs=None, m="w", k=None, c=chrom, a=num(r[5]), e=num(r[6]),
            p=num(r[3]), cm=None, lod=atlas.r6(r[4]), q=None, ef=None, ve=atlas.r6(r[1]), lv=None, af=None, cc=None,
            pa=None, anc=None, daf=None, up=-1, ups=None, dns=None, g=-1, gf=None, gt=None, pv=f"b13:S3#{rn}")
        orig[chrom][qid] = verbatim(r)

    nq = atlas.nlen(T["qtl"])
    add("datasets", id=DS_ID, sc=SC, label="Bloom 2013",
        cite="Bloom JS, Ehrenreich IM, Loo WT, Lite T-LV, Kruglyak L. 2013. Finding the sources of missing heritability in a yeast cross. Nature 494:234",
        doi="10.1038/nature11867", year=2013,
        design="segregant panel; one cross, BY × RM, segregants genotyped by sequencing (Supplementary Figure 1); end-point colony growth on agar plates for 46 traits",
        build="S288C R64-compatible: every SGD systematic name listed under an interval in Supplementary Table 3 lies within 5 kb of that interval on R64-5-1 coordinates. No transform.",
        method="Not in the bundle: the Online Methods are paywalled on the article page. Supplementary Table 3 reports, per QTL, the LOD, the peak, a 1.5-LOD confidence interval and the fraction of phenotypic variance explained; some LODs are below 1, so the column is read as the LOD in the final multiple-QTL model (most defensible reading).",
        fdr=None,
        interval="1.5-LOD confidence interval, left and right positions in bp, as given in Supplementary Table 3.",
        effect="Not signed in the bundle. Supplementary Table 3 gives the fraction of phenotypic variance explained by each QTL (qtl.ve); Supplementary Figure 3 shows absolute standardized allelic differences as an image only. qtl.ef is empty.",
        nq=nq,
        notes=f"1 cross (BY × RM), {len(tidx)} traits, {nq} QTL. crosses.n = {nmax}, the largest number of segregants with at least one measurement for any trait (Supplementary Table 1, column D); per trait it ranges from "
              f"{min(num(r[3]) for _, r in s1rows)} to {nmax}; the paper's own total is in the paywalled Methods. "
              "vexp: wve is 'Variance Explained Additive QTL' and h2 the narrow-sense heritability (average of replicates) from Supplementary Table 1; whether wve was cross-validated is not stated. "
              "No candidate genes: the table lists every gene under each interval, which the atlas keeps in the source row, not as candidates.")

    atlas.save_fragment(DS_ID, dict(T))

    for chrom in atlas.chrom_lengths():
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"),
                         {"cols": {"b13:S3": s3cols}, "rows": {"b13:S3": dict(orig.get(chrom, {}))}})
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"),
                     {"b13:S1": {"cols": s1cols, "rows": [verbatim(r) for _, r in s1rows]},
                      "b13:S2": {"cols": s2cols, "rows": [verbatim(r) for _, r in s2rows]}})

    data = atlas.assemble()
    print(f"{DS_ID}: 2 parents, 1 cross (n={nmax}), {len(tidx)} traits, {nq} QTL; vexp {atlas.nlen(T['vexp'])}", file=sys.stderr)
    print(f"  docs/data.json: {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
