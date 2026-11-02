"""Ingest Jerison et al. 2017 (eLife 6:e27167): genetic variation in adaptability and pleiotropy.

    python3 -I ingest/jerison2017.py          # from the repository root

230 BY x RM segregants (from the Bloom et al. 2013 panel) were measured for initial fitness and
evolved for 500 generations in two environments (OT: YPD 30 C; HT: synthetic complete 37 C).
QTL were mapped for initial fitness and for fitness gain (adaptability) in each environment, and
loci were tested as terms in combined fitness + QTL models for adaptability and pleiotropy.

Reads only from sources/ (the unpacked eLife bundle in sources/jerison2017/raw/, the PMC full
text sources/jerison2017/article_pmc5580887.xml, and sources/jerison2017/extracted/). Writes
build/fragments/jerison2017.json, docs/v/j17_*.json and, through ingest/atlas.py, docs/data.json.

Semantics (see OPEN_QUESTIONS.md, Jerison 2017 section):
- Supplementary file 5 gives each of the 8 loci as one position ("ChrXIV 376315") plus the genes
  inside its 1.5-LOD confidence interval; no interval bounds and no LOD are tabulated. QTL rows
  are therefore points at that position (a = e = p).
- One QTL row per (locus, trait) marked "y" in Supplementary file 5. The first four trait columns
  come from the genome scan (Methods, "Mapping QTLs"); the last four mark loci kept as terms in
  the combined fitness + QTL models (Figure 3/8 source data).
- Coefficients come from Figure 3-source data 1 and Figure 8-source data 1. Three loci there are
  quoted at positions that differ from Supplementary file 5 by 2.6-28.5 kb; they are matched to the
  nearest Supplementary file 5 locus on the same chromosome (within 50 kb), and the match is kept.
- Sign (most defensible reading): positive = the RM allele raises the trait. The text says the RM
  allele at the KRE33 locus is less fit and makes founders more adaptable, and the KRE33
  coefficient is positive in both adaptability QTL models.
"""
import csv
import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "jerison2017", "j17"
SRC = os.path.join(atlas.SOURCES, DS_ID)
RAW = os.path.join(SRC, "raw")
PFX = "elife_poa_e27167_"
F = {
    "SF5": PFX + "Supplementary_File_5.csv",
    "SF2": PFX + "Supplementary_File_2.csv",
    "SF3": PFX + "Supplementary_File_3.csv",
    "F3SD1": PFX + "Figure_3_Source_Data_1.xlsx",
    "F8SD1": PFX + "Figure_8_Source_Data_1.xlsx",
    "F6SD1": PFX + "Figure_6_Source_Data_1.xlsx",
}
XML = os.path.join(SRC, "article_pmc5580887.xml")
TF = os.path.join(SRC, "extracted", "text_facts.tsv")

# Supplementary file 5 trait columns -> atlas traits (raw, display name, category, ck, kind label)
TRAITS = [
    ("QTL for trait: fitness at OT", "Initial fitness, YPD 30 °C (OT)", "reference medium", "ypd", "genome scan"),
    ("QTL for trait: fitness at HT", "Initial fitness, SC 37 °C (HT)", "temperature", "sc 37c", "genome scan"),
    ("QTL for trait: adaptability at OT", "Adaptability (fitness gain), OT", "adaptability", "adaptability ypd", "genome scan"),
    ("QTL for trait: adaptability at HT", "Adaptability (fitness gain), HT", "adaptability", "adaptability sc 37c", "genome scan"),
    ("QTL for trait: adaptability at OT beyond fitness", "Adaptability beyond founder fitness, OT", "adaptability", None, "term in combined fitness + QTL model"),
    ("QTL for trait: adaptability at HT beyond fitness", "Adaptability beyond founder fitness, HT", "adaptability", None, "term in combined fitness + QTL model"),
    ("QTL for trait: pleiotropy at OT beyond fitness", "Pleiotropy beyond founder fitness: HT-evolved lines at OT", "adaptability", None, "term in combined fitness + QTL model"),
    ("QTL for trait: pleiotropy at HT beyond fitness", "Pleiotropy beyond founder fitness: OT-evolved lines at HT", "adaptability", None, "term in combined fitness + QTL model"),
]
# coefficient tables (sheet in the xlsx) -> Supplementary file 5 trait column they describe
COEF_SHEETS = {
    ("F3SD1", "QTL model, adaptability at OT"): "QTL for trait: adaptability at OT",
    ("F3SD1", "QTL model, adaptability at HT"): "QTL for trait: adaptability at HT",
    ("F3SD1", "Combined Fitness and QTL Model "): "QTL for trait: adaptability at OT beyond fitness",
    ("F3SD1", "Combined Fitness and QTL Model1"): "QTL for trait: adaptability at HT beyond fitness",
    ("F8SD1", "Combined fitness and QTL model "): "QTL for trait: pleiotropy at HT beyond fitness",
    ("F8SD1", "Combined fitness and QTL model1"): "QTL for trait: pleiotropy at OT beyond fitness",
}
ROMAN = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7, "VIII": 8, "IX": 9, "X": 10, "XI": 11, "XII": 12, "XIII": 13, "XIV": 14, "XV": 15, "XVI": 16}


def num(v):
    if isinstance(v, float) and v.is_integer() and abs(v) < 1e15:
        return int(v)
    return v


def locus(s):
    """'ChrXIV 376315' -> ('chrXIV', 376315)."""
    m = re.match(r"\s*Chr([IVX]+)\s+(\d+)\s*$", s)
    return ("chr" + m.group(1), int(m.group(2)))


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def xlsx_rows(path, sheet):
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = [list(r) for r in wb[sheet].iter_rows(values_only=True)]
    return [[None if v is None else num(v) for v in r] for r in rows]


def xml_tables():
    """Table 1 and Table 2 of the article, as lists of rows of cell text."""
    r = ET.parse(XML).getroot()
    out = {}
    for tw in r.findall(".//table-wrap"):
        rows = [[" ".join("".join(c.itertext()).split()) for c in tr] for tr in tw.findall(".//tr")]
        out[tw.get("id")] = rows
    return out


def first_num(s):
    m = re.match(r"\s*(-?[\d.]+)", s or "")
    return float(m.group(1)) if m else None


def main():
    sha = manifest()
    T = {k: {} for k in atlas.TABLES}

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

    srcs = [("j17:SF5", F["SF5"], "Supplementary file 5 (csv)", 1),
            ("j17:F3SD1", F["F3SD1"], "Figure 3-source data 1 (sheets: QTL and combined models)", 1),
            ("j17:F8SD1", F["F8SD1"], "Figure 8-source data 1 (sheets: combined models for pleiotropy)", 1),
            ("j17:SF2", F["SF2"], "Supplementary file 2 (csv)", 1),
            ("j17:SF3", F["SF3"], "Supplementary file 3 (csv)", 1),
            ("j17:T1", "../article_pmc5580887.xml", "Table 1 (article, PMC full text); row 1 = header, row 3 = h2, row 4 = r2 of QTL models; the trait is the column", 1),
            ("j17:TF", "../extracted/text_facts.tsv", "quoted statements from the article", 3)]
    for code, fname, sheet, hrow in srcs:
        add("sources", code=code, ds=0, file=fname, sheet=sheet, hrow=hrow,
            sha256=sha.get(fname))

    add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1="BY4716-derived (Bloom et al. 2013 panel)", s2=None, p18=None, org="laboratory strain (text)")
    add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1="RM11-1a-derived (Bloom et al. 2013 panel)", s2=None, p18=None, org="wine strain (text)")
    add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="RM × BY (230 MATa founders)", p1=1, p2=0, g1=None, n=230, dip=None, mm=None, pv="j17:TF#4")

    tidx = {}
    for col, name, cat, ck, kind in TRAITS:
        tidx[col] = add("traits", id=f"{SC}-" + re.sub(r"[^A-Za-z0-9]+", "-", col.replace("QTL for trait: ", "")).strip("-"),
                        ds=0, raw=col.replace("QTL for trait: ", ""), name=name, cond=None, bt=None, ck=ck, cat=cat,
                        catb="atlas assignment (Supplementary file 5 trait column)")

    # ---- Supplementary file 5 -----------------------------------------------------------------
    with open(os.path.join(RAW, F["SF5"]), newline="") as f:
        sf5 = list(csv.reader(f))
    head = sf5[0]
    loci = []
    for k, r in enumerate(sf5[1:]):
        c, p = locus(r[0])
        loci.append({"row": k + 2, "c": c, "p": p, "r": r, "label": r[0].strip()})

    # ---- coefficients from the model tables ---------------------------------------------------
    coef = {}       # (sf5 trait column, locus index) -> dict(coef, ve, term, pv, interaction, src_locus)
    matches = []    # (table locus, sf5 locus, offset)
    for (fc, sheet), col in COEF_SHEETS.items():
        rows = xlsx_rows(os.path.join(RAW, F[fc]), sheet)
        hdr_i = next(i for i, r in enumerate(rows) if "Term" in r)
        off = rows[hdr_i].index("Term")          # some sheets start in column B
        rows = [r[off:] for r in rows]
        hdr = rows[hdr_i]
        for i, r in enumerate(rows[hdr_i + 1:], start=hdr_i + 2):
            if not r or r[0] is None or r[1] in ("intercept",) or (isinstance(r[1], str) and r[1].startswith("founder fitness")):
                continue
            parts = [s.strip() for s in str(r[1]).split(",")]
            inter = len(parts) > 1
            for part in parts:
                c, p = locus(part)
                best = min((abs(L["p"] - p), j) for j, L in enumerate(loci) if L["c"] == c)
                assert best[0] <= 50000, (part, best)
                j = best[1]
                if best[0]:
                    matches.append((part, loci[j]["label"], best[0]))
                if inter and part == parts[0] and len(parts) == 2 and loci[j]["label"] == "ChrXIV 376315":
                    continue    # interaction terms are recorded on the partner locus
                key = (col, j)
                entry = {"coef": None if inter else r[2], "ve": r[3] if isinstance(r[3], (int, float)) else None,
                         "term": r[0], "pv": f"j17:{fc}#{sheet.strip()}!{i}", "interaction": inter, "src": part, "row": r, "hdr": hdr}
                if key not in coef or (coef[key]["interaction"] and not inter):
                    coef[key] = entry

    for (col, j) in coef:     # every model-table term must sit on a locus x trait cell flagged "y" in SF5
        assert loci[j]["r"][head.index(col)].strip().lower() == "y", (col, loci[j]["label"])

    # ---- QTL rows -----------------------------------------------------------------------------
    G = atlas.genes()
    gname = {}
    for i, n in enumerate(G["n"]):
        if n:
            gname.setdefault(n, i)
    kre33 = gname["KRE33"]
    orig = {}
    cols = {"j17:SF5": head}
    qn = 0
    for j, L in enumerate(loci):
        for ci, (col, name, cat, ck, kind) in enumerate(TRAITS):
            flag = L["r"][head.index(col)].strip().lower()
            if flag != "y":
                continue
            qid = f"{SC}-w-{L['row']}-{ci + 1}"
            e = coef.get((col, j))
            ef = e["coef"] if e else None
            k = kind
            if e and e["interaction"]:
                k = "term in combined fitness + QTL model (interaction with the ChrXIV 376315 / KRE33 locus)"
            is_kre = L["label"] == "ChrXIV 376315"
            add("qtl", id=qid, ds=0, t=tidx[col], x=0, m="w", c=L["c"], a=L["p"], e=L["p"], p=L["p"],
                lod=None, q=None, ef=atlas.r6(ef), ve=atlas.r6(e["ve"]) if e else None, lv=L["label"],
                up=(1 if ef > 0 else 0) if ef else -1, k=k,
                g=kre33 if is_kre else -1, gf=None, gt="c" if is_kre else None, pv=f"j17:SF5#{L['row']}")
            orig.setdefault(L["c"], {}).setdefault("j17:SF5", {})[qid] = [num(v) if v != "" else None for v in L["r"]]
            if e:
                cols.setdefault(e["pv"].split("#")[0], e["hdr"])
                orig[L["c"]].setdefault(e["pv"].split("#")[0], {})[qid] = e["row"]
            qn += 1

    # candidate gene: KRE33 for the ChrXIV 376315 locus, named in the text
    kre_q = [i for i, l in enumerate(T["qtl"]["lv"]) if l == "ChrXIV 376315"]
    add("causal", id=f"{SC}-c-1", ds=0, t=T["qtl"]["t"][kre_q[0]], g=kre33, orf=G["id"][kre33], name="KRE33",
        c=G["c"][kre33], a=G["a"][kre33], e=G["e"][kre33], pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None,
        q1=kre_q[0], q2=-1, tier="c", pv="j17:TF#14")
    add("cand", ds=0, g=kre33, tier="c", fdr=None, nt=len({T["qtl"]["t"][i] for i in kre_q}), n=1)

    # ---- heritability and variance explained (Table 1) ------------------------------------------
    tabs = xml_tables()
    t1 = tabs["table1"]
    t1cols = {1: "QTL for trait: fitness at OT", 2: "QTL for trait: fitness at HT", 3: "QTL for trait: adaptability at OT",
              4: "QTL for trait: adaptability at HT", 5: "QTL for trait: pleiotropy at OT beyond fitness", 6: "QTL for trait: pleiotropy at HT beyond fitness"}
    rowby = {r[0]: r for r in t1[1:]}
    for ci, col in t1cols.items():
        add("vexp", t=tidx[col], x=0, wve=first_num(rowby["r2, QTL models"][ci]) if rowby["r2, QTL models"][ci] != "NA" else None,
            jve=None, h2=first_num(rowby["h2"][ci]), h2se=None, re=None, rese=None, pv="j17:T1#3")

    add("datasets", id=DS_ID, sc=SC, label="Jerison 2017",
        cite="Jerison ER, Kryazhimskiy S, Mitchell JK, Bloom JS, Kruglyak L, Desai MM. 2017. Genetic variation in adaptability and pleiotropy in budding yeast. eLife 6:e27167",
        doi="10.7554/eLife.27167", year=2017,
        design="segregant panel: 230 BY × RM MATa segregants from Bloom et al. 2013, 8 replicate lines each evolved 500 generations in OT (YPD 30 °C) or HT (SC 37 °C); fitness by competition against a common reference",
        build="S288C R64-compatible: positions are those of the Bloom et al. 2013 markers (Supplementary file 1 names them cumulativePos_chr_pos_ref_alt); the genes listed under each locus lie near it on R64-5-1. No transform.",
        method="Iterative forward search (after Bloom et al. 2013): LOD = -n/2 log10(1 - r^2) between residual trait and genotype, weighted by replicate counts; each step's top locus kept if significant at 0.05 by permutation; 1.5-LOD-drop confidence intervals. Run for initial fitness and fitness gain in OT and HT. Loci then tested as terms in combined fitness + QTL models (Bonferroni-corrected F tests).",
        fdr=None,
        interval="1.5-LOD drop (Methods), but bounds are not tabulated: Supplementary file 5 gives one position per locus and the genes inside its interval. QTL are drawn at that position.",
        effect="Coefficients from the linear models in Figure 3/8 source data, in fitness units (selection coefficient per generation relative to the reference strain); the combined models code each locus +0.5/-0.5. positive = the RM allele raises the trait (most defensible reading: the text calls the RM allele at the KRE33 locus less fit and more adaptable, and its coefficient is positive in both adaptability QTL models). Initial-fitness QTL have no coefficients in the bundle.",
        nq=qn,
        notes="1 cross, 230 founders, 8 loci, 8 traits (initial fitness and adaptability by genome scan in two environments; adaptability and pleiotropy beyond founder fitness as model terms). Candidate gene: KRE33, named in the text for the ChrXIV 376315 locus (tier c); de novo mutations in evolved lines are kept verbatim, not drawn.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy/verbatim files ------------------------------------------------------------------------
    for chrom in atlas.chrom_lengths():
        out = {"cols": cols, "rows": orig.get(chrom, {})}
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"), out)
    extra = {}
    for code, key in (("j17:SF2", "SF2"), ("j17:SF3", "SF3")):
        with open(os.path.join(RAW, F[key]), newline="") as f:
            rows = list(csv.reader(f, delimiter=";" if key == "SF2" else ","))
        extra[code] = {"cols": rows[0], "rows": rows[1:]}
    for (fc, sheet) in list(COEF_SHEETS) :
        extra[f"j17:{fc}#{sheet.strip()}"] = {"cols": [], "rows": xlsx_rows(os.path.join(RAW, F[fc]), sheet)}
    extra["j17:F6SD1"] = {"cols": [], "rows": xlsx_rows(os.path.join(RAW, F["F6SD1"]), "Sheet 1 - Mutual Information wi")}
    extra["j17:T1"] = {"cols": tabs["table1"][0], "rows": tabs["table1"][1:]}
    extra["j17:T2"] = {"cols": tabs["table2"][0], "rows": tabs["table2"][1:]}
    extra["j17:locus_matches"] = {"cols": ["model-table locus", "Supplementary file 5 locus", "offset bp"], "rows": sorted(set(matches))}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: 8 loci, {len(TRAITS)} traits, {qn} QTL rows, coefficients on {sum(1 for v in T['qtl']['ef'] if v is not None)}; "
          f"locus matches {sorted(set(matches))}; data.json {len(data['datasets']['id'])} datasets", file=sys.stderr)


if __name__ == "__main__":
    main()
