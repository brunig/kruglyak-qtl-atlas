"""Ingest Foss et al. 2011 (PLoS Biol 9:e1001144): protein (mass spectrometry) and transcript
(microarray) linkages for 354 genes in 95 BY x RM segregants.

    python3 -I ingest/foss2011.py          # from the repository root

Reads only from sources/foss2011/ (and the shared SGD tables). Table S4 lists every linkage at
5% FDR as one marker per row, with no p-value, interval or direction. The 521 linkages go into
two lazily loaded per-gene tables, like Albert 2014's X-pQTL:
  - "mspqtl"  docs/v/f11_mspqtl.json  protein levels from mass spectrometry (179 linkages)
  - "eqtl"    docs/v/f11_eqtl.json    transcript levels from microarrays (342 linkages)
with per-50-kb density rows for the heatmap. No hotspot table is published, so there are no
qtl rows. Marker positions are used as given (consistent with R64 at marker resolution; see
OPEN_QUESTIONS). Deterministic: same inputs, same bytes out.
"""
import collections
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "foss2011", "f11"
SRC = os.path.join(atlas.SOURCES, DS_ID)
RAW = os.path.join(SRC, "raw")
P = "journal.pbio.1001144."
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
HB = atlas.HEAT_BIN
LOCAL_BP = 20000   # atlas rule: a linkage is local when its marker lies within 20 kb of the regulated gene (same chromosome)
TF = "../extracted/text_facts.tsv"
KIND = {"peptide": ("mspqtl", "protein"), "transcript": ("eqtl", "mRNA")}


def num(v):
    if isinstance(v, float) and v.is_integer() and abs(v) < 1e15:
        return int(v)
    return v


def verbatim(row):
    return [None if v == "" else num(v) for v in row]


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def xls(name):
    import xlrd
    sh = xlrd.open_workbook(os.path.join(RAW, P + name)).sheet_by_index(0)
    rows = [sh.row_values(i) for i in range(sh.nrows)]
    return sh.name, [str(v) for v in rows[0]], [(i + 1, r) for i, r in enumerate(rows) if i >= 1]


def facts():
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        lines = [l.rstrip("\n") for l in f if not l.startswith("#")]
    head = lines[0].split("\t")
    rows = [(i + 2, l.split("\t")) for i, l in enumerate(lines[1:])]   # spreadsheet-style row numbers, header row 1 after the comment
    return head, rows


def main():
    sha = manifest()
    sha[TF] = hashlib.sha256(open(os.path.join(RAW, TF), "rb").read()).hexdigest()
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

    s1sheet, s1cols, s1rows = xls("s006.xls")
    s3sheet, s3cols, s3rows = xls("s008.xls")
    s4sheet, s4cols, s4rows = xls("s009.xls")
    add("sources", code="f11:S1", ds=0, file=P + "s006.xls", sheet=s1sheet, hrow=1, sha256=sha[P + "s006.xls"])
    add("sources", code="f11:S3", ds=0, file=P + "s008.xls", sheet=s3sheet, hrow=1, sha256=sha[P + "s008.xls"])
    add("sources", code="f11:S4", ds=0, file=P + "s009.xls", sheet=s4sheet, hrow=1, sha256=sha[P + "s009.xls"])
    add("sources", code="f11:TF", ds=0, file=TF, sheet="text_facts", hrow=1, sha256=sha[TF])
    tfh, tfr = facts()
    design_row = [rn for rn, r in tfr if r[0] == "design"][0]

    # ---- parents, cross, traits ------------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1="BY (laboratory)", s2=None, p18=None, org="laboratory")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1="RM (vineyard)", s2=None, p18=None, org="vineyard")
    add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=95, dip=None, mm=None,
        pv=f"f11:TF#{design_row}")
    add("traits", id=f"{SC}-protein-abundance-mass-spec", ds=0, raw="protein abundance (mass spectrometry)",
        name="Protein abundance (mass spec)", cond=None, bt=None, ck="protein abundance", cat="expression",
        catb="atlas assignment: one peptide per protein quantified by LC-MS in each segregant (Foss et al. 2011)")
    add("traits", id=f"{SC}-transcript-abundance-microarray", ds=0, raw="transcript abundance (microarray)",
        name="Transcript abundance (microarray)", cond=None, bt=None, ck="transcript abundance", cat="expression",
        catb="atlas assignment: transcript levels of the same 354 genes from earlier microarray studies (Foss et al. 2011)")

    # ---- linkages (Table S4) -> two lazy tables + density -------------------------------------
    G = atlas.genes()
    gid = {}
    for i, s in enumerate(G["id"]):
        gid.setdefault(s, i)
    h4 = {c: i for i, c in enumerate(s4cols)}
    tables = {k: {c: [] for c in ["id", "g", "c", "a", "e", "p", "lod", "ef", "loc", "pv"]} for k, _ in KIND.values()}
    dens = collections.Counter()
    orig = collections.defaultdict(dict)
    for rn, r in s4rows:
        kind, label = KIND[r[h4["type_of_measurement"]]]
        orf = r[h4["y_name"]]
        c = "chr" + ROMAN[int(r[h4["chromosome"]]) - 1]
        p = int(r[h4["position"]])
        g = gid.get(orf, -1)
        loc = 1 if g >= 0 and G["c"][g] == c and max(0, G["a"][g] - p, p - G["e"][g]) <= LOCAL_BP else 0
        lid = f"{SC}-{'p' if kind == 'mspqtl' else 't'}-{rn}"
        L = tables[kind]
        for k, v in (("id", lid), ("g", g), ("c", c), ("a", p), ("e", p), ("p", p), ("lod", None), ("ef", None), ("loc", loc), ("pv", f"f11:S4#{rn}")):
            L[k].append(v)
        dens[(label + (", local" if loc else ", distant"), c, p // HB)] += 1
        orig[c][lid] = verbatim(r)
    chroms = list(atlas.chrom_lengths())
    for (k, c, b), n in sorted(dens.items(), key=lambda kv: (kv[0][0], chroms.index(kv[0][1]), kv[0][2])):
        add("dens", ds=0, k=k, c=c, b=b, n=n)
    np_, nt = len(tables["mspqtl"]["id"]), len(tables["eqtl"]["id"])
    add("lazy", ds=0, name="mspqtl", file=f"v/{SC}_mspqtl.json", n=np_,
        desc=f"The {np_} protein linkages of Table S4 (5% FDR; protein levels by mass spectrometry): regulated gene and linked marker. "
             "No p-value, interval or direction is published, so each is a point at the marker.")
    add("lazy", ds=0, name="eqtl", file=f"v/{SC}_eqtl.json", n=nt,
        desc=f"The {nt} transcript linkages of Table S4 (5% FDR; transcript levels by microarray, for the same 354 genes): regulated gene and linked marker. "
             "No p-value, interval or direction is published, so each is a point at the marker.")

    # ---- dataset row -------------------------------------------------------------------------
    nloc = {k: sum(tables[k]["loc"]) for k in tables}
    add("datasets", id=DS_ID, sc=SC, label="Foss 2011",
        cite="Foss EJ, Radulovic D, Shaffer SA, Goodlett DR, Kruglyak L, Bedalov A. 2011. Genetic variation shapes protein networks mainly through non-transcriptional mechanisms. PLoS Biol 9:e1001144",
        doi="10.1371/journal.pbio.1001144", year=2011,
        design="segregant panel BY × RM, 95 segregants; protein levels of 354 genes by LC-MS (new quantitation of previously published spectra) and transcript levels of the same genes from earlier microarray studies; linkage mapped for both",
        build="S288C; not stated. All marker positions fall inside R64 chromosomes, and the transcript and protein hotspot markers sit on the R64 positions of LEU2 (chrIII ~91 kb), AMN1 (chrII 555,787 vs gene 556,201) and near HAP1 (chrXII). Used as given, no transform; agreement is at marker resolution, not base-pair.",
        method="t test of each of 1,969 unique-pattern markers (2,955 genotyped) against each protein and transcript level; per-gene FDR from 100 permutations. Table S4 lists the 5% FDR linkages (the text also reports a 1% FDR subset, not separable in the table).",
        fdr=0.05,
        interval="None: each linkage is the single linked marker (a = e = p). Hotspots are described but not tabulated.",
        effect="Not reported per linkage; Table S4 gives no direction or size.",
        nq=atlas.nlen(T["qtl"]),
        notes=f"1 cross, 95 segregants, 354 genes measured at both levels; {np_} protein linkages (125 proteins) and {nt} transcript linkages (200 transcripts) at 5% FDR, loaded per gene on demand. "
              f"Local here = marker within {LOCAL_BP // 1000} kb of the regulated gene ({nloc['mspqtl']} protein, {nloc['eqtl']} transcript linkages); an atlas rule, the paper's cis list is not published.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy, verbatim and extra files ------------------------------------------------------------
    for k, L in tables.items():
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{k}.json"), L)
    for c in chroms:
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{c}.json"), {"cols": {"f11:S4": s4cols}, "rows": {"f11:S4": dict(orig.get(c, {}))}})
    extra = {"f11:S1": {"cols": s1cols, "rows": [verbatim(r) for _, r in s1rows]},
             "f11:S3": {"cols": s3cols, "rows": [verbatim(r) for _, r in s3rows]},
             "f11:TF": {"cols": tfh, "rows": [r for _, r in tfr]}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: {np_} protein + {nt} transcript linkages; local (≤{LOCAL_BP // 1000} kb) {nloc['mspqtl']} protein, {nloc['eqtl']} transcript; "
          f"{atlas.nlen(T['dens'])} density bins; data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
