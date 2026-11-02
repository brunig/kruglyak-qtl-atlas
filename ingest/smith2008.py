"""Ingest Smith & Kruglyak 2008 (PLoS Biol 6:e83): expression QTL in 109 BY x RM segregants grown in
glucose and in ethanol, and gene-environment interaction QTL (linkage of the ethanol - glucose difference).

    python3 -I ingest/smith2008.py          # from the repository root

Reads only from sources/smith2008/ (and the shared SGD tables).
  - Table S4 (st004.xls): every significant linkage, 11,392 rows, with LOD, 1.5-LOD support interval
    markers in bp, local/distant class, distant-peak name and BY / RM allele averages per condition.
    All rows go into the lazily loaded table "eqtl" (docs/v/s08_eqtl.json), sharing the page's mRNA eQTL
    track with Albert 2018. Extra columns: env (glucose / ethanol / gxe), gxe (1 for interaction QTL),
    up (parent whose allele raises the value), st (probe status: rows with a polymorphism in the
    microarray probe are kept but flagged, because the paper excluded them from all analyses; the
    text's counts are the polymorphism-free rows).
  - Table S1 (st001.xls): the 41 distant-linkage peaks (13 glucose, 13 ethanol, 15 gxe) as qtl rows,
    one trait per condition. The table gives cM ranges only; the atlas places each peak at the span of
    the peak markers (bp) of the polymorphism-free Table S4 linkages the paper assigned to it (column
    "distant peak"), with p at the most frequent peak marker. Candidate or verified genes -> causal rows.
Positions: marker bp positions as given. Markers named after ORFs sit inside their R64-5-1 genes
(99.9% within 2 kb), so no transform. Deterministic: same inputs, same bytes out.
"""
import collections
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "smith2008", "s08"
SRC = os.path.join(atlas.SOURCES, DS_ID)
RAW = os.path.join(SRC, "raw")
P = "journal.pbio.0060083."
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
HB = atlas.HEAT_BIN
TF = "../extracted/text_facts.tsv"
ENV = {"glucose": "glucose", "ethanol": "ethanol", "gxeQTL": "gxe"}
PEAK_ENV = {"glu": "glucose", "eth": "ethanol", "gxe": "gxe"}
# genes the paper reports as characterized at the gene level: AMN1, GPA1, MAT, HAP1 (cited) and IRA2,
# DIG1, MKT1 by allele replacement in this paper (text, Figure S3, Figure S6) -> tier x; others -> c
CONFIRMED = {"AMN1", "GPA1", "MAT", "HAP1", "IRA2", "DIG1", "MKT1"}


def num(v):
    if isinstance(v, float) and v.is_integer() and abs(v) < 1e15:
        return int(v)
    return v


def verbatim(row):
    return [None if v == "" else num(v) for v in row]


def fnum(v):
    return None if v in ("", None, "NA") else float(v)


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
    with open(os.path.join(SRC, "extracted", "text_facts.tsv"), encoding="utf-8") as f:
        lines = [l.rstrip("\n") for l in f]
    head = lines[1].split("\t")   # line 1 is a comment, line 2 the header
    rows = [(i + 3, l.split("\t")) for i, l in enumerate(lines[2:])]
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

    s1sheet, s1cols, s1rows = xls("st001.xls")
    s2sheet, s2cols, s2rows = xls("st002.xls")
    s3sheet, s3cols, s3rows = xls("st003.xls")
    s4sheet, s4cols, s4rows = xls("st004.xls")
    s1rows = [(rn, r) for rn, r in s1rows if isinstance(r[1], float)]    # peak rows only (drop the blank line and footnote)
    s1note = [r[0] for _, r in xls("st001.xls")[2] if str(r[0]).startswith("*")]
    for code, f, sh in (("s08:S1", "st001.xls", s1sheet), ("s08:S2", "st002.xls", s2sheet), ("s08:S3", "st003.xls", s3sheet), ("s08:S4", "st004.xls", s4sheet)):
        add("sources", code=code, ds=0, file=P + f, sheet=sh, hrow=1, sha256=sha[P + f])
    add("sources", code="s08:TF", ds=0, file=TF, sheet="text_facts", hrow=2, sha256=sha[TF])
    tfh, tfr = facts()
    fact_row = {r[0]: rn for rn, r in tfr}

    # ---- parents, cross, traits ------------------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY4716", s1="BY4716 (MATalpha lys2D0; Table S2)", s2=None, p18=None, org=None)
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM11-1a", s1="RM11-1a (MATa leu2D0 ura3D0; Table S2)", s2=None, p18=None, org=None)
    X = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=109, dip=None, mm=None,
            pv=f"s08:TF#{fact_row['design']}")
    tix = {}
    for env, label in (("glucose", "glucose"), ("ethanol", "ethanol"), ("gxe", "gene–environment interaction (ethanol − glucose)")):
        tix[env] = add("traits", id=f"{SC}-hotspots-{env}", ds=0, raw=f"distant-linkage peaks ({env})",
                       name="Expression hotspots, " + label, cond=env, bt=None, ck="expression hotspots " + env, cat="expression",
                       catb="atlas assignment: peaks of distant expression linkages, Table S1 (Smith & Kruglyak 2008)")

    # ---- Table S4: every linkage -> lazy eqtl table ------------------------------------------------
    G = atlas.genes()
    gid = {}
    for i, s in enumerate(G["id"]):
        gid.setdefault(s, i)
    name2g = {}
    for i, n in enumerate(G["n"]):
        if n:
            name2g.setdefault(n, i)
    h = {c: i for i, c in enumerate(s4cols)}
    cols = ["id", "g", "c", "a", "e", "p", "lod", "ef", "loc", "pv", "env", "gxe", "up", "st"]
    L = {c: [] for c in cols}
    dens = collections.Counter()
    orig = collections.defaultdict(dict)
    peak_rows = collections.defaultdict(list)
    for rn, r in s4rows:
        env = ENV[r[h["linkage type"]]]
        c = "chr" + ROMAN[int(r[h["chromosome of peak marker"]]) - 1]
        a = int(r[h["LOD support interval left marker position (bp)"]])
        p = int(r[h["marker closest to peak position (bp)"]])
        e = int(r[h["LOD support interval right marker position (bp)"]])
        if env == "gxe":
            byv, rmv = fnum(r[h["BY allele average difference between ethanol and glucose"]]), fnum(r[h["RM allele average difference between ethanol and glucose"]])
        else:
            byv, rmv = fnum(r[h[f"BY allele average effect in {env}"]]), fnum(r[h[f"RM allele average effect in {env}"]])
        ef = atlas.r6(rmv - byv) if byv is not None and rmv is not None else None
        up = "RM" if ef is not None and ef > 0 else "BY" if ef is not None and ef < 0 else None
        poly = r[h["polymorphism in probe?"]]
        st = None if poly == "no" else ("probe polymorphism: excluded by the paper" if poly == "yes" else "probe polymorphism not determined (NA)")
        lid = f"{SC}-e-{rn}"
        vals = {"id": lid, "g": gid.get(r[h["ORF"]], -1), "c": c, "a": a, "e": e, "p": p, "lod": atlas.r6(r[h["LOD score"]]), "ef": ef,
                "loc": 1 if r[h["local or distant"]] == "local" else 0, "pv": f"s08:S4#{rn}", "env": env, "gxe": 1 if env == "gxe" else 0, "up": up, "st": st}
        for k in cols:
            L[k].append(vals[k])
        orig[c][lid] = verbatim(r)
        if st is None:
            for b in range(a // HB, e // HB + 1):
                dens[(env, c, b)] += 1
            pk = str(r[h["distant peak"]])
            if pk[:3] in PEAK_ENV:
                peak_rows[pk].append((c, p))
    chroms = list(atlas.chrom_lengths())
    for (k, c, b), n in sorted(dens.items(), key=lambda kv: (kv[0][0], chroms.index(kv[0][1]), kv[0][2])):
        add("dens", ds=0, k=k, c=c, b=b, n=n)
    nl = len(L["id"])
    add("lazy", ds=0, name="eqtl", file=f"v/{SC}_eqtl.json", n=nl,
        desc=f"All {nl} linkages of Table S4 (glucose, ethanol and gene–environment interaction; 5% FDR) with 1.5-LOD support intervals; "
             "rows whose microarray probe carries a BY/RM polymorphism are flagged (the paper excluded them).")

    # ---- Table S1: distant-linkage peaks -> qtl rows; candidates -> causal ------------------------
    qkey = {}
    causal_rows = []
    h1 = {c: i for i, c in enumerate(s1cols)}
    for rn, r in s1rows:
        name = r[h1["Peak"]]
        env = PEAK_ENV[name[:3]]
        c = "chr" + ROMAN[int(r[h1["Chr"]]) - 1]
        pts = [pp for cc, pp in peak_rows.get(name, []) if cc == c]
        if pts:
            cnt = collections.Counter(pts)
            top = max(cnt.values())
            p = min(pp for pp, n in cnt.items() if n == top)
            a, e = min(pts), max(pts)
        else:
            a = e = p = None
        gene = str(r[h1["Candidate or verified gene"]]).strip()
        qid = f"{SC}-w-{rn}"
        q = add("qtl", id=qid, ds=0, t=tix[env], x=X, m="w", c=c, a=a, e=e, p=p, lod=None, q=None, ef=None, ve=None, lv=None,
                cc=len(peak_rows.get(name, [])), up=-1, g=-1, gf=None, gt=None,
                k=f"distant-linkage peak {name} ({r[h1['Position (cM)']]} cM); extent = span of its linkages' peak markers",
                pv=f"s08:S1#{rn}")
        qkey[name] = q
        orig[c][qid] = verbatim(r)
        if gene:
            star = gene.endswith("*")
            gname = gene.rstrip("*")
            g = name2g.get(gname, -1)
            tier = "x" if gname in CONFIRMED and not star else "c"
            causal_rows.append((rn, env, q, g, gname, tier, star))
    cz = []
    for rn, env, q, g, gname, tier, star in causal_rows:
        k = add("causal", id=f"{SC}-c-{rn}", ds=0, t=tix[env], g=g, orf=G["id"][g] if g >= 0 else None, name=gname,
                c=G["c"][g] if g >= 0 else T["qtl"]["c"][q], a=G["a"][g] if g >= 0 else None, e=G["e"][g] if g >= 0 else None,
                pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier=tier, pv=f"s08:S1#{rn}")
        cz.append(k)
        if g >= 0:
            T["qtl"]["g"][q] = g
            T["qtl"]["gt"][q] = tier
    rank = {"x": 0, "c": 1}
    per = collections.OrderedDict()
    for k in cz:
        g = T["causal"]["g"][k]
        if g >= 0:
            per.setdefault(g, []).append(k)
    for g, ks in per.items():
        add("cand", ds=0, g=g, tier=min((T["causal"]["tier"][k] for k in ks), key=rank.get), fdr=None,
            nt=len({T["causal"]["t"][k] for k in ks}), n=len(ks))

    # ---- dataset row --------------------------------------------------------------------------------
    clean = [i for i in range(nl) if L["st"][i] is None]
    cnt = collections.Counter(L["env"][i] for i in clean)
    add("datasets", id=DS_ID, sc=SC, label="Smith 2008",
        cite="Smith EN, Kruglyak L. 2008. Gene–Environment Interaction in Yeast Gene Expression. PLoS Biol 6:e83",
        doi="10.1371/journal.pbio.0060083", year=2008,
        design="segregant panel BY × RM, 109 segregants expression-profiled (Agilent two-colour arrays) in glucose and in ethanol; linkage of expression in each condition and of the ethanol − glucose difference (gene–environment interaction)",
        build="S288C; positions as printed in Table S4. Markers named after ORFs sit inside their R64-5-1 genes (99.9% within 2 kb, none on another chromosome), so positions are used as R64 with no transform.",
        method="Nonparametric interval mapping (R/qtl) with 2,894 markers; genome-wide FDR by permutation: LOD 3.25 (glucose), 3.3 (ethanol), 3.7 (difference) for about 5% FDR. Transcripts whose probe carries a BY/RM polymorphism were excluded from the paper's analyses.",
        fdr=0.05,
        interval="Linkages: 1.5-LOD support interval (left and right markers, bp, from Table S4). Hotspots (Table S1, cM only): span of the peak markers of the polymorphism-free linkages assigned to each peak (atlas placement).",
        effect="positive = the RM allele raises expression (RM minus BY allele average of the log2 ratio to the common reference, in that condition; for gene–environment rows, of the ethanol − glucose difference). Computed from Table S4's allele averages; per-row 'up' names the favoured parent.",
        nq=atlas.nlen(T["qtl"]),
        notes=f"1 cross, 109 segregants; {nl} Table S4 linkages ({cnt['glucose']} glucose, {cnt['ethanol']} ethanol, {cnt['gxe']} gene–environment without probe polymorphisms, as in the text; "
              f"{nl - len(clean)} more with probe polymorphisms kept and flagged), loaded per gene on demand. 41 distant-linkage peaks drawn as hotspots. "
              "Candidate genes from Table S1; AMN1, GPA1, HAP1, IRA2, DIG1 and MKT1 as experimentally confirmed (tier x; IRA2, DIG1, MKT1 by allele replacement in this paper), the rest tier c.")
    atlas.save_fragment(DS_ID, T)

    # ---- lazy, verbatim and extra files ---------------------------------------------------------------
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_eqtl.json"), L)
    for c in chroms:
        rows = orig.get(c, {})
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{c}.json"),
                         {"cols": {"s08:S4": s4cols, "s08:S1": s1cols},
                          "rows": {"s08:S4": {k: v for k, v in rows.items() if k.startswith(f"{SC}-e-")},
                                   "s08:S1": {k: v for k, v in rows.items() if k.startswith(f"{SC}-w-")}}})
    extra = {"s08:S1": {"cols": s1cols, "rows": [verbatim(r) for _, r in s1rows] + [[n] + [None] * (len(s1cols) - 1) for n in s1note]},
             "s08:S2": {"cols": s2cols, "rows": [verbatim(r) for _, r in s2rows]},
             "s08:S3": {"cols": s3cols, "rows": [verbatim(r) for _, r in s3rows]},
             "s08:TF": {"cols": tfh, "rows": [r for _, r in tfr]}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)
    data = atlas.assemble()
    print(f"{DS_ID}: {nl} linkages ({dict(cnt)} polymorphism-free); {atlas.nlen(T['qtl'])} peaks; {atlas.nlen(T['causal'])} candidate rows; "
          f"{atlas.nlen(T['dens'])} density bins; data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
