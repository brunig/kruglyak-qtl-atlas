"""Ingest Breunig et al. 2014 (PLOS Genetics 10:e1004142), metabolite QTL in BY x RM segregants.

    python3 -I ingest/breunig2014.py          # from the repository root

Reads only from sources/: the supporting-information files in sources/breunig2014/raw/
(Table S1 is parsed straight from its PDF), sources/breunig2014/extracted/text_facts.tsv
(per-linkage statements transcribed from the article text, each with its quote), and the
SGD tables in sources/sgd/. Writes build/fragments/breunig2014.json, docs/v/b14_*.json and,
through ingest/atlas.py, docs/data.json. Deterministic.

Genome build: none. Table S1 gives each linkage's closest marker as a chromosome and a
genetic position in cM; no physical position is reported and the marker map is not in the
bundle, so a, e and p are null (OPEN_QUESTIONS Q1). The cM value is kept in qtl.cm.
"""
import collections
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "breunig2014", "b14"
SRC = os.path.join(atlas.SOURCES, "breunig2014")
RAW = os.path.join(SRC, "raw")
S1 = "journal.pgen.1004142.s006.pdf"
S2 = "journal.pgen.1004142.s007.pdf"
S3 = "journal.pgen.1004142.s008.xlsx"
TF = "../extracted/text_facts.tsv"   # path relative to raw/, as sources.file is
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
S1_COLS = ["Name", "PubChem ID", "LOD score", "Chromosome", "Position (cM)", "Parent's P-value"]
S1_ROW = re.compile(r"^(.*\S)\s+(\S+)\s+(\d+\.\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\S+)$")


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def pdf_lines(name):
    import pdfplumber
    with pdfplumber.open(os.path.join(RAW, name)) as pdf:
        return [l for pg in pdf.pages for l in (pg.extract_text() or "").split("\n")]


def norm(name):
    """Compound name without the table's footnote marks (* same compound, two ionization modes; † see Q4)."""
    return re.sub(r"[*†]", "", name).strip()


def display(name):
    """The table's name without footnote marks; "_posi" (positive ion mode) spelled out."""
    return norm(name).replace("_posi", " (positive mode)")


def slug(s):
    return re.sub(r"-+", "-", re.sub(r"[^A-Za-z0-9]+", "-", s)).strip("-")


def main():
    import hashlib
    sha = manifest()
    tf_path = os.path.join(SRC, "extracted", "text_facts.tsv")
    sha[TF] = hashlib.sha256(open(tf_path, "rb").read()).hexdigest()
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

    srcs = {"b14:S1": (S1, "page 1 (text lines)", 1), "b14:S2": (S2, "page 1 (text lines)", None),
            "b14:S3": (S3, "ST3", 4), "b14:TF": (TF, "text_facts", 1)}
    for code, (f, sheet, hrow) in srcs.items():
        add("sources", code=code, ds=0, file=f, sheet=sheet, hrow=hrow, sha256=sha[f])

    # ---- parents and the one cross ----------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY4716", s1="BY4716 (MATa)", s2=None, p18=None, org=None)
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM11-1a", s1="RM11-1a", s2=None, p18=None, org=None)
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        design_row = [i + 1 for i, l in enumerate(f) if l.startswith("design\t")][0]
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=114,
            dip=None, mm=None, pv=f"b14:TF#{design_row}")
    pidx = {"BY": by, "RM": rm}

    # ---- Table S1: one row per linkage --------------------------------------------------
    lines = pdf_lines(S1)
    assert lines[0].startswith("Name PubChem ID LOD score Chromosome Position (cM)"), lines[0]
    s1rows = []
    for i, l in enumerate(lines[1:], start=2):
        m = S1_ROW.match(l)
        if m:
            s1rows.append((i, list(m.groups())))
        elif l.startswith("Sup. Table 1"):
            break
        else:
            raise SystemExit(f"Table S1 line {i} did not parse: {l!r}")
    caption = [l for l in lines if not S1_ROW.match(l)][1:]

    tidx = {}
    for name in sorted({norm(r[0]) for _, r in s1rows}, key=str.lower):
        raws = sorted({r[0] for _, r in s1rows if norm(r[0]) == name})
        tidx[name] = add("traits", id=f"{SC}-{slug(name)}", ds=0, raw=" / ".join(raws), name=display(name), cond=None,
                         bt=None, ck=name.lower(), cat="metabolite",
                         catb="atlas assignment: intracellular metabolite level (LC-MS/MS), Breunig et al. 2014")

    # text facts: direction and candidate genes stated per linkage in the article
    with open(tf_path) as f:
        tfr = [l.rstrip("\n").split("\t") for l in f]
    tf_head, tf = tfr[0], [(i + 2, dict(zip(tfr[0], r))) for i, r in enumerate(tfr[1:])]
    G = atlas.genes()
    gname = {}
    for i, n in enumerate(G["n"]):
        if n:
            gname.setdefault(n, i)

    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    qrows = {}
    for rn, r in s1rows:
        name, pub, lod, chrn, cm, pval = r
        chrom = "chr" + ROMAN[int(chrn) - 1]
        qid = f"{SC}-w-{rn}"
        dirs = [(tr, d) for tr, d in tf if d["kind"] == "direction" and d["metabolite"] == name and d["chromosome"] == chrn]
        up = pidx[dirs[0][1]["allele_up"]] if dirs else -1
        q = add("qtl", id=qid, ds=0, t=tidx[norm(name)], x=x, m="w", c=chrom, a=None, e=None, p=None,
                cm=float(cm), lod=float(lod), q=None, ef=None, ve=None, lv=None, af=None, cc=None, pa=None, anc=None,
                daf=None, up=up, ups=[up] if up >= 0 else None, dns=[rm if up == by else by] if up >= 0 else None,
                xs=None, k="ionization mode not stated (same compound measured twice)" if name.endswith("*") else None,
                g=-1, gf=None, gt=None, pv=f"b14:S1#{rn}")
        qrows.setdefault((name, chrn), []).append(q)
        orig[chrom]["b14:S1"][qid] = r

    # candidate genes named by the paper for specific linkages (tier "c")
    for tr, d in tf:
        if d["kind"] != "candidate":
            continue
        qs = qrows[(d["metabolite"], d["chromosome"])]
        g = gname[d["gene"]]
        for q in qs:
            cid = f"{SC}-c-{tr}-{T['qtl']['id'][q].rsplit('-', 1)[1]}"
            add("causal", id=cid, ds=0, t=T["qtl"]["t"][q], g=g, orf=G["id"][g], name=d["gene"], c=G["c"][g], a=G["a"][g],
                e=G["e"][g], pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier="c", pv=f"b14:TF#{tr}")
            if T["qtl"]["g"][q] < 0:
                T["qtl"]["g"][q] = g; T["qtl"]["gt"][q] = "c"
            orig[G["c"][g]]["b14:TF"][cid] = [d[k] for k in tf_head]
    cz = T["causal"]
    per_gene = C()
    for i in range(atlas.nlen(cz)):
        per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier="c", fdr=None, nt=len({cz["t"][i] for i in ii}), n=len(ii))

    nseg = 114
    add("datasets", id=DS_ID, sc=SC, label="Breunig 2014",
        cite="Breunig JS et al. 2014. Genetic basis of metabolome variation in yeast. PLOS Genetics 10:e1004142",
        doi="10.1371/journal.pgen.1004142", year=2014,
        design="segregant panel; BY4716 × RM11-1a, 114 segregants measured in duplicate; intracellular metabolite levels by targeted LC-MS/MS in exponential growth on glucose minimal medium",
        build="none: Table S1 gives the closest marker's chromosome and genetic position (cM) only; no physical coordinates are reported and the 2,820-marker map is not in the bundle, so these QTL are not placed on R64 (qtl.a, e, p are null; qtl.cm holds the cM position).",
        method="R/qtl single-QTL scan (normal model and nonparametric method) at 2,820 SNP markers; per-metabolite significance from 100 permutations, linkages in the top 10% (empirical FDR 10%) kept, thresholds LOD 3.14-3.58 (mean 3.35). Phenotype: log2(segregant / same-day RM) ion counts, duplicates averaged, OD-normalized.",
        fdr=0.10,
        interval="95% Bayes credible intervals (R/qtl bayesint) are drawn in Figures S2-S4 for chromosomes V, VIII and XV only; no table gives their bounds, so no interval is stored.",
        effect="Not reported per linkage. The text states which parent's allele raises the metabolite for the chrXV ira2 linkages (BY), orotate and orotidine-5'-phosphate at ura3 (RM, allele replacement) and SAM on chrVIII (RM); qtl.up holds those.",
        nq=atlas.nlen(T["qtl"]),
        notes=f"1 cross, {nseg} segregants, {len(tidx)} metabolites with at least one linkage (74 compounds tested). "
              "Candidate genes are those the paper names for a linkage (tier c: IRA2, URA3, SLT2, ERC1); not claimed causal beyond the paper.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy files: verbatim rows ------------------------------------------------------------
    tf_cols = tf_head
    for chrom in atlas.chrom_lengths():
        out = {"cols": {"b14:S1": S1_COLS, "b14:TF": tf_cols}, "rows": {"b14:S1": dict(orig[chrom].get("b14:S1", {})),
                                                                       "b14:TF": dict(orig[chrom].get("b14:TF", {}))}}
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"), out)
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(RAW, S3), read_only=True, data_only=True)
    s3 = [[None if v is None else v for v in r] for r in wb["ST3"].iter_rows(values_only=True)]
    extra = {"b14:S1caption": {"cols": ["line"], "rows": [[l] for l in caption]},
             "b14:S2": {"cols": ["line"], "rows": [[l] for l in pdf_lines(S2)]},
             "b14:S3": {"cols": [str(v) if v is not None else "" for v in s3[3]], "rows": s3}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: {atlas.nlen(T['qtl'])} linkages, {len(tidx)} metabolites, {atlas.nlen(T['causal'])} candidate rows "
          f"({atlas.nlen(T['cand'])} genes); {sum(1 for u in T['qtl']['up'] if u >= 0)} with a stated direction; "
          f"data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
