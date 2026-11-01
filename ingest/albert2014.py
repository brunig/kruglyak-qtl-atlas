"""Ingest Albert et al. 2014 (Nature 506:494), X-pQTL for 160 GFP-tagged proteins in BY x RM.

    python3 -I ingest/albert2014.py          # from the repository root

Reads only from sources/: the Nature supplementary files in sources/albert2014/raw/, the
transcriptions in sources/albert2014/extracted/ (Supplementary Tables 1-3, quoted facts, the
member list of Supplementary Data 3) and the SGD tables in sources/sgd/. Writes
build/fragments/albert2014.json, docs/v/a14_*.json and, through ingest/atlas.py, docs/data.json.
Deterministic.

Display design (as for Albert 2018): the 1,025 X-pQTL of Supplementary Data 2 go to the lazily
loaded docs/v/a14_pqtl.json (table `lazy`), with per-50-kb counts for the heatmap in table
`dens`. The paper defines hotspots, but the table that lists them (Extended Data Table 2) is not
in the bundle, so there are no hotspot rows in `qtl` and no candidate genes.

Semantics (Extended Data Figure 2-4 legends and the Supplementary Information; the main text is
paywalled): BY x RM segregant pools carrying one GFP-tagged gene were sorted by FACS into high-
and low-GFP tails within a cell-size range and sequenced. `alleleFrequencyDifference` is the BY
allele frequency in the high-GFP pool minus that in the low-GFP pool: positive = the BY allele
raises protein abundance. LOD from MULTIPOOL; genome-wide threshold LOD 4.5; intervals are 2-LOD
drops. In these experiments the tagged gene is fixed for the BY allele, so the genome-wide X-pQTL
are mostly distant; local effects were measured in separate 'local' experiments (Supplementary
Table 2). `loc` here is an atlas definition: 1 when the 2-LOD interval overlaps the tagged gene's
coordinates as given in Supplementary Data 1.
"""
import collections
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "albert2014", "a14"
SRC = os.path.join(atlas.SOURCES, "albert2014")
RAW = os.path.join(SRC, "raw")
P = "41586_2014_BFnature12904_"
HB = atlas.HEAT_BIN
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
EXT = {"a14:ST": "../extracted/supp_tables.tsv", "a14:TF": "../extracted/text_facts.tsv", "a14:SD3": "../extracted/sd3_members.tsv"}


def chrom(c):
    """'chr01' -> 'chrI'."""
    return "chr" + ROMAN[int(str(c)[3:]) - 1]


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def xlsx(name):
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(RAW, P + name), read_only=True, data_only=True)
    ws = wb.worksheets[0]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    width = max(len(r) for r in rows)
    rows = [r + [None] * (width - len(r)) for r in rows]
    return ws.title, ["" if v is None else str(v) for v in rows[0]], [(i + 1, r) for i, r in enumerate(rows) if i >= 1]


def tsv(name):
    with open(os.path.join(SRC, "extracted", name)) as f:
        lines = [l.rstrip("\n") for l in f if not l.startswith("#")]
    head = lines[0].split("\t")
    return head, [(i + 2, l.split("\t")) for i, l in enumerate(lines[1:])]


def main():
    sha = manifest()
    for code, rel in EXT.items():
        sha[rel] = hashlib.sha256(open(os.path.join(RAW, rel), "rb").read()).hexdigest()
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

    s1sheet, s1cols, s1rows = xlsx("MOESM49_ESM.xlsx")
    s2sheet, s2cols, s2rows = xlsx("MOESM50_ESM.xlsx")
    add("sources", code="a14:SD1", ds=0, file=P + "MOESM49_ESM.xlsx", sheet=s1sheet, hrow=1, sha256=sha[P + "MOESM49_ESM.xlsx"])
    add("sources", code="a14:SD2", ds=0, file=P + "MOESM50_ESM.xlsx", sheet=s2sheet, hrow=1, sha256=sha[P + "MOESM50_ESM.xlsx"])
    add("sources", code="a14:ST", ds=0, file=EXT["a14:ST"], sheet="Supplementary Tables 1-3 (from MOESM48_ESM.pdf)", hrow=1, sha256=sha[EXT["a14:ST"]])
    add("sources", code="a14:TF", ds=0, file=EXT["a14:TF"], sheet="text_facts", hrow=1, sha256=sha[EXT["a14:TF"]])
    add("sources", code="a14:SD3", ds=0, file=EXT["a14:SD3"], sheet="Supplementary Data 3 member list", hrow=1, sha256=sha[EXT["a14:SD3"]])

    tfh, tfr = tsv("text_facts.tsv")
    design_row = [rn for rn, r in tfr if r[0] == "design"][0]

    # ---- parents, cross, trait ----------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1="BY (S288C background; GFP-tagged strains from the GFP collection)", s2=None, p18=None, org=None)
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1="RM", s2=None, p18=None, org=None)
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=None,
            dip=None, mm=None, pv=f"a14:TF#{design_row}")
    t = add("traits", id=f"{SC}-protein-abundance-X-pQTL", ds=0, raw="protein abundance (X-pQTL)", name="Protein abundance (X-pQTL)",
            cond=None, bt=None, ck="protein abundance", cat="expression",
            catb="atlas assignment: single-cell GFP abundance of 160 tagged proteins (Albert et al. 2014)")

    # ---- genes studied (Supplementary Data 1) ----------------------------------------------
    h1 = {c: i for i, c in enumerate(s1cols)}
    gene_pos, gene_status = {}, {}
    analysed = set()
    for rn, r in s1rows:
        orf = r[h1["ORF"]]
        lo, hi = sorted((int(r[h1["startPosition"]]), int(r[h1["stopPosition"]])))
        gene_pos[orf] = (chrom(r[h1["chr"]]), lo, hi)
        st = r[h1["localExperiment?"]]
        gene_status[orf] = "analysed" if st in ("yes", "no") else str(st)
        if st in ("yes", "no"):
            analysed.add(orf)

    G = atlas.genes()
    gid = {}
    for i, s in enumerate(G["id"]):
        gid.setdefault(s, i)

    # ---- X-pQTL (Supplementary Data 2): lazy file + density ----------------------------------
    h2 = {c: i for i, c in enumerate(s2cols)}
    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    L = {k: [] for k in ["id", "g", "c", "a", "e", "p", "lod", "ef", "loc", "st", "pv"]}
    dens = collections.Counter()
    for rn, r in s2rows:
        orf = r[h2["gene"]]
        c = chrom(r[h2["chromosome"]])
        a, e, p = int(r[h2["2LODIntervalLeft"]]), int(r[h2["2LODIntervalRight"]]), int(r[h2["peakPosition"]])
        gc, glo, ghi = gene_pos[orf]
        loc = 1 if (gc == c and a <= ghi and e >= glo) else 0
        pid = f"{SC}-p-{rn}"
        for k, v in (("id", pid), ("g", gid.get(orf, -1)), ("c", c), ("a", a), ("e", e), ("p", p), ("lod", atlas.r6(r[h2["LOD"]])),
                     ("ef", atlas.r6(r[h2["alleleFrequencyDifference"]])), ("loc", loc), ("st", gene_status[orf]), ("pv", f"a14:SD2#{rn}")):
            L[k].append(v)
        if orf in analysed:
            for b in range(a // HB, e // HB + 1):
                dens[("local" if loc else "distant", c, b)] += 1
        orig[c]["a14:SD2"][pid] = list(r)
    chroms = list(atlas.chrom_lengths())
    for (k, c, b), n in sorted(dens.items(), key=lambda kv: (kv[0][0], chroms.index(kv[0][1]), kv[0][2])):
        add("dens", ds=0, k=k, c=c, b=b, n=n)
    n_an = sum(1 for g in L["st"] if g == "analysed")
    add("lazy", ds=0, name="pqtl", file=f"v/{SC}_pqtl.json", n=len(L["id"]),
        desc=f"All {len(L['id']):,} X-pQTL of Supplementary Data 2 (LOD ≥ 4.5, 2-LOD intervals): the tagged gene, peak, interval, LOD and the allele frequency difference "
             f"(positive = BY allele raises protein abundance). {n_an:,} belong to the 160 genes the paper analysed; the rest to genes it marked failed or excluded (column st).")

    # ---- dataset row --------------------------------------------------------------------------
    nloc = sum(1 for i in range(len(L["id"])) if L["loc"][i] == 1)
    add("datasets", id=DS_ID, sc=SC, label="Albert 2014",
        cite="Albert FW, Treusch S, Shockley AH, Bloom JS, Kruglyak L. 2014. Genetics of single-cell protein abundance variation in large yeast populations. Nature 506:494-497",
        doi="10.1038/nature12904", year=2014,
        design="X-pQTL (bulk segregant): for each of 160 GFP-tagged genes, a BY × RM segregant pool sorted by FACS into high- and low-GFP tails within a cell-size range, pools sequenced; 'local' experiments for 41 genes in which the tagged gene's own locus segregates",
        build="S288C; not stated. All peaks and intervals fall inside R64 chromosomes and the tagged genes' coordinates in Supplementary Data 1 lie within 1.9 kb of their R64-5-1 genes on the same chromosome. Used as R64, no transform.",
        method="Allele frequencies in the high- and low-GFP pools compared along the genome with MULTIPOOL; X-pQTL called at a genome-wide threshold of LOD 4.5 (main text with the full procedure is paywalled).",
        fdr=None,
        interval="2-LOD drop interval (Supplementary Data 2). An interval may start at position 0 (the chromosome start) as given.",
        effect="alleleFrequencyDifference: BY allele frequency in the high-GFP pool minus the low-GFP pool; positive = the BY allele raises protein abundance (Extended Data Figure 3 legend).",
        nq=atlas.nlen(T["qtl"]),
        notes=f"1 cross, pools of segregants (not counted), 160 genes analysed of 174 attempted; {len(L['id']):,} X-pQTL in Supplementary Data 2 "
              f"({n_an:,} for the 160 analysed genes). Local here = 2-LOD interval overlaps the tagged gene ({nloc} X-pQTL); the paper measured local effects in separate experiments (Supplementary Table 2, kept verbatim). "
              "Hotspots are described in the paper but listed only in Extended Data Table 2, which is not in the bundle, so none are drawn and no causative genes are assigned.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy and verbatim files -----------------------------------------------------------------
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_pqtl.json"), L)
    for c in chroms:
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{c}.json"),
                         {"cols": {"a14:SD2": s2cols}, "rows": {"a14:SD2": dict(orig[c].get("a14:SD2", {}))}})
    sth, str_ = tsv("supp_tables.tsv")
    sdh, sdr = tsv("sd3_members.tsv")
    extra = {"a14:SD1": {"cols": s1cols, "rows": [r for _, r in s1rows]},
             "a14:ST": {"cols": sth, "rows": [r for _, r in str_]},
             "a14:TF": {"cols": tfh, "rows": [r for _, r in tfr]},
             "a14:SD3": {"cols": sdh, "rows": [r for _, r in sdr]}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    size = os.path.getsize(os.path.join(atlas.DOCS, "v", f"{SC}_pqtl.json"))
    print(f"{DS_ID}: {len(L['id'])} X-pQTL ({n_an} for analysed genes; {nloc} local by overlap), "
          f"{len(set(r[h2['gene']] for _, r in s2rows))} genes with X-pQTL, {atlas.nlen(T['dens'])} density bins; "
          f"a14_pqtl.json {size:,} bytes; data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
