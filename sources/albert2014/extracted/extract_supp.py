"""One-off: transcribe Supplementary Tables 1-3 (SI PDF pages 3, 6, 7-8), quoted facts, and the SD3 zip member list."""
import pdfplumber, re, zipfile, os, sys
R = "sources/albert2014/raw/41586_2014_BFnature12904_MOESM"
OUT = "sources/albert2014/extracted"
pages = pdfplumber.open(R + "48_ESM.pdf").pages
txt = {i + 1: (p.extract_text() or "") for i, p in enumerate(pages)}
rows = [["table", "page", "line", "col1", "col2", "col3", "col4", "text"]]
# ST1: page 3
for ln, l in enumerate(txt[3].split("\n"), 1):
    m = re.match(r"^(\d+(?:\.\d+)?) (\d+) (\d+) \((\d+)%\) (\d+) \((\d+)%\)$", l)
    if m: rows.append(["ST1", "3", str(ln), m.group(1), m.group(2), m.group(3), m.group(5), l])
# ST2: page 6
for ln, l in enumerate(txt[6].split("\n"), 1):
    m = re.match(r"^(Y[A-P][LR]\d{3}[WC](?:-[A-Z])?) (-?\d+\.\d+) (\d+\.\d+)$", l)
    if m: rows.append(["ST2", "6", str(ln), m.group(1), m.group(2), m.group(3), "", l])
# ST3: pages 7-8 (gene, hotspot effect, function text, ribosome flag)
for pg in (7, 8):
    for ln, l in enumerate(txt[pg].split("\n"), 1):
        m = re.match(r"^([A-Z0-9]+) (-?\d\.\d\d) (.*?)( yes)?$", l)
        if m and not l.startswith("WWW"): rows.append(["ST3", str(pg), str(ln), m.group(1), m.group(2), m.group(3), "yes" if m.group(4) else "", l])
with open(os.path.join(OUT, "supp_tables.tsv"), "w") as f:
    f.write("# Supplementary Tables 1-3 transcribed line by line from sources/albert2014/raw/41586_2014_BFnature12904_MOESM48_ESM.pdf (pdfplumber text);\n")
    f.write("# ST1: col1 LOD threshold, col2 X-pQTL in detection set, col3 replicated at genome-wide significance, col4 replicated at relaxed criterion.\n")
    f.write("# ST2: col1 gene, col2 allele frequency difference, col3 local LOD. ST3: col1 gene name, col2 hotspot effect, col3 function, col4 translation/ribosome related.\n")
    f.write("# ST3 function text that wraps onto the previous PDF line is kept only in 'text' of that line's row when it is on the same line.\n")
    for r in rows: f.write("\t".join(r) + "\n")
n = {k: sum(1 for r in rows[1:] if r[0] == k) for k in ("ST1", "ST2", "ST3")}
print(n, file=sys.stderr)
assert n == {"ST1": 5, "ST2": 41, "ST3": 49}, n
# quoted facts, each verified verbatim against the source text
art = open(sys.argv[1]).read()
art1 = re.sub(r"\s+", " ", art)
si1 = re.sub(r"\s+", " ", "\n".join(txt.values()))
facts = [
    ("design", "abstract", "We applied this method to 160 genes and detected many more loci per gene than previous studies."),
    ("design", "abstract", "use pooled sequencing to compare allele frequencies across the genome in thousands of individuals with high versus low protein abundance"),
    ("hotspots", "abstract", "Most loci that we detected were clustered in ‘hotspots’ that influence multiple proteins, and some hotspots were found to influence more than half of the proteins that we examined."),
    ("sign", "Extended Data Figure 3 legend", "Positive difference values correspond to a higher frequency of the BY allele in the high GFP population, which we interpret as higher expression being caused by the BY allele at that locus."),
    ("threshold", "Extended Data Figure 3 legend", "Bottom panel: lod scores obtained from MULTIPOOL 16 . The red horizontal line is the genome-wide significance threshold (lod = 4.5)."),
    ("interval", "Extended Data Figure 3 legend", "The blue vertical boxes indicate positions of genome-wide X-pQTL, with the width representing the 2-lod drop interval."),
    ("fixation", "Extended Data Figure 3 legend", "Note the fixation of the BY allele in all segregants at the gene position and at the mating type locus on chromosome III, as well as the fixation of the RM allele at the synthetic genetic array marker integrated at the CAN1 locus on the left arm of chromosome V."),
    ("LEU1", "Extended Data Figure 3 legend", "For this gene, 14 X-pQTL are called."),
    ("local", "Extended Data Figure 4 legend", "note that YMR315W and GCN1 were ‘local’ experiments where peaks at the gene position are visible."),
    ("FACS", "Extended Data Figure 2 legend", "The low GFP (blue) and high GFP (red) gates sample extreme levels of GFP within a defined range of cell sizes."),
    ("reproducibility", "Supplementary Note 1", "We performed replication experiments for 27 genes."),
    ("reproducibility", "Supplementary Note 1", "Across the 27 detection experiments, we discovered 240 X-pQTL at genome-wide significance."),
    ("reproducibility", "Supplementary Note 1", "The direction of effect was concordant at 234 (97.5%) of the X-pQTL."),
    ("eQTL", "Supplementary Note 2", "The analyses were restricted to the 701 X-pQTL that did not overlap a significant eQTL, were not located on chromosomes II and III, and were not located on the chromosome on which the given gene is located."),
    ("eQTL", "Supplementary Note 2", "We obtained an estimate of π = 36%."),
]
out = [["kind", "where", "quote", "found_in"]]
for k, w, q in facts:
    qq = re.sub(r"\s+", " ", q)
    where = "article.html" if qq in art1 else "MOESM48_ESM.pdf" if qq in si1 else None
    if where is None and k == "eQTL" and "π" in q:   # pdf splits the subscript: 'π = 36%' appears as 'π = 36%' with the 1 on the next line
        where = "MOESM48_ESM.pdf" if "We obtained an estimate of π = 36%." in si1 else None
    assert where, q
    out.append([k, w, q, where])
with open(os.path.join(OUT, "text_facts.tsv"), "w") as f:
    for r in out: f.write("\t".join(r) + "\n")
# SD3 zip member list (not extracted)
z = zipfile.ZipFile(R + "51_ESM.zip")
with open(os.path.join(OUT, "sd3_members.tsv"), "w") as f:
    f.write("member\tbytes_uncompressed\tbytes_compressed\n")
    for i in z.infolist(): f.write(f"{i.filename}\t{i.file_size}\t{i.compress_size}\n")
print("ok", len(out) - 1, "facts", file=sys.stderr)
