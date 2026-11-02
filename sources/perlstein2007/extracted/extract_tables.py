"""Flatten Perlstein 2007 Supplementary Tables 1, 3, 4, 5 and 6 (PDFs) into TSVs, line for line.

    python3 -I sources/perlstein2007/extracted/extract_tables.py      (from the repo root)

Each output row carries the source page and line number so provenance points at the PDF text.
Rows that do not match the table's pattern are written to unparsed_lines.tsv rather than dropped.
Needs pdfplumber. Deterministic.
"""
import hashlib, os, re
import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(os.path.dirname(HERE), "raw")
P = "41588_2007_BFng1991_MOESM2{}_ESM.pdf"


def lines_of(n):
    path = os.path.join(RAW, P.format(n))
    out = []
    with pdfplumber.open(path) as pdf:
        for pi, pg in enumerate(pdf.pages, 1):
            for li, l in enumerate((pg.extract_text() or "").split("\n"), 1):
                out.append((pi, li, l.strip()))
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    return out, os.path.basename(path), sha


def write(name, src, sha, head, rows):
    with open(os.path.join(HERE, name), "w") as f:
        f.write(f"# extracted from raw/{src} (sha256 {sha}) by extract_tables.py\n")
        f.write("\t".join(["page", "line"] + head) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


unparsed = []
# Table 3 (all linkages) and Table 4 (collapsed): compound dose µM T=h marker chrom pos lod
RX = re.compile(r"^(?P<cmp>.+?)\s+(?P<dose>[\d.]+)\s*(?P<unit>µM|mM|nM|%|µg/ml|mg/ml|ug/ml)\s+T\s*=\s*(?P<t>[\d.]+)\s*h\s+(?P<mk>\S+)\s+(?P<chr>\d+)\s+(?P<pos>[\d.]+)\s+(?P<lod>[\d.]+)$")
for n, name in (("38", "table_s3.tsv"), ("39", "table_s4.tsv")):
    ls, src, sha = lines_of(n)
    rows = []
    for pi, li, l in ls:
        m = RX.match(l)
        if m:
            rows.append((pi, li, m["cmp"], m["dose"], m["unit"], m["t"], m["mk"], m["chr"], m["pos"], m["lod"]))
        elif l and not l.startswith(("Supplementary Table", "A B C", "Compound response", "Legend", "Column", "Page", "(in micro")):
            unparsed.append((src, pi, li, l))
    write(name, src, sha, ["compound", "dose", "unit", "time_h", "marker", "chrom", "pos", "lod"], rows)
# Table 5 (collapsed, basepair): trait marker chrom pos lod ci_start ci_end ci_length
ls, src, sha = lines_of("40")
rows = []
for pi, li, l in ls:
    t = l.split()
    if len(t) >= 8 and all(re.fullmatch(r"[\d.]+", x) for x in t[-6:]):
        rows.append((pi, li, " ".join(t[:-7]), t[-7], *t[-6:]))
    elif l and not l.startswith(("Supplementary Table", "A B C", "Compound response", "Column", "Page", "Legend")):
        unparsed.append((src, pi, li, l))
write("table_s5.tsv", src, sha, ["trait", "marker", "chrom", "pos_bp", "lod", "ci_start", "ci_end", "ci_length"], rows)
# Table 1: numbered SMP list
ls, src, sha = lines_of("36")
rows = []
for pi, li, l in ls:
    m = re.match(r"^(\d+)\s+(.+)$", l)
    if m:
        rows.append((pi, li, m[1], m[2]))
write("table_s1.tsv", src, sha, ["n", "smp_name"], rows)
# Table 6: hotspots; header line of 8 triples, then compound columns (kept as text lines)
ls, src, sha = lines_of("41")
trip = [l for pi, li, l in ls if re.fullmatch(r"(\d+ \d+ \d+ ?){8}", l)][0].split()
rows = [(1, [li for pi, li, l in ls if l.split() == trip][0], k + 1, trip[3 * k], trip[3 * k + 1], trip[3 * k + 2]) for k in range(8)]
write("table_s6_hotspots.tsv", src, sha, ["column", "chrom", "pos_bp", "n_smp"], rows)
body = [(pi, li, l) for pi, li, l in ls if l and l.split() != trip and not l.startswith(("QTL hot", "Legend", "the second", "the third", "Page"))]
write("table_s6_compound_lines.tsv", src, sha, ["text"], [(pi, li, l) for pi, li, l in body])
# Table 6 compound lists: assign each word to the hotspot column whose header number is nearest in x
with pdfplumber.open(os.path.join(RAW, P.format("41"))) as pdf:
    words = pdf.pages[0].extract_words()
hdr_top = min(w["top"] for w in words if w["text"] == trip[0])
hdr = sorted([w for w in words if abs(w["top"] - hdr_top) < 2], key=lambda w: w["x0"])
col_x = [(hdr[3 * k]["x0"] + hdr[3 * k + 2]["x1"]) / 2 for k in range(8)]
legend_top = min(w["top"] for w in words if w["text"].startswith("Legend"))
cells = []
for w in sorted(words, key=lambda w: (round(w["top"]), w["x0"])):
    if w["top"] <= hdr_top + 2 or w["top"] >= legend_top - 1:
        continue
    xm = (w["x0"] + w["x1"]) / 2
    k = min(range(8), key=lambda c: abs(col_x[c] - xm))
    cells.append((1, round(w["top"], 1), k + 1, w["text"]))
write("table_s6_hotspot_smps.tsv", src, sha, ["column", "smp"], [(pg, top, c, t) for pg, top, c, t in cells])
with open(os.path.join(HERE, "unparsed_lines.tsv"), "w") as f:
    f.write("file\tpage\tline\ttext\n")
    for r in unparsed:
        f.write("\t".join(str(x) for x in r) + "\n")
print("done")
