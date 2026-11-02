"""One-off extraction for Foss et al. 2007 (run from the repository root):
    python3 -I sources/foss2007/extracted/extract.py
Writes, next to this script:
  supp_table1.tsv  Supplementary Table 1 (page 4 of the Supplementary Information PDF), cell for cell
  text_facts.tsv   sentences quoted verbatim from the Supplementary Information PDF and the article page
Needs pdfplumber. The ingest reads only the TSVs.
"""
import hashlib, html, os, re
import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
PDF = os.path.join(SRC, "raw", "41588_2007_BFng200722_MOESM17_ESM.pdf")
ART = os.path.join(SRC, "article.html")
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

with pdfplumber.open(PDF) as pdf:
    pages = [(pg.extract_text() or "") for pg in pdf.pages]

# Supplementary Table 1
lines = pages[3].split("\n")
assert lines[0] == "peptide protein linkage marker lod score marker distance from middle of gene", lines[0]
rows = []
for k, l in enumerate(lines[1:], start=2):
    m = re.fullmatch(r"(\S+) (\S+) (\S+;) ([\d.]+) ([\d.]+)", l)
    if not m:
        break
    rows.append((k, m.groups()))
assert len(rows) == 6
with open(os.path.join(HERE, "supp_table1.tsv"), "w") as f:
    f.write(f"# Supplementary Table 1 of {os.path.basename(PDF)} (sha256 {sha(PDF)}), page 4, extracted with pdfplumber; pdf_line = line on that page\n")
    f.write("pdf_line\tpeptide\tprotein\tlinkage marker\tlod score\tmarker distance from middle of gene\n")
    for k, g in rows:
        f.write(f"{k}\t" + "\t".join(g) + "\n")

# Quotes
flat = lambda t: re.sub(r"\s+", " ", t)
si = flat(" ".join(pages))
art = flat(html.unescape(re.sub(r"<[^>]+>", " ", open(ART, encoding="utf-8", errors="replace").read())))
Q = [
    ("design", "art", "we apply this approach to elucidate the genetic basis of variation in protein abundance in a cross between two diverse strains of yeast."),
    ("finding", "art", "Loci that influenced protein abundance differed from those that influenced transcript levels, emphasizing the importance of direct analysis of the proteome."),
    ("st1_caption", "si", "Supplementary Table 1. Peptides showing local linkage. 6 peptides showed linkage to a marker that was less than 20 kilobases from the gene encoding that peptide."),
    ("snps", "si", "The average number of SNPs across the set of 278 genes encoding proteins that exactly matched peptides used in the linkage analysis was 3.4, vs. 7.4 in the total set of 6,229 genes."),
    ("corr", "si", "The correlations ranged from +0.64 to –0.3, with an average of 0.186"),
    ("parents", "si", "59 proteins differ in abundance between the parents and have transcripts that were different at P < 0.005 1, constituting 43% (59/137) of the proteins that differ between the parents."),
    ("peptide_links", "si", "71 of the 163 peptides had at least one linkage with LOD > 3, and there was a total of 93 such distinct linkages."),
    ("transcript_links", "si", "110 of the 163 transcripts had at least one linkage with LOD > 3, and there was a total of 140 such distinct linkages."),
    ("peptide_overlap", "si", "Of the 93 distinct loci linking to peptide abundance, 18 also linked to the abundance of the corresponding transcript with LOD > 3, 19 linked to the abundance of the corresponding peptide with 1 < LOD < 3, and 56 did not show linkage to the abundance of the corresponding peptide (LOD < 1)."),
    ("transcript_overlap", "si", "Of the 140 distinct loci linking to transcript abundance, 17 also linked to the abundance of the corresponding peptide with LOD > 3, 31 linked to the abundance of the corresponding transcript with 1 < LOD < 3, and 92 did not show linkage to the abundance of the corresponding transcript (LOD < 1)."),
    ("covariate", "si", "Of the 109 loci detected originally, 79 were also detected in this new analysis, along with 17 additional linkages."),
    ("hotspots", "si", "3 out of the 4 hotspots for protein regulation are reflected at the transcript level, something that could not have happened by chance (p < 0.0001)."),
    ("fig4", "art", "Figure 4: Linkage hot spots plotted against genome location."),
]
with open(os.path.join(HERE, "text_facts.tsv"), "w") as f:
    f.write(f"# Verbatim quotes. si = {os.path.basename(PDF)} (sha256 {sha(PDF)}), whitespace collapsed; art = article.html (sha256 {sha(ART)}), tags stripped, whitespace collapsed\n")
    f.write("key\tsource\tquote\n")
    for key, src, q in Q:
        assert q in (si if src == "si" else art), key
        f.write(f"{key}\t{src}\t{q}\n")
print("wrote supp_table1.tsv (6 rows) and text_facts.tsv (%d quotes)" % len(Q))
