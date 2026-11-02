"""One-off extraction for Gatbonton et al. 2006 (run from the repository root):
    python3 -I sources/gatbonton2006/extracted/extract.py
Writes, next to this file:
  text_facts.tsv  statements quoted verbatim from the article page (article.html) or the correction
                  page (correction.html); each quote is asserted to occur in the whitespace-normalised text.
  table_s1.tsv    Table S1 (journal.pgen.0020035.st001.xls, Sheet1) cell for cell, with its spreadsheet row.
Needs xlrd. Deterministic.
"""
import hashlib, html, os, re
import xlrd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)


def page_text(name):
    t = open(os.path.join(SRC, name), encoding="utf-8", errors="replace").read()
    m = re.search(r'<div class="article-text"[^>]*>(.*)', t, re.S)
    t = m.group(1) if m else t
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    t = html.unescape(re.sub(r"<[^>]+>", " ", t))
    return re.sub(r"\s+", " ", t)


ART, COR = page_text("article.html"), page_text("correction.html")
FACTS = [  # kind, source page, quote
    ("design", "article.html", "analyzed telomere length in 122 progeny"),
    ("design", "article.html", "The 122 segregants and their parent strains (S288C derivative BY4716 MATa lys2 ) and natural isolate derived from Bb32(2) (RM11-1a MATa leu2, ura3, HO::KanMx ) were described"),
    ("design", "article.html", "linkage analysis was performed with 3,312 genetic markers that were previously genotyped in all of the segregants using oligonucleotide microarrays"),
    ("method", "article.html", "Linkage analysis of telomere length in the BY and RM11 cross was carried out using publicly available R/qtl software"),
    ("method", "article.html", "To evaluate the significance of the mapping results, we carried out linkage analysis of 100 datasets with permutations of the segregant identities."),
    ("qtl", "article.html", "with significant linkage peaks occurring on Chromosome 12 and Chromosome 13 ( Table 2 )"),
    ("lod", "article.html", "The LOD scores of 7.4 and 4.3 at Chromosome 12 and Chromosome 13 loci were significant (genome-wide corrected p < 0.01), as we did not observe any scores higher than four in the analysis of 100 permuted datasets."),
    ("ve", "article.html", "The Chromosome 12 and Chromosome 13 loci explained only 25% and 13% of the telomere length-variation among the segregants, respectively"),
    ("direction", "article.html", "the Chromosome 12 RM11 allele conferred longer telomeres and the BY allele conferred shorter telomeres, even though the RM11 parent has shorter telomeres than the BY parent"),
    ("abstract", "article.html", "Genome-wide linkage analysis among these segregants revealed two loci that account for 30%–35% of telomere length-variation between the strains."),
    ("candidate", "article.html", "Several genes known to regulate telomere length reside in the mapped regions including EST1, VPS34, RPS28B, ARV1, and VPS9 [ 20 ]."),
    ("candidate", "article.html", "Polymorphisms that result in amino-acid substitutions in these genes ( Table S1 ) are likely candidates that control telomere length in the segregants."),
    ("phenotype", "article.html", "The average telomere length in the parents, obtained through three independent measurements, was 350 ± 7 and 256 ± 5 for BY and RM11, respectively."),
    ("phenotype", "article.html", "Average telomere length of the progeny is continuously distributed between 222 and 342 bp"),
    ("phenotype", "article.html", "a fraction (27%) have telomere length shorter than RM11"),
    ("survey", "article.html", "We found 72 genes that when deleted confer short telomeres, and 80 genes that confer long telomeres relative to those of wild-type yeast."),
    ("survey", "article.html", "Among identified genes, 88 have not been previously implicated in telomere length control."),
    ("survey", "article.html", "We also found 39 genes that have an important role in telomere maintenance or cell proliferation in the absence of telomerase"),
    ("survey", "article.html", "we investigated telomere lengths in 13 wild yeast strains"),
    ("correction", "correction.html", "Table 1 was incomplete. The full table is printed below."),
]


def main():
    out = ["kind\tsource\tquote"]
    for kind, src, q in FACTS:
        text = ART if src == "article.html" else COR
        assert q in text, (kind, q[:60])
        out.append(f"{kind}\t{src}\t{q}")
    open(os.path.join(HERE, "text_facts.tsv"), "w").write("\n".join(out) + "\n")
    xls = os.path.join(SRC, "raw", "journal.pgen.0020035.st001.xls")
    sh = xlrd.open_workbook(xls).sheet_by_name("Sheet1")
    rows = [f"# Table S1 of Gatbonton et al. 2006, from raw/journal.pgen.0020035.st001.xls (sha256 {hashlib.sha256(open(xls, 'rb').read()).hexdigest()}), Sheet1, cell for cell",
            "row\tcol_B\tcol_C\tcol_D\tcol_E"]
    for i in range(sh.nrows):
        v = sh.row_values(i)
        if not any(str(x).strip() for x in v):
            continue
        rows.append("\t".join([str(i + 1)] + [str(x).replace("\t", " ") for x in v[1:5]]))
    open(os.path.join(HERE, "table_s1.tsv"), "w").write("\n".join(rows) + "\n")
    print(len(FACTS), "facts;", len(rows) - 2, "table rows")


if __name__ == "__main__":
    main()
