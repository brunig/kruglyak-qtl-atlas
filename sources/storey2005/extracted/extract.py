"""One-off extraction for Storey, Akey & Kruglyak 2005 (PLoS Biol 3:e267).

    python3 -I sources/storey2005/extracted/extract.py      # from the repository root; needs pdfplumber

Writes, next to this script:
  s1_table.tsv   Table S1 (journal.pbio.0030267.st001.pdf), one row per trait, every printed column,
                 plus the PDF page and text line it came from
  s2_table.tsv   Table S2 (journal.pbio.0030267.st002.pdf), same treatment
  text_facts.tsv quoted sentences from the article page, each checked to occur verbatim
Each TSV starts with a comment line naming its source file and that file's SHA-256, so the ingest
can refuse a mismatch. The ingest reads these TSVs, not the PDFs, so the build needs no PDF library.
"""
import hashlib
import html
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
RAW = os.path.join(SRC, "raw")


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def pdf_lines(name):
    import pdfplumber
    out = []
    with pdfplumber.open(os.path.join(RAW, name)) as p:
        for i, pg in enumerate(p.pages):
            for j, line in enumerate((pg.extract_text() or "").split("\n")):
                out.append((i + 1, j + 1, line.strip()))
    return out


def write(name, src, head, rows):
    with open(os.path.join(HERE, name), "w") as f:
        f.write(f"# source: raw/{src} sha256={sha(os.path.join(RAW, src))}\n")
        f.write("\t".join(head) + "\n")
        for r in rows:
            f.write("\t".join(str(v) for v in r) + "\n")


S1 = re.compile(r"^(\S+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) ([\d.]+) ([012])$")
S2 = re.compile(r"^(\d+) (\S+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) ([\d.]+) ([012])(?: (.*))?$")

# Quotes from the article page; each must occur verbatim in its whitespace-normalized text.
QUOTES = [
    ("design", "Gene expression measurements were obtained for 6,216 open reading frames in 112 haploid segregants, and dichotomous genotypes were identified at 3,312 markers covering 99% of the genome"),
    ("design", "112 F1 segregants (one from each tetrad) were grown from a cross involving parental strains BY4716, isogenic to the lab strain S288C, and the wild isolate RM11-1a"),
    ("result", "Pairs of jointly linking quantitative trait loci are identified with high confidence for 170 gene expression traits, where it is expected that both loci are true positives for at least 153 traits"),
    ("result", "Of these 2,300 traits, 170 can be identified at a FDR of 10%"),
    ("result", "we estimate that 2,300 traits (approximately 37%) are jointly linked to two loci"),
    ("result", "Therefore, we expect that at most 17 of these joint linkages include a single false-positive locus, and about zero include two false-positive loci"),
    ("result", "In total, 58 traits demonstrate a cis linkage"),
    ("result", "we are able to show that epistatic interactions contribute to gene expression variation for at least 14% of all traits"),
    ("result", "At FDR cut-offs of 1% and 5%, there are 2,780 and 4,271 significant pair-wise linkages, respectively"),
    ("result", "At these same cut-offs, the 2D search yielded substantially fewer: 1,715 and 3,540 significant linkages"),
    ("cluster", "We observed 10 pair-wise bins with three or more traits"),
    ("cluster", "the 12 traits that jointly link to nearly identical positions on Chromosomes 3 and 8 are predominantly involved in the mating response"),
    ("cluster", "The linkage peak on Chromosome 3 maps to the precise location of the yeast mating type locus MAT"),
    ("cluster", "amino acid and mitochondria metabolism (13 traits defined by linkage to regions on Chromosomes 3 and 13), mitochondrial tricarboxylic acid cycle (ten traits defined by linkage to regions on opposite ends of Chromosome 15), and response to stress (five traits defined by linkage to Chromosomes 6 and 10)"),
    ("cluster", "of the 10 pair-wise bins with three or more linked traits, there are three that share a common linkage to the exact same position on Chromosome 15"),
    ("method", "The significance of each locus pair selected for linkage was computed using a standard permutation technique against the null hypothesis of no linkage to either locus"),
]


def article_text():
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    t = html.unescape(re.sub(r"<[^>]+>", " ", t))
    return " ".join(t.split())


def main():
    rows = []
    for pg, ln, line in pdf_lines("journal.pbio.0030267.st001.pdf"):
        m = S1.match(line)
        if m:
            rows.append(list(m.groups()) + [pg, ln])
    assert len(rows) == 170, len(rows)
    write("s1_table.tsv", "journal.pbio.0030267.st001.pdf",
          ["Trait", "Chr", "Start", "Stop", "Locus 1 Chr", "Locus 1 Pos", "Locus 2 Chr", "Locus 2 Pos", "q-value Joint", "Cis", "pdf_page", "pdf_line"], rows)
    rows = []
    for pg, ln, line in pdf_lines("journal.pbio.0030267.st002.pdf"):
        m = S2.match(line)
        if m:
            rows.append([g if g is not None else "" for g in m.groups()] + [pg, ln])
    assert len(rows) == 61, len(rows)
    write("s2_table.tsv", "journal.pbio.0030267.st002.pdf",
          ["Group", "Trait", "Chr", "Pos_start", "Pos_end", "Locus 1 Chr", "Locus 1 Pos", "Locus 2 Chr", "Locus 2 Pos", "Joint FDR", "Cis", "Putative Group Function", "pdf_page", "pdf_line"], rows)
    text = article_text()
    for _, q in QUOTES:
        assert " ".join(q.split()) in text, q[:60]
    with open(os.path.join(HERE, "text_facts.tsv"), "w") as f:
        f.write(f"# source: ../article.html sha256={sha(os.path.join(SRC, 'article.html'))}; each quote occurs verbatim in the page text\n")
        f.write("kind\tquote\n")
        for k, q in QUOTES:
            f.write(f"{k}\t{q}\n")
    print("wrote s1_table.tsv (170), s2_table.tsv (61), text_facts.tsv", len(QUOTES), file=sys.stderr)


if __name__ == "__main__":
    main()
