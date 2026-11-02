"""One-off: write text_facts.tsv, quotes from sources/ronald2005/article.html, each asserted verbatim
against the article text (tags stripped, whitespace collapsed, as below). Run from the repository root:
    python3 -I sources/ronald2005/extracted/extract_facts.py
"""
import html, os, re

HERE = os.path.dirname(os.path.abspath(__file__))


def article_text(path):
    t = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r'<div class="article-text"[^>]*>(.*)', t, re.S)
    t = m.group(1) if m else t
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    t = re.sub(r"<[^>]+>", " ", t)
    return re.sub(r"\s+", " ", html.unescape(t))


FACTS = [
    ("strains", "Results", "BY4716 [BY], isogenic to S288C, and RM11-1a [RM], a wild vineyard strain"),
    ("design", "Results", "using previously reported gene expression and genotype data for 112 segregants from the same cross"),
    ("marker", "Materials and Methods", "only a single marker per transcript, the one closest to the gene's start site, was tested"),
    ("fdr", "Materials and Methods", "A p -value of 0.012 corresponded to a false discovery rate less than 0.05 for the 5,727 transcripts tested"),
    ("n_self", "Materials and Methods", "Thus, of the 1,428 transcripts showing significant self-linkage, 1,357 were expected to be true positives"),
    ("effect_by", "Materials and Methods", "represents the mean expression level of the self-linking gene across all segregants bearing the BY allele at the marker closest to the gene"),
    ("effect_ci", "Materials and Methods", "We obtained confidence intervals for the effect sizes by bootstrap resampling from the 112 segregants"),
    ("ase_77", "Results", "we assayed 77 genes showing self-linkage for the presence of ASE in a diploid hybrid of the two parent strains, BY and RM"),
    ("ase_44", "Results", "Of the 77 assayed genes, 44 (57%) showed ASE at a nominal p -value of less than 0.05"),
    ("ase_opposite", "Results", "In only two of the 44 cases, ASE favored the allele associated by linkage analysis with lower expression"),
    ("ase_controls", "Results", "In this set of 16, we observed only two results with a nominal p -value of less than 0.05"),
    ("ase_78", "Results", "This procedure estimated a true rate of ASE of 78% in the 77 genes tested"),
    ("cis_range", "Abstract", "found that in 52%–78% of these genes, local regulatory variation acts directly in cis"),
    ("amn1", "Results", "Segregants that carry the BY allele of AMN1 show a 2.2-fold increase in its expression relative to segregants that carry the RM allele"),
    ("n_446", "Discussion", "all 446 genes with effect size greater than 1.2"),
    ("aligned", "Results", "These 5,182 genes included 1,233 of the 1,428 genes showing self-linkage"),
    ("chr7", "Results", "three genes with self-linkage and 88 without"),
    ("chr2", "Results", "46 genes with self-linkage and 52 without"),
    ("quarter", "Abstract", "nearly a quarter of all yeast genes contain local regulatory variation between two divergent strains"),
]


def main():
    text = article_text(os.path.join(HERE, "..", "article.html"))
    with open(os.path.join(HERE, "text_facts.tsv"), "w") as f:
        f.write("# quotes from sources/ronald2005/article.html, each verbatim in the tag-stripped, whitespace-collapsed article text\n")
        f.write("key\tsection\tquote\n")
        for k, sec, q in FACTS:
            assert q in text, (k, q)
            f.write(f"{k}\t{sec}\t{q}\n")
    print(len(FACTS), "facts")


if __name__ == "__main__":
    main()
