"""Quote, verbatim, the sentences of the article that the ingest relies on.

    python3 -I sources/wang2014/extracted/extract_facts.py      (from the repository root)

Reads sources/wang2014/article.html, reduces it to plain text (tags removed, whitespace
collapsed), splits it into sentences, and writes text_facts.tsv with one row per fact:
kind, the sentence quoted whole, and the section heading it sits under. A fact is found by a
substring that must occur in exactly one sentence. tests/test_wang2014.py checks each quote
occurs in the article text.
"""
import html
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ART = os.path.join(os.path.dirname(HERE), "article.html")
FACTS = [
    ("design", "we used a large panel of 1008 segregants from a cross between a laboratory strain BY4724"),
    ("design", "for growth at five different concentrations of haloperidol (40, 80, 120, 160, 200 µM) and identified nine distinct significant QTL"),
    ("direction", "At the major locus on the right arm of chromosome I, the allele from RM"),
    ("direction", "RM alleles at loci on chromosomes V, XII and XV also confer greater resistance"),
    ("interval", "Significant loci are shown through their 1.5-LOD drop confidence intervals"),
    ("effect", "Positive values represent better growth in the presence of the BY allele"),
    ("dose", "the loci on chromosomes VII, XII, and XIII were only detected at the two lower doses"),
    ("ve", "This model explained between 35.9% and 54.7% of the total phenotypic variance"),
    ("ve", "The locus on chromosome I (right arm) alone explained 21.4%, 26.7%, 21.4%, 11.6%"),
    ("candidate", "The confidence interval of this peak contained the gene SWH1"),
    ("candidate", "Thus, SWH1 RM is the resistant allele relative to SWH1 BY"),
    ("candidate", "Introducing the functional BY allele into the RM background slightly reduced resistance"),
    ("candidate", "resides within the locus on chromosome XV"),
    ("candidate", "The IRA2 RM allele in both BY and RM genetic backgrounds was more resistant to the drug"),
    ("candidate", "we confirmed MKT1 as the gene underlying the chromosome XIV locus"),
    ("candidate", "The observation that the BY allele of MKT1 led to greater growth in haloperidol"),
    ("interaction", "found 9, 10, 3, 5 and 5 significant pair-wise QTL interactions"),
    ("interaction", "explained an additional 11.0%, 6.6%, 1.7%, 7.8%, and 7.2% of phenotypic variance"),
    ("method", "Linkage was determined by calculating LOD scores for each genotype marker"),
    ("method", "QTL were called at a LOD cutoff of 3"),
    ("method", "Residuals were used for QTL mapping"),
    ("data", "All relevant data are within the paper and its Supporting Information files"),
]


def text():
    t = open(ART, encoding="utf-8", errors="replace").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    parts = re.split(r"(<h[1-6][^>]*>.*?</h[1-6]>)", t, flags=re.S)
    out, head = [], ""
    for p in parts:
        if re.match(r"<h[1-6]", p):
            head = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", p))).strip()
            continue
        s = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", p))).strip()
        for sent in re.split(r"(?<=[.!?])\s+(?=[A-Z(])", s):
            out.append((head, sent.strip()))
    return out


def main():
    sents = text()
    rows = []
    for kind, sub in FACTS:
        hits = [(h, s) for h, s in sents if sub in s]
        assert len(hits) >= 1, sub
        h, s = hits[0]
        rows.append((kind, s, h))
    with open(os.path.join(HERE, "text_facts.tsv"), "w", encoding="utf-8") as f:
        f.write("kind\tquote\tsection\n")
        for r in rows:
            f.write("\t".join(x.replace("\t", " ") for x in r) + "\n")


if __name__ == "__main__":
    main()
