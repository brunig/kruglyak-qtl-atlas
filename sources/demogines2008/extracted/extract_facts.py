"""Write text_facts.tsv: statements quoted verbatim from the article page, with their section.

    python3 -I sources/demogines2008/extracted/extract_facts.py      (from the repository root)

The article page (sources/demogines2008/article.html) is reduced to text: inline tags (italics,
links, sub/superscripts) are removed without adding space, block tags become line breaks, and
whitespace is collapsed. Each quote below must occur in that text exactly once, or the script
stops. Tables 1-3 of the article are images on the page and are not transcribed.
"""
import html
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAGE = os.path.join(os.path.dirname(HERE), "article.html")


def article_text(path):
    t = open(path, encoding="utf-8").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    t = re.sub(r"<h[1-6][^>]*>", "\n## ", t)
    t = re.sub(r"</?(p|div|h[1-6]|li|tr|table|br|section)[^>]*>", "\n", t)
    t = re.sub(r"<[^>]+>", "", t)
    t = html.unescape(t)
    t = re.sub(r"[ \t ]+", " ", t)
    return re.sub(r"\n\s*\n+", "\n", t)


# kind, quote
FACTS = [
    ("design", "we examined 123 haploid meiotic progeny of a BY/RM diploid"),
    ("markers", "2956 genetic markers were identified between the two parental strains"),
    ("doses", "onto YPD media containing 0.03% MMS (vol/vol) (Sigma), 0.25 µg/ml or 0.30 µg/ml 4-NQO (Sigma), 10 mM caffeine (Sigma), or 10 µg/ml bleomycin (Sigma)"),
    ("method", "Linkage analysis on 121 segregants with genotype information was performed in R/qtl"),
    ("method", "This corresponded to LOD scores of 3.9, 3.8, 3.8, and 3.7 for MMS-, 4-NQO-, bleomycin-, and caffeine-sensitivity, respectively"),
    ("method", "Support intervals defined by a 1.5 drop in LOD score from the peak were calculated using the lodint function in R/qtl."),
    ("qtl_mms", "MMS-sensitivity, which was seen in the RM parent and half of the hybrid segregants tested, showed strong linkage to a chromosome 12 locus near SMF3"),
    ("qtl_mms", "LOD score of 16.6 for YLR034C_1989"),
    ("qtl_4nqo", "Linkage analysis of 4-NQO-sensitivity identified the same region (LOD score 10.1 for YLR034C_1989, support intervals completely overlap)"),
    ("qtl_bleo", "The bleomycin-sensitive phenotype showed linkage to chromosome 2 (LOD score 5.2 for YBR138C_275"),
    ("qtl_caff", "The caffeine-sensitive phenotype also showed linkage to chromosome 2 (LOD score 5.0 for YBR161W_293"),
    ("qtl_chr2", "the linkages described here are not likely to be due to the same locus because their support intervals do not overlap"),
    ("fig2_legend", "Linkage analysis shows that the MMS- and 4-NQO-sensitive phenotypes observed in hybrid progeny are linked to a locus on chromosome 11."),
    ("interval", "Twenty-two genes lie within the support interval associated with the MMS-sensitive and 4-NQO sensitive phenotypes."),
    ("marker_dist", "A prime candidate is RAD5, located ∼2 kb away from the peak marker."),
    ("marker_dist", "It is located on the edge of the region, approximately 12.5 kb from the peak marker."),
    ("rad5_mms", "RM strains, which are sensitive to MMS, became resistant when the RAD5 gene in RM was replaced with the BY copy."),
    ("rad5_mms", "the MMS-sensitivity observed in RM strains could be created in the BY strain by replacing RAD5 with the RM copy"),
    ("rad5_4nqo", "consistent with the RM copy of RAD5 being associated with 4-NQO sensitivity"),
    ("bsa_design", "Genomic DNA was isolated from pooled resistant and sensitive segregants (40 from each set)."),
    ("bsa_regions", "These regions are on chromosome 7, between 433,546 bp and 551,683 bp (smoothed value of >1), with the peak at 463,019 bp, and on chromosome 14, between 437,935 bp and 473,640 bp (smoothed value of >1.5), with the peak at 467,209 bp"),
    ("bsa_regions", "approximately 118 kb and 36 kb in length for chromosome 7 and 14, respectively"),
    ("bsa_method", "Watson strand 25-mers tiled every 4 bp (BY reference genome)"),
    ("backcross", "At the end of backcrossing, the parental contribution of the final segregants was approximately 94% RM and 6% BY."),
    ("bsa_genes", "For the chromosome 7 linkage region, the peak includes KAP122"),
    ("bsa_genes", "For the chromosome 14 linkage region, the peak lies near MKT1"),
    ("negative", "None of the KAP122 and PDR1 plasmids suppressed the 4-NQO phenotype of these segregants (data not shown)."),
    ("mkt1", "indicating that the MKT1-BY allele contributed to 4-NQO-sensitivity"),
    ("mkt1", "we were unable to do so because the selectable marker used for integration interfered with MKT1 function"),
    ("mkt1", "Because we did not see a significant LOD score for the MKT1 locus in the original linkage analysis"),
    ("additive", "the RAD5-RM, MKT1-BY genotype is the most sensitive, while the RAD5-BY, MKT1-RM genotype is the most resistant"),
]


def main():
    text = article_text(PAGE)
    flat = re.sub(r"\s+", " ", text)
    section, sec_at = "", []
    for m in re.finditer(r"\n## ([^\n]+)", text):
        sec_at.append((len(re.sub(r"\s+", " ", text[:m.start()])), m.group(1).strip()))
    out = ["kind\tsection\tquote"]
    for kind, q in FACTS:
        n = flat.count(q)
        if n != 1:
            sys.exit(f"quote found {n} times: {q[:80]}")
        at = flat.index(q)
        section = ""
        for pos, name in sec_at:
            if pos <= at:
                section = name
        out.append(f"{kind}\t{section}\t{q}")
    with open(os.path.join(HERE, "text_facts.tsv"), "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    print(f"wrote {len(FACTS)} quotes", file=sys.stderr)


if __name__ == "__main__":
    main()
