"""Write text_facts.tsv: sentences quoted verbatim from the article page (sources/smith2008/article.html).

    python3 -I sources/smith2008/extracted/extract_facts.py      # from the repository root

The page text is normalized (tags stripped, HTML entities decoded, whitespace collapsed); each fact is
the sentence containing a fixed anchor. tests/test_smith2008.py re-normalizes the page the same way and
checks every quote is a substring. Deterministic.
"""
import html
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PAGE = os.path.join(HERE, "..", "article.html")

FACTS = [
    ("design", "109 segregant strains from a cross between the parental strains were grown in two conditions"),
    ("linkages_within", "3,997 and 3,489 linkages were observed in glucose and ethanol"),
    ("markers_transcripts", "subjected to linkage analysis with 2,894 markers"),
    ("gxe_total", "There were 1,382 transcripts that showed gene–environment interaction with at least one locus"),
    ("local_gxe", "172/1,555 (11%) of all gxe QTL being local"),
    ("distant_bins", "81% (1,126/1,383) of distant gxe QTL occurring in 31/563 bins"),
    ("peaks", "Significant bins that were located immediately next to each other were merged into a single peak to form 13, 13, and 15 peaks"),
    ("characterized", "Some of these regions coincided with loci that have been characterized at the gene level"),
    ("peaks_shared", "Thirteen peaks were found in each condition, but only seven were found in both."),
    ("gpa1", "In glucose, 49 transcripts show distant linkage here, whereas in ethanol, there are no distant linkages"),
    ("dig1", "Allele replacements of DIG1 in both backgrounds show that variants in DIG1 make a major contribution"),
    ("eth11", "A total of 386 transcripts show distant linkage to this locus in ethanol"),
    ("genomewide_overlap", "only 965/1,555 (62%) gxe QTL also reached genome-wide significance levels"),
    ("cond_specific", "All 1,555 gxe QTL showed an effect in at least one condition"),
    ("local_vs_distant_specific", "546/1,383 (40%) of distant linkages were significant in only one condition"),
    ("ira2_peak", "A total of 372 transcripts link here distantly"),
    ("ira2_replacement", "we generated allele replacement strains that carried the coding region of each parent strain in the background of the other parent"),
    ("ira2_within", "Many transcripts link to this region in glucose and ethanol (1,159 and 410, respectively)"),
    ("segregant_phenotypes", "Three phenotypes were used in linkage analysis"),
    ("markers_final", "Our final set of markers consisted of 2,894 markers"),
    ("thresholds", "LOD scores of 3.25, 3.3, and 3.7 corresponded to FDRs of 4.9%, 4.9%, and 4.8%"),
    ("ci", "We also calculated confidence intervals as LOD support regions"),
    ("local_def", "A linkage was classified as a local linkage if the confidence interval of the linkage overlapped"),
    ("local_region", "The two markers flanking the region from 500 bp upstream of the start of the ORF"),
    ("transcripts_used", "We ultimately used 4,342 transcripts for the parental analysis and 4,482 transcripts for the segregant analysis."),
    ("probe_exclusion", "We thus excluded them from all further analysis."),
    ("probe_reference", "The sequence of each probe was found in the S288c sequence"),
    ("mkt1", "The D30G Variant in MKT1 Is Responsible for Transcript Variation Linking to the MKT1 Region"),
    ("parental_anova", "Of the 4,342 transcripts with high-quality data, 2,037 (47%)"),
]


def page_text():
    t = open(PAGE, encoding="utf-8", errors="replace").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    t = html.unescape(re.sub(r"<[^>]+>", " ", t))
    return re.sub(r"\s+", " ", t)


def sentence(text, anchor):
    i = text.index(anchor)
    a = max(text.rfind(". ", 0, i) + 2, text.rfind(" ## ", 0, i) + 1, 0)
    j = text.find(". ", i + len(anchor) - 1)
    e = j + 1 if j >= 0 else len(text)
    return text[a:e].strip()


def main():
    text = page_text()
    out = ["# Sentences quoted verbatim from sources/smith2008/article.html (normalized: tags stripped, entities decoded, whitespace collapsed).",
           "id\tanchor\tquote"]
    for fid, anchor in FACTS:
        q = sentence(text, anchor)
        assert anchor in q and q in text, fid
        out.append("\t".join((fid, anchor, q.replace("\t", " "))))
    with open(os.path.join(HERE, "text_facts.tsv"), "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
