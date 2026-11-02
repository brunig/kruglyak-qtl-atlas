"""Write text_facts.tsv: statements quoted verbatim from the PMC author manuscript (raw/PMC1409747.html)
and the Nature page (../article.html). Each quote must occur exactly in the normalised page text, or this fails.
    python3 -I sources/brem2005nature/extracted/extract_facts.py
"""
import html, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)

def text(path):
    t = open(path, encoding="utf-8", errors="replace").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", t)))

PMC = text(os.path.join(SRC, "raw", "PMC1409747.html"))
NAT = text(os.path.join(SRC, "article.html"))
# (kind, section, source, quote)
F = [
 ("design", "main text", "PMC", "For this analysis, we used previously described genotype and gene expression data from 112 segregants"),
 ("design", "main text", "PMC", "We sought to identify the loci involved in such interactions in a cross between a lab strain of yeast, BY, and a wild strain, RM."),
 ("estimate", "abstract", "PMC", "Such locus pairs were estimated to play a role in the inheritance of 57% of transcripts; statistically significant pairs were identified for 225 transcripts."),
 ("estimate", "main text", "PMC", "it identified locus pairs for 225 transcripts at a false discovery rate of 5% ( i.e. fewer than 12 of these are expected to have either locus as a false positive)"),
 ("interaction", "main text", "PMC", "the latter term was significant at p < 0.05 for 65% of transcripts"),
 ("interaction", "main text", "PMC", "For comparison, we considered 547 transcripts for which two loci were identified as significant by their individual effects. Of these, only 13% had an interaction term significant at p < 0.05."),
 ("interaction", "main text", "PMC", "the interaction term represented a true positive for 91% of locus pairs mapped by the two-stage search, compared to 13% of pairs mapped independently"),
 ("secondary", "main text", "PMC", "found that 87% of these were significant at p < 0.05"),
 ("secondary", "main text", "PMC", "we found that only 74 (33%) of the 225 were, indicating that 67% of secondary loci would have been missed without the two-stage search"),
 ("bins", "Figure 1 legend", "PMC", "On each axis the genome is divided into 611 bins of 20 kb each, shown in chromosomal order."),
 ("pair", "main text", "PMC", "The largest number of transcripts linking to a single locus pair was 14, with four additional transcripts linking nearby"),
 ("pair", "main text", "PMC", "we eliminated the latter from further analysis, leaving 17 transcripts that linked to this locus pair"),
 ("pair", "main text", "PMC", "The primary QTL of the pair lay near MAT"),
 ("pair", "main text", "PMC", "The secondary locus contained the gene encoding the G-protein subunit Gpa1, in which a single polymorphism in the BY parent strain, S469I, has previously been shown to affect expression of pheromone response genes"),
 ("pair", "main text", "PMC", "Of the 17 transcripts whose linkage to the MAT and GPA1 loci was detected by the two-stage search, 11 were previously shown to be regulated by mating type, of which seven are also regulated by the pheromone response pathway; two have known function that does not involve mating, and four have no known function"),
 ("pair", "main text", "PMC", "For 14 of the transcripts, the MAT and GPA1 loci showed direct evidence for interaction in the regression model at p < 0.05; the independent effect of the GPA1 locus exceeded the genome-wide significance threshold for single-locus linkage for only one of these 17 transcripts."),
 ("engineered", "main text", "PMC", "we engineered isogenic yeast strains carrying each of the four possible combinations of alleles at the two loci"),
 ("engineered", "main text", "PMC", "For 7 of 11 transcripts previously known to be regulated by mating type, and 3 of 4 transcripts of unknown function, the pattern of effect of MAT - GPA1 genotype on expression was the same in engineered strains as in segregants"),
 ("engineered", "main text", "PMC", "YOR090C/PTC1 and YOR162C/YRR1 , which have known functions unrelated to mating, did not show agreement between engineered strains and segregants"),
 ("genotype", "Figure 2 legend", "PMC", "RM is MATa and BY is MATα"),
 ("method", "Methods", "PMC", "For each transcript and marker, a Wilcoxon rank-sum statistic was formed to quantify expression differences between the segregants grouped by inheritance at the locus."),
 ("method", "Methods", "PMC", "With a cut-off corresponding to a 5% false discovery rate, 225 expression traits were called significant for two-locus linkage using this method."),
 ("method", "Methods", "PMC", "the I469 and S469 alleles of GPA1 were each engineered into the S288c derivatives BY4709 (MATα URA3Δ0) and BY4724 (MAT a LYS2Δ0 URA3Δ0)"),
 ("supplement", "Supplementary information (Nature page)", "NAT", "Comparison of transcript levels in segregants and engineered strains for the 18 transcripts with significant linkage to the MAT-GPA1 locus pair under the two-stage linkage test"),
 ("supplement", "Supplementary information (Nature page)", "NAT", "This file contains details regarding an error in the Supplementary Figure at the time of publication."),
 ("abstract", "abstract (Nature page)", "NAT", "Such locus pairs were estimated to be involved in the inheritance of 57% of transcripts; statistically significant pairs were identified for 225 transcripts."),
]
rows = []
for kind, sec, src, q in F:
    page = PMC if src == "PMC" else NAT
    assert page.count(q) >= 1, (src, q[:60])
    rows.append((kind, sec, "raw/PMC1409747.html" if src == "PMC" else "article.html", q))
with open(os.path.join(HERE, "text_facts.tsv"), "w") as f:
    f.write("kind\tsection\tfile\tquote\n")
    for r in rows:
        f.write("\t".join(r) + "\n")
print(len(rows), "facts")
