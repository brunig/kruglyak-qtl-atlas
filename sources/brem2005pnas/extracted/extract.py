"""One-off: write extracted/text_facts.tsv and extracted/tables.tsv from ../article.html (the PMC
page for PMC547855). Every quote is checked to occur verbatim (whitespace-normalised) in the page.
    python3 -I sources/brem2005pnas/extracted/extract.py
"""
import html, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
t = open(os.path.join(HERE, "..", "article.html"), encoding="utf-8").read()
body = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
txt = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", body)))
FACTS = [
 ("design", "112 segregants, BY4716 x RM11-1a", "Materials and Methods", "We grew 6 independent cultures of BY, 12 of RM, and 1 of each of 112 segregants, isolated RNA, and hybridized cDNA to microarrays as described in ref. 6 ."),
 ("design", "parents BY4716 and RM11-1a", "Materials and Methods", "we used BY4716 (BY for short), isogenic to the lab strain S288C, and the wild isolate RM11-1a (RM for short)"),
 ("design", "6,216 ORFs on the array", "Materials and Methods", "Each array ( 20 ) assayed 6,216 yeast ORFs"),
 ("method", "2,957 markers; Wilcoxon linkage; permutations", "Materials and Methods", "We genotyped segregants at 2,957 markers, performed linkage calculations using the Wilcoxon test, and assessed significance via permutations, as described in ref. 17 ."),
 ("method", "FDR 0.05 = nominal P < 5.7e-5", "Materials and Methods", "the FDR = 0.05 cutoff corresponded to a nominal P < 5.7 × 10 -5 ."),
 ("result", "2,984 transcripts with at least one QTL at FDR 0.05", "Results", "We detected at least one QTL for each of 2,984 transcripts at a FDR (see Materials and Methods ) of 0.05."),
 ("result", "1,038 transcripts with unbiased estimates (average of 10)", "Results", "In 10 realizations of the method, on average, we made estimates for 1,038 transcripts."),
 ("result", "median variance explained 27%", "Results", "The median fraction of variance explained was 27%, with a range of QTL strengths from <10% to near 100%."),
 ("result", "strongest QTL >50% for 23% of mapped transcripts", "Results", "The strongest QTL explained >50% of genetic variance for only 23% of mapped transcripts ( Table 1 ) and <30% of genetic variance for more than half of mapped transcripts ( Table 1 )"),
 ("result", "3,546 highly heritable transcripts, heritability > 69%", "Results", "3,546 transcripts met this criterion, all with heritabilities >69%."),
 ("result", "2,091 (59%) heritable transcripts with a QTL", "Results", "Of these, only 2,091 transcripts (59%) showed linkage to at least one QTL."),
 ("result", "no QTL for 40% of highly heritable transcripts", "Abstract", "no QTLs were detected for 40% of highly heritable transcripts"),
 ("result", "3% single locus; 17-18% one or two loci; half > 5 loci", "Abstract", "only 3% of highly heritable transcripts are consistent with single-locus inheritance, 17–18% are consistent with control by one or two loci, and half require more than five loci under additive models."),
 ("result", "almost a third require n > 8", "Results", "Approximately half of the transcripts require models with n > 5, and almost a third require models with n > 8."),
 ("result", "best geometric model lambda 1.15, alpha 0.45", "Results", "The best fit to the distribution of linkage statistics was obtained for λ = 1.15 and α = 0.45."),
 ("result", "406 (11%) directional", "Results", "Of 3,546 highly heritable transcripts, 406 (11%) passed the directional test"),
 ("result", "2,093 (59%) transgressive", "Results", "Of 3,546 highly heritable transcripts, 2,093 (59%) passed our test at FDR = 0.05"),
 ("result", "583 (16%) epistatic", "Results", "Of 3,546 highly heritable transcripts, 583 (16%) passed the epistasis test"),
 ("data", "expression data in GEO GSE1990", "Footnotes", "The expression data reported in this paper have been deposited in the Gene Expression Omnibus database (accession no. GSE1990)."),
 ("supplement", "Supporting Text, Figs. 5-6 on the PNAS web site", "Results", "Fig. 6, which is published as supporting information on the PNAS web site"),
]
with open(os.path.join(HERE, "text_facts.tsv"), "w") as f:
    f.write("kind\tvalue\tlocation\tsource\tquote\n")
    for k, v, loc, q in FACTS:
        assert q in txt, q[:60]
        f.write(f"{k}\t{v}\t{loc}\t../article.html\t{q}\n")
# tables
out = ["table\trow\t" + "\t".join(f"c{i}" for i in range(1, 6))]
tabs = re.findall(r"<table.*?</table>", t, re.S)
ti = 0
for tb in tabs:
    rows = re.findall(r"<tr.*?</tr>", tb, re.S)
    cells = [[re.sub(r"\s+", " ", html.unescape(re.sub("<[^>]+>", "", c))).strip() for c in re.findall(r"<t[hd][^>]*>(.*?)</t[hd]>", r, re.S)] for r in rows]
    if len(cells) < 3:
        continue
    ti += 1
    for ri, c in enumerate(cells):
        out.append(f"Table {ti}\t{ri}\t" + "\t".join(c + [""] * (5 - len(c))))
open(os.path.join(HERE, "tables.tsv"), "w").write("\n".join(out) + "\n")
print(len(FACTS), "facts;", ti, "tables")
