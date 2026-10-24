# Open questions

Semantics the files and code did not settle, numbers that disagree with the paper, and
steps that were blocked. Each entry gives the reading the atlas uses and labels it.
Code references are to https://github.com/joshsbloom/yeast-16-parents at commit
`c913c9a` unless noted.

## Bloom et al. 2019

### Q1. Ecological origins of the 16 parents are not in the bundle
The brief expected origins in Figure 1-source data 1. The sheet `Crosses and Strains`
gives strain descriptions, the Peter et al. 2018 isolate ID for Parent 1 of each row, and
segregant counts; there is no origin column, and no other file in the bundle has one.
**Atlas:** `parents.org` is `null`; the Peter et al. 2018 ID (`parents.p18`) is shown so a
reader can look the isolate up. Origins could come from Peter et al. 2018 (Nature 556:339)
Supplementary Table 1 in a later step; they are not filled from memory.

### Q2. Which parent is genotype 1 in the within-cross model (the sign of `Beta`)
Within-cross genotypes are coded 1 / 2 by the order of the two parents' genotype columns
in each cross's VCF (`recode.as.parental`, `genotyping/code/segregants_hmm_fx.R`); that
order is not written anywhere in the code or bundle. `QTL_causality.R` (lines 377-383)
labels each shared-parent QTL pair `increase` / `decrease` from `sign(betas)` and whether
the shared parent's genotype column comes first. Inverting that rule on all 427 rows of
the `Causal_Genes` sheet gives, for every one of the 16 crosses, the **alphabetically first
code label** as genotype 1 (e.g. `BYa` in M22×BY, `273614xa` in PW5×273614), which is what a
VCF with samples sorted by name would give. Two rows in cross A point the other way
against 74 that agree. The within-cross `Beta` uses the same coding as the sheet's `betas`
(same sign on all 422 QTL present in both).
**Atlas (most defensible reading, inferred):** `crosses.g1` = alphabetically first parent;
`Beta > 0` means the other parent's allele increases colony radius. The QTL card says
"inferred" next to the favoured parent.

### Q3. Effect units
Settled by the code, recorded here because it is easy to misread. Both models regress the
cross-scaled phenotype on **standardized** genotypes (`scale(g)`), so effects are SD of
phenotype per SD of genotype. At a 1:1 locus the allele-class difference is about twice the
coefficient. The paper's Figure 3 axis ("average effect, SD units") plots the same numbers.

### Q4. Joint-model QTL have no interval, LOD or q-value in the bundle
`joint model (average effect)` gives a lead variant, its 1,011-panel frequency, cross
count, ancestral state and mean effect; the forward search's p and q are not exported.
**Atlas:** joint QTL are points (`a = e = p` = lead variant), `lod` and `q` are `null`, and
the dataset-level FDR (0.05) applies. The browser draws them as ticks and never invents a
width.

### Q5. 105 vs 102 candidate genes at FDR 5%
`Causal_Genes` has 102 unique ORFs with `localFDR < 0.05` (also 102 with `<= 0.05`), and
195 overall. The paper's 105 was not reproduced. The atlas keeps all 427 rows with their
FDR; the tier follows the sheet.

### Q6. 458 joint QTL "not detected within crosses" is not reproducible
Overlapping each joint lead variant with the 1.5-LOD intervals of the same trait, in the
crosses where it segregates (the criterion in `jointModelAF.R`), leaves 1,022 joint QTL
with no within-cross match; with any cross, 526. The code compared against `cross.peaks`,
which holds within-cross QTL up to q < 0.2 (`doTraitFDR(..., FDR_thresh=.2)` in
`mapping.R`), while the published sheet holds q < 0.05 only. Not fixable from the bundle.

### Q7. Supplementary file 3 vs Figure 3-source data 1
Page headers of Supplementary file 3 ("N total QTL | M joint QTL") match the tables for all
38 within-cross totals and 35 of 38 joint totals. EtOH + glucose (111 vs 110), glycerol
(116 vs 115) and lithium chloride (159 vs 158) have one more joint QTL in the figure than in
the table. The atlas uses the table.

### Q8. The 26 experimentally validated genes
The bundle has no list. `QTL_causality.R` builds its "previously identified" list from a
file not in the bundle (Table S1 of NIHMS544073) plus eight genes added by hand; only the
names that also appear among the candidate genes survive in the code (as a pasted R output
and the hand-added line). Those 25 names are in
`sources/bloom2019/code_tables/previously_identified.tsv`. 23 of them reach FDR 5%; HAP1
and RGA1 reach only FDR 20%. "24 of 26" cannot be checked. **Atlas:** tier `v` ("previously
identified, per the analysis code") is given to those 25 names; their FDR is shown beside it.

### Q9. Of the six cAMP-pathway genes, BCY1 is FDR 20%, not FDR 5%
IRA1, IRA2, CYR1, RAS1 and GPR1 have rows with local FDR < 0.05; BCY1's single row has
0.151. All six are in the causal table.

### Q10. The third field of trait names
Trait strings are `compound;condition;n` (`Caffeine;15mM;2`). The code calls `YPD;;2` and
`YPD;;3` "the two additional YPD replicates" and drops them, so the last field indexes an
experiment or batch. **Atlas (most defensible reading):** stored as `traits.bt`, shown as
"batch", not used for anything.

### Q11. Candidate-gene coordinates and QTL intervals
`Causal_Genes` start/end are the gene extended by half of each flanking intergenic gap
(`geneExtend.GR2`), not SGD bounds; the atlas stores them as given and uses SGD R64-5-1
bounds for the gene track. Genes were searched within 25 kb of each peak, which can reach
past the 1.5-LOD interval: 97 of 836 within-cross QTL with an assigned gene have that gene
(SGD bounds) outside the QTL's 1.5-LOD interval. All 836 lie inside the 25 kb search window.
For joint QTL, 152 of 385 resolved lead variants sit inside the assigned gene; all are within
50 kb (the variant and gene each lie within 25 kb of the resolved pair's peaks).

### Q12. Six causal-sheet QTL are not in the within-cross table
`QTL_causality.R` used every QTL in `cross.peaks` (q < 0.2). Six QTL in `Causal_Genes`
have q between 0.054 and 0.113 and so are absent from the published within-cross sheet
(rows 20, 142, 294, 376, 399, 415 of the sheet). Their `causal.q1` / `q2` is `-1`; the card
shows the pair from the sheet's own columns.

### Q13. Allele count among the 16 parents
The bundle has no parental genotypes. `crossCount` (number of crosses in which the lead
variant segregates, rounded up to even) is what the paper uses; `jointModelAF.R` treats
`crossCount / 2` as the number of parents carrying the minor allele, and `mapping_fx.R` calls
variants in fewer than 3 crosses rare ("seg.rare"). `cross.cnt` on per-cross rows is
sometimes odd (impossible in a perfect ring), so an exact carrier set cannot be derived.
**Atlas:** shows `crossCount`, `crossCount / 2` labelled "≈ parents with the minor allele",
and the crosses with an effect estimate.

### Q14. Which parents carry the alternate allele of a joint QTL
Joint effects are signed for the alternate vs reference allele. Without parental genotypes
the atlas cannot say which parents carry the alternate allele; the card says which allele
increases the trait, not which parent.

### Q15. Article title
The bundle does not contain the article title, and Crossref was unreachable from the build
environment (proxy refused `api.crossref.org`). The citation is stored as
"Bloom JS et al. 2019. eLife 8:e49212", DOI 10.7554/eLife.49212, from the brief.

### Q16. Genome build
The code loads `BSgenome.Scerevisiae.UCSC.sacCer3`, and its chromosome lengths
(`segregants_hmm_fx.R` line 74) equal the R64 lengths used with the SGD R64-5-1 annotation
(`microscale-xqtl/browser/scripts/export.R`). The atlas treats all coordinates as R64 and
applies no transform; the gene track is R64-5-1 (annotation updates since R64-1-1 do not
move the reference coordinates; this is an assumption, not checked against SGD here).

### Q17. Cross A segregant count
Figure 1-source data 1 lists 951 segregants for BY×RM (cross A), and the 13,950 total uses
it. `mapping.R` drops plate A11 from cross A for within-cross mapping "so statistical power
is similar to other crosses", so fewer segregants entered that model. The atlas reports 951
as the source does.
