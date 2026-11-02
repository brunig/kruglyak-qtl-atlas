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

## Treusch et al. 2015

Sources: the 14 supporting-information files in `sources/treusch2015/raw/` and the article
page `sources/treusch2015/article.html`. The paper's VCFs, mapping data and R scripts are on
Dryad (doi:10.5061/dryad.j20f0), which is not reachable from the build environment, so no
code could settle the questions below.

### T1. What the grouped QTL's LOD, peak and interval are
S4 gives one row per region per condition, grouped over crosses: "Combined LOD score of
grouped QTL", "Position with highest LOD score across grouped QTL", and a 2-LOD interval
built from "the mean of replicated QTL intervals" (Methods). How LODs were combined across
crosses, mating types and doses is not stated beyond "a conservative estimate"; values such as
41.2375 suggest averages over replicate selections. **Atlas:** stored as given in `qtl.lod`,
`qtl.p`, `qtl.a`/`qtl.e`; the region is one row with `x = -1` and `xs` = the crosses listed.
The per-cross LOD and interval of a grouped QTL are not in the bundle.

### T2. 155 per-cross QTL vs 118 in S4
The paper reports 155 QTL "identified across the twelve crosses and the different
conditions". Summing "# Crosses with QTL" over the 60 regions of S4 gives 118. The 155 count
QTL per cross before doses of the same condition were merged (0.5 M with 1 M NaCl, 15 mM with
20 mM caffeine); those per-cross rows are not published. Not reproducible from the bundle.

### T3. 13 vs 14 regions shared between salt and caffeine
The paper reports 13 regions recurring between conditions. Overlapping 2-LOD intervals in S4
give 14 (in both directions). The paper does not define recurrence; it may require more than
interval overlap.

### T4. S8 bold flags do not mark all multi-condition MAT QTL
The S8 caption says QTL found in multiple conditions are in bold; 11 rows are bold. The text
says 14 of the 19 stress MAT QTL were found at another dose or in the other stress, and
interval overlap within a cross reproduces exactly 14 (the bold rows plus Cross 4 chr12,
Cross 8 chr11 and Cross 10 chr13, which are not bold). **Atlas:** the bold flag is kept in
`qtl.k` and as an extra column of the verbatim row; the 14 is tested by overlap.

### T5. Candidate genes outside their interval
Two of the paper's named candidates fall outside the row's 2-LOD interval (SGD R64-5-1):
CYR1 for Caff.chr10.2 by 63 bp, and GPA1 for Cross 10 / 15 mM caffeine / chr08 by 9.5 kb
(the other five GPA1 MAT rows contain it). Kept as the paper assigned them.

### T6. Candidate-gene tokens that are not single SGD genes
"ENA locus" (NaCl.chr04.1) names the multi-copy ENA region, not one gene, and "STK2"
(NaCl.chr10.1) matches no SGD R64-5-1 standard or systematic name. Both are kept as text in
`causal.name` with `g = -1`; `qtl.g` takes the first listed token that is an SGD gene. "YJL107c &
YJL108c" is split into the two systematic names.

### T7. Evidence tiers for this dataset
Tier `x` = confirmed by allele replacement in this paper ("(confirmed)" in S4: TOR1, WHI2);
tier `c` = candidate gene named by the paper (S4 and S8 "Candidate QTG"). There is no FDR:
QTL were called at LOD 5 (`datasets.fdr` is null).

### T8. Direction of effects
X-QTL give no effect size. For grouped QTL, S7 lists the strains carrying beneficial and
deleterious alleles (`qtl.ups`, `qtl.dns`). At multi-allelic loci a strain can appear in both
lists (ENA: CLIB219 and M22), which is kept as given. For S3 and S8 rows no direction is
reported.

### T9. Strain names and identities
Parents use Bloom 2019's display names where the background is the same (273614X → `273614`);
`parents.s1` keeps the isolate label used in S1 (e.g. `M22_1b`, `273614X_b`). Whether these
are the same isolates as the Bloom 2019 parents (which list different strain IDs and markers)
is not established by either bundle. Ecological origin (`parents.org`) comes from the side
table in S7 ("Subgroup: Source"). BY and RM, used only in the replicate experiments (S3), have
no S1 or S7 entry.

### T10. Positions at 14,900
Many lower bounds (and some peaks) are 14,900 or 15,000, at the left end of chromosomes. This
looks like the first marker after masking the chromosome ends rather than a true interval
bound; it is kept as given.

### T11. MAT QTL traits
S8 LODs compare MATa with MATα pools under one condition, a different quantity from the
selection-vs-YPD LOD of S3/S4. They are stored as separate traits ("MATa vs MATα, <condition>").

### T12. Trait matching across datasets
Grouped QTL combine doses, so their condition key is `caffeine 15-20mM` / `nacl 0.5-1M`
rather than a single dose; S3's `15mM Caffeine` row is keyed `caffeine 15mM` like Bloom 2019.


## Ehrenreich et al. 2012

Files are the PLOS supporting-information files in `sources/ehrenreich2012/raw/`
(`journal.pgen.1002570.sNNN`); "the text" is the article page `sources/ehrenreich2012/article.html`.

### Q1. Genome build of the positions
The arrays were designed from "the BY genome sequence obtained from the Saccharomyces Genome
Database" (Methods); the release is not stated. The S288C reference was revised to R64 in
February 2011, and the paper was received in July 2011, so the design may predate R64. What the
files show: no position in File S10 or in the intensity files (S2-S7) exceeds an R64 chromosome
length; the loci the text places near CAN1, ENT3/RSF2/VPS70 and MKT1 sit at those genes' R64
positions; and every File S11 candidate gene is on the locus chromosome within 15 kb in R64.
**Atlas (most defensible reading):** positions are used as given, with no transform, and
`datasets.build` says the release is unknown. Any offset between the design release and R64 is far
smaller than the 200 kb grouping window and the 50 kb heatmap bins, but single-gene placement near
a locus should not be read more finely than that.

### Q2. What `position` is
File S10 gives one position per distinct locus. The Methods group peaks "starting with the most
strongly selected peak on each chromosome" within 200 kb. **Atlas (most defensible reading):** the
position is that strongest peak's position. Loci are stored as points (`a = e = p`); there is no
interval, LOD or effect size anywhere in the files, and the peak caller (Files S8/S9) never
computes one.

### Q3. Direction
`<cross>_sign` is +1 when the first-named (MATalpha) parent's allele is enriched, -1 for the MATa
parent's, 0 when no peak was detected (always consistent with the TRUE/FALSE columns). This
follows from the Methods (log10 MATalpha minus MATa probe intensity, selection minus controls) and
is confirmed by the data: the BY-resistant singleton pattern `1_1_1_0_0_0` has count 10 in both
File S10 and Table S1, and Copper sulfate in BYxRM has 18 BY-favoured and 2 RM-favoured loci, as in
the text. **Atlas:** `qtl.xs` lists the crosses with a peak and `qtl.ups` / `qtl.dns` hold, element
by element, the favoured and disfavoured parent in that cross. This pairing is a convention of this
dataset that the schema text does not yet state.

### Q4. 837 vs 836 peaks
The text reports 837 peaks; File S10 has 836 detected per-cross peaks. Not reproducible from the
files. The text's mean of 10.7 peaks per experiment holds either way.

### Q5. Pleiotropic 50-kb bins
The text names three bins with 8 or more distinct loci (chrV near CAN1, chrX near ENT3/RSF2/VPS70,
chrXIV near MKT1). Binning File S10 positions from each chromosome start gives chrV 0-50 kb,
chrX 650-700 kb, chrXII 650-700 kb and chrXV 600-650 kb; the MKT1 bin has fewer than 8. The
paper's bin boundaries (and whether it counted distinct loci or traits) are not given.

### Q6. Candidate genes and the two validated genes
File S11 ("candidate genes for bi-allelic loci") covers the 135 allelic-singleton loci only, not the
28 doubletons, with 1-14 genes per locus (mean 6.4), which matches the text's list restricted to
genes with nonsynonymous polymorphisms ("average ... 6 per locus"), not the broader list (mean 10,
range 1-18). **Atlas:** each gene is a `causal` row with tier `c` (candidate named by the paper);
no single gene is assigned to a locus (`qtl.g = -1`). The text validates RED1 (tunicamycin) and
HXT6 (rich medium) by allele replacement, but RED1 (YLR263W) appears in File S11 only for 4-NQO,
benomyl, copper sulfate, cycloheximide and LY83583, and HXT6 not at all; no row is therefore
given a "validated" tier. `causal.eff` (whether the singleton strain's allele increases
resistance) is derived from that locus's File S10 signs.

### Q7. YER109C
One File S11 gene, YER109C, is not in the SGD R64-5-1 gene table the atlas uses; its row is kept
with `causal.g = -1` and no coordinates.

### Q8. Tables S1 and S2 are pictures
Both are images inside Word files. The singleton/doubleton rules of Table S1 reproduce the text's
135 singletons and 28 doubletons from File S10, and its legible class totals (BY 26/16, YJM 7/11,
YPS 9/17, BY=RM 7/5, BY=YPS 7/3) were transcribed by eye for the tests. Table S2 (allelic series)
is not ingested; the 248 loci outside the two simple classes equal the text's 210 uninterpretable
loci plus 29 + 9 allelic series.

### Q9. Doses and condition keys
Each chemical was used at a cross-specific dose (File S1), so a trait has no single dose:
`traits.cond` lists the six doses, and `traits.ck` is the compound alone (`sds` for sodium dodecyl
sulfate, to match Bloom 2019's spelling). The SDS doses in File S1 are bare numbers (0.001,
0.00125), most likely fractions (0.1%); kept as written.

### Q10. Intensity files not in docs/
Files S2-S7 (6 MB each, about 20,000 probes per cross) are checksummed in `MANIFEST.tsv` but not
copied into `docs/`. Five values in them (one or two in each of S3-S6) are written in R's scientific notation (`1e+06`); a
reader must parse them as numbers.


## Breunig et al. 2014

Sources: `sources/breunig2014/raw/` (PLOS Genetics supporting information) and the article page
`sources/breunig2014/article.html`. Statements taken from the article text are transcribed, with
their quotes and sections, in `sources/breunig2014/extracted/text_facts.tsv`.

### Q1. No physical positions: linkages are not placed on R64
Table S1 (`journal.pgen.1004142.s006.pdf`) gives, for each linkage, the closest marker's
chromosome (as a number) and its position in **cM** (column `Position (cM)`). No base-pair
position appears in any table or in the text, and the 2,820-marker genotype map (Brem &
Kruglyak 2005, ref. [32]) is not in the bundle. Converting with the paper's average spacing
(4.3 kb per 1.5 cM) would invent positions.
**Atlas:** `qtl.a`, `e`, `p` are `null`, `qtl.c` is the chromosome and the new column `qtl.cm`
holds the cM position. The dataset's `build` says "none". The genome and region views
cannot draw these rows until a marker map with physical positions is added as a source
(`sources/breunig2014/` would need the BY×RM marker table from Brem & Kruglyak 2005 or
Smith & Kruglyak 2008; blocked here). The table view and the QTL card can show them.

### Q2. Confidence intervals are only in figures
95% Bayes credible intervals (R/qtl `bayesint`) are drawn in Figures S2-S4 (EPS images) for the
chrV, chrVIII and chrXV hot spots. They are not tabulated and are not read off the images.

### Q3. 53 rows in Table S1, not 52
The text and the table caption say "All 52 linkages"; the table has 53 data rows. By compound,
21 have one linkage and 9 have two (as stated), but two have three (ribose-phosphate, valine)
and two have four (hexose-phosphate, Sedoheptulose 7-phosphate), where the text says three and
one. The fourth hexose-phosphate row (chr13, 14.31 cM) is the only one of its rows with
parental p-value `N/A`, which may mean it came from a separate measurement (as the `*` rows for
S-adenosyl homocysteine did), but the table does not say so. **Atlas:** all 53 rows kept.

### Q4. Compounds and footnote marks
Names keep their footnote marks in `traits.raw`. `*` (S-adenosyl homocysteine): "same compound
but in different ionization modes"; its two rows (chrVIII 98.70 and 92.85 cM) are one trait
with two QTL, labelled in `qtl.k`. `†` (glutathione, glutathione disulfide): "considered same
compound"; they stay two traits, as Table S1 lists them under two names, which is how the count
of 34 compounds in the text comes out. `NAD+_posi`, `NADP+_posi` are read as positive-ion-mode
measurements (most defensible reading of the suffix).

### Q5. Parental differences: 28 of 34, not 29
The text says 29 of the 34 compounds with mQTL differed between parents at FDR 5%. With the
caption's cutoff (p = 0.0898) applied to column `Parent's P-value`, 28 compounds pass on any of
their rows; the others are `N/A` (not detected in both parents) or larger.

### Q6. Direction and candidate genes come from the text
Effect sizes are not reported. Direction (`qtl.up`) is set only where the text states it for a
specific linkage: BY allele higher for the five chrXV (ira2) metabolites; RM allele (the
engineered ura3 deletion) higher for orotate and orotidine-5'-phosphate (allele replacement);
RM allele higher for SAM on chrVIII (Figure 5 title). Candidate genes (tier `c`, "candidate
named by the paper") are IRA2 for the five chrXV linkages, URA3 for the five chrV linkages, and
SLT2 and ERC1 for SAM and SAH on chrVIII. Thiamine shares the chrVIII hot spot but the paper
does not tie it to slt2/erc1, so it has no candidate. The candidate rows have no FDR.

### Q7. Hot spots
The text defines 8 hot spots (m1-m8) by binning linkage peaks into 10 kb bins. Without physical
positions the bins cannot be rebuilt, and hot-spot labels are not in Table S1, so they are not
stored.

### Q8. Compounds tested without linkage
74 compounds were tested; only the 34 with a linkage are named in a table. The other 40 are not
in the atlas.

### Q9. Text-derived provenance
Provenance keys for candidate genes, directions and the cross point to rows of
`sources/breunig2014/extracted/text_facts.tsv` (`b14:TF#<row>`), a transcription from the
article page rather than a file of the bundle. `sources.file` for it is
`../extracted/text_facts.tsv` (relative to `raw/`), with its own SHA-256.

## Bloom et al. 2013

Sources: the four supplementary files in `sources/bloom2013/raw/` and the article page
`sources/bloom2013/article.html`. The article is paywalled: the page carries the abstract,
the supplement captions, references and change history, but not the main text or the Online
Methods. The authors' data and code site named in the Supplementary Information
(http://genomics-pubs.princeton.edu/YeastCross_BYxRM/) is not reachable from the build
environment.

### Q1. No methods in the bundle
How QTL were detected (model, threshold, FDR or permutation level) and how intervals and
variance explained were computed are in the paywalled Online Methods. **Atlas:** reports
what Supplementary Table 3 states (1.5-LOD confidence interval, LOD, fraction of variance
explained) and leaves `datasets.fdr` null.

### Q2. The LOD column is not a detection LOD
34 of 591 QTL have LOD < 2.5 and 5 have LOD < 1 (minimum 0.0001), yet the table calls them
detected QTL and Table S1's per-trait QTL counts match its rows exactly. **Atlas (most
defensible reading):** the column is each QTL's LOD in the trait's final multiple-QTL model,
conditional on its other QTL, not the LOD at which it was detected. Stored as given.

### Q3. No effect direction
The caption of Supplementary Table 3 says it lists "effect sizes", but the only effect column
is "Fraction of Phenotypic Variance Explained", which has no sign. Figure S3 plots absolute
standardized differences in allelic means, as an image. **Atlas:** `qtl.ve` holds the
fraction; `qtl.ef` and `qtl.up` are empty, so the QTL card cannot say which parent's allele
increases the trait.

### Q4. Doses stored as fractions with a percent format
Table S1 stores eleven doses as numbers with the Excel format `0%` (0.02 displays as "2%").
**Atlas:** uses the displayed value, so the carbon sources, ethanol and formamide read 2%.
SDS is stored as 0.1 and therefore displays as "10%", which is unusually high for SDS
(Bloom 2019 used 0.075%); the atlas shows "10%" as the table does and flags it here rather
than guessing the intended dose.

### Q5. Segregant count
The bundle does not state the panel size (the abstract says "a large cross"). Table S1 gives,
per trait, the number of segregants with at least one measurement (599 to 1,007).
**Atlas:** `crosses.n` = 1,007, the largest of these, labelled as such in the dataset notes.

### Q6. Genome build
Not stated in the bundle. The coordinates agree with S288C R64: each of the 29,205 SGD
systematic names listed under an interval lies within 5 kb of that interval on R64-5-1
coordinates (28,189 lie entirely inside it); 4,018 listed names (tRNAs, Ty and LTR
elements, centromeres, snoRNAs and similar) are not in the atlas's gene table. No transform.

### Q7. Genes under each interval are not candidate genes
Supplementary Table 3 lists every gene under each 1.5-LOD interval (up to hundreds per QTL).
**Atlas:** keeps them in the verbatim source row and does not assign `qtl.g`; no `causal`
rows for this dataset.

### Q8. Trait categories and condition keys
The bundle has no media table, so categories follow the Bloom 2019 assignment for the same
compounds (carbon sources; sorbitol as osmotic stress) and are labelled "atlas assignment by
analogy with Bloom 2019". Whether the 2% ethanol plates also contained glucose is not stated.
Condition keys follow Bloom 2019's style, with "/ml" written "/mL"; 16 conditions then match
Bloom 2019 exactly (e.g. 6-azauracil 50 µg/mL, zeocin 25 µg/mL, YPD 37 °C, YNB pH 8).
Whether these are the same experiments is not established.

### Q9. Interactions not listed
Table S1 counts 78 pairs of interacting QTL and the additional variance they explain per
trait, but the pairs themselves are not in the bundle. Not ingested as QTL.

### Q10. "Nearly the entire additive contribution"
The abstract's headline claim. From Table S1, the variance explained by the additive QTL
divided by narrow-sense heritability (average of replicates) has median 0.88 across the 46
traits (range 0.67 to 1.21; 38 traits at 0.8 or above). Reported here, not tested: the
abstract gives no number to compare against.

## Ehrenreich et al. 2010

Sources: the Supplementary Information PDF in `sources/ehrenreich2010/raw/` and the article
page `sources/ehrenreich2010/article.html`. The article body is behind Nature's paywall, so
only the abstract, editorial summary and citation metadata are readable there. Statements the
ingest uses are quoted verbatim, with location, in `sources/ehrenreich2010/extracted/text_facts.tsv`.

### Q1. Almost none of the detected loci are tabulated
The abstract reports loci for all 17 chemicals and mitochondrial function, "some traits
influenced by one major locus and others by at least 20 loci", and the Figure S4 legend counts
14 loci for 4-NQO alone. Figures S5A-S5Q draw the loci per chemical at a global FDR of 0.05,
but no table gives their positions, LODs, intervals or directions, and the main text (which may
hold a table) is not readable. Table S2 is the only table with coordinates: the two known 4-NQO
loci, RAD5 and MKT1. **Atlas:** those two loci are ingested; nothing is read off the figures.
The strict xfail `test_14_4nqo_loci` records the gap. Getting the rest needs the main text, or
data the authors deposited elsewhere (none is linked from the files).

### Q2. What the Table S2 coordinates are
Columns: `Gene`, `Chromosome`, `Position (bp)` (a range), `Estimated Peak Location (bp)`,
`Distance from Causal Gene (bp)`. The distance equals the gap from the peak to the nearest end of
the gene range in both rows (1553 = 204992 - 203439; 1304 = 470929 - 469625), so `Position` is the
causal gene's span and `Estimated Peak Location` the X-QTL summit. **Atlas:** each locus is a point
at the peak (`a = e = p`; no interval is reported); the gene span is stored on the `causal` row.

### Q3. Genome build
Not stated. Both gene spans lie inside the R64-5-1 gene features of RAD5 (204970-208580) and MKT1
(466813-469930), and both peaks are inside R64 chromosome lengths. The source spans are narrower
than the R64-5-1 `gene` features (the SGD feature may include more than the ORF), so the exact
release cannot be pinned from two genes. **Atlas (most defensible reading):** positions used as
given, no transform; any offset is far below locus resolution.

### Q4. Direction of effect
Table S2 does not say which parent's allele confers resistance at RAD5 or MKT1, and Figure S4's
allele effects are a plot. **Atlas:** `up` is unset; nothing read off the figure.

### Q5. HAP1 and mitochondrial function
Figure S6 compares BY4716 with "BY4716 with a functional HAP1": the strain carrying a functional
HAP1 stains more with Mitotracker Red (p = 6 x 10^-7). The files give no locus row for
mitochondrial function. **Atlas (most defensible reading):** a `causal` row for HAP1 on the trait
"Mitochondrial function", tier `x` (confirmed by allele replacement, as the Figure S6 test swaps a
nonfunctional for a functional HAP1 in BY), with SGD R64-5-1 gene bounds since the source gives no
coordinate, and no QTL row. The trait gets the new category `cellular phenotype`, since none of the
existing categories fits.

### Q6. Tier for RAD5 and MKT1
Table S2 calls them the "Causal Gene" of "known loci for 4-NQO resistance". **Atlas:** tier `v`
(previously identified), the closest existing tier; the paper's files do not cite where they were
first identified.

### Q7. Doses that look misprinted
Table S1 prints benomyl at "600 mg/ml" and cobalt chloride at "3 M", both implausible as selective
doses (600 ug/ml and 3 mM would be ordinary). **Atlas:** kept exactly as printed, in `traits.cond`,
the display name and `traits.ck`; not corrected.

### Q8. 4-NQO naming
Table S1 lists "4-nitroquinoline"; Figures S4-S5A call it 4-NQO (4-nitroquinoline 1-oxide by the
usual abbreviation, which the files do not spell out). **Atlas:** display name "4-NQO 5 µg/mL",
condition key `4-nqo 5ug/mL`, source name kept in `traits.raw`.

### Q9. Strains and pool size
Only BY4716 is named (Figure S6 legend); the RM isolate is not identified in the readable files.
Pools held "more than 10^5 haploid individuals" (SI page 3; the superscript is lost in PDF text
extraction, read as "105"). **Atlas:** `crosses.n` is null; RM's `code` is just "RM".

### Q10. Effects in the 96 validation segregants
Figure S4 reports, as a plot, the effects of the 14 4-NQO loci in 96 individually genotyped
segregants; the legend says nine were significant and RAD5 explained ~40% of the variance. These
are a different experiment from the X-QTL, so the atlas does not put them on the QTL rows; the
legend text is in `text_facts.tsv` and the full SI text in `docs/v/e10_extra.json`.

## Bloom et al. 2015

Sources: Supplementary Information and Supplementary Data 1-5 in `sources/bloom2015/raw/`, the article
page `sources/bloom2015/article.html`, and two extractions in `sources/bloom2015/extracted/`
(`markers.tsv`, the genotype-matrix column names from Supplementary Data 5, and
`supp_data_5_contents.tsv`). No analysis code is published with the paper ("custom R code").

### Q1. Marker count: 28,220, 28,820 or 28,880
The Methods say 28,220 SNPs were retained; the Introduction says 28,820 variant sites; the
Supplementary Data 5 caption says a 4390x28880 genotype matrix. The matrix in the zip
(`genotypes.RData`, object `gdata`) has 4,390 rows and 28,220 columns. **Atlas:** 28,220, the file.

### Q2. Sign and units of `lm_coefficient`
The paper defines the effect sizes of Fig. 4 by a per-trait multiple regression on all
significant peak markers, and Supplementary Data 5 states genotypes are coded -1 (BY) and 1 (RM).
**Most defensible reading:** `lm_coefficient` > 0 means the RM allele increases growth (larger
colony radius); `qtl.up` follows it. The phenotypes in Supplementary Data 5 are centred near 0
but not scaled to unit variance (raw SD 0.5-1.0 by trait), so the coefficient is in normalized
growth units per genotype step (half the BY-RM allele difference), not SD units. Not checked
against an independent statement (no example direction is given in the text).

### Q3. Which peak, and what `max.lod` is
Supplementary Data 3 has two peaks: `pos` (with `index`, the forward-search peak marker) and
`pos.peak_refineqtl`, plus the 1.5-LOD interval from R/qtl. They differ for 460 QTL, and `pos`
lies outside the interval for 11. **Atlas:** `qtl.p` = the refineqtl peak (inside the interval),
`qtl.lv` = the forward-search marker named from `markers.tsv`. `max.lod` is read as the maximum of
the refineqtl LOD profile used for the interval (most defensible reading); 125 of 797 values are
below 3 and the smallest is about 0.00001, far below what a 5% FWER scan would call significant, so
it is not the forward-search statistic, which is not reported. The browser's LOD stems for this
dataset therefore show this refineqtl value, not the detection strength.

### Q4. 31 vs 25 single-gene QTL
The text says 31 QTL were refined to a single gene (Supplementary Data 3). The `systematic name`
column lists exactly one gene for 25 QTL (13 distinct genes: COQ5, DPS1, HAP1, IMH1, IRA2, PMR1,
RIM15, RPN3, SKY1, SNN1, SWH1, WHI2, ZPS1). **Atlas:** those 25, tier `c`.

### Q5. "172 interactions detected by both approaches"
No matching rule between the full and marginal scans is given. Same trait with both loci within
50 kb matches 170 full-scan rows, within 100 kb 173, identical marker pairs 13, identical closest
additive QTL pair 165. Not reproducible; kept as a strict xfail.

### Q6. Variance-component ranges
Supplementary Data 2 (y=Za+Zi+Zp+e) reproduces the additive range 8.6-70.4% (median 43.3%, mean
43.7%) and the interaction median 9.2%, but gives an interaction range of 2.29-21.26% where the
text says 2.2-21.2% (truncated rather than rounded), and a minimum residual repeatability of
0.56% where the text says 0.05%.

### Q7. Trait doses
The 20 trait names carry no doses, and the article says phenotyping was "as described previously"
(Bloom et al. 2013). Bloom 2013's Supplementary Table 1 (`sources/bloom2013/raw/41586_2013_BFnature11867_MOESM86_ESM.xls`)
lists all 20 traits with doses, e.g. cobalt chloride 1.875 mM, copper 6 mM, diamide 1.5 mM, E6-berbamine
100 uM, ethanol 0.02, formamide 0.02, hydroxyurea 200 mM, indoleacetic acid 100 uM, magnesium chloride
100 mM, manganese sulfate 10 mM, menadione 100 uM, neomycin 5 mg/ml, zeocin 25 ug/ml, and 0.02 for the
carbon sources. Whether the 2015 phenotyping used the same doses is not stated. **Atlas:** `traits.cond`
is empty; `traits.ck` is the compound name only, which links to other datasets by compound.

### Q8. Ethanol's category
"Ethanol" could be ethanol replacing glucose (carbon source) or added to glucose (stress); Bloom 2019
has both. Neither this paper nor Bloom 2013's table (dose 0.02) says which. **Atlas:** `carbon source`,
labelled as an ambiguous atlas assignment.

### Q9. QTL-QTL interactions have no home in the schema
The 205 full-scan and 266 marginal-scan interactions (Supplementary Data 4) are pairs of loci with a
LOD, coefficient and variance explained. The atlas tables have no pair type, so they are kept verbatim in
`docs/v/b15_extra.json` and not drawn. A top-level `pairs` table (dataset, trait, two loci with
positions and marker names, LOD, effect, VE, scan, provenance) would let the browser show them.

### Q10. Genome-wide LOD profiles could be recomputed, but are not
Supplementary Data 5 holds `gdata` (4,390 x 28,220, -1/1, R64 positions in the column names) and
`pheno_raw` (20 traits, two replicates per segregant). These would allow genome-wide single-marker LOD
scans for a future track clearly labelled "recomputed by the atlas", not reported by the paper. The
paper's own scans used BLUP residuals from a mixed model, so a simple scan would not reproduce its
statistics exactly.

## Albert et al. 2018

Sources: the 36 files of `elife-35471-supp-v1.zip` (unpacked in `sources/albert2018/raw/`; the zip's
checksum is in `sources/albert2018/SOURCE_ZIP.tsv`), the article page, the eLife figures page
(`extracted/figures_page.html`, for file captions) and the analysis code
[joshsbloom/eQTL_BYxRM](https://github.com/joshsbloom/eQTL_BYxRM) at commit `ad83786`.

### Q1. Display design: hotspots drawn, eQTL on demand
Chosen by the user. The 102 hotspots are ordinary `qtl` rows (trait "Expression hotspots
(trans-eQTL)", category `expression`); the 36,498 eQTL are in `docs/v/a18_eqtl.json` (table
`lazy`), with per-50-kb counts in `dens` (`local`, `distant`). There is no trait per gene.

### Q2. Sign of r — settled
Genotypes are coded −1 = BY, +1 = RM (code README; Source data 3 caption on the figures page), and
`r` is the Pearson correlation between residual expression and that genotype at the peak marker
(`fasterLOD` in `code/eQTL_BYxRM1000_fx_stranded.R`). **r > 0: the RM allele raises expression.**
`r` and `lm.coeff` agree in sign for 36,494 of 36,498 rows. `crosses.g1` is set to BY so the page's
generic effect orientation matches.

### Q3. Two definitions of "local"
Column `cis` (used for `loc`) is the paper's local eQTL: interval overlaps the gene padded 1 kb
upstream and 200 bp downstream, and, in the code, peak within 40 kb of the gene (2,969 eQTL in 2,884
genes, as the Methods state). Column `gene.to.CI.dist` measures distance to the *unpadded* gene
(`NA` on another chromosome). The sentence "33,529 of 36,498 did not overlap the genes they
influenced" matches 36,498 − 2,969 (cis rule); counting `gene.to.CI.dist` NA or ≥ 1 gives 33,588.
The "5,606 genes with a distant eQTL" matches the unpadded rule (the cis rule gives 5,605). Strict
xfail for the 33,529 under the unpadded rule; the atlas uses `cis`.

### Q4. Hotspot rows carry no LOD or effect
The hotspot overview gives bootstrap peak (median bootstrap position), 95% bootstrap interval and
counts. `lod` and `ef` are null; `qtl.cc` holds `numberEQTLInHotspot` (eQTL whose interval overlaps
the hotspot; an eQTL can count for two overlapping hotspots, so these sum to 55,049). The number of
genes affected (`numberNonzeroEffects`, the paper's "affected N genes") and per-gene effects (Source
data 9) are in the verbatim rows and `docs/v/a18_hotspot_effects.json`.

### Q5. Candidate genes for hotspots
Only genes the text names: MKT1, HAP1, IRA2, GPA1 get tier `x` because the paper says it
"previously identified and experimentally confirmed" them (in earlier papers, not by allele
replacement in this one; the page's `x` label "Confirmed by allele replacement" is slightly
narrower than that). GIS1, STB5, MOT3 (single-gene hotspots), KRE33 ("likely caused by"), ERC1, and
the six transcription factors with damaging variants located in hotspots (GAT1, HMS1, PUT3, RFX1,
SRD1, TBS1) get tier `c`. Each is placed on the hotspot whose interval overlaps its SGD bounds, and
each is also in that hotspot's `allGenesInInterval`. Not assigned: the mating-type locus (confirmed,
but not one gene), AMN1 and SSY1 (cited as known causal genes of yeast eQTL hotspots, but neither
lies in any of the 102 intervals), and GCR2/YNL198C/WHI3/SLZ1 (the text says the chrXIV 267 kb
hotspot "spans" them, without naming a candidate).

### Q6. Genes not in the R64-5-1 table
59 genes with eQTL (371 eQTL; e.g. the mitochondrial Q0140 and Ty-element ORFs such as YAR009C)
are not `gene` features of the R64-5-1 table; their eQTL keep `g = -1` in the lazy file, with the
systematic name in the verbatim row.

### Q7. Markers
The genotype file (Source data 3) holds all 42,052 markers; the 11,530-marker collapsed map used
for mapping is not listed separately. Peak and interval markers in Source data 4 are named
`chr:pos_BY/RM` on sacCer3 and all lie inside R64 chromosomes.

### Q8. Large files kept by member list only
Source data 1 (expression, 97 MB unpacked), 3 (genotypes, 107 MB) and 10 (all GO results,
1.07 million rows) are inner zips kept untouched; nothing from them is copied to `docs/`.
Source data 9 (hotspot effects) and 11 (TFBS enrichment) are kept verbatim in their own lazy files
(`a18_hotspot_effects.json`, `a18_tfbs.json`) to keep `a18_extra.json` small.

## Gou et al. 2019

Sources: the article PDF and the Figshare supplement (doi:10.6084/m9.figshare.7343036), both supplied
by the user, and the authors' notebooks from github.com/gouliangke/Mutation-rate at commit `bcaf987`.
Statements used are quoted in `sources/gou2019/extracted/text_facts.tsv`.

### Q1. No QTL coordinates anywhere
The four mutation-rate QTL are named by chromosome in the Results text (XII, XIV, I, V) with their
direction and variance explained. Their positions, 1.5-LOD intervals and LODs are only drawn: Figure 1B
plots LOD against the genetic map, as do Figure 4 (4NQO, MMS) and Figure S5 (H2O2). The supplement has
no QTL table, and the QTL notebook computes the scan and plots it without printing peaks or intervals.
**Atlas:** QTL rows carry the chromosome only (`a`, `e`, `p`, `cm` null). The page lists them in the
table and card but does not draw them. Peaks could be recomputed from the repository's genotype and
phenotype files, but that would be re-mapping, so it is not done.

### Q2. Which QTL the mutagen traits have
The text says "two large-effect QTL for 4NQO and MMS sensitivity ... overlapped with the major QTL for
mutation rate", and the next sentence refers to "The QTL on chromosome XII and XIV". **Atlas (most
defensible reading):** one 4NQO and one MMS QTL on each of chrXII and chrXIV. For H2O2 the chrXII and
chrXIV QTL "were still observed ... but they had small effects", and its large-effect QTL on other
chromosomes are not named. Whether the chrXII/XIV H2O2 signals pass the 3.64 threshold is not stated,
so H2O2 gets no QTL rows.

### Q3. Candidate tiers
RAD5 (chrXII) and MKT1 (chrXIV) are confirmed for mutation rate by allele-replacement strains (tier
`x`). For 4NQO the paper cites Demogines et al. 2008 for RAD5 and MKT1 as "the underlying causal genes for
4NQO sensitivity", so those rows are tier `c` (named by the paper, from earlier work). MMS gets no
gene: the text does not assign one.

### Q4. Effect direction for the mutagen traits
The Discussion cites Demogines et al. 2008 for the RM RAD5 allele increasing 4NQO and MMS sensitivity.
That is not a measurement in this paper's QTL scan, so `up` is left empty for the mutagen QTL.

### Q5. Engineered-strain mean differences
File S2 reproduces three of the five stated differences at one decimal (I791S − BY 0.9, E783D − BY 0.3,
I791S − E783D 0.6, × 10⁻⁷). BY::RAD5-RM − BY gives 2.84 × 10⁻⁷ against the stated 2.9, and RM::MKT1-BY −
RM gives 6.17 × 10⁻⁷ against the stated 6.1 (strict xfails).

### Q6. RM's mutation rate differs between Table S1 and File S2
Table S1 gives RM11-1a 5.8 × 10⁻⁷ (mean of 96 replicates). File S2's 96 RM replicates average 2.1 × 10⁻⁷.
BY agrees (1.7 vs 1.78 × 10⁻⁷). The two are presumably separate experiments (Table S1 with the seven
isolates, File S2 with the allele replacements, Figure 3), but the paper does not say so. Recorded, not
resolved.

### Q7. Variance explained adds to 20.6%, not 20.7%
8.8 + 6.1 + 3.1 + 2.6 = 20.6, and the text gives 20.7% for the four together. This is consistent with
rounding of the parts and is not marked as a discrepancy.

### Q8. Category
No category is specific to mutation rate. **Atlas:** `cellular phenotype`, as for Ehrenreich 2010's
mitochondrial function. The condition keys for 4NQO (`4nqo 0.1mg/mL`), MMS and H2O2 (`hydrogen peroxide
4mM`) match no other dataset's dose, so they do not link to Bloom 2013 or Ehrenreich 2010's 4NQO and
H2O2 traits.

### Q9. Segregant count
1040 segregants were assayed and 843 were mapped after removing 197 with unreliable cell counts.
`crosses.n` is 843, the mapped number (and File S1's row count).

### Q10. The 4NQO dose may be a unit misprint
The paper prints 4NQO at "0.1 mg/ml" (Materials and Methods and the Figure 4 legend). Bloom 2013 used
4NQO at 0.125 µg/mL and Ehrenreich 2010 at 5 µg/ml, so 0.1 mg/ml would be 800-fold higher than the
earlier BY × RM work. It may be meant as µg/ml. **Atlas:** kept as printed.

## Albert et al. 2014

Sources: the four supplementary files in `sources/albert2014/raw/` (Supplementary Information
PDF, Supplementary Data 1-3) and the Nature article page `sources/albert2014/article.html`. The
main text and Methods are paywalled; semantics come from the abstract, the Extended Data figure
legends and the Supplementary Information, quoted in `sources/albert2014/extracted/text_facts.tsv`.

### Q1. No hotspot table, so no hotspot rows and no causative genes
The abstract says most X-pQTL "were clustered in 'hotspots'", and Extended Data Figure 8 refers
to hotspots "labelled with the names of the causative genes", but the hotspots and their genes are
listed only in Extended Data Table 2 ("Hotspot regulators of protein expression"), which is on the
paywalled article page and not in the bundle. **Atlas:** no rows in `qtl`, no candidate genes; the
hotspots are visible as dense bins in the heatmap rows built from the X-pQTL.

### Q2. Sign of the effect
Settled by the Extended Data Figure 3 legend: "Positive difference values correspond to a higher
frequency of the BY allele in the high GFP population, which we interpret as higher expression
being caused by the BY allele at that locus." `ef` is `alleleFrequencyDifference` as given, so
positive = BY allele raises protein abundance. This is the opposite orientation to Albert 2018's
`r` (positive = RM allele raises expression); the dataset rows say so.

### Q3. Local versus distant
In the genome-wide experiments the tagged gene is fixed for the BY allele (Extended Data Figure 3
legend), so its own locus cannot produce an X-pQTL there; local effects were measured in separate
"local" experiments for 41 genes (Supplementary Table 2, kept verbatim in the extra file, with an
allele frequency difference and a local LOD but no interval). Supplementary Data 2 has no
local/distant column. **Atlas (definition, labelled):** `loc = 1` when the 2-LOD interval overlaps
the tagged gene's coordinates from Supplementary Data 1. That gives 5 of 1,025; only one of them
(YDL126C) is among the five genes whose genome-wide data came from a local experiment (Supplementary
Data 1 comment "distant experiment failed, local experiment used in analyses across all genes").

### Q4. X-pQTL for genes the paper excluded
Supplementary Data 2 has 66 X-pQTL for 9 genes that Supplementary Data 1 marks "failed" or
"exclude due to …" (nonsynonymous GFP mutation, possible nonfixation of the GFP cassette). The
paper's 160 genes exclude them. **Atlas:** all 1,025 rows are kept in the lazy file with the gene's
status in column `st`; the heatmap density counts only the 959 X-pQTL of the 160 analysed genes.

### Q5. What Supplementary Data 1's eQTL/pQTL columns describe
`number_of_trans_pQTL` matches the number of Supplementary Data 2 rows for only 4 of the 91 analysed
genes that have a value, and `local_pQTL_LOD` is missing for 72 of 174 genes. **Most defensible
reading:** these columns (and the eQTL columns) report earlier eQTL and mass-spectrometry pQTL
studies used for comparison, not this study's X-pQTL. Kept verbatim; not used.

### Q6. Supplementary Table 3 effects versus Supplementary Data 2
The 49 proteins affected by the chromosome II hotspot: 48 have a chrII X-pQTL of the same sign in
Supplementary Data 2 (SUR4 has no SGD standard-name match in the R64-5-1 table), but only 24 equal
the nearest chrII peak's allele frequency difference at two decimals. Probably the table's effect
is read at the hotspot position, not at each gene's own peak; not settled without the Methods.

### Q7. Intervals starting at 0
69 of 1,025 intervals have `2LODIntervalLeft = 0`, the chromosome start. Kept as given (the
atlas's coordinates are otherwise 1-based).

### Q8. Genome build
Not stated. All peaks and interval ends lie within R64 chromosome lengths, and the tagged genes'
coordinates in Supplementary Data 1 lie within 1.9 kb of their R64-5-1 genes (162 of 174 inside the
SGD gene span, consistent with coding-sequence rather than gene-feature coordinates). Used as R64,
no transform.

### Q9. Numbers not checkable from the files
Supplementary Note 1 (27 replicated genes, 240 X-pQTL in detection sets, 234 concordant) and
Supplementary Table 1 come from the replicate count data in Supplementary Data 3, which the atlas
does not re-map. Supplementary Note 2's 701 X-pQTL also depend on overlap with published eQTL whose
positions are not in the bundle; the other two filters (not chrII/III, not the gene's own chromosome)
leave at least 701, which is consistent.

### Q10. Pool sizes
No segregant count is given for the sorted pools (the abstract says "thousands of individuals with
high versus low protein abundance"). `crosses.n` is null.

## Torabi and Kruglyak 2012

Sources: the nine supporting-information files in `sources/torabi2012/raw/`, the article page
`sources/torabi2012/article.html`, and `sources/torabi2012/extracted/` (Tables S1-S3 written from
the .doc files, and verbatim quotes, by `extract_sources.py`). No code repository is cited; the raw
arrays are in NCBI GEO (GSE33817), which is not in the bundle.

### Q1. The readthrough-dependent loci are not tabulated
The paper reports 18 readthrough-dependent locus-condition pairs (ten distinct loci; one to six per
condition) and, per cross, one to 20 X-QTL loci per condition. All of these are drawn only (Figure 2,
Figure S2 A-I). **Atlas:** only the locus the text names is ingested: chromosome XIII, right arm,
for diamide and for hydrogen peroxide (2 rows). The text says it affected five conditions; the other
three are visible only in Figure 2B and are not ingested. A strict xfail records the gap.

### Q2. No coordinates for the chromosome XIII locus
Figure S3 draws the 50 kb region around the locus, and the text names SKY1 and MRE11 in it, but no
interval or peak is given. **Atlas:** the rows carry `c = chrXIII` and null `a`, `e`, `p` (listed,
not drawn), as for Gou 2019. They are not placed at SKY1, because the gene is the causal candidate, not
the reported QTL position.

### Q3. Which cross the QTL row belongs to
A readthrough-dependent locus is defined by the difference between the wildtype and sup35 pools. The
skew that defines it is in the sup35 pool ("a frequency skew in favor of the RM allele in the sup35
pool but not in the wildtype pool"). **Atlas (most defensible reading):** `x` = the sup35 cross;
`k` says the locus is the pool difference.

### Q4. Direction on hydrogen peroxide
The text states the RM-favouring skew for diamide only. For hydrogen peroxide the direction comes
from "sup35 C653R and the BY allele of SKY1 act together to lower growth rate on diamide and
hydrogen peroxide" and from Table S2 (SKY1 (BY) estimate −0.0708). **Atlas:** `up = RM` on both rows.

### Q5. Units in Table S1
Table S1 prints chlorpromazine 20 mM, E6-berbamine 40 mM and tunicamycin 3 mM, and CoCl2 3 M, as
X-QTL concentrations. Chlorpromazine, E6-berbamine and tunicamycin are usually used in µM, and 3 M
CoCl2 is implausible on agar (Ehrenreich 2010's Table S1 prints the same "3 M"). **Atlas:** shown as
printed; the condition keys therefore do not link to other datasets' doses.

### Q6. SKY1 tier
SKY1 is confirmed by allele replacement in the X-QTL (the skew disappears when both parents carry the
BY allele) and by growth assays of swapped strains. **Atlas:** tier `x`. MRE11 was tested the same
way and ruled out; it is recorded in the quotes, not as a candidate.

### Q7. Genome build
No coordinates are given, so none applies. SKY1's position is from SGD R64-5-1.

## Foss et al. 2011

Sources: the ten supporting-information files in `sources/foss2011/raw/` and the article page
`sources/foss2011/article.html` (PLoS Biol 9:e1001144). No code or data repository is cited;
the paper says mapping and quantitation code are available from the authors on request.

### Q1. Table S4 gives no p-value, LOD, interval or direction
Each row is one linkage: measurement type (peptide or transcript), regulated gene, peptide,
linked marker, chromosome and position. The text reports linkages at 5% and 1% FDR; the table's
179 + 342 rows equal the 5% FDR counts, so it is the 5% set and the 1% subset cannot be
separated. **Atlas:** every linkage is a point at its marker (`a = e = p`), `lod` and `ef` are
null, and the dataset FDR is 0.05.

### Q2. "49 and 97 loci" are not reproducible
The text says the 5% FDR linkages come from 49 protein and 97 transcript loci. Table S4 has 73
and 156 distinct markers. The rule that merged neighbouring markers into loci is not stated.
Recorded as a strict xfail; the atlas keeps markers as given and does not merge them.

### Q3. Local (cis) linkages are an atlas definition
The paper reports 22 cis transcripts and 3 cis proteins at 1% FDR but does not list them.
**Atlas (labelled rule):** a linkage is local when its marker lies on the regulated gene's
chromosome within 20 kb of the gene (R64-5-1 bounds). That gives 25 transcript linkages
(24 genes) and 3 protein linkages (3 genes) in the 5% set, consistent with (not equal to) the
1% counts. At 10 kb the transcript count (19 genes) would fall below the paper's 22, so 10 kb
is too strict for this marker density.

### Q4. Genome build
Not stated. Markers come from the earlier array genotyping of these segregants (names like
`6964_s_at_x08_2_555787`, ending in chromosome number and position). All positions fall inside
R64 chromosomes, and the hotspot markers sit at the R64 positions of LEU2 (chrIII 90,412-92,013
vs gene 91,179-92,594), AMN1 (chrII 555,787 vs 556,201) and near HAP1 (chrXII 659,357 vs
646,415-650,923). **Atlas:** positions used as given, no transform; agreement is checked at
marker resolution only, not to the base pair.

### Q5. No direction of effect
Table S4 has no sign. The text gives one example (ACS2 protein higher with the RM allele at
chrXII 662,627) in a figure legend; a single example is not ingested as data.

### Q6. Hotspots are not tabulated
Protein and transcript hotspots are described and plotted (Figure 3B, Figure 5) but not listed,
so there are no hotspot rows. The density rows show them.

### Q7. Two measurement kinds in one dataset
Protein linkages go to lazy table `mspqtl` and transcript linkages to `eqtl`. The transcript
data are earlier microarray measurements (references 17, 19, 21 of the paper) re-mapped here
for the same 354 genes; they overlap in kind with Albert 2018's RNA-seq eQTL but are a separate
measurement and are kept separate.

### Q8. Word-format supplements not parsed
Table S2 (community membership) and Text S1 (quantitation algorithm) are legacy Word files; they
are inventoried and checksummed, not parsed.

## Wang and Kruglyak 2014

Sources: the seven supporting-information files in `sources/wang2014/raw/` (fetched from
journals.plos.org), the article page `sources/wang2014/article.html`, the verbatim flattenings
of the Word tables in `sources/wang2014/extracted/S*_Table.tsv`, and the quoted sentences in
`extracted/text_facts.tsv`. The Data Availability statement says all data are in the paper and
its supporting information; no code or data repository is cited.

### Q1. Not X-QTL: linkage in 1008 segregants
The brief anticipated a bulk-segregant (X-QTL) design. The paper maps haloperidol resistance by
linkage in the 1008 individually genotyped BY × RM segregants of Bloom et al. 2013 (R/qtl,
Haley-Knott and non-parametric, LOD 3, 1000 permutations). **Atlas:** `m = "w"`, one cross,
n = 1008.

### Q2. No QTL positions or intervals are tabulated
S1 Table names each locus by chromosome and arm (1L, 1R, 2, 7, 12, 13, 14, 15-1, 15-2) and gives
its drop-one LOD and % variance in the per-dose additive model. The 1.5-LOD intervals and peaks
are drawn in Figure 3A only; main-text Tables 1 and 2 are images. **Atlas:** 32 locus-dose rows
with a chromosome and no coordinates (`a`, `e`, `p` null), listed in the table and card but not
drawn; the arm or locus number is in `qtl.k`. Positions could only come from re-mapping the
Bloom 2013 panel, which the atlas does not do.

### Q3. The direction sentence names a chrV locus that S1 Table does not have
Text: "RM alleles at loci on chromosomes V, XII and XV also confer greater resistance, while BY
alleles confer higher resistance at the remaining five loci." S1 Table has no chromosome V locus
at any dose. **Atlas (most defensible reading):** RM raises resistance at 1R, 12, 15-1 and 15-2
(four loci; "XV" covering both chrXV loci) and BY at 1L, 2, 7, 13 and 14 (five, matching "the
remaining five"); `qtl.up` carries this. "V" may be a slip for the second chrXV locus; the files
cannot settle it. Strict xfail `test_a_chrV_locus_exists`.

### Q4. Which chrXV locus is IRA2's
The text says IRA2 "resides within the locus on chromosome XV" and that the chrXIV and chrXV loci
"became the major QTL at higher doses". Locus 15-1 is present at 80-200 µM with 1-8% variance;
15-2 only at 40-80 µM below 1%. **Atlas (reading):** IRA2 assigned to 15-1.

### Q5. LOD values are drop-one LODs in the additive model
S1 Table's LOD column is the drop-one statistic in the multiple-QTL model at each dose, not the
genome-scan peak LOD used for calling (cutoff 3). Two rows are below 3: 15-2 at 80 µM (1.629) and
13 at 160 µM (2.313). Kept as given; stems in LOD mode would understate scan strength.

### Q6. chrXIII at 160 µM
Text: the chrVII, XII and XIII loci "were only detected at the two lower doses". S1 Table lists
chrXIII at 160 µM as well (drop-one LOD 2.313, 0.69%). Strict xfail
`test_chrXIII_only_at_the_two_lower_doses`.

### Q7. Variance explained by the chrI right-arm locus
Text: 21.4%, 26.7%, 21.4%, 11.6% at 40-160 µM. S1 Table: 21.3477, 26.6095, 21.390, 11.6005, so the
first two round to 21.3 and 26.6. Strict xfail `test_chrI_right_arm_variance_as_stated`.

### Q8. Not checkable from the files
- The additional variance explained by interactions (11.0, 6.6, 1.7, 7.8, 7.2%): S3 Table gives
  per-term variances but no model totals.
- Broad- versus narrow-sense heritability differences (14.6%, 12.9%): Table 1 is an image.
- Signed effect sizes per locus and dose: Figure 3B only.

## Blocked steps

### B1. GitHub Pages
Done: Pages serves branch `claude/main`, folder `/docs`, at https://brunig.github.io/kruglyak-qtl-atlas/.
Work from this session is pushed to `claude/focused-albattani-bvuv83` and reaches the site
when merged into `claude/main`.

### B2. Papers whose supplements could not be fetched
From the build environment (network policy as of 2026-11-01):
- **Nature / Nature Communications**: resolved once `idp.nature.com` was allowed; Bloom 2013,
  Bloom 2015 and Ehrenreich 2010 supplements fetched. Their main texts and Online Methods are
  paywalled, so semantics come from tables, captions and supplements only.
- **eLife** (Albert 2018): resolved; `elifesciences.org` and its CDN now answer plain requests.
- **Genetics / OUP** (Gou 2019): `academic.oup.com` answers 403 and figshare is not allowed;
  resolved by the user supplying the article PDF and the supplement zip.
- **Dryad** (Treusch 2015 code and mapping data): not an allowed host.
Fix: add `idp.nature.com`, `doi.org`, `*.figshare.com`, `datadryad.org` to the environment's
allowed domains, or supply the files directly (upload, or a shared folder).
