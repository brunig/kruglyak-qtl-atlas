# Datasets

## Ingested

| id | paper | design | build | QTL | traits | crosses / segregants | candidate genes | ingest |
|---|---|---|---|---|---|---|---|---|
| `bloom2019` (`b19`) | Bloom JS et al. 2019. eLife 8:e49212. doi:[10.7554/eLife.49212](https://doi.org/10.7554/eLife.49212) | 16-parent round-robin segregant panel; within-cross and joint (all crosses) models | S288C R64 (sacCer3), no transform | 4,552 joint (FDR 5%, lead variant only) + 7,751 within-cross (q < 0.05, 1.5-LOD intervals); 18,758 per-cross joint effects | 38 | 16 crosses of 16 parents / 13,950 | 427 gene-to-QTL-pair rows, 195 genes (102 at local FDR 5%) | `ingest/bloom2019.py` |
| `treusch2015` (`t15`) | Treusch S, Albert FW, Bloom JS, Kotenko IE, Kruglyak L. 2015. PLoS Genet 11:e1004913. doi:[10.1371/journal.pgen.1004913](https://doi.org/10.1371/journal.pgen.1004913) | X-QTL (bulk segregant) in a 12-strain round-robin cross; MATa/MATα pools; BY × RM replicate selections | S288C R64 (sacCer3), no transform | 97 X-QTL: 60 regions grouped across crosses (LOD > 5, 2-LOD intervals), 12 replicate-confirmed BY/RM QTL, 25 MAT-dependent QTL | 10 | 12 round-robin crosses + BY × RM / pools (not counted) | 60 rows, 37 genes (TOR1, WHI2 confirmed by allele replacement) | `ingest/treusch2015.py` |
| `ehrenreich2012` (`e12`) | Ehrenreich IM, Bloom J, Torabi N, Wang X, Jia Y, Kruglyak L. 2012. PLoS Genet 8:e1002570. doi:[10.1371/journal.pgen.1002570](https://doi.org/10.1371/journal.pgen.1002570) | X-QTL in all 6 pairwise crosses of BY, RM, YJM789, YPS163; 78 selections | S288C; array SGD release not stated, consistent with R64 at locus resolution; no transform | 411 loci (1% FDR; points, favoured parent by cross; 836 per-cross peaks) | 13 chemicals (doses differ by cross) | 6 crosses of 4 parents / pools (not counted) | 864 rows, 680 genes, candidates for the 135 singleton loci (tier `c`) | `ingest/ehrenreich2012.py` |
| `breunig2014` (`b14`) | Breunig JS, Hackett SR, Rabinowitz JD, Kruglyak L. 2014. PLoS Genet 10:e1004142. doi:[10.1371/journal.pgen.1004142](https://doi.org/10.1371/journal.pgen.1004142) | segregant panel BY × RM; intracellular metabolites by LC-MS/MS; R/qtl | none: Table S1 gives cM only and the marker map is not in the bundle, so QTL are listed, not drawn | 53 linkages (text says 52), LOD only | 34 metabolites | 1 cross / 114 segregants | 4 genes named in the text (IRA2, URA3, SLT2, ERC1), tier `c` | `ingest/breunig2014.py` |
| `bloom2013` (`b13`) | Bloom JS, Ehrenreich IM, Loo WT, Lite T-LV, Kruglyak L. 2013. Finding the sources of missing heritability in a yeast cross. Nature 494:234. doi:[10.1038/nature11867](https://doi.org/10.1038/nature11867) | segregant panel, BY × RM, genotyped by sequencing; end-point colony growth | R64-compatible (listed genes sit under their intervals on R64-5-1); no transform | 591 QTL (1.5-LOD intervals, LOD, variance explained; no signed effects) | 46 | 1 cross / up to 1,007 segregants per trait | none (genes under intervals kept in source rows) | `ingest/bloom2013.py` |
| `ehrenreich2010` (`e10`) | Ehrenreich IM, Torabi N, Jia Y, Kent J, Martis S, Shapiro JA, Gresham D, Caudy AA, Kruglyak L. 2010. Dissection of genetically complex traits with extremely large pools of yeast segregants. Nature 464:1039. doi:[10.1038/nature08923](https://doi.org/10.1038/nature08923) | X-QTL in very large BY × RM pools (>10^5 segregants), 17 chemicals plus mitochondrial function | not stated; Table S2 spans consistent with R64; no transform | 2 tabulated loci (RAD5, MKT1 for 4-NQO); all other loci only in figures (not read off) | 18 | 1 cross / pools (not counted) | RAD5, MKT1 (tier `v`), HAP1 (tier `x`) | `ingest/ehrenreich2010.py` |
| `bloom2015` (`b15`) | Bloom JS, Kotenko I, Sadhu MJ, Treusch S, Albert FW, Kruglyak L. 2015. Genetic interactions contribute less than additive effects to quantitative trait variation in yeast. Nat Commun 6:8712. doi:[10.1038/ncomms9712](https://doi.org/10.1038/ncomms9712) | segregant panel BY × RM, 4,390 segregants, 20 end-point growth traits; additive QTL (5% FWER per chromosome) and QTL–QTL interactions (FDR 10%) | S288C R64 (sacCer3), no transform | 797 additive QTL (1.5-LOD intervals, signed effects); 205 + 266 interactions kept verbatim, not drawn | 20 | 1 cross / 4,390 (genotypes and phenotypes published) | 25 single-gene QTL, 13 genes (tier `c`) | `ingest/bloom2015.py` |

Source files, checksums and sheet-by-sheet contents: `sources/<id>/INVENTORY.md`. Validation
results and discrepancies with each paper: `tests/test_<id>.py` and `OPEN_QUESTIONS.md`.

## Next up

Not ingested yet. Each row is a stub: the inputs, design and build still have to be
confirmed from the paper's own files when it is ingested. Blocked hosts are listed in
`OPEN_QUESTIONS.md` (B2): Nature papers need `idp.nature.com` allowed; eLife pages refuse
non-browser requests; OUP/figshare are blocked. A PLOS Genetics search for the lab
(2026-11-01) also lists Albert et al. 2014, "Genetic Influences on Translation in Yeast"
(PLoS Genet 10:e1004692), reachable now; whether it reports QTL to ingest is unchecked.

| id (proposed) | paper | design | what to ingest | notes |
|---|---|---|---|---|
| `albert2014` | Albert et al. 2014, eLife | bulk eQTL-style (protein levels by X-pQTL) | pQTL per gene | trait = gene; needs a gene-trait type |
| `albert2018` | Albert et al. 2018, eLife | eQTL in ~1,000 BY×RM segregants | local and distant eQTL, hotspots | large: lazy-load by chromosome |
| `gou2019` | Gou et al. 2019, Genetics | mutation-rate QTL | QTL table | |
| `bruni_microscale` | Bruni et al., microscale X-QTL (`brunig/microscale-xqtl`) | 96-well X-QTL, 36 isolates × BBQ, 4 conditions | fitness and stress-contrast peaks, candidate genes | data already in that repo's `docs/data.json`; reuse its export |
