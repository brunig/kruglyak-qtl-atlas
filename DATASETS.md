# Datasets

## Ingested

| id | paper | design | build | QTL | traits | crosses / segregants | candidate genes | ingest |
|---|---|---|---|---|---|---|---|---|
| `bloom2019` (`b19`) | Bloom JS et al. 2019. eLife 8:e49212. doi:[10.7554/eLife.49212](https://doi.org/10.7554/eLife.49212) | 16-parent round-robin segregant panel; within-cross and joint (all crosses) models | S288C R64 (sacCer3), no transform | 4,552 joint (FDR 5%, lead variant only) + 7,751 within-cross (q < 0.05, 1.5-LOD intervals); 18,758 per-cross joint effects | 38 | 16 crosses of 16 parents / 13,950 | 427 gene-to-QTL-pair rows, 195 genes (102 at local FDR 5%) | `ingest/bloom2019.py` |

Source files, checksums and sheet-by-sheet contents: `sources/bloom2019/INVENTORY.md`.
Validation results and discrepancies with the paper: `tests/test_bloom2019.py`,
`OPEN_QUESTIONS.md`.

## Next up

Not ingested yet. Each row is a stub: the inputs, design and build still have to be
confirmed from the paper's own files when it is ingested.

| id (proposed) | paper | design | what to ingest | notes |
|---|---|---|---|---|
| `bloom2013` | Bloom et al. 2013, Nature | BY×RM segregant panel (~1,000 segregants), 46 traits | QTL table with intervals and effects | same cross as `bloom2019` cross A; link by BY/RM parents |
| `bloom2015` | Bloom et al. 2015, Nat Commun | BY×RM panel, ~4,400 segregants, 20 traits | QTL table, QTL-QTL interactions | interactions need a new table (pairs) |
| `ehrenreich2010` | Ehrenreich et al. 2010, Nature | X-QTL (bulk segregant), BY×RM | X-QTL peaks with intervals | point/interval rows with `m` = `x`; no per-segregant LOD |
| `ehrenreich2012` | Ehrenreich et al. 2012, PLoS Genet | X-QTL, multiple crosses | X-QTL peaks | |
| `treusch2015` | Treusch et al. 2015, Genetics | X-QTL across 12 crosses | Table S4 2-LOD intervals | already used as published loci by `microscale-xqtl` |
| `albert2014` | Albert et al. 2014, eLife | bulk eQTL-style (protein levels by X-pQTL) | pQTL per gene | trait = gene; needs a gene-trait type |
| `albert2018` | Albert et al. 2018, eLife | eQTL in ~1,000 BY×RM segregants | local and distant eQTL, hotspots | large: lazy-load by chromosome |
| `breunig2014` | Breunig et al. 2014, PLoS Genet | metabolite QTL, BY×RM | mQTL table | trait = metabolite |
| `gou2019` | Gou et al. 2019, Genetics | mutation-rate QTL | QTL table | |
| `bruni_microscale` | Bruni et al., microscale X-QTL (`brunig/microscale-xqtl`) | 96-well X-QTL, 36 isolates × BBQ, 4 conditions | fitness and stress-contrast peaks, candidate genes | data already in that repo's `docs/data.json`; reuse its export |
