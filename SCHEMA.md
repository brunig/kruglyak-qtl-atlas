# Atlas schema

The atlas is one JSON file, `docs/data.json`, plus lazily loaded files in `docs/v/`.
Tables are **columnar**: each table is an object of equal-length arrays under short
keys, as in the Microscale X-QTL browser (`toJSON(..., dataframe="columns")`). Row `i`
of a table is `{k: table[k][i] for k in table}`. Cross-table references are **row
indices** into the referenced table (`-1` for none). Missing values are `null`.

Every dataset adds rows to the same tables; nothing in the page is specific to one paper.
The `ds` column on each table is the row index of the dataset in `datasets`.

Coordinates are 1-based S288C reference positions, build R64 (UCSC sacCer3); see
`datasets.build`. Chromosomes are named `chrI` … `chrXVI`.

## Top level

| key | type | meaning |
|---|---|---|
| `schema` | string | schema version, `"atlas-2"` |
| `build` | string | genome build all coordinates use |
| `chroms` | object | chromosome → length (bp), R64 |
| `hb` | int | heatmap bin width (bp), 50,000 |
| `datasets`, `sources`, `parents`, `crosses`, `traits`, `qtl`, `jeff`, `causal`, `cand`, `genes`, `vexp` | table | below |

## `datasets`

| key | meaning |
|---|---|
| `id` | atlas id, e.g. `bloom2019` |
| `sc` | short code used in row ids and URLs, e.g. `b19` |
| `label` | short label for chips, e.g. `Bloom 2019` |
| `cite` | citation string |
| `doi` | DOI |
| `year` | publication year |
| `design` | study design (`segregant panel`, `X-QTL`, `joint model`, or a combination) |
| `build` | genome build of the source coordinates, and how it was established |
| `method` | mapping method, in one sentence per model |
| `fdr` | FDR threshold at which the QTL in this dataset were called; `null` when the paper used a LOD threshold instead (stated in `method`) |
| `interval` | what `qtl.a`–`qtl.e` mean for this dataset |
| `effect` | units and sign convention of `qtl.ef` |
| `nq` | number of QTL rows the dataset contributes |
| `notes` | anything a reader must know to interpret the rows |

## `sources`

One row per source table a value can come from. `qtl.pv` and every other provenance
key point here.

| key | meaning |
|---|---|
| `code` | provenance code, e.g. `F3D1:w` |
| `ds` | dataset index |
| `file` | file name under `sources/<dataset>/raw/` |
| `sheet` | sheet name |
| `hrow` | 1-based spreadsheet row holding the column names |
| `sha256` | checksum of the file |

A provenance key is `<code>#<row>`, where `<row>` is the 1-based spreadsheet row as
Excel shows it. `F3D1:w#2` is the first data row of the within-cross sheet.

## `parents`

| key | meaning |
|---|---|
| `id` | atlas id, `<sc>-<name>` |
| `ds` | dataset index |
| `name` | display name (`BY`, `RM`, `CLIB413`, …) |
| `code` | label used by the analysis code (`BYa`, `RMx`, …) |
| `s1` | strain description when listed as Parent 1 in the source |
| `s2` | strain description when listed as Parent 2 in the source (can differ: another isolate of the same background or the other mating type) |
| `p18` | isolate ID in Peter et al. 2018, if the source gives one |
| `org` | ecological origin, if the source gives one; `null` otherwise |

## `crosses`

| key | meaning |
|---|---|
| `id` | atlas id, `<sc>-<code>` |
| `ds` | dataset index |
| `code` | cross name in the source and code (`375`, `A`, …) |
| `p1`, `p2` | parent indices, in source order (Parent 1, Parent 2) |
| `name` | cross name as the source writes it (`Cross 2`), may be `null` |
| `g1` | parent index whose allele is genotype 1 in the within-cross model (the other parent is genotype 2). Inferred; see OPEN_QUESTIONS Q2 |
| `n` | number of segregants analysed; `null` for bulk-segregant pools that were not counted |
| `dip` | diploid parent strain |
| `mm` | magic-marker plasmid |
| `pv` | provenance key |

## `traits`

| key | meaning |
|---|---|
| `id` | atlas id, `<sc>-<slug>` (URL safe) |
| `ds` | dataset index |
| `raw` | trait string exactly as in the source (`Cadmium_Chloride;75uM;2`) |
| `name` | display name (`Cadmium chloride 75 µM`) |
| `cond` | second field of `raw` (dose, temperature or pH), may be empty |
| `bt` | third field of `raw`. Most defensible reading: experiment batch / replicate index (OPEN_QUESTIONS Q10) |
| `ck` | condition key for matching the same condition across datasets: lower-case compound plus dose as the source gives it (`caffeine 15mM`); `null` if unsure |
| `cat` | category: `drug/chemical`, `carbon source`, `temperature`, `pH`, `osmotic stress`, `reference medium`, `metabolite` |
| `catb` | where the category comes from (`Figure 1-source data 1, Phenotypes sheet: <section>` or `atlas assignment`) |

## `qtl`

One row per reported QTL. Within-cross QTL and joint-model QTL share the table and are
told apart by `m`.

| key | meaning |
|---|---|
| `id` | atlas id, `<sc>-<m>-<source row>` |
| `ds` | dataset index |
| `t` | trait index |
| `x` | cross index; `-1` for a QTL not tied to one cross (joint model, or a locus grouped over several crosses) |
| `xs` | list of cross indices in which a grouped QTL was detected (`null` when `x` is set or the QTL is a joint-model QTL) |
| `k` | short subtype label from the source, e.g. `grouped round-robin QTL`, `replicate-confirmed BY/RM QTL`, `MAT-linked`; `null` if the dataset has one kind |
| `m` | model: `w` linkage in one segregant cross (individual genotyped segregants), `j` joint model over several crosses, `x` bulk-segregant / X-QTL peak (pooled selection, allele frequencies) |
| `c` | chromosome |
| `a`, `e` | interval start and end. Within-cross: positions of the left and right markers of the 1.5-LOD drop interval. Joint: both equal the lead variant position (the source reports no interval) |
| `p` | peak position. Within-cross: the relocalized peak marker (equals the source `pos` column). Joint: the lead variant |
| `cm` | genetic position (cM) when the source gives no physical position; `a`, `e`, `p` are then `null` and the QTL is listed but not drawn on the genome |
| `lod` | LOD at the peak (within-cross only) |
| `q` | FDR q-value of the forward-search step that added the QTL (within-cross only; joint QTL passed q < 0.05 but their q is not reported) |
| `ef` | standardized effect, SD units (see sign convention below) |
| `ve` | fraction of phenotypic variance explained within the cross (within-cross only; joint per-cross values are in `jeff`) |
| `lv` | lead variant `chr_pos_ref_alt` (the source marker without the cross-specific index suffix) |
| `af` | minor allele frequency of the lead variant in the 1,011-isolate panel (Peter et al. 2018); joint only |
| `cc` | number of crosses in which the lead variant segregates (`crossCount`, rounded up to even in the source); joint only |
| `pa` | approximate number of the 16 parents carrying the minor allele, `cc / 2`, as the analysis code uses it (`jointModelAF.R`); joint only |
| `anc` | source `ancestral` flag: 1 = reference allele is ancestral (alternate allele derived), 0 = reference allele derived; joint only |
| `daf` | derived allele frequency in the 1,011 panel (`unfolded allele frequency`); joint only |
| `up` | within-cross: index of the parent whose allele increases the trait; joint: `-1` (the sign of `ef` says whether the alternate or reference allele increases it) |
| `ups`, `dns` | lists of parent indices the source names as carrying the beneficial (trait-increasing) or deleterious allele, when it reports direction by strain rather than by effect size; `null` otherwise |
|  | When a QTL has `xs` and `ups`/`dns` of the same length, entry *k* belongs to cross `xs[k]` (the favoured and disfavoured parent in that cross, as Ehrenreich 2012 reports); otherwise they are unordered strain lists (Treusch 2015). |
| `g` | candidate gene: gene index, `-1` if the paper assigned none |
| `gf` | local FDR of that assignment |
| `gt` | evidence tier of that gene (see `cand`) |
| `pv` | provenance key of the source row |

**Effect sign.** The phenotype is end-point colony radius (`norm.s.radius.mean` in
`phenotyping/code/process_images.R`; on condition plates, the residual after regressing
on the same segregant's control-plate radius), scaled to mean 0, SD 1 within each cross.
Positive effects mean larger colonies.

- Within-cross `ef` (`Beta`): coefficient on the standardized genotype (1 / 2 coding,
  scaled) in a multiple regression of the scaled phenotype on all QTL for that trait in
  that cross. `ef > 0` means the allele of the genotype-2 parent increases colony radius;
  `ef < 0` means the genotype-1 parent's allele does. `crosses.g1` says which parent is
  genotype 1; `up` resolves it per QTL.
- Joint `ef` (`MeanBeta`): mean over the crosses in which the lead variant segregates of
  the per-cross coefficient on the standardized reference/alternate genotype.
  `ef > 0` means the **alternate** allele (non-S288C) increases colony radius
  (`analysis/jointModelAF.R`: "direction of effect is - if reference allele increases,
  + if alternate allele increases").
- Because genotypes are standardized, `ef` is per SD of genotype. At a locus segregating
  1:1, the difference between the two allele classes is about `2 × ef` SD.

## `jeff` — joint-model effects per cross

One row per (joint QTL, cross) from `joint model (effect per cross)`.

| key | meaning |
|---|---|
| `q` | qtl index of the joint QTL |
| `x` | cross index |
| `b` | `beta` in that cross (same convention as joint `ef`) |
| `ve` | variance explained within that cross |
| `cn` | `cross.cnt` as given on that row (not rounded) |
| `pv` | provenance key |

## `genes`

SGD S288C annotation R64-5-1, in the Microscale X-QTL browser's gene-track format,
read from `sources/sgd/genes_R64-5-1.tsv`.

| key | meaning |
|---|---|
| `c`, `a`, `e`, `s` | chromosome, start, end, strand |
| `id` | systematic name |
| `n` | standard name (empty if none) |
| `k` | ORF class: `V` verified, `U` uncharacterized, `D` dubious, empty otherwise |

## `causal` — the paper's candidate-gene rows

One row per row of `Causal_Genes` (Figure 4-source data 1). Each row is one gene
assigned to a pair of within-cross QTL for the same trait in two crosses that share a
parent, with the probability mass of the pair's joint posterior summed over the gene.

| key | meaning |
|---|---|
| `id` | atlas id, `<sc>-c-<source row>` |
| `ds`, `t` | dataset and trait indices |
| `g` | gene index into `genes` (matched on ORF) |
| `orf`, `name` | as in the source |
| `c`, `a`, `e` | the gene's region as in the source: the gene extended by half of each flanking intergenic gap (`geneExtend.GR2` in `QTL_causality.R`), not the SGD gene bounds |
| `pcs` | `pCausalSum`: summed joint posterior probability of causality over the region |
| `sp` | shared parent index |
| `eff` | `increase` / `decrease`: whether the shared parent's allele increases or decreases the trait in both crosses |
| `fdr` | `localFDR`: cumulative FDR over rows ranked by `pcs` |
| `jppc`, `jv` | maximum joint PPC and the variant that has it |
| `q1`, `q2` | qtl indices of the two within-cross QTL (`-1` if the QTL is not in the within-cross table; OPEN_QUESTIONS Q12) |
| `tier` | evidence tier of this row |
| `pv` | provenance key |

## `cand` — candidate genes, one row per gene per dataset

| key | meaning |
|---|---|
| `ds`, `g` | dataset and gene indices |
| `tier` | best evidence tier over the gene's rows: `v` previously identified (the analysis code's list, `sources/bloom2019/code_tables/previously_identified.tsv`), `5` local FDR < 0.05, `20` local FDR < 0.2 |
|  | Further codes from later datasets: `x` confirmed by allele replacement, `c` candidate gene named by the paper without a statistical tier. Datasets without fine-mapping leave `causal.pcs`, `sp`, `eff`, `fdr`, `jppc`, `jv` null and `q2 = -1`. |
| `fdr` | lowest local FDR |
| `nt` | number of traits the gene is assigned for |
| `n` | number of rows in `causal` |

A gene in tier `v` is also in tier 5 or 20 by its FDR; `fdr` keeps that. Tiers are the
paper's evidence, not a claim of causality beyond it.

## `vexp` — variance explained per trait and cross

From Figure 2-source data 1, `Additive QTL Model`.

| key | meaning |
|---|---|
| `t`, `x` | trait and cross indices |
| `wve` | within-cross variance explained by QTL, cross-validated |
| `jve` | joint-model variance explained by QTL, cross-validated |
| `h2`, `h2se` | additive heritability (whole genome) and its standard error |
| `re`, `rese` | residual error and its standard error |
| `pv` | provenance key |

## `dens` — precomputed QTL density for datasets too large to load up front

One row per (dataset, kind, chromosome, 50 kb bin) with a non-zero count. Used for the
heatmap row of a dataset whose individual QTL live in a lazy file (Albert 2018 eQTL).

| key | meaning |
|---|---|
| `ds` | dataset index |
| `k` | kind, e.g. `local`, `distant` |
| `c`, `b` | chromosome and bin index (bin = floor(position / `hb`)) |
| `n` | number of QTL whose interval touches the bin |

## `lazy` — per-dataset QTL tables loaded on demand

| key | meaning |
|---|---|
| `ds` | dataset index |
| `name` | what the table holds, e.g. `eqtl` |
| `file` | path under `docs/`, e.g. `v/a18_eqtl.json` |
| `n` | number of rows |
| `desc` | one sentence for the page |

A lazy QTL file is columnar like `qtl`: `id`, `g` (index into `genes`: the gene whose
expression is the trait), `c`, `a`, `e`, `p`, `lod`, `ef` (as the source gives it; sign
convention in `datasets.effect`), `loc` (1 local / cis, 0 distant / trans), `pv`.

## Lazy files in `docs/v/`

- `<sc>_<chrom>.json` — the **original source rows**, every column verbatim, for all
  QTL and causal rows on that chromosome. Shape:
  `{"cols": {code: [column names]}, "rows": {code: {row id: [values]}}}`, where `code`
  is a `sources.code` and `row id` is the atlas id (`qtl.id`, `causal.id`), or
  `qtl.id` for `jeff` rows (a list of rows per joint QTL). The QTL card loads it on
  demand.
- `<sc>_extra.json` — source tables the views do not draw but that the atlas keeps
  verbatim: GO enrichments, variance-component models, allele-frequency variance
  components, the media table and the key-resources table. Shape:
  `{code: {"cols": [...], "rows": [[...], ...]}}`.

## Adding a dataset

Write `ingest/<id>.py` on the pattern of `ingest/bloom2019.py`: read only from
`sources/<id>/`, build a *fragment* (the same tables, with indices local to the
dataset; gene indices are global) and hand it to `ingest/atlas.py`, which rebases
indices and merges all fragments into `docs/data.json`. Keep the
original rows in `docs/v/<sc>_*.json`, fill `datasets` completely (build, method,
FDR, interval, effect convention), and add a `tests/test_<id>.py` that checks the
paper's reported numbers. A combined build (`build.sh`) runs each ingest in turn;
see README.
