# Kruglyak Lab QTL Atlas

A static browser for every QTL the Kruglyak lab has mapped in yeast, one dataset at a time,
in the visual language of the
[Microscale X-QTL browser](https://github.com/brunig/microscale-xqtl). 26 datasets so far: 14,679 QTL on the maps and 58,636 more
(eQTL, pQTL, single-cell eQTL) loaded per gene on demand, across 319 traits.

| dataset | design | QTL drawn or listed | per-gene QTL on demand |
|---|---|---:|---:|
| Boocock 2025 | one-pot single-cell eQTL: pools of MATa segregants from three crosses (A = BY × RM, 27,744 | 52 | 6,183 |
| Bloom 2019 | segregant panel | 12,303 | 0 |
| Gou 2019 | segregant panel, BY × RM (1040 segregants genotyped in Bloom et al. 2013 | 8 | 0 |
| Albert 2018 | segregant panel | 102 | 36,498 |
| Jerison 2017 | segregant panel: 230 BY × RM MATa segregants from Bloom et al. 2013, 8 replicate lines eac | 30 | 0 |
| Bloom 2015 | segregant panel | 797 | 0 |
| Treusch 2015 | X-QTL (bulk segregant) in a 12-strain round-robin cross (each strain in two crosses) | 97 | 0 |
| Albert 2014 | X-pQTL (bulk segregant): for each of 160 GFP-tagged genes, a BY × RM segregant pool sorted | 0 | 1,025 |
| Breunig 2014 | segregant panel | 53 | 0 |
| Wang 2014 | segregant panel, BY4724 × RM11-1a, 1008 segregants (the Bloom et al. 2013 panel) | 32 | 0 |
| Bloom 2013 | segregant panel | 591 | 0 |
| Ehrenreich 2012 | X-QTL (bulk segregant analysis of millions of segregants selected for resistance) in all 6 | 411 | 0 |
| Torabi 2012 | X-QTL (bulk segregant) in BY × RM: for each of nine conditions, a wildtype pool and a pool | 2 | 0 |
| Foss 2011 | segregant panel BY × RM, 95 segregants | 0 | 521 |
| Ehrenreich 2010 | X-QTL (bulk segregant): very large BY × RM segregant pools (more than 10^5 haploids) selec | 2 | 0 |
| Dimitrov 2009 | segregant panel, BY4716 × RM11-1a (122 segregants, 2957 markers from Brem and Kruglyak 200 | 3 | 0 |
| Demogines 2008 | segregant panel, BY × RM (123 segregants phenotyped, 121 genotyped at 2,956 markers) on MM | 6 | 0 |
| Smith 2008 | segregant panel BY × RM, 109 segregants expression-profiled (Agilent two-colour arrays) in | 41 | 11,392 |
| Foss 2007 | segregant panel BY × RM | 0 | 6 |
| Perlstein 2007 | segregant panel, BY × RM, 104 genotyped segregants grown with a panel of small molecules ( | 132 | 0 |
| Gatbonton 2006 | segregant panel, BY4716 × RM11-1a, 122 segregants genotyped at 3,312 microarray markers | 2 | 0 |
| Brem 2005 (Nature) | segregant panel, BY × RM (112 segregants | 2 | 0 |
| Brem 2005 (PNAS) | segregant panel, BY4716 × RM11-1a, 112 segregants | 0 | 0 |
| Ronald 2005 | segregant panel BY4716 × RM11-1a, 112 segregants (expression and genotypes of Brem & Krugl | 0 | 1,406 |
| Storey 2005 | segregant panel BY4716 × RM11-1a, 112 segregants | 0 | 340 |
| Yvert 2003 | segregant panel BY × RM, 86 segregants plus parents | 13 | 1,265 |

Details, counts and the queue of papers still to ingest: `DATASETS.md`.

**Open it:** https://brunig.github.io/kruglyak-qtl-atlas/ (once GitHub Pages is enabled; see
below). Locally: `python3 -m http.server --directory docs` and open http://localhost:8000.

## What it shows

- **Genome view.** A heatmap of QTL per 50 kb bin for every trait, so hotspots stand out,
  and one lollipop track per cross for the selected traits, with the joint model on its own
  track: a stem at each reported peak over a strip for its 1.5-LOD interval. Stem height is
  the LOD (log scale) or the signed effect, switchable; up means the second-named parent's
  allele (joint track: the alternate allele) makes larger colonies. These are the paper's
  reported QTL, not LOD curves: no paper so far published genome-wide scans. With several
  datasets selected, tracks are grouped by paper; a QTL a paper grouped across crosses is
  drawn in each of those crosses, and one reported without LOD or without an effect size is a
  hollow marker. QTL with a genetic position only (Breunig 2014) are listed in the table and
  card, not drawn.
- **Region view.** The same per-cross tracks for the selected traits, then every QTL from
  every trait and cross overlapping the region (or only the selected traits), the paper's candidate genes coloured by evidence tier, and the SGD gene
  track.
- **QTL card.** Click a QTL: trait, cross with parents, position and interval width, LOD,
  q, effect with its sign explained (which parent's or which allele increases colony size),
  lead variant with its frequency among the 1,011 sequenced isolates and the number of
  crosses it segregates in, per-cross effects for joint QTL, genes in the interval, the
  candidate gene if the paper assigned one, overlapping QTL from other traits and crosses,
  and the verbatim source row with file, sheet and row number.
- **Table view.** The full QTL table, filterable by dataset, model, trait, cross,
  chromosome, LOD and |effect|, sortable, with a download of the filtered rows as TSV.
- **eQTL and pQTL on demand.** Search `TDH3 expression` (or open a gene and press "Show eQTL
  and pQTL for its expression") to draw that gene's mRNA eQTL (Albert 2018) and protein pQTL
  (Albert 2014) genome-wide on adjacent tracks, local and distant, with the gene marked and the
  parent whose allele raises the level named for each. The heatmap carries eQTL density rows for the whole set.
- **Search.** Genes (`IRA2`, `YOL081W`), coordinates (`chrIV:509600`, `IV:400k-600k`),
  traits (`caffeine`), crosses (`375`, `M22 x BY`), parent strains (`YJM454`) and datasets
  (`Bloom 2019`, `round robin`).
- **Dataset chips**, a model switch (all / within-cross / joint / X-QTL), trait and cross filters,
  and a URL hash that carries the full state, with copy-link buttons.

The atlas shows what the papers reported. Candidate genes carry the paper's evidence tier
(previously identified, local FDR 5%, local FDR 20%) and are not claimed as causal beyond
it. Nothing is re-mapped or re-estimated.

## Layout

    docs/                    the site, served by GitHub Pages from /docs (no build step;
      index.html             no dependencies beyond Google Fonts)
      data.json              all datasets' normalized tables, columnar, short keys (SCHEMA.md)
      v/<sc>_<chrom>.json    every original source row, loaded when a QTL card asks for it
      v/<sc>_extra.json      source tables kept verbatim but not drawn (GO, variance components)
    sources/<dataset>/       data of record: the raw supplementary files, untouched, with
      raw/                   INVENTORY.md and MANIFEST.tsv (SHA-256)
      code_tables/           small tables transcribed from the authors' code, with line refs
    sources/sgd/             SGD R64-5-1 genes and R64 chromosome lengths
    ingest/                  inventory.py, atlas.py (shared merge), <dataset>.py
    tests/                   test_<dataset>.py (data validation), render.mjs (browser),
      screenshots/           the render test's screenshots
    SCHEMA.md                tables, keys, units, sign conventions
    DATASETS.md              ingested datasets and the queue
    OPEN_QUESTIONS.md        what the files did not settle, and every discrepancy with a paper

## Rebuilding

    pip install xlrd openpyxl pdfplumber pytest
    ./build.sh               # inventory, every ingest, docs/data.json and docs/v/
    ./build.sh test          # the same, then the data validation tests
    node tests/render.mjs    # browser test (Playwright with Chromium)

The build is deterministic: the same sources give byte-identical output. Each ingest runs
with `python3 -I` and reads only from `sources/`.

## Data of record

| dataset | source | build |
|---|---|---|
| Bloom et al. 2019 | eLife bundle `elife-49212-supp-v1.zip` | S288C R64 (sacCer3), no transform |
| Treusch et al. 2015 | PLOS supporting information S1-S14 and the article page | R64, no transform |
| Ehrenreich et al. 2012 | PLOS supporting information S1-S17 and the article page | consistent with R64 at locus resolution; no transform |
| Breunig et al. 2014 | PLOS supporting information S1-S8 and the article page | none (cM only) |
| Bloom et al. 2013 | Nature supplementary information (Tables 1-3, SI PDF) and the article page | R64-compatible, no transform |
| Bloom et al. 2015 | Nature Communications supplementary information and Data 1-5 (including genotypes and phenotypes) | R64, no transform |
| Ehrenreich et al. 2010 | Nature supplementary information PDF and the article page | consistent with R64, no transform |
| Albert et al. 2018 | eLife bundle `elife-35471-supp-v1.zip`, unpacked (zip checksum in `SOURCE_ZIP.tsv`) | R64, no transform |
| Albert et al. 2014 | Nature supplementary information and Data 1-3 and the article page | consistent with R64, no transform |
| Gou et al. 2019 | article PDF and supplement zip supplied by the user; QTL notebooks from `gouliangke/Mutation-rate` @ `bcaf987` | none (chromosome only) |
| genes | SGD S288C annotation R64-5-1, via `brunig/microscale-xqtl` `docs/data.json` | R64 |

Each dataset's raw files sit untouched in `sources/<id>/raw/` with checksums in
`sources/<id>/MANIFEST.tsv`.

For Bloom 2019, column semantics (interval method, effect units and signs, FDR procedure, how candidate
genes were assigned) come from the authors' analysis code,
[joshsbloom/yeast-16-parents](https://github.com/joshsbloom/yeast-16-parents) at commit
`c913c9a`; for Ehrenreich 2012, from the peak caller published with the paper; for the others,
from the paper's methods. They are written out in `SCHEMA.md` and each ingest. Where nothing
settled a question the atlas uses the most defensible reading and says so in `OPEN_QUESTIONS.md`.

## Validation

Each dataset has `tests/test_<id>.py`, which checks the paper's own reported numbers
against the ingested tables, coordinate sanity (inside R64 chromosomes, start ≤ peak ≤ end),
provenance, and that every source row is kept verbatim. Run them all with `./build.sh test`:
466 pass, 5 are skipped (not checkable from the files), and 45 are strict expected failures,
each a number in a paper that its own supplementary files do not reproduce, named with file
and column. `OPEN_QUESTIONS.md` explains each one. `node tests/render.mjs` checks the browser
(37 checks).

## Adding a dataset

1. Put the paper's files, untouched, in `sources/<id>/raw/`, and run or extend
   `ingest/inventory.py` to write `INVENTORY.md` and `MANIFEST.tsv`.
2. Read the authors' code for column semantics. Record anything unsettled in
   `OPEN_QUESTIONS.md`.
3. Copy `ingest/bloom2019.py` to `ingest/<id>.py`. Build a fragment with the same tables
   (`SCHEMA.md`): a `datasets` row filled completely (build, method, FDR, interval and effect
   conventions), its parents, crosses, traits and QTL, every row with a provenance key, and
   the original rows written to `docs/v/<sc>_*.json`. `atlas.save_fragment()` then
   `atlas.assemble()` merge it with the others.
4. Add it to the loop in `build.sh`, write `tests/test_<id>.py` against the paper's reported
   numbers, and move its row in `DATASETS.md` from the queue to the ingested table.

The page needs no change: datasets, traits, crosses and models are all read from
`data.json`. New QTL types (X-QTL peaks, eQTL) fit the same `qtl` table with a new `m` code;
the page draws any row with an interval as a bar and any point row as a tick.

## Publishing

GitHub Pages serves `/docs` (`docs/.nojekyll` is present). In the repository settings,
Pages → Build and deployment → Deploy from a branch → `main`, folder `/docs`.
