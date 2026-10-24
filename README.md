# Kruglyak Lab QTL Atlas

A static browser for every QTL the Kruglyak lab has mapped in yeast, one dataset at a time,
in the visual language of the
[Microscale X-QTL browser](https://github.com/brunig/microscale-xqtl). The first dataset is
the 16-parent round-robin segregant panel of Bloom et al. 2019 (eLife 8:e49212): 4,552
joint-model QTL and 7,751 within-cross QTL for 38 traits in 16 crosses (13,950 segregants).

**Open it:** https://brunig.github.io/kruglyak-qtl-atlas/ (once GitHub Pages is enabled; see
below). Locally: `python3 -m http.server --directory docs` and open http://localhost:8000.

## What it shows

- **Genome view.** A heatmap of QTL per 50 kb bin for every trait, so hotspots stand out,
  and the selected traits' QTL by cross, with the joint model's lead variants on their own
  row. Within-cross intervals are shaded by LOD; joint lead variants are sized by |effect|.
- **Region view.** Every QTL from every trait and cross overlapping a region (or only the
  selected traits), the paper's candidate genes coloured by evidence tier, and the SGD gene
  track.
- **QTL card.** Click a QTL: trait, cross with parents, position and interval width, LOD,
  q, effect with its sign explained (which parent's or which allele increases colony size),
  lead variant with its frequency among the 1,011 sequenced isolates and the number of
  crosses it segregates in, per-cross effects for joint QTL, genes in the interval, the
  candidate gene if the paper assigned one, overlapping QTL from other traits and crosses,
  and the verbatim source row with file, sheet and row number.
- **Table view.** The full QTL table, filterable by dataset, model, trait, cross,
  chromosome, LOD and |effect|, sortable, with a download of the filtered rows as TSV.
- **Search.** Genes (`IRA2`, `YOL081W`), coordinates (`chrIV:509600`, `IV:400k-600k`),
  traits (`caffeine`), crosses (`375`, `M22 x BY`), parent strains (`YJM454`) and datasets
  (`Bloom 2019`, `round robin`).
- **Dataset chips**, a model switch (joint / within-cross / both), trait and cross filters,
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
| Bloom et al. 2019 | eLife supplementary bundle `elife-49212-supp-v1.zip` (Figure 1-4 source data, Supplementary files 1-3), checksums in `sources/bloom2019/MANIFEST.tsv` | S288C R64 (sacCer3), no transform |
| genes | SGD S288C annotation R64-5-1, via `brunig/microscale-xqtl` `docs/data.json` | R64 |

Column semantics (interval method, effect units and signs, FDR procedure, how candidate
genes were assigned) come from the authors' analysis code,
[joshsbloom/yeast-16-parents](https://github.com/joshsbloom/yeast-16-parents) at commit
`c913c9a`, and are written out in `SCHEMA.md`. Where the code did not settle a question the
atlas uses the most defensible reading and says so in `OPEN_QUESTIONS.md`.

## Validation

`tests/test_bloom2019.py` checks the paper's numbers against the ingested tables: 16
parents in 16 crosses, 13,950 segregants (mean 872), 38 traits, 4,552 joint QTL (mean 120 per
trait, range 52-195), a median of 12 within-cross QTL per trait and cross, numbers the
authors recorded in their code, coordinates inside R64 chromosomes, start ≤ peak ≤ end,
candidate-gene placement, the cAMP-pathway genes, Supplementary file 3's per-trait counts,
and that every row keeps its source row. Six paper numbers the data do not reproduce are
strict expected failures that name the file and column (105 genes at FDR 5%, 458
joint-only QTL, 24 of 26 validated genes, and others); see `OPEN_QUESTIONS.md`.

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
