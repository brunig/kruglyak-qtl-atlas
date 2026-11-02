"""Ingest Brem and Kruglyak 2005 (PNAS 102:1572), the landscape of genetic complexity across 5,700
gene expression traits in BY x RM.

    python3 -I ingest/brem2005pnas.py          # from the repository root

Reads only from sources/: sources/brem2005pnas/extracted/text_facts.tsv and tables.tsv (verbatim
extractions from the PMC article page, sources/brem2005pnas/article.html, made by extracted/extract.py),
the E-utilities record in raw/, and the SGD tables in sources/sgd/. Writes
build/fragments/brem2005pnas.json, docs/v/b05p_extra.json and, through ingest/atlas.py, docs/data.json.
Deterministic.

What the sources hold: summaries only. The article reports how many transcripts have a QTL, the
distribution of variance explained (Table 1), bounds on the number of loci per transcript (Table 2)
and counts from directional, transgression and epistasis tests. No per-transcript or per-QTL linkage
result is tabulated, and the PMC supplementary files sit behind a proof-of-work page that was not
circumvented. So the dataset contributes its design, one summary trait and no QTL rows. Nothing is
invented to fill the gap.
"""
import collections
import csv
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "brem2005pnas", "b05p"
SRC = os.path.join(atlas.SOURCES, "brem2005pnas")
TF = "../extracted/text_facts.tsv"   # relative to raw/, as sources.file is
TB = "../extracted/tables.tsv"
JATS = "PMC547855_efetch_jats.xml"


def sha_of(rel):
    return hashlib.sha256(open(os.path.join(SRC, "raw", rel), "rb").read()).hexdigest()


def tsv(rel):
    with open(os.path.join(SRC, "raw", rel), encoding="utf-8") as f:
        r = list(csv.reader(f, delimiter="\t"))
    return r[0], r[1:]


def main():
    C = collections.OrderedDict
    T = {k: C() for k in atlas.TABLES}

    def add(table, **kw):
        tab = T[table]
        n = atlas.nlen(tab)
        for k, v in kw.items():
            if k not in tab:
                tab[k] = [None] * n
            tab[k].append(v)
        for k in tab:
            if k not in kw:
                tab[k].append(None)
        return n

    for code, f, sheet in (("b05p:TF", TF, "text_facts"), ("b05p:T", TB, "tables"), ("b05p:JATS", JATS, "efetch record")):
        add("sources", code=code, ds=0, file=f, sheet=sheet, hrow=1, sha256=sha_of(f))

    th, trows = tsv(TF)
    facts = [dict(zip(th, r)) for r in trows]

    def fact(value_start):
        for i, d in enumerate(facts):
            if d["value"].startswith(value_start):
                return i + 2   # 1-based row in the TSV, header on row 1
        raise SystemExit("missing fact " + value_start)

    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY4716", s1="BY4716, isogenic to S288C", s2=None, p18=None, org="laboratory")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM11-1a", s1=None, s2="RM11-1a, wild isolate", p18=None, org=None)
    add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=112, dip=None, mm=None,
        pv=f"b05p:TF#{fact('112 segregants')}")
    add("traits", id=f"{SC}-transcript-levels", ds=0, raw="transcript levels (summary of ~5,700 expression traits)",
        name="Transcript levels (summary only)", cond="", bt="", ck="transcript levels", cat="expression",
        catb="atlas assignment (expression traits; per-transcript results not tabulated)")

    add("datasets", id=DS_ID, sc=SC, label="Brem 2005 (PNAS)",
        cite="Brem RB, Kruglyak L. 2005. The landscape of genetic complexity across 5,700 gene expression traits in yeast. PNAS 102:1572–1577",
        doi="10.1073/pnas.0408709102", year=2005,
        design="segregant panel, BY4716 × RM11-1a, 112 segregants; transcript levels of ~5,700 genes on two-colour microarrays (6 BY and 12 RM parental cultures as reference)",
        build="none: no marker or QTL positions are given",
        method="Linkage at 2,957 markers with a Wilcoxon test, significance by permutation; FDR 0.05 (nominal P < 5.7e-5). Variance explained estimated in a held-out half of the segregants; heritability, directional, transgression and epistasis tests on parent and segregant distributions.",
        fdr=0.05,
        interval="none reported",
        effect="No effect sizes or directions are tabulated; Table 1 gives the distribution of variance explained by detected QTLs.",
        nq=0,
        notes="Summary-only dataset: 2,984 transcripts with at least one QTL at FDR 0.05; of 3,546 highly heritable transcripts, 2,091 (59%) have a QTL, 406 (11%) pass the directional test, 2,093 (59%) the transgression test and 583 (16%) the epistasis test; median variance explained 27%. Per-transcript results are not tabulated in the article, and its PMC supplementary files could not be fetched (proof-of-work page). Tables 1 and 2 and every quoted number are kept in docs/v/b05p_extra.json.")

    atlas.save_fragment(DS_ID, T)

    bh, brows = tsv(TB)
    extra = {"b05p:TF": {"cols": th, "rows": trows}, "b05p:T": {"cols": bh, "rows": brows}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)
    data = atlas.assemble()
    print(f"{DS_ID}: 0 QTL rows (summary only); {len(trows)} quoted facts, {len(brows)} table rows; "
          f"docs/data.json holds {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
