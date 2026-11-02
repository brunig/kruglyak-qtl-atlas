"""One-off extraction for Torabi & Kruglyak 2012 (run from the repository root):

    pip install olefile
    python3 -I sources/torabi2012/extracted/extract_sources.py

1. Tables S1-S3 are legacy Word (.doc) files. Their text is read through the Word 97 piece table
   (olefile) and written cell for cell to table_s1.tsv, table_s2.tsv, table_s3.tsv. Nothing is
   edited: values, units and footnote-free headers are as the .doc holds them.
2. text_facts.tsv: statements quoted verbatim from the article page (sources/torabi2012/article.html),
   each with its section heading. Every quote is checked to occur in the page text.
"""
import hashlib, html, os, re, struct, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
SRC = os.path.join(ROOT, "sources", "torabi2012")
OUT = os.path.join(SRC, "extracted")


def doc_text(path):
    import olefile
    ole = olefile.OleFileIO(path); w = ole.openstream("WordDocument").read()
    flags = struct.unpack("<H", w[0x0A:0x0C])[0]
    tbl = ole.openstream("1Table" if flags & 0x0200 else "0Table").read()
    csw = struct.unpack("<H", w[32:34])[0]; off = 34 + csw * 2
    cslw = struct.unpack("<H", w[off:off + 2])[0]; off += 2 + cslw * 4 + 2
    fcClx, lcbClx = struct.unpack("<II", w[off + 33 * 8:off + 33 * 8 + 8])
    clx = tbl[fcClx:fcClx + lcbClx]; i = 0
    while clx[i] == 1:
        i += 1 + 2 + struct.unpack("<H", clx[i + 1:i + 3])[0]
    lcb = struct.unpack("<I", clx[i + 1:i + 5])[0]; pt = clx[i + 5:i + 5 + lcb]; n = (lcb - 4) // 12
    cps = struct.unpack("<%dI" % (n + 1), pt[:4 * (n + 1)]); out = ""
    for k in range(n):
        fc = struct.unpack("<I", pt[4 * (n + 1) + 8 * k + 2:4 * (n + 1) + 8 * k + 6])[0]; ln = cps[k + 1] - cps[k]
        if fc & 0x40000000:
            s = (fc & 0x3FFFFFFF) // 2; out += w[s:s + ln].decode("cp1252", "replace")
        else:
            out += w[fc:fc + 2 * ln].decode("utf-16le", "replace")
    return out


def doc_table(path):
    """Caption line, then rows: cells end with \\x07 and each row ends with an extra \\x07."""
    t = doc_text(path)
    cap, body = t.split("\r", 1)
    rows = [r.split("\x07") for r in body.replace("\r", "").split("\x07\x07")]
    rows = [[c.strip() for c in r if c.strip() != "" or len(r) > 1] for r in rows]
    rows = [r for r in rows if any(c for c in r)]
    width = len(rows[0])
    return re.sub(r"\s*SEQ Table \\\* ARABIC\s*", "", cap).strip(), [r[:width] + [""] * (width - len(r)) for r in rows]


def page_text():
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace").read()
    m = re.search(r'<div class="article-text"[^>]*>(.*)', t, re.S); t = m.group(1) if m else t
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    secs = []
    for part in re.split(r"<h[1-6][^>]*>", t)[1:]:
        head, _, body = part.partition("</h")
        head = html.unescape(re.sub(r"<[^>]+>", " ", head)); body = html.unescape(re.sub(r"<[^>]+>", " ", body))
        secs.append((re.sub(r"\s+", " ", head).strip(), re.sub(r"\s+", " ", body)))
    return secs


FACTS = [
    ("design", "strains", "BY4724 (a laboratory strain) and RM11_1a (a wine strain)"),
    ("design", "pools", "a wildtype pool from a cross between wildtype BY and RM, and a sup35 pool from a cross between BY and RM both carrying sup35 C653R"),
    ("design", "replicates", "two biological replicates for the wildtype BY×RM cross and two biological replicates for the sup35 cross"),
    ("design", "array", "assay ∼18,000 single nucleotide polymorphisms (SNPs) between BY and RM"),
    ("design", "geo", "GEO Series accession number GSE33817"),
    ("count", "loci per condition", "The number of loci detected for growth rates in the nine conditions at an FDR of 5% ranged from one to 20 in both wildtype and sup35 crosses"),
    ("count", "readthrough-dependent loci", "18 loci showed significant differences between these pools"),
    ("count", "distinct loci", "Each readthrough-dependent locus had an effect in one to five conditions, with a total of ten distinct loci detected"),
    ("locus", "chrXIII", "we focused on the locus on the right arm of chromosome XIII, which affected growth in five conditions"),
    ("direction", "diamide", "At this locus, we detected a frequency skew in favor of the RM allele in the sup35 pool but not in the wildtype pool"),
    ("locus", "hydrogen peroxide", "SKY1 variation underlies the locus on Chromosome XIII for growth on hydrogen peroxide"),
    ("candidate", "SKY1", "These results demonstrate that polymorphisms in SKY1 are causal for the effects of this locus"),
    ("candidate", "MRE11", "ruling it out as the causal gene for this locus"),
    ("direction", "anova", "SKY1 allelic status also had a significant effect, with the RM allele increasing growth"),
    ("direction", "both conditions", "sup35 C653R and the BY allele of SKY1 act together to lower growth rate on diamide and hydrogen peroxide"),
]


def main():
    for k, name in (("s1", "s007"), ("s2", "s008"), ("s3", "s009")):
        f = os.path.join(SRC, "raw", f"journal.pgen.1002546.{name}.doc")
        cap, rows = doc_table(f)
        sha = hashlib.sha256(open(f, "rb").read()).hexdigest()
        with open(os.path.join(OUT, f"table_{k}.tsv"), "w") as fh:
            fh.write(f"# {cap}\n# extracted from raw/{os.path.basename(f)} (sha256 {sha}) by extract_sources.py\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        print(k, len(rows) - 1, "data rows", file=sys.stderr)
    secs = page_text()
    with open(os.path.join(OUT, "text_facts.tsv"), "w") as fh:
        fh.write("kind\ttopic\tsection\tquote\n")
        for kind, topic, q in FACTS:
            sec = next((h for h, b in secs if q in b), None)
            assert sec, q
            fh.write(f"{kind}\t{topic}\t{sec}\t{q}\n")
    print(len(FACTS), "facts", file=sys.stderr)


if __name__ == "__main__":
    main()
