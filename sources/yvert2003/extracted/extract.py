"""One-off extraction of Yvert et al. 2003 supplementary tables (PDF text) into TSVs.

    python3 -I sources/yvert2003/extracted/extract.py      (needs pdfplumber)

Writes, next to this script:
  groups_summary.tsv  Supplementary Table 2 summary (13 trans-acting groups), as printed
  group_members.tsv   Supplementary Table 2 member messages per group (ORF, rest of the printed line)
  clusters.tsv        Supplementary Table 1: 593 co-expression clusters (id, n genes, average/min correlation, member ORFs)
  text_facts.tsv      quoted sentences used by the ingest/tests, with their source
Each output starts with a '# source: <file> sha256 <hash>' line; the header is on line 2.
"""
import hashlib, os, re, sys
import pdfplumber
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "raw")
SI1 = os.path.join(RAW, "41588_2003_BFng1222_MOESM1_ESM.pdf")
SI2 = os.path.join(RAW, "41588_2003_BFng1222_MOESM2_ESM.pdf")
ORF = r"(Y[A-P][LR]\d{3}[WC](?:-[A-Z])?|Q\d{4}|R\d{4}[WC]?)"

def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

def text(p):
    with pdfplumber.open(p) as pdf:
        return [(i + 1, pg.extract_text() or "") for i, pg in enumerate(pdf.pages)]

def write(name, src, head, rows):
    with open(os.path.join(HERE, name), "w") as f:
        f.write(f"# source: {os.path.basename(src)} sha256 {sha(src)}\n")
        f.write("\t".join(head) + "\n")
        for r in rows:
            f.write("\t".join(str(x).replace("\t", " ") for x in r) + "\n")

def main():
    pages2 = text(SI2)
    # summary table: lines "N msgs function bin regulator pct" on the page that starts with "Summary of groups"
    summ = []
    for pno, t in pages2:
        if not t.startswith("Summary of groups"):
            continue
        for line in t.split("\n"):
            m = re.match(r"^(\d+) (\d+) (.+?) ([IVX]+):(\d+) (\S+) (\d+|N/A)$", line)
            if m:
                summ.append([m.group(1), m.group(2), m.group(3), m.group(4), m.group(5), m.group(6), m.group(7), pno, line])
            elif summ and re.match(r"^[A-Za-z].*", line) and not line.startswith(("Group", "function", "Summary")):
                summ.append(["", "", "", "", "", "", "", pno, line])   # continuation line (second function / regulator), kept verbatim
    write("groups_summary.tsv", SI2, ["group", "n_msgs", "common_function", "chrom_roman", "position", "putative_regulator", "pct_recapitulated_brem2002", "page", "printed_line"], summ)
    # members
    full = "\n".join(t for _, t in pages2)
    pageof = []
    for pno, t in pages2:
        pageof += [pno] * (len(t.split("\n")))
    lines = full.split("\n")
    mem, grp = [], None
    for k, line in enumerate(lines):
        g = re.match(r"^Group (\d+)$", line)
        if g:
            grp = g.group(1); continue
        if grp is None:
            continue
        m = re.match("^" + ORF + r"\b ?(.*)$", line)
        if m:
            mem.append([grp, m.group(1), m.group(2), pageof[k] if k < len(pageof) else ""])
    write("group_members.tsv", SI2, ["group", "orf", "printed_rest_of_line", "page"], mem)
    # clusters from SI1
    pages1 = text(SI1)
    cl, cur = [], None
    for pno, t in pages1:
        for line in t.split("\n"):
            m = re.search(r"cluster (\d+) \((\d+) genes\) average_corr: ([\d.]+) min_corr: ([\d.]+)", line)
            if m:
                cur = [m.group(1), m.group(2), m.group(3), m.group(4), pno, []]; cl.append(cur); continue
            o = re.match("^" + ORF + r"\b", line)
            if o and cur is not None:
                cur[5].append(o.group(1))
    write("clusters.tsv", SI1, ["cluster", "n_genes", "average_corr", "min_corr", "page", "member_orfs"], [c[:5] + [",".join(c[5])] for c in cl])
    # text facts
    art = open(os.path.join(HERE, "..", "article.html"), encoding="utf-8", errors="replace").read()
    import html as H
    abst = re.sub(r"\s+", " ", H.unescape(re.sub(r"<[^>]+>", " ", re.search(r'id="Abs1-content"[^>]*>(.*?)</div>', art, re.S).group(1)))).strip()
    s1 = pages1[0][1].replace("\n", " "); s2 = full.replace("\n", " ")
    facts = []
    def q(src, quote, hay):
        assert quote in hay, quote
        facts.append([src, quote])
    q("article.html abstract", "Most expression differences mapped to trans -acting loci.", abst)
    q("article.html abstract", "Positional cloning and functional assays showed that polymorphisms in GPA1 and AMN1 affect expression of genes involved in pheromone response and daughter cell separation, respectively.", abst)
    q("article.html abstract", "Notably, transcription factors showed no enrichment", abst)
    q("Supplementary Table 1 p1", "List of 593 clusters of genes showing correlated expression in the cross.", s1)
    q("Supplementary Table 1 p1", "86 segregants: 1 culture each", s1)
    q("Supplementary Table 2 p1", "we divided the genome into 20-kb bins and counted the number of messages linking to each bin (after removing cis-linkages).", s2)
    q("Supplementary Table 2 p1", "We empirically computed a significance threshold of 12 messages per bin", s2)
    q("Supplementary Table 2 p1", "Twenty-two bins fit this criterion.", s2)
    q("Supplementary Table 2 p1", "Group 2 contains linkages from four bins, Group 11 from three bins, and Groups 3, 7, 8, and 12 from two bins each.", s2)
    q("Supplementary Table 2 p1", "Group 2 contains the cluster of 12 genes shown in Figure 1 of the main text, whose regulator we confirmed to be Amn1.", s2)
    q("Supplementary Table 2 p3", "Group 7 contains the cluster of 7 pheromone response genes shown in Figure 1 of the main text, whose regulator we confirmed to be GPA1.", s2)
    q("Supplementary Table 2 p2", "Multiple additional regulators may control subsets of the 197 Group 2 genes", s2)
    q("Supplementary Table 2 p2", "Group 4 contains 100 messages linking to a region near the gene LEU2 (an auxotrophic marker deleted in the RM strain).", s2)
    q("Supplementary Table 2 p2", "The 36 messages of Group 5 map to the MAT locus.", s2)
    q("Supplementary Table 2 p2", "Group 6 contains 19 messages linking to the URA3 locus (another auxotrophic marker deleted in the RM strain).", s2)
    q("Supplementary Table 2 p3", "We showed previously that our BY strain carries a Ty insertion in HAP1", s2)
    q("Supplementary Table 2 p3", "The messages in Group 9 link to a region on chromosome XII containing SIR3.", s2)
    q("Supplementary Table 2 p3", "A good candidate in the linkage region is CAT5", s2)
    q("Supplementary Table 2 p2", "A good candidate is the RNA helicase MAK5", s2)
    q("Supplementary Table 2 p4", "Group 11 contains 380 messages", s2)
    write("text_facts.tsv", SI2, ["source", "quote"], facts)
    print(len(summ), len(mem), len(cl), len(facts), file=sys.stderr)

if __name__ == "__main__":
    main()
