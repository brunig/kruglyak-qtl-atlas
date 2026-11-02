"""Flatten the Word supplementary tables of Wang & Kruglyak 2014 into TSV, cell text verbatim.

    python3 -I sources/wang2014/extracted/extract_tables.py      (from the repository root)

Each output row carries the source file, the table number inside that .docx (1-based), the row
number inside that table (1-based, header = 1), and the paragraph that heads the table (the dose,
e.g. "40 μM haloperidol:"). tests/test_wang2014.py re-parses the .docx files and checks these
TSVs are identical, so they cannot drift from the sources.
"""
import os
import zipfile
import xml.etree.ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(os.path.dirname(HERE), "raw")
FILES = {"S1_Table": "journal.pgen.1004894.s002.docx", "S2_Table": "journal.pgen.1004894.s003.docx",
         "S3_Table": "journal.pgen.1004894.s004.docx", "S4_Table": "journal.pgen.1004894.s005.docx",
         "S5_Table": "journal.pgen.1004894.s006.docx", "S6_Table": "journal.pgen.1004894.s007.docx"}


def flatten(path):
    root = ET.fromstring(zipfile.ZipFile(path).read("word/document.xml"))
    out, last_para, ti = [], "", 0
    for el in root.find(W + "body"):
        if el.tag == W + "p":
            t = "".join(x.text or "" for x in el.iter(W + "t")).strip()
            if t:
                last_para = t
        elif el.tag == W + "tbl":
            ti += 1
            for ri, tr in enumerate(el.iter(W + "tr")):
                cells = ["".join(x.text or "" for x in tc.iter(W + "t")).strip() for tc in tr.findall(W + "tc")]
                out.append([os.path.basename(path), str(ti), str(ri + 1), last_para] + cells)
    return out


def main():
    for name, f in FILES.items():
        rows = flatten(os.path.join(RAW, f))
        width = max(len(r) for r in rows)
        head = ["file", "table", "row", "heading"] + [f"c{i}" for i in range(1, width - 3)]
        with open(os.path.join(HERE, name + ".tsv"), "w", encoding="utf-8") as fh:
            fh.write("\t".join(head) + "\n")
            for r in rows:
                fh.write("\t".join(r + [""] * (width - len(r))) + "\n")


if __name__ == "__main__":
    main()
