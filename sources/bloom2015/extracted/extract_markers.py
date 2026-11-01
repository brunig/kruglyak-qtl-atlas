"""One-off extraction: marker names (genotype-matrix column names) from Supplementary Data 5.

    pip install rdata
    python3 -I sources/bloom2015/extracted/extract_markers.py <dir with genotypes.RData>

Unzip raw/41467_2015_BFncomms9712_MOESM729_ESM.zip, then Supplementary_Data_5/genotypes.zip, into
an empty scratch directory first (the archive is untrusted data; nothing is executed from it).
genotypes.RData holds one object, `gdata`: 4,390 segregants x 28,220 markers, coded -1 (BY) / 1 (RM).
Writes markers.tsv next to this script: the 1-based column index (the `index` / `Q1_index` /
`Q2_index` columns of Supplementary Data 3 and 4 refer to it) and the column name, split into its
parts. Only names are extracted; no genotype or phenotype value is copied.
"""
import os
import sys

import rdata

src = os.path.join(sys.argv[1], "genotypes.RData")
v = rdata.conversion.convert(rdata.parser.parse_file(src))["gdata"]
cols = [str(x) for x in v.coords[v.dims[1]].values]
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "markers.tsv")
with open(out, "w") as f:
    f.write("# Column names of gdata in Supplementary Data 5 (genotypes.RData), Bloom et al. 2015.\n")
    f.write(f"# {v.shape[0]} segregants x {v.shape[1]} markers; written by extract_markers.py. name = cumulative_chr_pos_ref_alt\n")
    f.write("index\tname\tchrom\tpos\tref\talt\n")
    for i, c in enumerate(cols, 1):
        cum, chrom, pos, ref, alt = c.split("_", 4)
        f.write(f"{i}\t{c}\t{chrom}\t{pos}\t{ref}\t{alt}\n")
print("wrote", out, len(cols), "markers", file=sys.stderr)
