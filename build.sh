#!/usr/bin/env bash
# Rebuild the atlas from sources/. One command, deterministic.
#   ./build.sh            regenerate the inventory, run every ingest, write docs/data.json and docs/v/
#   ./build.sh test       also run the data validation tests
# Needs Python 3 with xlrd, openpyxl and pdfplumber (pip install xlrd openpyxl pdfplumber pytest).
set -euo pipefail
cd "$(dirname "$0")"
rm -rf build/fragments
python3 -I ingest/inventory.py
for f in ingest/albert2014.py ingest/albert2018.py ingest/bloom2013.py ingest/bloom2015.py ingest/bloom2019.py ingest/boocock2025.py ingest/brem2005pnas.py ingest/breunig2014.py ingest/dimitrov2009.py ingest/ehrenreich2010.py ingest/ehrenreich2012.py ingest/foss2011.py ingest/gou2019.py ingest/jerison2017.py ingest/torabi2012.py ingest/treusch2015.py ingest/wang2014.py; do   # each merges into docs/data.json
  python3 -I "$f"
done
if [[ "${1:-}" == "test" ]]; then
  python3 -I -m pytest -q tests/
fi
