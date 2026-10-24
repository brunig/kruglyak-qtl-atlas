#!/usr/bin/env bash
# Rebuild the atlas from sources/. One command, deterministic.
#   ./build.sh            regenerate the inventory, run every ingest, write docs/data.json and docs/v/
#   ./build.sh test       also run the data validation tests
# Needs Python 3 with xlrd, openpyxl and pdfplumber (pip install xlrd openpyxl pdfplumber pytest).
set -euo pipefail
cd "$(dirname "$0")"
rm -rf build/fragments
python3 -I ingest/inventory.py
for f in ingest/bloom2019.py; do   # add new dataset ingests here; each merges into docs/data.json
  python3 -I "$f"
done
if [[ "${1:-}" == "test" ]]; then
  python3 -I -m pytest -q tests/test_bloom2019.py
fi
