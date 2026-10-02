#!/bin/sh
# Build an arXiv source bundle (main.tex + figures) from the AIMS manuscript.
# The arXiv version leaves out the journal label and the two declaration
# sections that the authors complete for the journal submission.
set -e
OUT=$(realpath -m "${1:-$(dirname "$0")/arxiv_source.zip}")
cd "$(dirname "$0")"
TMP=$(mktemp -d)
mkdir -p "$TMP/figures"
python3 - "$TMP/main.tex" <<'PY'
import re, sys
t = open("main.tex").read()
t = t.replace("\\noindent{\\small\\itshape Research article}\n\n\\vspace{6mm}\n", "")
t = re.sub(r"\\section\*\{Author contributions\}\n\n\[To be completed by the authors\.\]\n\n", "", t)
t = re.sub(r"\\section\*\{Use of AI tools declaration\}\n\n\[To be completed by the authors\.\]\n\n", "", t)
assert "To be completed" not in t and "Research article" not in t
open(sys.argv[1], "w").write(t)
PY
cp figures/fig1_cube.pdf figures/fig2_torus.pdf figures/fig3_exceptions.pdf figures/fig4_rho_sum.pdf "$TMP/figures/"
(cd "$TMP" && for i in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error main.tex > /dev/null; done)
grep -q "undefined" "$TMP/main.log" && { echo "undefined references"; exit 1; }
rm -f "$OUT"
(cd "$TMP" && zip -q -r "$OUT" main.tex figures)
cp "$TMP/main.pdf" "${OUT%.zip}_preview.pdf"
rm -rf "$TMP"
echo "wrote $OUT"
