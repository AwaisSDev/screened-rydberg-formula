"""Compile docs/pra/paper_pra.tex to docs/paper_pra.pdf with Tectonic (a self-contained LaTeX engine).

    set TECTONIC=path\to\tectonic.exe     (https://github.com/tectonic-typesetting/tectonic/releases)
    python tools/compile_pra.py

Compiles in a temporary copy (the first run downloads the TeX packages it needs), prints the error and overfull-box
counts from the log, and copies the PDF to docs/paper_pra.pdf. Overleaf (pdfLaTeX) is the reference compiler for
submission; this is a local check with the same REVTeX 4.2 class.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
exe = os.environ.get("TECTONIC") or shutil.which("tectonic")
if not exe:
    sys.exit("set the TECTONIC environment variable to tectonic.exe")
with tempfile.TemporaryDirectory() as tmp:
    for f in os.listdir(os.path.join(ROOT, "docs", "pra")):
        shutil.copy(os.path.join(ROOT, "docs", "pra", f), tmp)
    r = subprocess.run([exe, "-X", "compile", "paper_pra.tex", "--keep-logs"], cwd=tmp, capture_output=True, text=True)
    log = open(os.path.join(tmp, "paper_pra.log"), encoding="utf-8", errors="ignore").read()
    print("errors:", log.count("\n!"), "| overfull boxes:", log.count("Overfull"), "| undefined:", log.lower().count("undefined"))
    if not os.path.exists(os.path.join(tmp, "paper_pra.pdf")):
        sys.exit(r.stderr[-800:])
    shutil.copy(os.path.join(tmp, "paper_pra.pdf"), os.path.join(ROOT, "docs", "paper_pra.pdf"))
    print("wrote docs/paper_pra.pdf")
