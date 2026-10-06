"""Build docs/paper_draft.pdf from docs/paper_draft.md.

    py -3.13 tools/build_paper_pdf.py

Markdown is rendered in the browser (marked) with LaTeX math typeset by KaTeX, then printed to PDF by
headless Chrome. Needs internet access for the two CDN libraries. The intermediate HTML is kept next to
the Markdown (docs/paper_draft.html) so relative figure paths (../results/figures/...) resolve.
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "docs" / "paper_draft.md"
HTML = ROOT / "docs" / "paper_draft.html"
PDF = ROOT / "docs" / "paper_draft.pdf"
CHROME = [
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
]

CSS = r"""
@page { size: A4; margin: 20mm 18mm 20mm 18mm; }
html { font-size: 10.5pt; }
body { font-family: "Cambria", "Georgia", "Times New Roman", serif; color: #111; line-height: 1.42;
       text-align: justify; hyphens: auto; max-width: 100%; margin: 0; }
h1 { font-size: 19pt; text-align: center; line-height: 1.25; margin: 0 0 10pt; }
h2 { font-size: 13pt; margin: 16pt 0 6pt; border-bottom: 0.6pt solid #888; padding-bottom: 2pt;
     page-break-after: avoid; }
h3 { font-size: 11.5pt; margin: 12pt 0 4pt; page-break-after: avoid; }
h4 { font-size: 10.5pt; margin: 10pt 0 3pt; page-break-after: avoid; }
p { margin: 0 0 6pt; }
hr { border: 0; border-top: 0.5pt solid #bbb; margin: 10pt 0; }
code { font-family: "Consolas", "Courier New", monospace; font-size: 8.6pt; background: #f3f3f3;
       padding: 0 2px; border-radius: 2px; word-break: break-word; }
pre { background: #f6f6f6; padding: 6pt; font-size: 8.4pt; white-space: pre-wrap; page-break-inside: avoid; }
pre code { background: none; padding: 0; }
table { border-collapse: collapse; width: 100%; margin: 6pt 0 10pt; font-size: 8.2pt; line-height: 1.25;
        page-break-inside: auto; text-align: left; }
thead { display: table-header-group; }
tr { page-break-inside: avoid; }
th { border-top: 1.1pt solid #000; border-bottom: 0.7pt solid #000; padding: 3pt 4pt; font-weight: bold;
     background: #fafafa; }
td { border-bottom: 0.3pt solid #ccc; padding: 2.5pt 4pt; vertical-align: top; }
tbody tr:last-child td { border-bottom: 1.1pt solid #000; }
img { display: block; max-width: 92%; max-height: 120mm; margin: 8pt auto 4pt; page-break-inside: avoid; }
blockquote { margin: 6pt 12pt; color: #333; border-left: 2pt solid #ccc; padding-left: 8pt; }
ul, ol { margin: 0 0 6pt 16pt; padding: 0; }
li { margin: 1pt 0; }
.katex { font-size: 1.04em; }
.katex-display { margin: 6pt 0; overflow: hidden; page-break-inside: avoid; }
.draftnote { text-align: center; font-size: 8.5pt; color: #8a1c1c; border: 0.6pt solid #8a1c1c;
             padding: 3pt; margin-bottom: 10pt; }
"""

JS = r"""
(function () {
  const src = JSON.parse(document.getElementById('md-src').textContent);
  const math = [];
  const esc = s => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  // protect display math, then inline math, from the Markdown parser
  let t = src.replace(/\$\$([\s\S]+?)\$\$/g, (m) => { math.push(m); return '@@MATH' + (math.length - 1) + '@@'; });
  t = t.replace(/(^|[^\\$])\$([^\$\n]+?)\$/g, (m, pre, body) => { math.push('$' + body + '$'); return pre + '@@MATH' + (math.length - 1) + '@@'; });
  let html = marked.parse(t, { gfm: true, breaks: false });
  html = html.replace(/@@MATH(\d+)@@/g, (m, i) => esc(math[+i]));
  document.getElementById('content').innerHTML = html;
  renderMathInElement(document.getElementById('content'), {
    delimiters: [{ left: '$$', right: '$$', display: true }, { left: '$', right: '$', display: false }],
    throwOnError: false
  });
  // KaTeX does not line-break display math: once its fonts are loaded (widths are final),
  // shrink any equation wider than the text block
  document.fonts.ready.then(() => {
    document.querySelectorAll('.katex-display').forEach(el => {
      let s = 1.0;
      while (el.scrollWidth > el.clientWidth + 1 && s > 0.45) { s -= 0.03; el.style.fontSize = s + 'em'; }
    });
    document.title = 'rendered';
  });
})();
"""


def build_html():
    md = MD.read_text(encoding="utf-8")
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>A screened Rydberg formula (draft)</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
<script src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/contrib/auto-render.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/marked@12.0.2/marked.min.js"></script>
<style>{CSS}</style></head>
<body>
<div class="draftnote">MANUSCRIPT DRAFT &mdash; not for distribution &mdash; items marked [TODO]/[VERIFY] must be resolved before submission</div>
<div id="content"></div>
<script id="md-src" type="application/json">{json.dumps(md).replace("</", "<\\/")}</script>
<script>{JS}</script>
</body></html>"""
    HTML.write_text(page, encoding="utf-8")


def print_pdf():
    exe = next((c for c in CHROME if os.path.exists(c)), None)
    if exe is None:
        sys.exit("no Chrome/Edge found")
    url = "file:///" + quote(str(HTML).replace("\\", "/"), safe="/:")
    if PDF.exists():
        PDF.unlink()
    cmd = [exe, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", "--no-first-run",
           "--run-all-compositor-stages-before-draw", "--virtual-time-budget=30000",
           f"--print-to-pdf={PDF}", url]
    subprocess.run(cmd, check=True, timeout=180, capture_output=True)
    for _ in range(30):
        if PDF.exists() and PDF.stat().st_size > 0:
            break
        time.sleep(0.5)
    print(f"wrote {PDF} ({PDF.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    build_html()
    print_pdf()
