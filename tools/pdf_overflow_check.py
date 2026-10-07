"""Find content that extends past the page margins in a compiled PDF.   python tools/pdf_overflow_check.py file.pdf
Needs PyMuPDF (fitz). Reports, per page, text/image/drawing boxes whose right edge is within `tol` pt of the page edge
or beyond the typical text block, plus tiny font sizes."""
import sys
import fitz

d = fitz.open(sys.argv[1])
print("pages:", len(d), "| size:", d[0].rect.width, "x", d[0].rect.height)
xs0, xs1 = [], []
for p in d:
    for b in p.get_text("blocks"):
        if len(b[4].strip()) > 40:
            xs0.append(b[0]); xs1.append(b[2])
xs0.sort(); xs1.sort()
left, right = xs0[len(xs0) // 10], xs1[9 * len(xs1) // 10]
print(f"typical text block: x from {left:.0f} to {right:.0f}")
sizes = {}
for i, p in enumerate(d):
    issues = []
    for b in p.get_text("blocks"):
        if b[2] > right + 3 or b[0] < left - 3:
            issues.append(f"text x {b[0]:.0f}-{b[2]:.0f}: {b[4].strip()[:50]!r}")
    for im in p.get_image_info():
        r = im["bbox"]
        if r[2] > right + 3 or r[0] < left - 3:
            issues.append(f"image x {r[0]:.0f}-{r[2]:.0f}")
    for dr in p.get_drawings():
        r = dr["rect"]
        if r.width > 5 and (r.x1 > right + 3 or r.x0 < left - 3) and r.width > 100:
            issues.append(f"rule/box x {r.x0:.0f}-{r.x1:.0f}")
            break
    for blk in p.get_text("dict")["blocks"]:
        for ln in blk.get("lines", []):
            for sp in ln["spans"]:
                if sp["text"].strip():
                    sizes[round(sp["size"], 1)] = sizes.get(round(sp["size"], 1), 0) + len(sp["text"])
    if issues:
        print(f"page {i + 1}: {len(issues)} issue(s)")
        for x in issues[:6]:
            print("   ", x)
print("font sizes (pt: characters):", dict(sorted(sizes.items())))
