"""Clean-up after the second compiled PDF (srf-v6): stray spaces before full stops and a stale equation number."""
import io
import re

P = "docs/paper_draft.md"
s = io.open(P, encoding="utf-8").read()
n_before = len(re.findall(r"(?<=\S) \.(?=\s|$)", s)) + len(re.findall(r"(?<=\S) ,(?=\s)", s))
s = re.sub(r"(?<=\S) \.(?=\s|$)", ".", s)          # "Z ≤ 118 ." -> "Z ≤ 118."
s = re.sub(r"(?<=\S) ,(?=\s)", ",", s)
assert s.count("Same\nequation (2.2) and inputs") + s.count("Same equation (2.2) and inputs") + \
    len(re.findall(r"Same\s+equation \(2\.2\) and inputs", s)) >= 1
s = re.sub(r"Same(\s+)equation \(2\.2\) and inputs", r"Same\1equation and inputs", s)
io.open(P, "w", encoding="utf-8", newline="\n").write(s)
print("fixed", n_before, "stray spaces;", "equation (2.2) left:", len(re.findall(r"equation \(2\.2\)", s)))
