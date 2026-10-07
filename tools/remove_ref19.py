"""Remove the Di Rocco 1992 reference (could not be located) and renumber references 20-34 -> 19-33."""
import io
import re

P = "docs/paper_draft.md"
s = io.open(P, encoding="utf-8").read()
i = s.index("## References")
body, refs = s[:i], s[i:]

# the only sentence that cited [19]: support it with Pomarico et al. (old [20]), who state where the expressions are
pat = re.compile(r"\s+".join(re.escape(w) for w in
                 "because the original gives the coefficients as explicit expressions in a paper [19] that we could not access.".split()))
assert len(pat.findall(body)) == 1
body = pat.sub("because the original model gives the coefficients as explicit expressions in an earlier paper that we "
               "could not access [20].", body)
assert "[19]" not in body

renum = lambda m: "[%d]" % (int(m.group(1)) - 1) if 20 <= int(m.group(1)) <= 34 else m.group(0)
body = re.sub(r"\[(\d{2})\]", renum, body)

lines = []
for ln in refs.split("\n"):
    m = re.match(r"^\[(\d+)\] ", ln)
    if m:
        n = int(m.group(1))
        if n == 19:
            assert "Di Rocco" in ln and "Braz. J. Phys." in ln
            continue
        if n >= 20:
            ln = "[%d] " % (n - 1) + ln[m.end():]
    lines.append(ln)
refs = "\n".join(lines)
# remaining "as cited in" notes are gone with the entry; confirm
assert "as cited in" not in refs
io.open(P, "w", encoding="utf-8", newline="\n").write(body + refs)

T = "tools/to_pra.py"
t = io.open(T, encoding="utf-8").read()
a = t[t.index('    "H. O. Di Rocco, *Braz": '):t.index('    "R. D. Cowan":')]
t = t.replace(a, "")
io.open(T, "w", encoding="utf-8", newline="\n").write(t)

A = "docs/review/reference_audit.md"
r = io.open(A, encoding="utf-8").read().rstrip("\n")
r += """

## Decision: Di Rocco 1992 removed from the paper

The 1992 paper could not be located, so it was removed from the reference list. The one sentence that cited it now cites
Pomarico et al. 2005, who state that the explicit coefficient expressions are in the earlier paper. References 20-34 were
renumbered to 19-33; the paper now has 33 references. (The numbers printed by `tools/verify_refs.py` refer to the
earlier numbering.)
"""
io.open(A, "w", encoding="utf-8", newline="\n").write(r + "\n")
print("done")
