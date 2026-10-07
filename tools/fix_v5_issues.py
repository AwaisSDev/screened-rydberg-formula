"""One-off edits of docs/paper_draft.md after the first compiled PDF (srf-v5) was reviewed."""
import io
import re

P = "docs/paper_draft.md"
s = io.open(P, encoding="utf-8").read()


def rep(old, new, count=1):
    global s
    n = len(re.findall(r"\s+".join(re.escape(w) for w in old.split()), s))
    assert n == count, (n, old[:60])
    s = re.sub(r"\s+".join(re.escape(w) for w in old.split()), lambda m: new, s)


# Mendoza's own tables must not be renumbered as ours
rep("from Tables 1 and 2 of the authors' open-access deposit", "from the two constants tables of the authors' open-access deposit")
rep("the 84 IEs of its Table 3 to ≤ 0.023 % and its Tables 4–8 to ≤ 0.16 %",
    "the 84 IEs of its third table to ≤ 0.023 % and its fourth to eighth tables to ≤ 0.16 %")
# Table A1: the long class row is redundant (Table A2 lists every delta-tau_c)
rep("leaves the 19 of Table A1.", "leaves the 19 listed in Table A2.")
rep("the class definitions are in Table A2.", "the class definitions and their fitted deviations δτ_c are in Table A2.")
rep("(1.8027 as used in the code, which applies \\|κ\\| + 0.05)", "(used as |κ| + 0.05 = 1.8027)")
rep("δτ_c is from Table A1, and τ_g + δτ_c", "δτ_c is the fitted deviation, and τ_g + δτ_c")
lines = s.split("\n")
out = []
for ln in lines:
    if ln.startswith("| δτ_c (19 classes, defined in Table A2)"):
        out.append("| δτ_c (19 classes) | listed per class in Table A2 |")
    else:
        out.append(ln)
s = "\n".join(out)
io.open(P, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
