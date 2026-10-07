"""Structural sanity check of docs/pra/paper_pra.tex (no LaTeX needed): braces, environments, equations, tables."""
import collections
import re

t = open("docs/pra/paper_pra.tex", encoding="utf-8").read()
body = re.sub(r"(?<!\\)%.*", "", t)
op, cl = body.count("{") - body.count(r"\{"), body.count("}") - body.count(r"\}")
print("braces open/close:", op, cl, "OK" if op == cl else "MISMATCH")
env = collections.Counter(re.findall(r"\\(begin|end)\{(\w+\*?)\}", body))
names = {k[1] for k in env}
bad = {e: (env[("begin", e)], env[("end", e)]) for e in names if env[("begin", e)] != env[("end", e)]}
print("environment mismatches:", bad or "none")
print("equations:", len(re.findall(r"\\begin\{equation\}", body)), "| aligned:", len(re.findall(r"\\begin\{aligned\}", body)))
print("tables:", len(re.findall(r"\\begin\{table\*?\}", body)), "| adjustbox blocks:", len(re.findall(r"\\begin\{adjustbox\}", body)))
for m in re.finditer(r"\\begin\{equation\}(.*?)\\end\{equation\}", body, re.S):
    e = m.group(1).strip()
    print("  EQ lines:", e.count("\n") + 1, "| starts:", e[:36].replace("\n", " "), "| ends:", e[-26:].replace("\n", " "))
txt = re.sub(r"\\begin\{equation\}.*?\\end\{equation\}", "", body, flags=re.S)
print("code-like tokens left:", re.findall(r"\\texttt\{[^}]*\}|[A-Za-z]+\\_[A-Za-z]+", txt)[:15])
longest = max((len(l), i + 1) for i, l in enumerate(txt.split("\n")))
print("longest source line:", longest)
