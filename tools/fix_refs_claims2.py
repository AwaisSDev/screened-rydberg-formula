"""Second round of reference corrections after searching for the sources (2026-10-07)."""
import io
import re

P = "docs/paper_draft.md"
s = io.open(P, encoding="utf-8").read()


def rep(old, new):
    global s
    pat = re.compile(r"\s+".join(re.escape(w) for w in old.split()))
    assert len(pat.findall(s)) == 1, old[:60]
    s = pat.sub(lambda m: new, s)


rep("The regular variation of atomic energies along isoelectronic sequences is a standard tool of atomic spectroscopy [8], and screening theory has been applied to transition energies of highly charged ions [9].",
    "Edlén's review of atomic spectra [8] describes the regularities of energy levels along isoelectronic sequences, and screening theory combined with many-body perturbation theory has been applied to transition energies of highly charged ions [9].")
rep("[19] H. O. Di Rocco, *Braz. J. Phys.*, vol. 22, 1992, as cited in [20] (p. 227) and [22] (pp. 1–10); the original paper could not be located.",
    "[19] H. O. Di Rocco, *Braz. J. Phys.*, vol. 22, 1992, as cited in [20] and in H. O. Di Rocco, *Il Nuovo Cimento D*, vol. 20, pp. 131–140, 1998 (both give p. 227) and in [22] (pp. 1–10); the original paper could not be located.")
io.open(P, "w", encoding="utf-8", newline="\n").write(s)

T = "tools/to_pra.py"
t = io.open(T, encoding="utf-8").read()
a = t[t.index('    "H. O. Di Rocco, *Braz": '):t.index('    "R. D. Cowan":')]
b = ('    "H. O. Di Rocco, *Braz": "H. O. Di Rocco, Braz. J. Phys. **22** (1992), as cited in Ref. 20 and in H. O. Di Rocco, '
     'Il Nuovo Cimento D **20**, 131 (1998) (both give p. 227) and in Ref. 22 (pp. 1–10); the original paper could not '
     'be located.",\n')
t = t.replace(a, b)
io.open(T, "w", encoding="utf-8", newline="\n").write(t)

A = "docs/review/reference_audit.md"
r = io.open(A, encoding="utf-8").read()
r = r.rstrip("\n") + """

## Follow-up searches (same day)

* **Safronova 1993 [9]**: abstract found in the OSTI/ETDE record (via web search): the paper applies many-body perturbation
  theory together with the screening method to radiative transition energies of highly charged ions, and compares with
  Hartree–Fock–Pauli and relativistic model-potential results. This supports the sentence now in the paper.
* **Edlén 1964 [8]**: the Springer chapter description says the chapter treats the gross structure of atomic energy levels
  and the regularities along isoelectronic sequences. The paper's sentence now says exactly that.
* **Di Rocco 1992 [19]**: a third citing paper was read through Crossref's reference list, H. O. Di Rocco, *Il Nuovo Cimento D*
  **20**, 131 (1998), doi:10.1007/BF03036007, which prints "Braz. J. Phys., 22 (1992) 227". So p. 227 appears in two of
  the three citing papers and 1–10 in the third. The 1992 paper itself was not found (Crossref has no record, OSTI refused
  the connection, Springer required a login handshake that was not followed).
"""
io.open(A, "w", encoding="utf-8", newline="\n").write(r)
print("ok")
