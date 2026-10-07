"""Corrections after the reference audit (2026-10-07): make every sentence say only what the cited source supports."""
import io
import re

P = "docs/paper_draft.md"
s = io.open(P, encoding="utf-8").read()


def rep(old, new, count=1):
    global s
    pat = re.compile(r"\s+".join(re.escape(w) for w in old.split()))
    n = len(pat.findall(s))
    assert n == count, (n, old[:70])
    s = pat.sub(lambda m: new, s)


# abstract / contributions / section titles: "Edlen-type" implied a functional form we did not verify in the source
rep("depends on the ion charge in an Edlén-type form.", "depends on the ion charge and varies smoothly along isoelectronic sequences.")
rep("plus a bounded Edlén-type remainder with 33 global parameters", "plus a bounded remainder with 33 global parameters")
rep("### 2.3 Charge-dependent penetration as an Edlén-type term", "### 2.3 Charge-dependent penetration along isoelectronic sequences")
rep("the isoelectronic regularity that Edlén formalised [8], and it is consistent",
    "the kind of smooth isoelectronic regularity that atomic spectroscopy uses routinely [8], and it is consistent")
rep("The 1/(Z_a + κ) remainder is Edlén-type isoelectronic behaviour [8].",
    "The smooth variation of the 1/(Z_a + κ) remainder along isoelectronic sequences is a regularity of the kind long used in atomic spectroscopy [8].")
# introduction
rep("Higher orders need continuum sums [6], [7]. Along isoelectronic sequences, Edlén [8] described the smooth variation of screening with ion charge by expansions in 1/(ζ + s). Z-expansion codes with relativistic corrections were developed for selected sequences [9].",
    "Higher orders have been computed numerically for two-electron ions [6], [7]. The regular variation of atomic energies along isoelectronic sequences is a standard tool of atomic spectroscopy [8], and screening theory has been applied to transition energies of highly charged ions [9].")
rep("derived from analytical potentials [14], [15],", "based on analytical potentials [14], [15],")
rep("Average-atom and kinetics codes use SHMs because they are fast and cover all ions [23]–[25].",
    "Average-atom models [23] and kinetics and radiation-hydrodynamics codes [24], [25] need atomic data for many ions quickly; SpK, for example, uses screened hydrogenic atoms [25].")
rep("Dirac–Fock total energies for all ground configurations up to Z = 118 [30].",
    "Dirac–Fock total energies of ions with 3 to 105 electrons and Z up to 118 [30].")
# Z-expansion statements
rep("complex value E₁ = 1.5592742 agrees with Layzer's.", "complex value is E₁ = 1.5592742.")
rep("Z-expansion work [4], [5], [9] computed", "Z-expansion work [4], [5] computed")
rep("since SHMs and Dirac–Fock tables already cover all ions,", "since SHMs and Dirac–Fock tables already cover most ions,")
rep("E(Z) is analytic in 1/Z only up to a critical charge; for He, 1/Z_c ≈ 1.0975. Beyond first order",
    "The 1/Z series has a finite radius of convergence. Beyond first order")
rep("and is known only numerically [6].", "and has been evaluated only numerically [6].")
# Kregar/Di Rocco: Pomarico et al. state that no empirical adjustments are made and that the coefficient expressions are in [19]
rep("because the original uses fitted closed forms whose coefficients we could not obtain [19].",
    "because the original gives the coefficients as explicit expressions in a paper [19] that we could not access.")
# reference list entries
rep("[19] H. O. Di Rocco, *Braz. J. Phys.*, vol. 22, p. 227, 1992 (as cited in [20]; not independently verified).",
    "[19] H. O. Di Rocco, *Braz. J. Phys.*, vol. 22, 1992, as cited in [20] (p. 227) and [22] (pp. 1–10); the original paper could not be located.")
rep("*High Energy Density Phys.*, 2023, doi: 10.1016/j.hedp.2023.101053.",
    "*High Energy Density Phys.*, vol. 48, 101053, 2023, doi: 10.1016/j.hedp.2023.101053.")
io.open(P, "w", encoding="utf-8", newline="\n").write(s)

# to_pra.py: Physical Review style entries that are written out by hand
T = "tools/to_pra.py"
t = io.open(T, encoding="utf-8").read()
a = t[t.index('    "H. O. Di Rocco, *Braz": '):t.index('    "R. D. Cowan":')]
b = ('    "H. O. Di Rocco, *Braz": "H. O. Di Rocco, Braz. J. Phys. **22** (1992), as cited in Ref. 20 (p. 227) and Ref. 22 '
     '(pp. 1–10); the original paper could not be located.",\n')
t = t.replace(a, b)
io.open(T, "w", encoding="utf-8", newline="\n").write(t)
print("ok")
