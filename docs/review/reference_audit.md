# Reference audit (2026-10-07)

How each of the 34 references was checked. Scripts: `tools/verify_refs.py` (Crossref: first author, journal, volume, first
page, year, title), `tools/verify_refs2.py` (details Crossref search did not settle; Semantic Scholar abstracts),
`tools/verify_doi_resolve.py` (does every printed DOI resolve at doi.org). Web searches were used for the items Crossref
does not hold.

## Result

| refs | check | outcome |
|---|---|---|
| 1–7, 9, 11–18, 20–24, 26–31 (all journal articles) | Crossref record matches author, journal, volume, first page and year | **26 of 26 match.** Titles match too; [31] also checked for article number 033103, vol. 44, issue 3 |
| 4–9, 11, 12, 14, 15–18, 20–25, 30 (all printed DOIs, 14 distinct) | DOI resolves at doi.org | **14 of 14 resolve** (Royal Society and IOP answered the script with a block page but redirected to the right article) |
| 8 Edlén | Crossref chapter record (pp. 80–220, parent book "Spectroscopy I", ed. Flügge); library catalogues list Handbuch der Physik **Bd. XXVII**, Spektroskopie I, Springer 1964, 507 pp. | exists, vol. 27 confirmed |
| 10 Mayer LA-647 | Los Alamos report archive listing (sgp.fas.org/othergov/doe/lanl/lib-www/la-pubs/00407999) | exists, "Methods of opacity calculations", Oct. 1947 |
| 32 NIST ASD | NIST version-history page read directly | ASD **5.12**, 7 Nov 2024, DOI 10.18434/T4W30F |
| 33 Cowan | several library catalogues | exists; University of California Press, 1981, ISBN 0520038215, 731 pp. |
| 34 Froese Fischer, Brage, Jönsson | several library catalogues | exists; Institute of Physics Publishing, Bristol/Philadelphia, 1997, ISBN 0750303743 |
| 19 Di Rocco 1992 | **not in Crossref; only attested by citations** | see below |

## Errors found and corrected

1. **[19] page number was wrong.** The paper printed "22, 227 (1992)" from an earlier note. The two papers that cite it
   disagree: Pomarico et al. 2005 [20] print "22, 227 (1992)", Di Rocco and Lanzini 2016 [22] print "22, 1–10 (1992)".
   Crossref has no record of the 1992 paper itself. The entry now gives both and says the original could not be located.
2. **"The original uses fitted closed forms"** (Sec. IV E 1) was unsupported. Pomarico et al. say "no empirical
   adjustments are made" and that the explicit expressions for the coefficients are in Di Rocco 1992. The sentence now says
   the coefficients are given as explicit expressions in a paper we could not access.
3. **[25] Crilly** lacked volume and article number: now vol. 48, 101053 (2023), as Crossref gives.
4. **[30] Rodrigues 2004** covers Dirac–Fock energies for ions with 3 to 105 electrons, with Z up to 118 (its abstract),
   not "all ground configurations". Wording corrected, and "cover all ions" became "cover most ions".
5. **"1/Z_c ≈ 1.0975" for He** had no source (it came from an internal note). The literature value of the critical
   charge is Z_c = 0.9110282, which gives 1/Z_c ≈ 1.0977. The number was removed.
6. **"The Be-like value agrees with Layzer's"**: never checked against Layzer's paper. Removed.
7. **Edlén [8]**: the chapter exists, but the details I attributed to it ("expansions in 1/(ζ + s)", "formalised")
   could not be verified. The paper now says only that the regular variation of energies along isoelectronic sequences is
   a standard tool of atomic spectroscopy, and "Edlén-type" was dropped from the abstract, the contributions and the
   section title.
8. **Safronova [9]**: no abstract found. The claim "Z-expansion codes with relativistic corrections" was replaced by what
   the title supports (screening theory applied to transition energies of highly charged ions), and [9] was removed from
   the list of Z-expansion first-order tabulations in Sec. V A.
9. **Martel [14]** is about analytic expressions for charge-distribution momenta from analytical potentials (its
   abstract). "Derived from analytical potentials" became "based on analytical potentials".
10. **[6], [7]** (helium-sequence perturbation calculations): "higher orders need continuum sums" became "have been
    computed numerically for two-electron ions".
11. **[23]–[25]**: Rozsnyai is an average-atom Hartree–Fock–Slater model, not a screened hydrogenic one. The sentence now
    separates them, and only SpK [25] is said to use screened hydrogenic atoms (its abstract).

## Still not verified

* The content of Di Rocco 1992 [19], Edlén [8] and Safronova [9], which could not be read. They are cited only for
  their existence and general subject.
* Abstracts were not available for most older papers, so "cited for" was judged from titles, and from abstracts only for
  [16], [25], [14] and [30].

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

## Decision: Di Rocco 1992 removed from the paper

The 1992 paper could not be located, so it was removed from the reference list. The one sentence that cited it now cites
Pomarico et al. 2005, who state that the explicit coefficient expressions are in the earlier paper. References 20-34 were
renumbered to 19-33; the paper now has 33 references. (The numbers printed by `tools/verify_refs.py` refer to the
earlier numbering.)

