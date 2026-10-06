"""Convert docs/paper_draft.md citations [Key; Key] to IEEE numeric style and rewrite the reference list.

    py -3.11 tools/ieee_refs.py      (idempotent only on the key-style draft; run once)

Reference data were checked on 2026-10-06 against Crossref (DOIs given where Crossref confirmed them).
Layzer1967 could not be found in Crossref and is no longer cited; DiRocco1992 is not in Crossref and is labelled
"as cited in" Pomarico2005.
"""
import io
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = os.path.join(ROOT, "docs", "paper_draft.md")

R = {
    "Slater1930": "J. C. Slater, “Atomic shielding constants,” *Phys. Rev.*, vol. 36, pp. 57–64, 1930.",
    "ClementiRaimondi1963": "E. Clementi and D. L. Raimondi, “Atomic screening constants from SCF functions,” *J. Chem. Phys.*, vol. 38, pp. 2686–2689, 1963.",
    "ClementiRaimondiReinhardt1967": "E. Clementi, D. L. Raimondi, and W. P. Reinhardt, “Atomic screening constants from SCF functions. II. Atoms with 37 to 86 electrons,” *J. Chem. Phys.*, vol. 47, pp. 1300–1307, 1967.",
    "Layzer1959": "D. Layzer, “On a screening theory of atomic spectra,” *Ann. Phys. (N.Y.)*, vol. 8, pp. 271–296, 1959, doi: 10.1016/0003-4916(59)90023-5.",
    "LayzerEtAl1964": "D. Layzer, Z. Horák, M. N. Lewis, and D. P. Thompson, “Second-order Z-dependent theory of many-electron atoms,” *Ann. Phys. (N.Y.)*, vol. 29, no. 1, pp. 101–124, 1964, doi: 10.1016/0003-4916(64)90192-7.",
    "ScherrKnight1963": "C. W. Scherr and R. E. Knight, “Two-electron atoms III. A sixth-order perturbation study of the 1¹S ground state,” *Rev. Mod. Phys.*, vol. 35, no. 3, pp. 436–442, 1963, doi: 10.1103/RevModPhys.35.436.",
    "DalgarnoStewart1958": "A. Dalgarno and A. L. Stewart, “A perturbation calculation of properties of the helium iso-electronic sequence,” *Proc. R. Soc. Lond. A*, vol. 247, no. 1249, pp. 245–259, 1958, doi: 10.1098/rspa.1958.0182.",
    "Edlen1964": "B. Edlén, “Atomic spectra,” in *Handbuch der Physik / Encyclopedia of Physics*, vol. 27, *Spectroscopy I*, S. Flügge, Ed. Berlin, Germany: Springer, 1964, pp. 80–220, doi: 10.1007/978-3-662-35391-2_2.",
    "Safronova1993": "U. I. Safronova, I. Yu. Tolstikhina, R. Bruch, T. Tanaka, F. Hao, and D. Schneider, “Screening theory for transition energies of highly charged ions,” *Phys. Scr.*, vol. 47, no. 3, pp. 364–382, 1993, doi: 10.1088/0031-8949/47/3/007.",
    "Mayer1947": "H. Mayer, “Methods of opacity calculations,” Los Alamos Sci. Lab., Los Alamos, NM, USA, Rep. LA-647, 1947.",
    "More1982": "R. M. More, “Electronic energy levels in dense plasmas,” *J. Quant. Spectrosc. Radiat. Transf.*, vol. 27, no. 3, pp. 345–357, 1982, doi: 10.1016/0022-4073(82)90127-3.",
    "Faussurier1997": "G. Faussurier, C. Blancard, and A. Decoster, “New screening coefficients for the hydrogenic ion model including l-splitting for fast calculations of atomic structure in plasmas,” *J. Quant. Spectrosc. Radiat. Transf.*, vol. 58, no. 2, pp. 233–260, 1997, doi: 10.1016/S0022-4073(97)00018-6.",
    "Faussurier2008": "G. Faussurier, C. Blancard, and P. Renaudin, “Equation of state of dense plasmas using a screened-hydrogenic model with l-splitting,” *High Energy Density Phys.*, vol. 4, pp. 114–123, 2008.",
    "Martel1998": "P. Martel, J. G. Rubiano, J. M. Gil, L. Doreste, and E. Mínguez, “Analytical expressions for the n-order momenta of charge distribution for ions,” *J. Quant. Spectrosc. Radiat. Transf.*, vol. 60, no. 4, pp. 623–633, 1998, doi: 10.1016/S0022-4073(97)00226-4.",
    "Rubiano2002": "J. G. Rubiano, R. Rodríguez, J. M. Gil, F. H. Ruano, P. Martel, and E. Mínguez, “A screened hydrogenic model using analytical potentials,” *J. Quant. Spectrosc. Radiat. Transf.*, vol. 72, no. 5, pp. 575–588, 2002, doi: 10.1016/S0022-4073(01)00142-X.",
    "Mendoza2011": "M. A. Mendoza, J. G. Rubiano, J. M. Gil, R. Rodríguez, R. Florido, P. Martel, and E. Mínguez, “A new set of relativistic screening constants for the screened hydrogenic model,” *High Energy Density Phys.*, vol. 7, no. 3, pp. 169–179, 2011, doi: 10.1016/j.hedp.2011.04.006.",
    "Kregar1984": "M. Kregar, “The virial and the independent particle models of the atom,” *Phys. Scr.*, vol. 29, no. 5, pp. 438–447, 1984, doi: 10.1088/0031-8949/29/5/005.",
    "Kregar1985": "M. Kregar, “The virial as the atomic model potential energy operator,” *Phys. Scr.*, vol. 31, no. 4, pp. 246–254, 1985, doi: 10.1088/0031-8949/31/4/005.",
    "DiRocco1992": "H. O. Di Rocco, *Braz. J. Phys.*, vol. 22, p. 227, 1992 (as cited in @Pomarico2005@; not independently verified).",
    "Pomarico2005": "J. Pomarico, D. I. Iriarte, and H. O. Di Rocco, “An efficient screening approach to be used in plasma modeling and ion-surface collision experiments,” *Braz. J. Phys.*, vol. 35, no. 1, pp. 130–135, 2005, doi: 10.1590/S0103-97332005000100008.",
    "LanziniDiRocco2015": "F. Lanzini and H. O. Di Rocco, “Screening parameters for the relativistic hydrogenic model,” *High Energy Density Phys.*, vol. 17, pp. 240–247, 2015, doi: 10.1016/j.hedp.2015.08.002.",
    "DiRoccoLanzini2016": "H. O. Di Rocco and F. Lanzini, “Breit and quantum electrodynamics energy contributions in multielectron atoms from the relativistic screened hydrogenic model,” *Braz. J. Phys.*, vol. 46, pp. 175–183, 2016, doi: 10.1007/s13538-015-0397-9.",
    "Rozsnyai1972": "B. F. Rozsnyai, “Relativistic Hartree-Fock-Slater calculations for arbitrary temperature and matter density,” *Phys. Rev. A*, vol. 5, no. 3, pp. 1137–1149, 1972, doi: 10.1103/PhysRevA.5.1137.",
    "Chung2005": "H.-K. Chung, M. H. Chen, W. L. Morgan, Yu. Ralchenko, and R. W. Lee, “FLYCHK: Generalized population kinetics and spectral model for rapid spectroscopic analysis for all elements,” *High Energy Density Phys.*, vol. 1, pp. 3–12, 2005.",
    "Crilly2023": "A. J. Crilly *et al.*, “SpK: A fast atomic and microphysics code for the high-energy-density regime,” *High Energy Density Phys.*, 2023, doi: 10.1016/j.hedp.2023.101053.",
    "Koopmans1934": "T. Koopmans, “Über die Zuordnung von Wellenfunktionen und Eigenwerten zu den einzelnen Elektronen eines Atoms,” *Physica*, vol. 1, pp. 104–113, 1934.",
    "KohnSham1965": "W. Kohn and L. J. Sham, “Self-consistent equations including exchange and correlation effects,” *Phys. Rev.*, vol. 140, pp. A1133–A1138, 1965.",
    "Kotochigova1997": "S. Kotochigova, Z. H. Levine, E. L. Shirley, M. D. Stiles, and C. W. Clark, “Local-density-functional calculations of the energy of atoms,” *Phys. Rev. A*, vol. 55, pp. 191–199, 1997.",
    "Chakravorty1993": "S. J. Chakravorty, S. R. Gwaltney, E. R. Davidson, F. A. Parpia, and C. Froese Fischer, “Ground-state correlation energies for atomic ions with 3 to 18 electrons,” *Phys. Rev. A*, vol. 47, pp. 3649–3670, 1993.",
    "Rodrigues2004": "G. C. Rodrigues, P. Indelicato, J. P. Santos, P. Patté, and F. Parente, “Systematic calculation of total atomic energies of ground state configurations,” *At. Data Nucl. Data Tables*, vol. 86, pp. 117–233, 2004, doi: 10.1016/j.adt.2003.11.005.",
    "Cowan1981": "R. D. Cowan, *The Theory of Atomic Structure and Spectra*. Berkeley, CA, USA: Univ. California Press, 1981.",
    "FroeseFischer1997": "C. Froese Fischer, T. Brage, and P. Jönsson, *Computational Atomic Structure: An MCHF Approach*. Bristol, U.K.: Inst. Phys. Publ., 1997.",
    "YerokhinShabaev2015": "V. A. Yerokhin and V. M. Shabaev, “Lamb shift of n = 1 and n = 2 states of hydrogen-like atoms, 1 ≤ Z ≤ 110,” *J. Phys. Chem. Ref. Data*, vol. 44, 033103, 2015.",
    "NISTASD": "A. Kramida, Yu. Ralchenko, J. Reader, and NIST ASD Team, *NIST Atomic Spectra Database* (version 5.12). Gaithersburg, MD, USA: Nat. Inst. Standards Technol., 2024. [Online]. Available: https://physics.nist.gov/asd (accessed on or before Oct. 5, 2026), doi: 10.18434/T4W30F.",
}


def main():
    s = io.open(P, encoding="utf-8").read()
    i_ref, i_app = s.index("## References"), s.index("## Appendix A")
    body, tail = s[:i_ref], s[i_app:]
    body = body.replace("[Layzer1959; LayzerEtAl1964; Layzer1967]", "[Layzer1959; LayzerEtAl1964]")
    assert "Layzer1967" not in body
    kp = "|".join(sorted(R, key=len, reverse=True))
    pat = re.compile(r"\[((?:%s)(?:;\s*(?:%s))*)\]" % (kp, kp))
    order = []
    for m in pat.finditer(body):
        for k in re.split(r";\s*", m.group(1)):
            if k not in order:
                order.append(k)
    num = {k: i + 1 for i, k in enumerate(order)}

    def fmt(ks):
        ns = sorted(num[k] for k in ks)
        out, i = [], 0
        while i < len(ns):
            j = i
            while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
                j += 1
            out.append(f"[{ns[i]}]–[{ns[j]}]" if j - i >= 2 else ", ".join(f"[{x}]" for x in ns[i:j + 1]))
            i = j + 1
        return ", ".join(out)
    body = pat.sub(lambda m: fmt(re.split(r";\s*", m.group(1))), body)
    left = re.findall(r"\[[A-Z][A-Za-z]+\d{4}[^\]]*\]", body)
    assert not left, left
    refs = "## References\n\n" + "\n".join(
        f"[{num[k]}] " + R[k].replace("@Pomarico2005@", f"[{num['Pomarico2005']}]") + "\n" for k in order) + "\n---\n\n"
    s = body + refs + tail
    s = s.replace("## Abstract\n\n", "**Abstract—**", 1).replace("**Keywords:**", "**Index Terms—**", 1)
    io.open(P, "w", encoding="utf-8", newline="\n").write(s)
    print(f"cited {len(order)} references; not cited (dropped): {[k for k in R if k not in num]}")


if __name__ == "__main__":
    main()
