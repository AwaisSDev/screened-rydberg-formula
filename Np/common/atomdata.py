"""Shared data access + constants for all ionization-energy models.

Every model in this project is scored against the same table, data/nist_ie.csv
(NIST ASD, 5847 successive ionization energies, Z = 1..110), through evaluate.py.

Conventions
-----------
Z  : nuclear charge
N  : number of electrons in the ion BEFORE ionization (N = Z - ion_charge)
IE(Z, N) = E_total(Z, N-1) - E_total(Z, N)   [eV], always > 0
A "subshell" is a tuple (n, l, occ) with l an int (s=0, p=1, d=2, f=3).
"""
import csv
import os
from functools import lru_cache

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_CSV = os.path.join(ROOT, "data", "nist_ie.csv")
RESULTS_DIR = os.path.join(ROOT, "results")

# CODATA 2018
HARTREE_EV = 27.211386245988
RYDBERG_EV = 13.605693122994          # R_inf * h c
ALPHA = 7.2973525693e-3               # fine-structure constant
ELECTRON_MASS_U = 5.48579909065e-4    # electron mass in u

L_LETTER = "spdfghik"


def parse_shells(s):
    """'1s2.2s2.2p1' -> [(1,0,2),(2,0,2),(2,1,1)]"""
    out = []
    for tok in s.split("."):
        tok = tok.strip()
        if not tok:
            continue
        i = 0
        while tok[i].isdigit():
            i += 1
        n = int(tok[:i])
        l = L_LETTER.index(tok[i])
        occ = int(tok[i + 1:] or 1)
        out.append((n, l, occ))
    return out


def shells_to_str(shells):
    return ".".join(f"{n}{L_LETTER[l]}{q}" for n, l, q in shells)


@lru_cache(maxsize=1)
def load_records():
    """List of dicts, one per (Z, N), with parsed fields. Cached."""
    recs = []
    with open(DATA_CSV, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            Z, q = int(r["Z"]), int(r["ion_charge"])
            recs.append({
                "Z": Z,
                "N": Z - q,
                "ion_charge": q,
                "symbol": r["symbol"],
                "isoelectronic_seq": r["isoelectronic_seq"],
                "shells": parse_shells(r["shells_expanded"]),
                "ground_level": r["ground_level"],
                "IE_eV": float(r["IE_eV"]),
                "unc_eV": float(r["unc_eV"]) if r["unc_eV"] else None,
                "status": r["status"],          # experimental | semi-empirical | theoretical
            })
    by_key = {(x["Z"], x["N"]): x for x in recs}
    for x in recs:
        x["removed"], x["rearranged"] = _removed_subshell(x, by_key.get((x["Z"], x["N"] - 1)))
    return tuple(recs)


def _removed_subshell(rec, final):
    """Which (n,l) loses an electron going from the N-electron ground config to the
    (N-1)-electron ground config. If the ground configurations differ by more than one
    electron removal (e.g. Ni I 3d8 4s2 -> Ni II 3d9), rearranged=True and the subshell
    returned is the one that lost the most electrons (ties -> outermost)."""
    init = {(n, l): q for n, l, q in rec["shells"]}
    fin = {} if final is None else {(n, l): q for n, l, q in final["shells"]}
    diff = {k: init.get(k, 0) - fin.get(k, 0) for k in set(init) | set(fin)}
    lost = {k: d for k, d in diff.items() if d > 0}
    gained = {k: -d for k, d in diff.items() if d < 0}
    rearranged = bool(gained) or sum(lost.values()) != 1
    if not lost:
        return None, True
    k = max(lost, key=lambda k: (lost[k], k[0], k[1]))
    return k, rearranged


def get(Z, N):
    for r in load_records():
        if r["Z"] == Z and r["N"] == N:
            return r
    raise KeyError((Z, N))


def madelung_shells(N):
    """Aufbau / Madelung (n+l, then n) filling of N electrons. Used as the configuration
    fallback when a model is asked about an (Z, N) that is not in the NIST table."""
    order = sorted(((n, l) for n in range(1, 9) for l in range(min(n, 4))),
                   key=lambda nl: (nl[0] + nl[1], nl[0]))
    out, left = [], N
    for n, l in order:
        if left <= 0:
            break
        q = min(left, 2 * (2 * l + 1))
        out.append((n, l, q))
        left -= q
    return sorted(out)


def ground_shells(Z, N):
    """NIST ground configuration of the N-electron ion of element Z, else Madelung."""
    try:
        return get(Z, N)["shells"]
    except KeyError:
        return madelung_shells(N)


if __name__ == "__main__":
    recs = load_records()
    print(len(recs), "records")
    print("rearranged:", sum(r["rearranged"] for r in recs))
    for Z, N in [(2, 2), (26, 26), (28, 28), (29, 29), (64, 64), (92, 1)]:
        r = get(Z, N)
        print(Z, N, r["symbol"], shells_to_str(r["shells"]), "removed", r["removed"],
              "rearranged", r["rearranged"], r["IE_eV"], r["status"])
