"""Configuration helpers shared by all Track-B models.

A configuration ("shells") is a list of (n, l, occ) tuples, as produced by
common.atomdata.parse_shells / ground_shells.
"""
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from common.atomdata import madelung_shells, load_records  # noqa: E402

_TABLE = {(r["Z"], r["N"]): r["shells"] for r in load_records()}


def table_shells(Z, N):
    """NIST ground configuration if (Z, N) is in the table, else None (no IE values used)."""
    return _TABLE.get((Z, N))

CAP = {0: 2, 1: 6, 2: 10, 3: 14}


def normalize(shells):
    """Merge duplicates, drop empty subshells, sort by (n, l)."""
    d = {}
    for n, l, q in shells:
        d[(n, l)] = d.get((n, l), 0) + q
    return sorted((n, l, q) for (n, l), q in d.items() if q > 0)


def removed_subshell(shells_N, shells_Nm1):
    """Same rule as common.atomdata._removed_subshell (diff of the two configurations).
    Returns ((n,l), rearranged)."""
    init = {(n, l): q for n, l, q in shells_N}
    fin = {(n, l): q for n, l, q in shells_Nm1}
    diff = {k: init.get(k, 0) - fin.get(k, 0) for k in set(init) | set(fin)}
    lost = {k: d for k, d in diff.items() if d > 0}
    gained = {k: -d for k, d in diff.items() if d < 0}
    rearranged = bool(gained) or sum(lost.values()) != 1
    if not lost:
        return None, True
    k = max(lost, key=lambda k: (lost[k], k[0], k[1]))
    return k, rearranged


def remove_one(shells, nl):
    out = []
    for n, l, q in shells:
        if (n, l) == nl:
            q -= 1
        if q > 0:
            out.append((n, l, q))
    return out


def initial_final(Z, N, shells=None):
    """Return (shells_N, shells_Nm1, removed(n,l), rearranged).

    * shells is None  -> NIST ground configurations of the N- and (N-1)-electron ions;
      removed subshell by diffing them. If the N-electron ion is not tabulated its
      configuration is the Madelung one; if only the (N-1)-electron ion is missing, the
      electron is removed from the outermost subshell (largest n, then l).
    * shells given    -> the final configuration is shells minus one electron from the
      'outermost' subshell, using the simple configuration-only rule: largest n, ties ->
      largest l. (This reproduces the NIST removed subshell for 5522/5847 = 94.4% of the
      table; the exceptions are nf / (n+1)d vs outer s,p competition, e.g. 4f^k 5s2 5p6
      ions, where NIST removes 4f.)
    """
    if shells is None:
        t = table_shells(Z, N)
        sN = normalize(t if t is not None else madelung_shells(N))
        if N == 1:
            return sN, [], (sN[0][0], sN[0][1]), False
        tF = table_shells(Z, N - 1)
        if tF is None:   # final ion not tabulated: do not trust Madelung, remove outermost
            rem = max(((n, l) for n, l, q in sN), key=lambda k: (k[0], k[1]))
            return sN, remove_one(sN, rem), rem, False
        sF = normalize(tF)
        rem, rearr = removed_subshell(sN, sF)
        if rem is None:  # pathological table inconsistency -> fall back to outermost
            rem = max(((n, l) for n, l, q in sN), key=lambda k: (k[0], k[1]))
            sF = remove_one(sN, rem)
            rearr = True
        return sN, sF, rem, rearr
    sN = normalize(shells)
    rem = max(((n, l) for n, l, q in sN), key=lambda k: (k[0], k[1]))
    return sN, remove_one(sN, rem), rem, False


def hund_pairs(l, k):
    """Number of parallel-spin electron pairs in l^k under Hund's first rule."""
    m = 2 * l + 1
    up = min(k, m)
    dn = k - up
    return up * (up - 1) // 2 + dn * (dn - 1) // 2


def exchange_kink(l, k):
    """Centred change in Hund parallel pairs when one electron is removed from l^k:
    [P(k) - P(k-1)] - <same for a statistical (configuration-average) distribution>.
    For p: k=1..6 -> 0, .6, 1.2, -1.2, -.6, 0. Positive = extra exchange stabilization
    lost on ionization (IE raised), negative = IE lowered (p4, d6, f8 ...)."""
    if l == 0 or k < 1:
        return 0.0
    hund = hund_pairs(l, k) - hund_pairs(l, k - 1)
    avg = 2.0 * (k - 1) * l / (4 * l + 1)
    return hund - avg


def jj_j(l, k):
    """j of the electron removed from l^k when filled in jj coupling
    (l-1/2 first: 2l electrons; then l+1/2). Returns (j, kappa_sign)."""
    if l == 0:
        return 0.5
    return (l - 0.5) if k <= 2 * l else (l + 0.5)
