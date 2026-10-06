"""Delta-SCF (LSDA) and Slater transition-state ionization energies with ks_atom.py.

Rows computed: all ions with Z <= ZMAX_ALL plus neutral-atom first IEs for Z <= ZMAX_NEUTRAL.
Writes results/fp_dft_predictions.csv incrementally (columns Z,N,IE_pred_eV,IE_TS_eV,IE_HOMO_eV).
Usage: py -3.13 run_dft.py ZMAX_ALL ZMAX_NEUTRAL
"""
import csv
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
from common.atomdata import load_records, ground_shells, shells_to_str, HARTREE_EV  # noqa: E402
import ks_atom as K  # noqa: E402

H_STEP, TOL = 0.009, 1e-7
_E = {}


def energy(Z, shells):
    key = (Z, shells_to_str(shells) if shells else "-")
    if key not in _E:
        if not shells:
            _E[key] = (0.0, None)
        else:
            occ = K.occupations(shells, spin=True)
            a = K.Atom(Z, occ, spin=True, grid=K.Grid(Z, h=H_STEP))
            E = a.scf(tol=TOL)
            _E[key] = (E, a)
    return _E[key]


def transition_state(Z, shells, removed):
    n, l = removed
    q = dict(((a, b), c) for a, b, c in shells)[(n, l)]
    s = 1 if q > 2 * l + 1 else 0
    occ = K.occupations(shells, spin=True, remove=(n, l), remove_amount=0.5)
    a = K.Atom(Z, occ, spin=True, grid=K.Grid(Z, h=H_STEP))
    a.scf(tol=TOL)
    return -a.orbs[(n, l, s)][0]


def main(zall, zneu):
    recs = [r for r in load_records() if r["Z"] <= zall or (r["N"] == r["Z"] and r["Z"] <= zneu)]
    recs.sort(key=lambda r: (r["Z"] > zall, r["Z"], -r["N"]))
    out = os.path.join(ROOT, "results", "fp_dft_predictions.csv")
    f = open(out, "w", newline="", encoding="utf-8")
    w = csv.writer(f)
    w.writerow(["Z", "N", "IE_pred_eV", "IE_TS_eV", "IE_HOMO_eV", "converged"])
    t0 = time.time()
    for r in recs:
        Z, N = r["Z"], r["N"]
        sh = r["shells"]
        shi = ground_shells(Z, N - 1) if N > 1 else []
        try:
            E1, a1 = energy(Z, sh)
            E0, a0 = energy(Z, shi)
            ie = (E0 - E1) * HARTREE_EV
            n, l = r["removed"]
            try:
                ts = transition_state(Z, sh, (n, l)) * HARTREE_EV
            except Exception:
                ts = float("nan")
            q = dict(((x, y), c) for x, y, c in sh)[(n, l)]
            s = 1 if q > 2 * l + 1 else 0
            homo = -a1.orbs[(n, l, s)][0] * HARTREE_EV
            conv = int(a1.converged and (a0 is None or a0.converged))
        except Exception as e:  # noqa: BLE001
            print("FAIL", Z, N, e, flush=True)
            continue
        w.writerow([Z, N, f"{ie:.6f}", f"{ts:.6f}", f"{homo:.6f}", conv])
        f.flush()
        print(Z, N, f"{ie:.4f} ref {r['IE_eV']:.4f} TS {ts:.4f}  t={time.time() - t0:.0f}s", flush=True)
    f.close()


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
