"""Check the Kregar/Di Rocco SHM implementation against numbers PRINTED in Pomarico, Iriarte & Di Rocco,
Braz. J. Phys. 35, 130 (2005): Table 6 (Z -> infinity screening parameters, Kr-like), Table 1 (total
energies, hartree, incl. relativistic corr.), Tables 2-3 (Ar isonuclear / isoelectronic IEs, eV).
Writes results/lit_kregar_validation.json."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kregar_shm as K  # noqa: E402

KR = [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 6), (3, 2, 10), (4, 0, 2), (4, 1, 6)]
LAB = ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"]
# Table 6 of Pomarico et al. 2005 (f_ji above diagonal, g_ij below, k_ii diagonal), transcribed from the PDF
T6 = np.array([
    [0.3125, 0.6924, 0.8776, 0.8095, 0.9299, 0.9967, 0.8636, 0.9539],
    [0.0258, 0.3008, 0.2230, 0.6055, 0.6466, 0.7756, 0.7290, 0.7630],
    [0.0149, 0.3668, 0.3492, 0.6617, 0.6821, 0.8452, 0.7650, 0.7887],
    [0.0068, 0.0528, 0.0471, 0.2988, 0.2435, 0.2164, 0.5512, 0.5737],
    [0.0034, 0.0674, 0.0635, 0.3121, 0.3104, 0.2561, 0.5861, 0.5914],
    [0.0002, 0.0562, 0.0368, 0.4213, 0.4066, 0.3765, 0.6281, 0.6391],
    [0.0026, 0.0189, 0.0167, 0.0776, 0.0746, 0.0678, 0.2982, 0.2545],
    [0.0011, 0.0238, 0.0223, 0.0900, 0.0855, 0.0784, 0.2889, 0.2987]])

TAB1_HT = {"He": (2, 2.85), "Be": (4, 14.59), "Ne": (10, 128.44), "Mg": (12, 200.06), "Ar": (18, 528.78),
           "Ca": (20, 680.24), "Zn": (30, 1793.04), "Kr": (36, 2787.80), "Sr": (38, 3177.06),
           "Cd": (48, 5580.65), "Xe": (54, 7423.77)}
TAB2_AR = [14.72, 27.34, 41.20, 56.22, 72.34, 89.55, 127.75, 147.15, 418.49, 480.70, 545.47, 612.81,
           682.70, 755.17, 850.34, 918.27, 4121.25, 4427.30]      # Z-N+1 = 1..18
TAB3_ISO = {18: 14.73, 19: 30.74, 20: 50.16, 21: 72.86, 24: 159.93, 25: 195.16}


def main():
    out = {}
    one = np.ones(len(KR))
    S, F0 = K._monopole(KR, one, K._orbitals(KR, one))
    eps = K.exchange_fraction(KR, one)
    M = S * (1 - eps)
    # Our M[i,j] = screening of i by one electron of j.  Table: above diagonal f_ji (row=j inner, col=i outer)
    # -> our M[col,row]; below diagonal g_ij likewise M[col,row] (checked: 1s-2s f=0.69, g=0.026).
    ours = M.T          # Table6[row, col] = screening of electron "col" by one electron of "row"
    dev = np.abs(ours - T6)
    out["table6_max_abs_dev"] = float(dev.max())
    out["table6_mean_abs_dev"] = float(dev.mean())
    out["table6_worst"] = [LAB[i] + "," + LAB[j] for i, j in zip(*np.unravel_index(np.argsort(-dev, axis=None)[:5], dev.shape))]
    print("Table 6: max |dev| %.4f  mean %.4f" % (dev.max(), dev.mean()), out["table6_worst"])
    np.set_printoptions(precision=4, suppress=True, linewidth=150)
    print(ours)

    import evaluate  # noqa: F401  (only for path sanity)
    t1 = {}
    for el, (Z, ref) in TAB1_HT.items():
        e = K.energies(Z, Z)
        t1[el] = {"ref_Ht": ref, "nr": -e[0], "pauli": -e[1], "dirac": -e[2]}
        print("Tab1 %-2s ref %9.2f  nr %9.2f pauli %9.2f dirac %9.2f" % (el, ref, -e[0], -e[1], -e[2]))
    out["table1_total_energy_Ht"] = t1
    t2 = []
    for zn1, ref in enumerate(TAB2_AR, start=1):
        N = 18 - zn1 + 1
        t2.append({"N": N, "ref_eV": ref, "pauli": K.predict(18, N, "pauli"), "nr": K.predict(18, N, "nr")})
    out["table2_Ar_isonuclear"] = t2
    t3 = [{"Z": Z, "ref_eV": ref, "pauli": K.predict(Z, 18, "pauli"), "nr": K.predict(Z, 18, "nr")}
          for Z, ref in TAB3_ISO.items()]
    out["table3_Ar_isoelectronic"] = t3
    for row in t2 + t3:
        print(row)
    rel2 = [abs(r["pauli"] - r["ref_eV"]) / r["ref_eV"] * 100 for r in t2]
    rel3 = [abs(r["pauli"] - r["ref_eV"]) / r["ref_eV"] * 100 for r in t3]
    out["table2_mean_rel_dev_%"] = float(np.mean(rel2))
    out["table3_mean_rel_dev_%"] = float(np.mean(rel3))
    print("mean rel dev vs printed values: Tab2 %.2f%%  Tab3 %.2f%%" % (np.mean(rel2), np.mean(rel3)))
    root = os.path.dirname(os.path.dirname(HERE))
    with open(os.path.join(root, "results", "lit_kregar_validation.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    main()
