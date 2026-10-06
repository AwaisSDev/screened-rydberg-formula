"""Reproduce numbers PRINTED in Mendoza et al., HEDP 7 (2011) 169 with model.py.

Printed values were transcribed from the author-deposited scan https://oa.upm.es/11165/2/INVE_MEM_2011_102088.pdf
(Tables 3, 4, 5, 6, 7, 8; NRSHM columns only). All energies in the paper are in atomic units (hartree), except
Table 8 (eV).

    py -3.13 models/benchmarks/mendoza2011/validate_paper.py   -> prints a report, writes
    results/bench_mendoza2011_validation.json
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
import model as M  # noqa: E402
from common.atomdata import HARTREE_EV, parse_shells, ground_shells  # noqa: E402

# ---- Table 3: NRSHM ionization energies (a.u.) of C-, Ne-, Al-, Ar-, Fe-like ions -------------------
T3 = {
    "C-like": (6, {6: 0.3638, 7: 1.0450, 8: 1.9763, 9: 3.1581, 10: 4.5905, 11: 6.2738, 12: 8.2084, 13: 10.3946,
                   14: 12.8330, 15: 15.5239, 16: 18.4680, 17: 21.6657, 18: 25.1179, 19: 28.8251, 20: 32.7882,
                   21: 37.0079, 22: 41.4851, 23: 46.2208, 24: 51.2159, 25: 56.4715, 26: 61.9886, 27: 67.7684,
                   28: 73.8122, 29: 80.1212, 30: 86.6968, 31: 93.5404, 32: 100.6534}),
    "Ne-like": (10, {10: 0.5594, 11: 1.5454, 12: 2.7810, 13: 4.2664, 14: 6.0015, 15: 7.9863, 16: 10.2208,
                     17: 12.7052, 18: 15.4394, 19: 18.4234, 20: 21.6575, 21: 25.1415, 22: 28.8757, 23: 32.8601,
                     24: 37.0948, 25: 41.5799, 26: 46.3156, 27: 51.3020, 28: 56.5392, 29: 62.0273, 30: 67.7666,
                     31: 73.7571, 32: 79.9992}),
    "Al-like": (13, {13: 0.1810, 14: 0.6108, 15: 1.1519, 16: 1.8043, 17: 2.5681, 18: 3.4435, 19: 4.4307,
                     20: 5.5298, 21: 6.7409, 22: 8.0643, 23: 9.5002, 24: 11.0489, 25: 12.7105, 26: 14.4854,
                     27: 16.3739, 28: 18.3763, 29: 20.4928, 30: 22.7240, 31: 25.0701, 32: 27.5315}),
    "Ar-like": (18, {18: 0.6332, 19: 1.2652, 20: 2.0083, 21: 2.8623, 22: 3.8275, 23: 4.9038, 24: 6.0912,
                     25: 7.3898, 26: 8.7997, 27: 10.3209, 28: 11.9535, 29: 13.6975, 30: 15.5530, 31: 17.5201,
                     32: 19.5989}),
    "Fe-like": (26, {26: 0.1616, 27: 0.2609, 28: 1.0733, 29: 1.9964, 30: 3.0302, 31: 4.1747, 32: 5.4298}),
}

# ---- Table 4: configuration-average energies (a.u.), NRSHM column, jj configurations as printed ---------
T4 = [  # (label, Z, {nlj: occupation}, E_NRSHM)
    ("N I", 7, "1s2.2s2.2p3", -54.45), ("N II", 7, "1s2.2s2.2p2", -54.03), ("N III", 7, "1s2.2s2.2p1", -52.98),
    ("N IV", 7, "1s2.2s2", -51.22), ("N V", 7, "1s2.2s1", -48.38), ("N VI", 7, "1s2", -44.78),
    ("O I", 8, "1s2.2s2.2p4", -74.94), ("O II", 8, "1s2.2s2.2p3", -74.46), ("O III", 8, "1s2.2s2.2p2", -73.28),
    ("O IV", 8, "1s2.2s2.2p1", -71.31), ("O V", 8, "1s2.2s2", -68.44), ("O VI", 8, "1s2.2s1", -64.26),
    ("O VII", 8, "1s2", -59.18),
]
# Printed jj configs (N I 2p1/2^2 2p3/2^1, O I 2p1/2^2 2p3/2^2, ...) equal the 'low' j-split of these.

# ---- Table 5: total energies of Na-like ions (a.u.), NRSHM column --------------------------------------
T5 = {18: -513.3, 20: -651.2, 23: -889.7, 26: -1166.4, 30: -1595.3, 34: -2093.9, 36: -2369.7, 47: -4211.7,
      50: -4812.8, 54: -5683.0, 55: -5913.1, 67: -9086.5, 79: -13097.0}

# ---- Table 6: Ne-like electron binding energies (a.u.), NRSHM X-alpha / HF / derivation-in-1/2 ----------
T6 = {16: {"1s1/2": (91.66, 96.53, 92.65), "2s1/2": (11.96, 12.84, 13.84), "2p1/2": (9.34, 10.37, 11.68),
           "2p3/2": (9.29, 10.31, 11.70)},
      26: {"1s1/2": (276.17, 284.25, 274.21), "2s1/2": (49.86, 51.45, 48.91), "2p1/2": (44.93, 46.84, 44.43),
           "2p3/2": (44.52, 46.41, 44.43)},
      66: {"1s1/2": (2128.97, 2152.13, 2139.7), "2s1/2": (481.99, 486.71, 494.98),
           "2p1/2": (466.20, 472.03, 482.68), "2p3/2": (438.28, 443.78, 453.26)},
      91: {"1s1/2": (4434.03, 4471.47, 4451.87), "2s1/2": (1049.37, 1056.65, 1071.31),
           "2p1/2": (1023.35, 1032.36, 1051.16), "2p3/2": (896.20, 904.13, 917.81)}}

# ---- Table 7: Au+25 (N=54, 4f8) level energies (a.u.), NRSHM column ------------------------------------
T7 = {"1s1/2": 2984.4, "2s1/2": 542.7, "2p1/2": 526.1, "2p3/2": 461.9, "3s1/2": 149.0, "3p1/2": 139.3,
      "3p3/2": 129.9, "3d3/2": 112.6, "3d5/2": 109.7, "4s1/2": 53.7, "4p1/2": 48.6, "4p3/2": 46.1,
      "4d3/2": 39.1, "4d5/2": 38.6, "4f5/2": 30.3, "4f7/2": 29.7}

# ---- Table 8: iron ions, binding energies (eV), NRSHM X-alpha and HF columns ----------------------------
T8 = {"Fe+2 3d6": (24, {"1s1/2": (6989.28, 7209.17), "2s1/2": (850.70, 893.45), "2p1/2": (737.21, 788.47),
                        "2p3/2": (704.61, 754.83), "3s1/2": (127.51, 141.43), "3p1/2": (86.48, 100.01),
                        "3p3/2": (83.62, 96.59)}),
      "Fe+6 3d2": (20, {"1s1/2": (7092.70, 7312.59), "2s1/2": (959.02, 1001.78), "2p1/2": (843.66, 894.96),
                        "2p3/2": (812.33, 862.62), "3s1/2": (219.63, 234.15), "3p1/2": (176.59, 191.32),
                        "3p3/2": (173.83, 188.04)}),
      "Fe+9 3p5": (17, {"1s1/2": (7195.87, 7415.76), "2s1/2": (1064.47, 1107.25), "2p1/2": (944.59, 995.98),
                        "2p3/2": (915.62, 966.05), "3s1/2": (298.66, 313.68)}),
      "Fe+12 3p2": (14, {"1s1/2": (7320.14, 7540.03), "2s1/2": (1187.43, 1230.29), "2p1/2": (1055.50, 1107.06),
                         "2p3/2": (1033.26, 1083.98), "3s1/2": (380.68, 396.34)})}


def madelung_shells_for(N):
    from common.atomdata import madelung_shells
    return madelung_shells(N)


def rel(a, b):
    return 100.0 * (a - b) / b


def hf_binding(Z, P, k):
    """App. B: eps_HF = eps_Xa - sigma_kk / r_k, r_k = nbar*sqrt(nbar^2 + a^2 Q^2)/Q (B.2-B.3), a.u.
    (alpha written 'a' in B.2). nbar = n - |kappa| + s, s = sqrt(kappa^2 - alpha^2 Q^2)."""
    import math
    n, l, j = M.KEYS[k]
    Q = M.charges(Z, P)[k]
    kap = j + 0.5  # |kappa|
    s = math.sqrt(kap * kap - (M.ALPHA * Q) ** 2)
    nbar = n - kap + s
    r = nbar * math.sqrt(nbar ** 2 + (M.ALPHA * Q) ** 2) / Q
    return M.binding_xalpha(Z, P, k) - M.SIG[k][k] / r


def main():
    report = {}
    # Table 3 -- primary IE prescription: total-energy difference with ground configurations
    t3 = {}
    allerr = []
    for seq, (N, vals) in T3.items():
        rows = []
        for Z, printed in vals.items():
            # configuration choice: (a) NIST ground configs for N and N-1
            ia = M.predict(Z, N) / HARTREE_EV
            # (b) Madelung (standard filling) configs for N and N-1
            ib = M.predict(Z, N, shells=madelung_shells_for(N), shells_final=madelung_shells_for(N - 1)) / HARTREE_EV
            rows.append({"Z": Z, "printed_au": printed, "nist_cfg_au": round(ia, 5), "madelung_cfg_au": round(ib, 5),
                         "dev_nist_%": round(rel(ia, printed), 4), "dev_madelung_%": round(rel(ib, printed), 4)})
            allerr.append((seq, Z, min(abs(rel(ia, printed)), abs(rel(ib, printed)))))
        t3[seq] = rows
    report["table3_IE_delta_E"] = t3
    # Table 4
    t4 = []
    for lab, Z, cfg, printed in T4:
        E = M.total_energy(Z, M.jsplit(parse_shells(cfg), "low"))
        t4.append({"ion": lab, "printed_au": printed, "model_au": round(E, 4), "dev_%": round(rel(E, printed), 4)})
    report["table4_total_energy"] = t4
    # Table 5
    t5 = []
    for Z, printed in T5.items():
        E = M.total_energy(Z, M.jsplit(parse_shells("1s2.2s2.2p6.3s1"), "low"))
        t5.append({"Z": Z, "printed_au": printed, "model_au": round(E, 3), "dev_%": round(rel(E, printed), 4)})
    report["table5_total_energy_Na_like"] = t5
    # Table 6
    t6 = []
    P = M.jsplit(parse_shells("1s2.2s2.2p6"), "low")
    for Z, d in T6.items():
        for lab, (xa, hf, half) in d.items():
            k = M.LABELS.index(lab)
            mxa = -M.binding_xalpha(Z, P, k)
            mhf = -hf_binding(Z, P, k)
            Ph = list(P); Ph[k] -= 0.5
            mhalf = -M.binding_xalpha(Z, Ph, k)
            t6.append({"Z": Z, "nlj": lab, "printed_xa": xa, "model_xa": round(mxa, 3), "dev_xa_%": round(rel(mxa, xa), 4),
                       "printed_hf": hf, "model_hf": round(mhf, 3), "dev_hf_%": round(rel(mhf, hf), 4),
                       "printed_half": half, "model_half": round(mhalf, 3), "dev_half_%": round(rel(mhalf, half), 4)})
    report["table6_binding_Ne_like"] = t6
    # Table 7 (Au+25, 4f8)
    t7 = []
    P = M.jsplit(ground_shells(79, 54), "low")
    for lab, printed in T7.items():
        k = M.LABELS.index(lab)
        mxa = -M.binding_xalpha(79, P, k)
        mhf = -hf_binding(79, P, k)
        t7.append({"nlj": lab, "printed": printed, "model_xa": round(mxa, 2), "dev_xa_%": round(rel(mxa, printed), 3),
                   "model_hf": round(mhf, 2), "dev_hf_%": round(rel(mhf, printed), 3)})
    report["table7_Au25"] = t7
    # Table 8 (eV)
    t8 = []
    for ion, (N, d) in T8.items():
        P = M.jsplit(ground_shells(26, N), "low")
        for lab, (xa, hf) in d.items():
            k = M.LABELS.index(lab)
            mxa = -M.binding_xalpha(26, P, k) * HARTREE_EV
            mhf = -hf_binding(26, P, k) * HARTREE_EV
            t8.append({"ion": ion, "nlj": lab, "printed_xa_eV": xa, "model_xa_eV": round(mxa, 2),
                       "dev_xa_%": round(rel(mxa, xa), 3), "printed_hf_eV": hf, "model_hf_eV": round(mhf, 2),
                       "dev_hf_%": round(rel(mhf, hf), 3)})
    report["table8_Fe_ions"] = t8

    # summaries
    def mx(rows, key):
        return max(abs(r[key]) for r in rows)
    summ = {}
    for seq, rows in t3.items():
        summ[f"T3 {seq} max|dev| NIST cfg %"] = mx(rows, "dev_nist_%")
        summ[f"T3 {seq} max|dev| Madelung cfg %"] = mx(rows, "dev_madelung_%")
    summ["T4 max|dev| %"] = mx(t4, "dev_%")
    summ["T5 max|dev| %"] = mx(t5, "dev_%")
    summ["T6 max|dev| Xalpha %"] = mx(t6, "dev_xa_%")
    summ["T6 max|dev| HF %"] = mx(t6, "dev_hf_%")
    summ["T6 max|dev| half %"] = mx(t6, "dev_half_%")
    summ["T7 max|dev| Xalpha %"] = mx(t7, "dev_xa_%")
    summ["T7 max|dev| HF %"] = mx(t7, "dev_hf_%")
    summ["T8 max|dev| Xalpha %"] = mx(t8, "dev_xa_%")
    summ["T8 max|dev| HF %"] = mx(t8, "dev_hf_%")
    report["summary"] = summ
    for k, v in summ.items():
        print(f"{k:40s} {v:9.4f}")
    with open(os.path.join(ROOT, "results", "bench_mendoza2011_validation.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
    return report


if __name__ == "__main__":
    main()
