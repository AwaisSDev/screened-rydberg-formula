"""Run the Kregar/Di Rocco SHM on all 5847 NIST rows (no IE values read).
usage: py -3.13 run_lit.py part <k> <nparts>   -> results/lit_kregar_part<k>.json
       py -3.13 run_lit.py merge <nparts>      -> results/lit_kregar_{nr,pauli,dirac}_predictions.csv"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
import kregar_shm as K  # noqa: E402
from common.atomdata import load_records, HARTREE_EV  # noqa: E402

RES = os.path.join(ROOT, "results")


def part(k, nparts):
    Zs = sorted({r["Z"] for r in load_records()}, reverse=True)
    mine = Zs[k::nparts]
    out = {}
    for Z in mine:
        e = K._iso(Z)
        out[Z] = {N: list(v) for N, v in e.items()}
        print(Z, flush=True)
    with open(os.path.join(RES, f"lit_kregar_part{k}.json"), "w") as f:
        json.dump(out, f)


def merge(nparts):
    E = {}
    for k in range(nparts):
        with open(os.path.join(RES, f"lit_kregar_part{k}.json")) as f:
            for Z, d in json.load(f).items():
                E[int(Z)] = {int(N): v for N, v in d.items()}
    for name, idx in K.VARIANTS.items():
        path = os.path.join(RES, f"lit_kregar_{name}_predictions.csv")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Z", "N", "IE_pred_eV"])
            for r in load_records():
                Z, N = r["Z"], r["N"]
                w.writerow([Z, N, (E[Z][N - 1][idx] - E[Z][N][idx]) * HARTREE_EV])
        print("wrote", path)


if __name__ == "__main__":
    if sys.argv[1] == "part":
        part(int(sys.argv[2]), int(sys.argv[3]))
    else:
        merge(int(sys.argv[2]))
