"""Coverage test of the public API (ionization.py) for every Z = 1..118 and every N = 1..Z.

    py -3.13 tests/test_coverage.py [--model pa_hier_rel] [--json results/uni_coverage.json]

Asserts every predicted IE is finite and > 0 (exit code 1 on any crash or non-finite / non-positive value).
Counts and lists monotonicity violations: removing one more electron should cost more energy, i.e.
IE(Z, N-1) > IE(Z, N); a violation is IE(Z, N-1) <= IE(Z, N). Violations are reported, not failures.
Also cross-checks that the scalar API ionization_energy(Z, N) equals the vectorised
successive_ionization_energies(Z) on a sample.
"""
import argparse
import json
import math
import os
import sys
import time
import traceback

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
from ionization import ionization_energy, successive_ionization_energies, SYMBOLS  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=None)
    ap.add_argument("--zmax", type=int, default=118)
    ap.add_argument("--json", default=os.path.join(_ROOT, "results", "uni_coverage.json"))
    args = ap.parse_args()
    t0 = time.time()
    crashes, bad, viol, n_ok = [], [], [], 0
    table = {}
    for Z in range(1, args.zmax + 1):
        try:
            ies = successive_ionization_energies(Z, model=args.model)   # [IE(N=Z), IE(N=Z-1), ..., IE(N=1)]
        except Exception as e:                                          # noqa: BLE001
            crashes.append({"Z": Z, "error": repr(e), "trace": traceback.format_exc(limit=3)})
            continue
        if len(ies) != Z:
            crashes.append({"Z": Z, "error": f"expected {Z} values, got {len(ies)}"})
            continue
        byN = {Z - i: v for i, v in enumerate(ies)}
        table[Z] = byN
        for N, v in byN.items():
            if not (isinstance(v, float) and math.isfinite(v) and v > 0):
                bad.append({"Z": Z, "N": N, "IE": repr(v)})
            else:
                n_ok += 1
        for N in range(2, Z + 1):
            a, b = byN[N - 1], byN[N]                                   # IE(Z, N-1) must exceed IE(Z, N)
            if math.isfinite(a) and math.isfinite(b) and a <= b:
                viol.append({"Z": Z, "sym": SYMBOLS[Z - 1], "N": N, "IE_N": b, "IE_N_minus_1": a,
                             "ratio": a / b})
    # scalar API cross-check on a sample
    mism = []
    for Z in (1, 2, 6, 26, 54, 55, 79, 92, 110, 111, 118):
        if Z > args.zmax or Z not in table:
            continue
        for N in sorted(n for n in {1, 2, max(1, Z // 2), Z} if n <= Z):
            try:
                v = ionization_energy(Z, N, model=args.model)
            except Exception as e:                                      # noqa: BLE001
                crashes.append({"Z": Z, "N": N, "error": "scalar API: " + repr(e)})
                continue
            if abs(v / table[Z][N] - 1) > 1e-9:
                mism.append({"Z": Z, "N": N, "scalar": v, "vector": table[Z][N]})
    total = args.zmax * (args.zmax + 1) // 2
    print(f"model: {args.model or 'default (pa_hier_rel)'}; Z = 1..{args.zmax}; {total} (Z, N) cases; {time.time()-t0:.1f}s")
    print(f"finite and > 0: {n_ok}/{total}; crashes: {len(crashes)}; non-finite/non-positive: {len(bad)}; "
          f"scalar/vector mismatches: {len(mism)}")
    print(f"monotonicity violations (IE(Z,N-1) <= IE(Z,N)): {len(viol)}")
    for v in viol:
        print(f"  {v['sym']:>2s} Z={v['Z']:3d}  IE(N={v['N']-1:3d}) = {v['IE_N_minus_1']:12.4f} <= IE(N={v['N']:3d}) = "
              f"{v['IE_N']:12.4f} eV  (ratio {v['ratio']:.4f}){'   [Z>110: outside NIST table]' if v['Z'] > 110 else ''}")
    for c in crashes:
        print("CRASH", c)
    for b in bad[:50]:
        print("BAD", b)
    for m in mism:
        print("MISMATCH", m)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump({"model": args.model or "pa_hier_rel", "zmax": args.zmax, "n_cases": total, "n_ok": n_ok,
                       "crashes": crashes, "bad": bad, "scalar_vector_mismatches": mism,
                       "n_monotonicity_violations": len(viol), "monotonicity_violations": viol}, f, indent=1)
    ok = not crashes and not bad and n_ok == total
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
