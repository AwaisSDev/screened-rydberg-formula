"""Clean the raw NIST ASD ionization-energy export into a tidy CSV.

Source: NIST Atomic Spectra Database, Ionization Energies form
(https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html), all spectra H-Ds, units eV.

NIST notation: a value in ( ) is theoretical, in [ ] is semi-empirical/interpolated,
bare is experimental.
"""
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "nist_ie_raw.csv")
OUT = os.path.join(HERE, "nist_ie.csv")


def unwrap(s):
    # NIST CSV wraps every field as ="value"
    s = s.strip()
    m = re.fullmatch(r'="(.*)"', s)
    return m.group(1) if m else s


def parse_config(cfg):
    """'1s2.2s2.2p6.3s' -> list of (n, l_letter, occupancy). Handles '[Ne].3s' style cores."""
    cores = {
        "He": "1s2", "Ne": "1s2.2s2.2p6", "Ar": "1s2.2s2.2p6.3s2.3p6",
        "Kr": "1s2.2s2.2p6.3s2.3p6.3d10.4s2.4p6",
        "Xe": "1s2.2s2.2p6.3s2.3p6.3d10.4s2.4p6.4d10.5s2.5p6",
        "Rn": "1s2.2s2.2p6.3s2.3p6.3d10.4s2.4p6.4d10.5s2.5p6.4f14.5d10.6s2.6p6",
        "Cd": "1s2.2s2.2p6.3s2.3p6.3d10.4s2.4p6.4d10.5s2",
        "Hg": "1s2.2s2.2p6.3s2.3p6.3d10.4s2.4p6.4d10.5s2.5p6.4f14.5d10.6s2",
    }
    cfg = cfg.strip()
    for k, v in cores.items():
        cfg = cfg.replace(f"[{k}]", v)
    out = []
    for tok in cfg.split("."):
        tok = tok.strip()
        if not tok:
            continue
        m = re.fullmatch(r"(\d+)([spdfghik])(\d*)", tok)
        if not m:
            return None
        out.append((int(m.group(1)), m.group(2), int(m.group(3) or 1)))
    return out


def main():
    rows = []
    with open(RAW, newline="", encoding="utf-8") as f:
        rd = csv.reader(f)
        header = [unwrap(h) for h in next(rd)]
        for rec in rd:
            if not rec or not rec[0].strip():
                continue
            d = dict(zip(header, [unwrap(x) for x in rec]))
            ie = d.get("Ionization Energy (b) (eV)", "").strip()
            if not ie:
                continue
            ie_clean = re.sub(r"[^\d.eE+-]", "", ie)
            try:
                ie_val = float(ie_clean)
            except ValueError:
                continue
            prefix = d.get("Prefix", "").strip()
            status = {"(": "theoretical", "[": "semi-empirical"}.get(prefix, "experimental")
            Z = int(d["At. num"])
            q = int(d["Ion Charge"].replace("+", ""))
            N = Z - q
            unc = d.get("Uncertainty (c) (eV)", "").strip()
            shells = d.get("Ground Shells (a)", "").strip()
            conf = d.get("Ground Config.", "").strip()
            parsed = parse_config(shells) if shells else None
            if parsed and sum(p[2] for p in parsed) != N:
                raise ValueError(f"electron count mismatch Z={Z} q={q}: {shells}")
            if parsed:
                n_out = max(p[0] for p in parsed)
                last = parsed[-1]
                last_sub = f"{last[0]}{last[1]}{last[2]}"
            else:
                n_out, last_sub = "", ""
            rows.append({
                "Z": Z, "symbol": d["Sp. Name"].split()[0], "ion_charge": q, "N": N,
                "isoelectronic_seq": d.get("Isoel. Seq.", ""),
                "ground_shells": shells, "ground_config": conf,
                "shells_expanded": ".".join(f"{a}{b}{c}" for a, b, c in parsed) if parsed else "",
                "ground_level": d.get("Ground Level", ""),
                "ionized_level": d.get("Ionized Level", ""),
                "IE_eV": ie_val, "unc_eV": unc, "status": status,
                "n_max": n_out, "last_subshell": last_sub,
            })
    rows.sort(key=lambda r: (r["Z"], r["ion_charge"]))
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    n_exp = sum(r["status"] == "experimental" for r in rows)
    print(f"wrote {len(rows)} rows ({n_exp} experimental) -> {OUT}")


if __name__ == "__main__":
    main()
