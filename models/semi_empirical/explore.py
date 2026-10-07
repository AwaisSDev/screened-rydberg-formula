"""Pattern-recognition / data-exploration figures for Track B.

Run:  py -3.13 models/semi_empirical/explore.py
Writes results/figures/se_explore_*.png and results/se_moseley_fits.csv
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _ROOT)
sys.path.insert(0, _HERE)
from common.atomdata import load_records, RYDBERG_EV, ALPHA, L_LETTER  # noqa: E402

FIG = os.path.join(_ROOT, "results", "figures")
RY = RYDBERG_EV
R = load_records()
D = {(r["Z"], r["N"]): r for r in R}


def fig_ie_vs_N():
    fig, ax = plt.subplots(figsize=(8, 5))
    for Z in (18, 36, 54, 80, 103):
        Ns = sorted(N for (z, N) in D if z == Z)
        ax.semilogy(Ns, [D[(Z, N)]["IE_eV"] for N in Ns], ".-", label=f"Z={Z}")
    for N in (2, 10, 18, 28, 36, 46, 54, 68, 78, 86):
        ax.axvline(N + 0.5, color="0.85", lw=0.8, zorder=0)
    ax.set_xlabel("N (electrons before ionization)")
    ax.set_ylabel("IE (eV)")
    ax.set_title("Successive IEs at fixed Z: jumps at closed shells N=2,10,18,28,36,...")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "se_explore_ie_vs_N.png"), dpi=130)
    plt.close(fig)


def moseley_fits():
    """y = n*sqrt(IE/Ry) vs Z along each isoelectronic sequence N (removed n from data).
    Linear fit y = a Z + b on the 'non-relativistic, ionized' window q>=2, Z<=min(N+30,60)."""
    rows = []
    for N in range(1, 111):
        pts = [(r["Z"], r) for r in R if r["N"] == N]
        if len(pts) < 4:
            continue
        sel = [(Z, r) for Z, r in pts if Z - N >= 2 and Z <= max(N + 30, 30) and Z <= 60]
        if len(sel) < 4:
            continue
        Zs = np.array([Z for Z, _ in sel], float)
        n = np.array([r["removed"][0] for _, r in sel], float)
        y = n * np.sqrt(np.array([r["IE_eV"] for _, r in sel]) / RY)
        a, b = np.polyfit(Zs, y, 1)
        resid = y - (a * Zs + b)
        rem = sel[0][1]["removed"]
        rows.append(dict(N=N, n=rem[0], l=rem[1], slope=a, intercept=b, sigma=-b / a,
                         n_eff=rem[0] / a, rms=float(np.sqrt(np.mean(resid ** 2))), npts=len(sel)))
    import csv
    with open(os.path.join(_ROOT, "results", "se_moseley_fits.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def fig_moseley(rows):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    ax = axs[0]
    for N in (1, 2, 3, 10, 11, 18, 19, 29, 36, 37, 55):
        pts = sorted((r["Z"], r) for r in R if r["N"] == N)
        Zs = np.array([p[0] for p in pts])
        y = np.array([p[1]["removed"][0] * np.sqrt(p[1]["IE_eV"] / RY) for p in pts])
        ax.plot(Zs, y, "-", lw=1.2, label=f"N={N}")
    ax.plot([0, 110], [0, 110], "k:", lw=0.8, label="y = Z")
    ax.set_xlabel("Z")
    ax.set_ylabel(r"$n\sqrt{IE/\mathrm{Ry}}$  (= $Z_{\rm eff}$)")
    ax.set_title("Moseley plot along isoelectronic sequences")
    ax.legend(fontsize=7, ncol=2)
    Ns = np.array([r["N"] for r in rows])
    col = [f"C{r['l']}" for r in rows]
    axs[1].scatter(Ns, [r["sigma"] for r in rows], c=col, s=14)
    axs[1].plot([0, 110], [0, 110], "k:", lw=0.8)
    axs[1].set_xlabel("N")
    axs[1].set_ylabel(r"screening $\sigma_N$ = -intercept/slope")
    axs[1].set_title(r"$\sigma_N \approx N-1-p$ (colour = l of removed e$^-$: s,p,d,f)")
    axs[2].scatter(Ns, [r["slope"] for r in rows], c=col, s=14)
    axs[2].axhline(1, color="k", lw=0.8, ls=":")
    axs[2].set_xlabel("N")
    axs[2].set_ylabel("slope a = n / n_eff")
    axs[2].set_title("Slope ~1: hydrogenic n is the right asymptotic n")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "se_explore_moseley.png"), dpi=130)
    plt.close(fig)
    # penetration per electron of the inner shells: p = N-1-sigma
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter(Ns, [r["N"] - 1 - r["sigma"] for r in rows], c=col, s=14)
    ax.set_xlabel("N")
    ax.set_ylabel(r"$p_N = N-1-\sigma_N$ (unscreened charge)")
    ax.set_title("Penetration charge p_N; sawtooth = same-shell screening < 1")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "se_explore_penetration.png"), dpi=130)
    plt.close(fig)


def fig_kinks():
    """Exchange kinks: IE/(q+1)^2 across p, d, f filling, at several ionization stages."""
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
    blocks = [("2p", 5, 10, axs[0]), ("3d", 19, 28, axs[1]), ("4f", 57, 70, axs[2])]
    for name, N0, N1, ax in blocks:
        for q in (0, 1, 2, 3, 5, 8):
            xs, ys = [], []
            for N in range(N0, N1 + 1):
                r = D.get((N + q, N))
                if r is None:
                    continue
                xs.append(N)
                ys.append(r["IE_eV"] / (q + 1) ** 2)
            ax.plot(xs, ys, "o-", ms=3, label=f"q={q}")
        ax.set_xlabel("N")
        ax.set_ylabel("IE / (q+1)^2  (eV)")
        ax.set_title(f"{name} filling: half-filled kink")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "se_explore_kinks.png"), dpi=130)
    plt.close(fig)


def fig_neutral():
    Zs = sorted(Z for (Z, N) in D if Z == N)
    fig, ax = plt.subplots(figsize=(10, 4.3))
    ax.plot(Zs, [D[(Z, Z)]["IE_eV"] for Z in Zs], "k.-")
    for Z in Zs:
        r = D[(Z, Z)]
        if Z in (2, 3, 7, 8, 10, 11, 15, 16, 18, 19, 24, 29, 30, 36, 37, 54, 55, 57, 79, 80, 86, 87):
            ax.annotate(r["symbol"], (Z, r["IE_eV"]), fontsize=7, xytext=(2, 3),
                        textcoords="offset points")
    ax.set_xlabel("Z")
    ax.set_ylabel("first IE (eV)")
    ax.set_title("Neutral-atom first IEs: alkali minima, noble-gas maxima, N/O & P/S kinks")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "se_explore_neutral.png"), dpi=130)
    plt.close(fig)


def fig_relativity(rows):
    """Ratio of NIST IE to the non-relativistic Moseley line, vs (Z alpha)^2."""
    fit = {r["N"]: r for r in rows}
    fig, ax = plt.subplots(figsize=(7, 4.8))
    Zs = np.arange(1, 111)
    for N in (1, 2, 3, 10, 11, 28, 29, 46, 47):
        if N not in fit:
            continue
        a, b = fit[N]["slope"], fit[N]["intercept"]
        pts = sorted((r["Z"], r) for r in R if r["N"] == N and r["Z"] - N >= 2)
        x = np.array([(p[0] * ALPHA) ** 2 for p in pts])
        n = pts[0][1]["removed"][0]
        y = np.array([p[1]["IE_eV"] / (RY * ((a * p[0] + b) / n) ** 2) for p in pts])
        ax.plot(x, y, "-", label=f"N={N} ({n}{L_LETTER[pts[0][1]['removed'][1]]})")
    x = (Zs * ALPHA) ** 2
    ax.plot(x, 1 / np.sqrt(1 - x) * 0 + (1 - np.sqrt(1 - x)) * 2 / x, "k:", label="Dirac 1s factor")
    ax.set_xlabel(r"$(Z\alpha)^2$")
    ax.set_ylabel("IE(NIST) / non-rel. Moseley line")
    ax.set_title("Relativistic growth along isoelectronic sequences")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "se_explore_relativity.png"), dpi=130)
    plt.close(fig)


def fig_removal_map():
    fig, ax = plt.subplots(figsize=(8, 7))
    for l in range(4):
        pts = [(r["Z"], r["N"]) for r in R if r["removed"] and r["removed"][1] == l]
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=3, c=f"C{l}",
                   label=f"removed l={L_LETTER[l]}")
    pts = [(r["Z"], r["N"]) for r in R if r["rearranged"]]
    ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=18, facecolors="none",
               edgecolors="k", label="rearranged (63)")
    ax.set_xlabel("Z")
    ax.set_ylabel("N")
    ax.set_title("Which subshell is ionized (4s/3d, 6s/5d/4f competition near N=Z)")
    ax.legend(markerscale=3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "se_explore_removal_map.png"), dpi=130)
    plt.close(fig)


def neutral_quantum_defects():
    """Effective quantum number n* = sqrt(Ry/IE) of neutral atoms (asymptotic charge 1)."""
    out = []
    for Z in range(1, 111):
        r = D.get((Z, Z))
        if r:
            out.append((Z, r["symbol"], r["removed"], np.sqrt(RY / r["IE_eV"]),
                        r["removed"][0] * np.sqrt(r["IE_eV"] / RY)))
    return out


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_ie_vs_N()
    rows = moseley_fits()
    fig_moseley(rows)
    fig_kinks()
    fig_neutral()
    fig_relativity(rows)
    fig_removal_map()
    print(f"{'N':>3} {'nl':>3} {'slope':>7} {'sigma':>7} {'p=N-1-s':>8} {'rms':>6}")
    for r in rows:
        if r["N"] in (1, 2, 3, 4, 5, 8, 10, 11, 12, 18, 19, 20, 21, 25, 26, 28, 29, 30, 31, 36, 37, 46,
                      47, 54, 55, 56, 60):
            print(f"{r['N']:3d} {r['n']}{L_LETTER[r['l']]} {r['slope']:7.4f} {r['sigma']:7.3f} "
                  f"{r['N'] - 1 - r['sigma']:8.3f} {r['rms']:6.3f}")
    print("\nneutral atoms: Z sym removed n*=sqrt(Ry/IE)  Zeff(true n)")
    for Z, s, rem, ns, ze in neutral_quantum_defects():
        if Z in (1, 2, 3, 4, 5, 10, 11, 18, 19, 36, 37, 54, 55, 86, 87):
            print(Z, s, rem, round(ns, 3), round(ze, 3))
