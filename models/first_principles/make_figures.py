"""Figures for docs/first_principles.md (results/figures/fp_*.png)."""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zexp import ROOT  # noqa: E402
import relativity as rel  # noqa: E402
from common.atomdata import load_records, HARTREE_EV  # noqa: E402

FIG = os.path.join(ROOT, "results", "figures")
RES = os.path.join(ROOT, "results")
R = load_records()
REF = {(r["Z"], r["N"]): r["IE_eV"] for r in R}
C = ["#2a6fdb", "#e0672b", "#2e9e5b", "#9b4dca", "#c23b4b", "#7a7a7a", "#c9a227", "#1fa2a8"]


def load_pred(name, col="IE_pred_eV"):
    p = {}
    path = os.path.join(RES, f"fp_{name}_predictions.csv")
    if not os.path.exists(path):
        return p
    for r in csv.DictReader(open(path, encoding="utf-8")):
        try:
            v = float(r[col])
        except (ValueError, KeyError):
            continue
        if np.isfinite(v):
            p[(int(r["Z"]), int(r["N"]))] = v
    return p


def ape(p, k):
    return abs(p[k] / REF[k] - 1) * 100


def fig_hlike():
    Z = np.arange(1, 111)
    ref = np.array([REF[(z, 1)] for z in Z])
    curves = {
        "Bohr  Z²Ry": [z * z * 0.5 * HARTREE_EV for z in Z],
        "Dirac (point nucleus)": [rel.hydrogenic_binding(z, recoil=False, fns=False, qed=False) * HARTREE_EV for z in Z],
        "+ recoil + finite size": [rel.hydrogenic_binding(z, qed=False) * HARTREE_EV for z in Z],
        "+ QED (SE + VP)": [rel.hydrogenic_binding(z) * HARTREE_EV for z in Z],
    }
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for i, (k, v) in enumerate(curves.items()):
        ax.semilogy(Z, np.abs(np.array(v) / ref - 1) + 1e-12, "-", color=C[i], lw=1.8, label=k)
    ax.set_xlabel("Z (H-like ion, N = 1)")
    ax.set_ylabel("|relative error| vs NIST")
    ax.set_ylim(1e-11, 1)
    ax.axhline(1e-4, color="k", ls=":", lw=0.8)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title("Hydrogen-like ions: Bohr → Dirac → nuclear → QED")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fp_hlike_error.png"), dpi=140)
    plt.close(fig)


def fig_isoelectronic():
    preds = {"2-term exact (Z²ΔE0+ZΔE1)+rel": load_pred("zexp2_rel"),
             "completed square (Z−σ₁)² + rel": load_pred("zexp"),
             "+ dE2 per sequence (fitted)": load_pred("zexp3_seq")}
    Ns = [2, 3, 4, 10, 11, 18, 28, 36, 46, 54]
    fig, axes = plt.subplots(2, 5, figsize=(15, 6), sharey=True)
    for ax, N in zip(axes.flat, Ns):
        Zs = sorted(z for (z, n) in REF if n == N)
        for i, (lab, p) in enumerate(preds.items()):
            e = [ape(p, (z, N)) if (z, N) in p else np.nan for z in Zs]
            ax.semilogy(Zs, np.maximum(e, 1e-4), "-", color=C[i], lw=1.4, label=lab)
        ax.set_title(f"N = {N}", fontsize=10)
        ax.axvline(N, color="k", ls=":", lw=0.7)
        ax.set_ylim(1e-4, 1e4)
        ax.grid(alpha=0.25)
    for ax in axes[1]:
        ax.set_xlabel("Z")
    for ax in axes[:, 0]:
        ax.set_ylabel("APE %")
    axes[0, 0].legend(frameon=False, fontsize=7.5, loc="upper right")
    fig.suptitle("1/Z expansion along isoelectronic sequences (dotted line: neutral atom Z = N)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fp_zexp_isoelectronic.png"), dpi=130)
    plt.close(fig)


def fig_vs_ratio():
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    bins = np.linspace(0, 1, 21)
    for i, (lab, name) in enumerate((("Bohr (1 term)", "bohr"), ("2 terms + rel", "zexp2_rel"),
                                     ("completed square + rel", "zexp"), ("3 terms, dE2 per sequence (fit)", "zexp3_seq"))):
        p = load_pred(name)
        x = np.array([k[1] / k[0] for k in p])
        y = np.array([ape(p, k) for k in p])
        med = [np.median(y[(x > a) & (x <= b)]) if ((x > a) & (x <= b)).any() else np.nan
               for a, b in zip(bins[:-1], bins[1:])]
        ax.semilogy(0.5 * (bins[1:] + bins[:-1]), med, "o-", color=C[i], lw=1.6, ms=4, label=lab)
    ax.set_xlabel("N / Z  (0 = bare nucleus limit, 1 = neutral atom)")
    ax.set_ylabel("median APE % in bin")
    ax.axhline(1, color="k", ls=":", lw=0.8)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    ax.set_title("The 1/Z series is asymptotic in N/Z")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fp_zexp_error_vs_N_over_Z.png"), dpi=140)
    plt.close(fig)


def fig_sigma():
    rows = [r for r in csv.DictReader(open(os.path.join(RES, "fp_zexp_coefficients.csv"))) if r["which"] == "neutral"]
    N = np.array([int(r["N"]) for r in rows])
    sa = np.array([float(r["sigma_abinitio"]) for r in rows])
    ss = np.array([float(r["sigma_slater"]) for r in rows])
    se = []
    for r in rows:
        n = int(r["removed"][0])
        ie = REF.get((int(r["N"]), int(r["N"])))
        se.append(int(r["N"]) - n * np.sqrt(ie / 13.605693) if ie else np.nan)
    m = N <= 86
    fig, ax = plt.subplots(figsize=(7.5, 4.3))
    ax.plot(N[m], sa[m], "o-", color=C[0], ms=3, lw=1.2, label="ab initio σ₁ = −ΔE1/(2ΔE0)  (exact, first order)")
    ax.plot(N[m], ss[m], "s-", color=C[1], ms=3, lw=1.2, label="Slater (1930) σ")
    ax.plot(N[m], np.array(se)[m], "^-", color=C[2], ms=3, lw=1.2, label="'experimental' σ = Z − n√(IE/Ry)")
    ax.set_xlabel("Z = N (neutral atom, first IE)")
    ax.set_ylabel("screening constant of removed electron")
    ax.legend(frameon=False, fontsize=8.5)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fp_sigma_abinitio_vs_slater.png"), dpi=140)
    plt.close(fig)


def fig_dft():
    pth = os.path.join(RES, "fp_dft_predictions.csv")
    if not os.path.exists(pth):
        return
    rows = list(csv.DictReader(open(pth, encoding="utf-8")))
    neu = [r for r in rows if r["Z"] == r["N"]]
    Z = np.array([int(r["Z"]) for r in neu])
    ref = np.array([REF[(int(r["Z"]), int(r["N"]))] for r in neu])
    fig, ax = plt.subplots(figsize=(8, 4.3))
    ax.plot(Z, ref, "k-", lw=2, label="NIST")
    for i, (col, lab) in enumerate((("IE_pred_eV", "LSDA ΔSCF"), ("IE_TS_eV", "Slater transition state"),
                                    ("IE_HOMO_eV", "−ε_HOMO (LSDA)"))):
        ax.plot(Z, [float(r[col]) for r in neu], "o-", color=C[i], ms=3, lw=1.1, label=lab)
    ax.set_xlabel("Z")
    ax.set_ylabel("first ionization energy (eV)")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    ax.set_title("Kohn-Sham LSDA (own solver): neutral-atom first IEs")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fp_dft_neutral_first_IE.png"), dpi=140)
    plt.close(fig)
    # all ions Z<=18: APE vs N/Z
    fig, ax = plt.subplots(figsize=(7, 4))
    for i, (col, lab) in enumerate((("IE_pred_eV", "ΔSCF"), ("IE_TS_eV", "transition state"))):
        x = [int(r["N"]) / int(r["Z"]) for r in rows if int(r["Z"]) <= 18]
        y = [abs(float(r[col]) / REF[(int(r["Z"]), int(r["N"]))] - 1) * 100 for r in rows if int(r["Z"]) <= 18]
        ax.semilogy(x, y, "o", color=C[i], ms=3, alpha=0.7, label=lab)
    ax.set_xlabel("N / Z")
    ax.set_ylabel("APE %")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    ax.set_title("LSDA ΔSCF, all ions with Z ≤ 18")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fp_dft_ions_Zle18.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    for f in (fig_hlike, fig_isoelectronic, fig_vs_ratio, fig_sigma, fig_dft):
        try:
            f()
            print("ok", f.__name__)
        except Exception as e:  # noqa: BLE001
            print("FAIL", f.__name__, e)
