"""Render the 33-parameter Screened Rydberg formula with its fitted values as an image.

    py -3.11 tools/formula_card.py   -> results/figures/formula_33param.png
Values are read from results/pa_params.json (candidate pa_hier_rel); nothing is typed by hand.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
c = json.load(open(os.path.join(ROOT, "results", "pa_params.json"), encoding="utf-8"))["candidates"]["pa_hier_rel"]
P = dict(zip(c["names"], c["values"]))
f4 = lambda v: f"{v:+.4f}".replace("+", " ") if v >= 0 else f"{v:.4f}"

fig = plt.figure(figsize=(12, 13.2), dpi=160)
fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
y = 0.965
def t(s, size=13, dy=0.04, x=0.04, **kw):
    global y
    ax.text(x, y, s, fontsize=size, va="top", ha="left", **kw); y -= dy

t("Screened Rydberg formula — final model (pa_hier_rel), 33 fitted parameters", 18, 0.05, weight="bold")
t("Ionization energy of the N-electron ion of element Z (eV); removed electron in subshell (n, l), occupancy k",
  11, 0.045, color="#444")
t(r"$\mathrm{IE}=\mu(Z)\left[\mathrm{Ry}\,\frac{Z_{\rm eff}^2}{n^2}\,F_{n,j}(Z_{\rm eff})\,"
  r"\left(1+r_c\,\frac{(Z\alpha)^2}{n}\left(\frac{Z_{\rm eff}}{Z_a}-1\right)\right)"
  r"+\mathrm{Ry}\,\frac{x_l\,K_l(k)}{n^2}\right]-\Delta E_{\rm QED+FNS}\left(\frac{Z_{\rm eff}}{Z}\right)^2$", 17, 0.085)
t(r"$Z_{\rm eff}=Z-\sigma_1-D,\qquad D=\frac{T}{Z_a+\kappa+|T|/h},\qquad Z_a=Z-N+1$", 17, 0.065)
t(r"$T=\sum_{g}\tau_g\,\nu_g+\sum_{c}\delta\tau_c\,\nu_c,\qquad h=(N-1)-\sigma_1\ \ (T\geq0),\quad h=\sigma_1\ \ (T<0)$",
  17, 0.06)
t(r"$\sigma_1$ = exact first-order (1/Z) screening constant (no parameters);  $F_{n,j}$ = exact Dirac/Schrödinger "
  r"ratio;  $\mu$ = reduced mass;", 11.5, 0.03, color="#333")
t(r"$K_l(k)$ = Hund kink ($K_p$ = 0, 0.6, 1.2, −1.2, −0.6, 0);  QED+FNS only for 1s/2s removal.  "
  r"Bound: $Z_a \leq Z_{\rm eff} \leq Z$ for any parameter values.", 11.5, 0.05, color="#333")

def block(title, items, ncol):
    global y
    t(title, 13.5, 0.035, weight="bold", color="#1f4e79")
    rows = [items[i:i + ncol] for i in range(0, len(items), ncol)]
    for r in rows:
        for j, (k, v) in enumerate(r):
            ax.text(0.06 + j * (0.92 / ncol), y, k, fontsize=12, va="top", ha="left")
            ax.text(0.06 + j * (0.92 / ncol) + 0.92 / ncol * 0.62, y, f4(v), fontsize=12, va="top", ha="left",
                    family="monospace")
        y -= 0.03
    y -= 0.012

G = ["same", "in", "core", "df", "out"]
block("Group screening  τ_g  (5)", [(fr"$\tau_{{\rm {g}}}$", P["tau_" + g]) for g in G], 5)
block("Saturation  κ  (1)   — used as |κ| + 0.05", [(r"$\kappa$", P["kappa"])], 5)
dev = [n for n in c["names"] if n.startswith("dtau_")]
block("Class deviations  δτ_c  (19, ridge-shrunk; classes defined in paper Table A2)",
      [(n[5:], P[n]) for n in dev], 4)
block("Relativistic penetration  r_c  (5)",
      [(fr"$r_{{{k}}}$", P["rel_" + k]) for k in ["s", "p1", "p3", "d", "f"]], 5)
block("Hund amplitudes  x_l  (3)", [(fr"$x_{{{k}}}$", P["x_" + k]) for k in ["p", "d", "f"]], 5)

t(f"Total fitted parameters: {len(c['values'])}.   Accuracy on 5847 NIST IEs: 1.87 % mean (0.94 % median); "
  "neutral atoms 7.5 %; blind Z≥55 test 6.60 %.", 11.5, 0.03, color="#444")
t("Source: results/pa_params.json; paper Table A1.  M. Awais, Screened Rydberg formula (2026).", 10, 0.0, color="#777")
out = os.path.join(ROOT, "results", "figures", "formula_33param.png")
fig.savefig(out, dpi=160, facecolor="white")
print("wrote", out, len(c["values"]), "params")
