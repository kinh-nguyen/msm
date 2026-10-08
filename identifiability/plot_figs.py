import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#9a9993"

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.edgecolor": MUTED, "axes.linewidth": .6,
    "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "grid.color": "#e6e5e0", "grid.linewidth": .6,
})

R = json.load(open("audit_f.json"))
fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.5))

# ---- (a) profile likelihood in the level of the dissolution intensities
d = np.array(R["profile_delta"])
c = np.array(R["profile_correct"]); c = c - c.max()
i = np.array(R["profile_implemented"]); i = i - i.max()
a0 = ax[0]
a0.axhline(0, color=MUTED, lw=.6, zorder=1)
a0.plot(d, c / 1e3, color=BLUE, lw=2, zorder=3)
a0.plot(d, i / 1e3, color=ORANGE, lw=2, zorder=3)
a0.set_ylim(-60, 8)
a0.axvline(0, color=MUTED, lw=.6, ls=(0, (3, 3)), zorder=1)
a0.plot([d[i.argmax()]], [0], "o", ms=5, color=ORANGE, mec="white", mew=1.2, zorder=4)
a0.plot([0], [0], "o", ms=5, color=BLUE, mec="white", mew=1.2, zorder=4)
a0.annotate("correct\nlikelihood", (-0.55, -32), color=BLUE, ha="center", fontsize=7.5)
a0.annotate("as implemented", (1.95, -18), color=ORANGE, ha="center", fontsize=7.5)
a0.annotate(r"$\hat\delta=%.2f$" % d[i.argmax()], (d[i.argmax()], 4),
            color=ORANGE, ha="center", fontsize=7.5)
a0.set_xlabel(r"$\delta$: common shift of $\log\alpha_D,\log\alpha_W$")
a0.set_ylabel(r"population log-lik. $-$ max ($\times 10^3$)")
a0.set_title("(a) level of the dissolution rates", loc="left", color=INK)
a0.grid(axis="y")

# ---- (b) sensitivity vs duration
t = np.array(R["duration_order"]["t"])
bp = R["duration_order"]["by_param"]
a1 = ax[1]
for nm, col, lab in [("c_MD", BLUE, r"$\alpha_D$  (slope %.2f)"),
                     ("c_DR", AQUA, r"$\rho_D$  (slope %.2f)"),
                     ("c_RD", VIOLET, r"$\alpha'_D$  (slope %.2f)")]:
    s = np.array(bp[nm]["sensitivity"])
    a1.loglog(t, s, color=col, lw=2, marker="o", ms=3.5, mec="white", mew=.8)
    a1.annotate(lab % bp[nm]["slope"], (t[-1] * 1.15, s[-1]), color=col,
                fontsize=7.5, va="center")
a1.set_xlim(.9, 130)
a1.set_xlabel("years since first marriage, $t$")
a1.set_ylabel(r"$\max_c\,|\partial \pi_c/\partial\theta|$")
a1.set_title("(b) order of information in $t$", loc="left", color=INK)
a1.grid(which="major", axis="y")

# ---- (c) standard errors vs augmentation
af = R["augmentation_fraction"]
fs = sorted(float(k) for k in af)
a2 = ax[2]
for nm, col, lab in [("c_WR", BLUE, r"$\rho_W$ (remarry after widowhood)"),
                     ("c_MW", ORANGE, r"$\alpha_W$ (widowhood, 1st union)"),
                     ("c_RW", VIOLET, r"$\alpha'_W$ (widowhood, later unions)")]:
    y = [af[("%g" % f) if ("%g" % f) in af else str(f)]["se"][nm] for f in fs]
    a2.plot(np.array(fs) * 100, y, color=col, lw=2, marker="o", ms=3.5,
            mec="white", mew=.8)
    a2.annotate(lab.split()[0], (fs[-1] * 100 + 2, y[-1]), color=col, fontsize=7.5,
                va="center")
a2.set_xlabel("% of sample with exact dissolution date")
a2.set_ylabel("asymptotic s.e. (log scale)")
a2.set_title("(c) value of the missing variable", loc="left", color=INK)
a2.set_xlim(-3, 118)
a2.grid(axis="y")

fig.tight_layout()
fig.savefig("fig_ident.pdf", bbox_inches="tight")
print("wrote fig_ident.pdf")
