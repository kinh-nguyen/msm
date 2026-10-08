import sys, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
r = pd.read_csv(sys.argv[1]); ccs = sorted(r.cc.unique())
fig, axs = plt.subplots(7, 6, figsize=(13, 14), sharex=True)
col = {"deb": "#2a6fbb", "mar": "#c4572a"}
for ax, cc in zip(axs.flat, ccs):
    for ev, g in r[r.cc == cc].groupby("event"):
        x = (g.coh.values - 1975) / 10; y = g["median"].values
        ax.plot(g.coh, y, "o", ms=3, color=col[ev])
        if len(g) >= 3:
            b = np.polyfit(x, np.log(y), 1, w=np.sqrt(g.n.values))
            xx = np.linspace(x.min(), x.max(), 20)
            ax.plot(1975 + 10 * xx, np.exp(np.polyval(b, xx)), "-", lw=1, color=col[ev])
    ax.set_title(cc, fontsize=9); ax.tick_params(labelsize=7)
for ax in axs.flat[len(ccs):]: ax.axis("off")
fig.text(0.5, 0.005, "Birth cohort (5-year)", ha="center")
fig.text(0.005, 0.5, "Median age (Kaplan-Meier)", va="center", rotation=90)
fig.legend(handles=[plt.Line2D([], [], color=col["deb"], marker="o", label="Sexual debut"),
                    plt.Line2D([], [], color=col["mar"], marker="o", label="First marriage")],
           loc="lower right", bbox_to_anchor=(0.98, 0.04), fontsize=10)
fig.suptitle("Empirical median age at debut and first marriage by birth cohort, women, with linear fit on log scale", fontsize=11)
fig.tight_layout(rect=(0.01, 0.01, 1, 0.98)); fig.savefig(sys.argv[2], dpi=110)
