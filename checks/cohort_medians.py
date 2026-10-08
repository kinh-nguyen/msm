import sys, numpy as np, pandas as pd, pyreadr
d = pyreadr.read_r(sys.argv[1])[None]
d = d[(d.sex == "female") & (d.weights > 0)].copy()
d["cc"] = d.lab.str[:2]
d["w"] = d.weights / d.groupby("lab").weights.transform("mean")
d["agex"] = np.floor((d.doi - d.dob) / 12)                  # exact completed age
d["yob"] = 1900 + np.floor((d.dob - 1) / 12)
d["coh"] = (d.yob // 5) * 5
d["deb"] = d.afs.where(d.afs > 0)                            # NaN = not yet
d["mar"] = d.marriage_age.where(d.marriage_age > 0)
def km_median(t_ev, cens, w):
    ages = np.arange(5, 50); S = 1.0
    for a in ages:
        cz = np.isnan(t_ev)
        risk = w[t_ev >= a].sum() + w[cz & (cens > a)].sum() + 0.5 * w[cz & (cens == a)].sum()
        if risk < 30: return np.nan
        ev = w[t_ev == a].sum()
        if risk <= 0: return np.nan
        S2 = S * (1 - ev / risk)
        if S2 <= 0.5:   # linear interpolation within the year
            return a + (S - 0.5) / (S - S2)
        S = S2
    return np.nan
rows = []
for (cc, coh), g in d.groupby(["cc", "coh"]):
    if len(g) < 200: continue
    for ev in ["deb", "mar"]:
        m = km_median(g[ev].values, g.agex.values, g.w.values)
        if m >= np.quantile(g.agex.values, 0.5): m = np.nan   # median must lie below the cohort's median age
        rows.append(dict(cc=cc, coh=coh, event=ev, n=len(g), median=m))
r = pd.DataFrame(rows).dropna()
r.to_csv(sys.argv[2], index=False)
# linear vs quadratic in log median, per country and event
out = []
for (cc, ev), g in r.groupby(["cc", "event"]):
    if g.coh.nunique() < 4: continue
    x = (g.coh.values - 1975) / 10; y = np.log(g["median"].values); w = g.n.values
    b1 = np.polyfit(x, y, 1, w=np.sqrt(w)); b2 = np.polyfit(x, y, 2, w=np.sqrt(w))
    r1 = y - np.polyval(b1, x); r2 = y - np.polyval(b2, x)
    out.append(dict(cc=cc, event=ev, ncoh=len(g), coh0=g.coh.min(), coh1=g.coh.max(),
        med_first=g["median"].iloc[0], med_last=g["median"].iloc[-1],
        slope_pct_per_decade=100*b1[0], rmse_lin_pct=100*np.sqrt(np.average(r1**2, weights=w)),
        rmse_quad_pct=100*np.sqrt(np.average(r2**2, weights=w)), quad_pct=100*b2[0]))
o = pd.DataFrame(out); o.to_csv(sys.argv[3], index=False)
pd.set_option("display.width", 200)
print(o.round(2).to_string())
