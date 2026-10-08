"""Driver: runs the full identifiability audit and dumps results to JSON."""
import json, sys, time
import numpy as np
from scipy.optimize import minimize
import msm_ident as M

t0 = time.time()
RES = {}
SEX = sys.argv[1] if len(sys.argv) > 1 else "f"

# ---------------------------------------------------------------- 1. real data
counts, totals = M.load_real(sex=SEX)
design_counts = M.counts_to_design(counts)
design_w = [(a, b, sum(c.values())) for a, b, c in design_counts]
N_total = sum(w for _, _, w in design_w)
RES["cell_totals"] = {k: float(v) for k, v in totals.items()}
RES["cell_share"] = {k: float(v) / sum(totals.values()) for k, v in totals.items()}
RES["n_ever_married"] = float(N_total)
RES["n_design_points"] = len(design_w)
print("cell totals", RES["cell_totals"], flush=True)

# --------------------------------------- 2. fit the constrained model properly
def nll(tp):
    if not np.all(np.isfinite(tp)):
        return 1e12
    try:
        return M.negll_real(tp, design_counts)
    except Exception:
        return 1e12

start = np.array([-6.0, 0.02, -7.5, 0.05, -3.0, -0.02, 0.0])
fit = minimize(nll, start, method="Nelder-Mead",
               options=dict(maxiter=4000, maxfev=8000, xatol=1e-6, fatol=1e-4))
fit = minimize(nll, fit.x, method="Nelder-Mead",
               options=dict(maxiter=4000, maxfev=8000, xatol=1e-8, fatol=1e-6))
tp_hat = fit.x
RES["fit_correct_likelihood"] = dict(zip(M.PAR_PAPER, tp_hat.tolist()))
RES["fit_nll"] = float(fit.fun)
print("theta_hat", tp_hat, "nll", fit.fun, f"[{time.time()-t0:.0f}s]", flush=True)

th0 = M.expand(tp_hat)

# implied annual rates at reference ages, for the write-up
ref = {}
for afm, age in [(18, 20), (18, 25), (18, 30), (18, 40)]:
    r = M._rates(float(age), float(afm), th0)
    ref[f"afm{afm}_age{age}"] = dict(zip(["MD", "MW", "DR", "WR", "RD", "RW"],
                                        [float(x) for x in r]))
RES["implied_rates"] = ref

# ------------------------------------------- 3. Fisher spectra, three models
print("Fisher: constrained ...", flush=True)
I7 = M.fisher(design_w, tp_hat, scheme="dhs", constrained=True)
RES["fisher_constrained"] = M.spectrum(I7, M.PAR_PAPER)

print("Fisher: unconstrained ...", flush=True)
I13 = M.fisher(design_w, th0, scheme="dhs", constrained=False)
RES["fisher_unconstrained"] = M.spectrum(I13, M.PAR_FULL)

for scheme, label in [("v538", "unconstrained_plus_v538"),
                      ("v511a", "unconstrained_plus_v511a"),
                      ("exact", "unconstrained_plus_exact_dissolution")]:
    print("Fisher:", scheme, f"[{time.time()-t0:.0f}s]", flush=True)
    I = M.fisher(design_w, th0, scheme=scheme, constrained=False)
    RES["fisher_" + label] = M.spectrum(I, M.PAR_FULL)

# ------------------------- 4. profile of the population log-likelihood in the
#                             overall level of the dissolution rates
print("profiles ...", flush=True)
deltas = np.round(np.arange(-1.5, 3.01, 0.25), 3)
prof_correct, prof_impl = [], []
for d in deltas:
    tp = tp_hat.copy()
    tp[0] += d          # c_MD
    tp[2] += d          # c_MW  (scale both dissolution intensities by e^delta)
    prof_correct.append(M.expected_ll(tp, design_w, th0, implemented=False))
    prof_impl.append(M.expected_ll(tp, design_w, th0, implemented=True))
RES["profile_delta"] = deltas.tolist()
RES["profile_correct"] = [float(x) for x in prof_correct]
RES["profile_implemented"] = [float(x) for x in prof_impl]

# ------------------- 5. profile along the two constraints the paper imposes
#  (a) second-union dissolution relative to first: c_RD, c_RW shifted by kappa
#  (b) remarriage after widowhood vs divorce: c_WR shifted by psi
prof2 = {}
grid = np.round(np.arange(-2.0, 2.01, 0.2), 3)
for name, idx in [("kappa_second_union_dissolution", (8, 10)),
                  ("psi_remarriage_after_widowhood", (6,))]:
    vals = []
    for g in grid:
        th = th0.copy()
        for i in idx:
            th[i] += g
        vals.append(M.expected_ll(th, design_w, th0, implemented=False))
    prof2[name] = dict(grid=grid.tolist(), ll=[float(v) for v in vals])
RES["profiles_constraints"] = prof2

# ---------------------------- 6. how the information scales with duration:
#    leading power of t at which each parameter first perturbs the cells
hom = np.array([-6.0, 0.0, -7.0, 0.0, -3.0, 0.0, -3.0, 0.0,
                -6.0, 0.0, -7.0, 0.0, 0.0])
ts = np.array([1, 2, 3, 4, 6, 8, 12, 16, 24, 32])
orders = {}
for j, nm in enumerate(M.PAR_FULL):
    if nm.startswith("g") or nm == "b_tm":
        continue
    h = 1e-4
    tp, tm = hom.copy(), hom.copy()
    tp[j] += h; tm[j] -= h
    sp, sm = M.sweep(0, tp, amax=40), M.sweep(0, tm, amax=40)
    d = [float(np.abs(sp[1 + t] - sm[1 + t]).max() / (2 * h)) for t in ts]
    d = np.array(d)
    ok = d > 1e-14
    slope = np.polyfit(np.log(ts[ok][:5]), np.log(d[ok][:5]), 1)[0] if ok.sum() >= 3 else np.nan
    orders[nm] = dict(sensitivity=[float(x) for x in d], slope=float(slope))
RES["duration_order"] = dict(t=ts.tolist(), by_param=orders)

# ---------------------------- 7. observability (Krylov) rank of the chain
Q = np.zeros((6, 6))
aD, aW, rD, rW, bD, bW = .02, .006, .05, .04, .02, .006
Q[0, 1] = aD; Q[0, 2] = aW; Q[1, 3] = rD; Q[2, 3] = rW
Q[3, 4] = bD; Q[3, 5] = bW; Q[4, 3] = rD; Q[5, 3] = rW
np.fill_diagonal(Q, -Q.sum(1))
K = np.vstack([np.eye(6)[0] @ np.linalg.matrix_power(Q, k) for k in range(6)])
s = np.linalg.svd(K, compute_uv=False)
RES["krylov"] = dict(rank=int(np.linalg.matrix_rank(K, tol=1e-12)),
                     singular_values=[float(x) for x in s],
                     cond=float(s[0] / s[-1]))

# ---------------------------- 8. minimal augmentation: what fraction of the
#    sample needs the dissolution date before the unconstrained model is usable
frac_res = {}
I_dhs = M.fisher(design_w, th0, scheme="dhs", constrained=False)
I_ex = M.fisher(design_w, th0, scheme="exact", constrained=False)
for f in [0.0, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 1.0]:
    I = (1 - f) * I_dhs + f * I_ex
    sp = M.spectrum(I, M.PAR_FULL)
    frac_res[f] = dict(cond=sp["cond"], se=sp["se"])
RES["augmentation_fraction"] = {str(k): v for k, v in frac_res.items()}

json.dump(RES, open(f"audit_{SEX}.json", "w"), indent=1)
print("done", time.time() - t0, flush=True)
