"""
Identifiability audit for the marital-transition block of the DHS multistate
model in MultistageSurv/code/msm.cpp.

Everything after first marriage is reproduced exactly as implemented, then
analysed for structural and practical identifiability.
"""
from __future__ import annotations

import json
import os
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

AMAX = 60  # oldest age modelled

# ---------------------------------------------------------------------------
# States of the expanded ("observation-complete") chain
#   0 M   first union intact
#   1 D1  first union ended by divorce/separation, not remarried
#   2 W1  first union ended by widowhood, not remarried
#   3 R   currently in a union of order >= 2
#   4 D2  a union of order >= 2 ended by divorce/separation
#   5 W2  a union of order >= 2 ended by widowhood
# DHS (v501, v503) observes exactly this partition.
CELLS = ["M", "D1", "W1", "R", "D2", "W2"]

# Unconstrained parameter vector (13)
PAR_FULL = ["c_MD", "g_MD", "c_MW", "g_MW", "c_DR", "g_DR",
            "c_WR", "g_WR", "c_RD", "g_RD", "c_RW", "g_RW", "b_tm"]
# Constrained model of the paper (7): c_RD=c_MD, g_RD=g_MD, c_RW=c_MW,
# g_RW=g_MW, c_WR=c_DR, g_WR=g_DR
PAR_PAPER = ["c_MD", "g_MD", "c_MW", "g_MW", "c_DR", "g_DR", "b_tm"]


def expand(tp):
    c_MD, g_MD, c_MW, g_MW, c_DR, g_DR, b = tp
    return np.array([c_MD, g_MD, c_MW, g_MW, c_DR, g_DR,
                     c_DR, g_DR, c_MD, g_MD, c_MW, g_MW, b])


def softplus(x):
    return np.log1p(np.exp(-abs(x))) + max(x, 0.0)


def _rates(age, afm, th):
    tm = softplus(age - afm)
    c = th
    return (np.exp(c[0] + c[1] * age + c[12] * tm),   # MD
            np.exp(c[2] + c[3] * age + c[12] * tm),   # MW
            np.exp(c[4] + c[5] * age + c[12] * tm),   # DR
            np.exp(c[6] + c[7] * age + c[12] * tm),   # WR
            np.exp(c[8] + c[9] * age + c[12] * tm),   # RD
            np.exp(c[10] + c[11] * age + c[12] * tm))  # RW


def _P_year(age, afm, th):
    MD, MW, DR, WR, RD, RW = _rates(age, afm, th)
    Q = np.zeros((6, 6))
    Q[0, 1] = MD; Q[0, 2] = MW
    Q[1, 3] = DR
    Q[2, 3] = WR
    Q[3, 4] = RD; Q[3, 5] = RW
    Q[4, 3] = DR; Q[5, 3] = WR       # 3rd+ unions merged with the 2nd
    np.fill_diagonal(Q, -Q.sum(1))
    return expm(Q)


def sweep(afm, th, amax=AMAX):
    """Forward sweep: returns dict aai -> length-6 cell probability vector.
    One sweep serves every age at interview, as in msm.cpp (annual expm)."""
    out = {}
    v = np.zeros(6); v[0] = 1.0
    a = int(afm) + 1
    out[a] = v.copy()
    while a < amax:
        v = v @ _P_year(float(a), float(afm), th)
        a += 1
        out[a] = v.copy()
    return out


# ---------------------------------------------------------------------------
# The likelihood as ACTUALLY implemented (submatrix truncation by `fit`)
# ---------------------------------------------------------------------------
def _P_year_paper(age, afm, th, size):
    """makeQ() of msm.cpp restricted to post-marriage states {M,D,W,R}.
    size=2 reproduces fit=5 (R column absent -> D,W absorbing);
    size=3 reproduces fit=6 (full chain)."""
    MD, MW, DR, WR, RD, RW = _rates(age, afm, th)
    Q = np.zeros((4, 4))
    Q[0, 1] = MD; Q[0, 2] = MW
    if size == 3:
        Q[1, 3] = DR; Q[2, 3] = WR
        Q[3, 1] = RD; Q[3, 2] = RW
    np.fill_diagonal(Q, -Q.sum(1))
    return expm(Q)


def sweep_implemented(afm, th, amax=AMAX):
    """dict aai -> the six quantities the current code evaluates.
    Entry 0 is identically 1: with fit=3 the state M is absorbing, so a
    currently-married-once respondent contributes log(1)=0."""
    out = {}
    v5 = np.zeros(4); v5[0] = 1.0
    v6 = np.zeros(4); v6[0] = 1.0
    a = int(afm) + 1
    out[a] = np.array([1.0, v5[1], v5[2], v6[3], v6[1], v6[2]])
    while a < amax:
        v5 = v5 @ _P_year_paper(float(a), float(afm), th, 2)
        v6 = v6 @ _P_year_paper(float(a), float(afm), th, 3)
        a += 1
        out[a] = np.array([1.0, v5[1], v5[2], v6[3], v6[1], v6[2]])
    return out


# ---------------------------------------------------------------------------
# Augmented observation schemes
# ---------------------------------------------------------------------------
def sweep_v538(afm, th, amax=AMAX):
    """DHS v538 ('how did your previous union end') observed for respondents in
    a union of order >= 2: splits cell R into R|divorced and R|widowed.
    Seven cells: M, D1, W1, R_D, R_W, D2, W2."""
    out = {}
    v = np.zeros(7); v[0] = 1.0
    a = int(afm) + 1
    out[a] = v.copy()
    while a < amax:
        MD, MW, DR, WR, RD, RW = _rates(float(a), float(afm), th)
        Q = np.zeros((7, 7))
        # 0 M, 1 D1, 2 W1, 3 R_D, 4 R_W, 5 D2, 6 W2
        Q[0, 1] = MD; Q[0, 2] = MW
        Q[1, 3] = DR
        Q[2, 4] = WR
        Q[3, 5] = RD; Q[3, 6] = RW
        Q[4, 5] = RD; Q[4, 6] = RW
        Q[5, 3] = DR
        Q[6, 4] = WR
        np.fill_diagonal(Q, -Q.sum(1))
        v = v @ expm(Q)
        a += 1
        out[a] = v.copy()
    return out


def sweep_remarriage_age(afm, th, amax=AMAX):
    """DHS-8 v511a (age at start of CURRENT union) observed for respondents in
    a union of order >= 2.  Cells: M, D1, W1, (R, r) for each remarriage year
    r, D2, W2."""
    out = {}
    pre = np.zeros(3); pre[0] = 1.0          # M, D1, W1
    post = {}                                 # r -> length-3 vector (R, D2, W2)
    a = int(afm) + 1
    def snapshot():
        d = {"M": pre[0], "D1": pre[1], "W1": pre[2],
             "D2": sum(v[1] for v in post.values()),
             "W2": sum(v[2] for v in post.values())}
        for r, v in post.items():
            d[("R", r)] = v[0]
        return d
    out[a] = snapshot()
    while a < amax:
        MD, MW, DR, WR, RD, RW = _rates(float(a), float(afm), th)
        # pre-remarriage block with exit to R
        Qp = np.zeros((4, 4))                # M, D1, W1, ->R
        Qp[0, 1] = MD; Qp[0, 2] = MW
        Qp[1, 3] = DR; Qp[2, 3] = WR
        np.fill_diagonal(Qp, -Qp.sum(1))
        Pp = expm(Qp)
        newpre = np.array([pre @ Pp[:3, 0], pre @ Pp[:3, 1], pre @ Pp[:3, 2]])
        flux = float(pre @ Pp[:3, 3])
        # post-remarriage block, cycles among R, D2, W2
        Qq = np.zeros((3, 3))
        Qq[0, 1] = RD; Qq[0, 2] = RW
        Qq[1, 0] = DR; Qq[2, 0] = WR
        np.fill_diagonal(Qq, -Qq.sum(1))
        Pq = expm(Qq)
        post = {r: v @ Pq for r, v in post.items()}
        pre = newpre
        a += 1
        if flux > 0:
            post[a] = post.get(a, np.zeros(3)) + np.array([flux, 0.0, 0.0])
        out[a] = snapshot()
    return out


def sweep_dissolution_age(afm, th, amax=AMAX):
    """Age (and type) at the end of the FIRST union observed exactly -- the
    single variable DHS does not collect.  Cells: M, (j, u, current) with
    j in {D,W}, u the exit year."""
    out = {}
    surv = 1.0                                # still in M
    cohorts = {}                              # (j, u) -> length-5 vector over D1,W1,R,D2,W2
    a = int(afm) + 1
    def snapshot():
        d = {"M": surv}
        for (j, u), v in cohorts.items():
            d[(j, u, "notR")] = v[0] + v[1]
            d[(j, u, "R")] = v[2]
            d[(j, u, "D2")] = v[3]
            d[(j, u, "W2")] = v[4]
        return d
    out[a] = snapshot()
    while a < amax:
        MD, MW, DR, WR, RD, RW = _rates(float(a), float(afm), th)
        Qm = np.array([[-(MD + MW), MD, MW], [0, 0, 0], [0, 0, 0]])
        Pm = expm(Qm)
        newsurv = surv * Pm[0, 0]
        fD, fW = surv * Pm[0, 1], surv * Pm[0, 2]
        Qc = np.zeros((5, 5))                 # D1, W1, R, D2, W2
        Qc[0, 2] = DR; Qc[1, 2] = WR
        Qc[2, 3] = RD; Qc[2, 4] = RW
        Qc[3, 2] = DR; Qc[4, 2] = WR
        np.fill_diagonal(Qc, -Qc.sum(1))
        Pc = expm(Qc)
        cohorts = {k: v @ Pc for k, v in cohorts.items()}
        surv = newsurv
        a += 1
        if fD > 0:
            cohorts[("D", a)] = cohorts.get(("D", a), np.zeros(5)) + np.array([fD, 0, 0, 0, 0])
        if fW > 0:
            cohorts[("W", a)] = cohorts.get(("W", a), np.zeros(5)) + np.array([0, fW, 0, 0, 0])
        out[a] = snapshot()
    return out


SWEEPS = {
    "dhs": sweep,                       # current status only (what DHS gives)
    "v538": sweep_v538,                 # + manner of previous dissolution
    "v511a": sweep_remarriage_age,      # + age at start of current union
    "exact": sweep_dissolution_age,     # + exact age/type of first dissolution
}


# ---------------------------------------------------------------------------
# Design and real cell counts
# ---------------------------------------------------------------------------
def load_real(path=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "cc6.csv.bz2"),
              sex="f", country=None):
    import pandas as pd
    df = pd.read_csv(path)
    df = df[df.A == "M"].copy()
    if sex is not None:
        df = df[df.sex == sex]
    if country is not None:
        df = df[df.cc == country]

    def cell(r):
        if r.Z == "M":
            return "M"
        if r.Z == "R":
            return "R"
        return ({"D": "D1", "W": "W1"} if r.fit == 6 else {"D": "D2", "W": "W2"})[r.Z]

    df["cell"] = df.apply(cell, axis=1)
    df = df[(df.afm >= 10) & (df.afm <= 45) & (df.aai > df.afm + 1) & (df.aai <= AMAX)]
    counts = (df.groupby(["afm", "aai", "cell"])["n"].sum().reset_index())
    return counts, df.groupby("cell")["n"].sum().to_dict()


def counts_to_design(counts):
    """[(afm, aai, {cell: n})] sorted, plus total N."""
    d = {}
    for _, r in counts.iterrows():
        d.setdefault((int(r.afm), int(r.aai)), {})[r.cell] = float(r.n)
    return [(k[0], k[1], v) for k, v in sorted(d.items())]


def synthetic_design(n=200000.0):
    out = {}
    for aai in range(17, 50):
        for afm in range(12, min(aai - 1, 46)):
            w = np.exp(-0.5 * ((np.log(afm) - np.log(18.5)) / 0.22) ** 2)
            out[(afm, aai)] = w
    s = sum(out.values())
    return [(a, b, {"_N": n * w / s}) for (a, b), w in out.items()]


# ---------------------------------------------------------------------------
# Likelihoods
# ---------------------------------------------------------------------------
def negll_real(tp, design):
    """Correct multinomial log-likelihood on the real cell counts."""
    th = expand(tp)
    tot = 0.0
    by_afm = {}
    for afm, aai, cells in design:
        by_afm.setdefault(afm, []).append((aai, cells))
    for afm, items in by_afm.items():
        sw = sweep(afm, th)
        for aai, cells in items:
            p = sw.get(aai)
            if p is None:
                continue
            for c, n in cells.items():
                tot += n * np.log(max(p[CELLS.index(c)], 1e-300))
    return -tot


def expected_ll(tp, design_w, th0, implemented=False):
    """Population log-likelihood E_{th0} log L(th): its maximiser is the
    probability limit of the estimator; a flat profile is exactly structural
    non-identifiability."""
    th = expand(tp) if len(tp) == 7 else np.asarray(tp)
    tot = 0.0
    by_afm = {}
    for afm, aai, N in design_w:
        by_afm.setdefault(afm, []).append((aai, N))
    for afm, items in by_afm.items():
        sw0 = sweep(afm, th0)
        sw = sweep_implemented(afm, th) if implemented else sweep(afm, th)
        for aai, N in items:
            p0, p = sw0.get(aai), sw.get(aai)
            if p0 is None or p is None:
                continue
            tot += N * float(np.sum(p0 * np.log(np.maximum(p, 1e-300))))
    return tot


# ---------------------------------------------------------------------------
# Fisher information
# ---------------------------------------------------------------------------
def _probvecs(design_w, th, scheme="dhs"):
    """[(N, dict cell->prob)] for one parameter value."""
    fn = SWEEPS[scheme]
    by_afm = {}
    for afm, aai, N in design_w:
        by_afm.setdefault(afm, []).append((aai, N))
    out = []
    for afm, items in by_afm.items():
        sw = fn(afm, th)
        for aai, N in items:
            p = sw.get(aai)
            if p is None:
                continue
            if scheme == "dhs":
                p = {c: p[i] for i, c in enumerate(CELLS)}
            elif scheme == "v538":
                p = {c: p[i] for i, c in enumerate(
                    ["M", "D1", "W1", "R_D", "R_W", "D2", "W2"])}
            out.append((N, p))
    return out


def fisher(design_w, th, free=None, scheme="dhs", h=1e-4, constrained=False):
    th = np.asarray(th, float)
    k = len(th)
    free = list(range(k)) if free is None else free
    base = _probvecs(design_w, expand(th) if constrained else th, scheme)
    derivs = []
    for j in free:
        tp = th.copy(); tp[j] += h
        tm = th.copy(); tm[j] -= h
        pp = _probvecs(design_w, expand(tp) if constrained else tp, scheme)
        pm = _probvecs(design_w, expand(tm) if constrained else tm, scheme)
        derivs.append([{c: (a[1][c] - b[1][c]) / (2 * h) for c in a[1]}
                       for a, b in zip(pp, pm)])
    m = len(free)
    I = np.zeros((m, m))
    for i in range(m):
        for j in range(i, m):
            s = 0.0
            for d, (N, p) in enumerate(base):
                for c in p:
                    if p[c] > 1e-12:
                        s += N * derivs[i][d][c] * derivs[j][d][c] / p[c]
            I[i, j] = I[j, i] = s
    return I


def spectrum(I, names):
    w, V = np.linalg.eigh(I)
    o = np.argsort(w)[::-1]
    w, V = w[o], V[:, o]
    Iinv = np.linalg.pinv(I, rcond=1e-14)
    se = np.sqrt(np.clip(np.diag(Iinv), 0, None))
    return dict(eigenvalues=[float(x) for x in w],
                cond=float(w[0] / max(w[-1], 1e-300)),
                se={n: float(s) for n, s in zip(names, se)},
                weakest={n: float(v) for n, v in zip(names, V[:, -1])})
