import sys, json, numpy as np, pandas as pd
sys.path.insert(0, sys.argv[1]); import msm_ident as M
d0 = json.load(open(sys.argv[1] + "/audit_f.json"))["fit_correct_likelihood"]
cg = [d0[k] for k in ["c_MD","g_MD","c_MW","g_MW","c_DR","g_DR"]]; b0 = d0["b_tm"]

def rates(age, afm, th):          # th: 12 c/g (MD,MW,DR,WR,RD,RW) + bD,bW,bR
    tm = M.softplus(age - afm); c = th; bD, bW, bR = th[12], th[13], th[14]
    return (np.exp(c[0]+c[1]*age+bD*tm), np.exp(c[2]+c[3]*age+bW*tm),
            np.exp(c[4]+c[5]*age+bR*tm), np.exp(c[6]+c[7]*age+bR*tm),
            np.exp(c[8]+c[9]*age+bD*tm), np.exp(c[10]+c[11]*age+bW*tm))
M._rates = rates
VAR = {  # name -> (b names, function b-vector -> (bD,bW,bR))
 "A shared":            (["b"],              lambda b: (b[0], b[0], b[0])),
 "B D,W,R separate":    (["bD","bW","bR"],   lambda b: (b[0], b[1], b[2])),
 "C D,R; W none":       (["bD","bR"],        lambda b: (b[0], 0.0,  b[1])),
 "D D=W, R":            (["bDiss","bR"],     lambda b: (b[0], b[0], b[1])),
}
def make_expand(f):
    def ex(tp):
        c = tp[:6]; bD, bW, bR = f(tp[6:])
        return np.array([c[0],c[1],c[2],c[3],c[4],c[5],c[4],c[5],c[0],c[1],c[2],c[3],bD,bW,bR])
    return ex
df = pd.read_csv(sys.argv[2]); df = df[(df.A=="M")&(df.sex=="f")]
ccs = sorted(df.cc.unique()); out = []
for cc in ccs + ["ALL"]:
    counts, tot = M.load_real(path=sys.argv[2], country=None if cc=="ALL" else cc)
    dw = [(a, b, sum(c.values())) for a, b, c in M.counts_to_design(counts)]
    N = sum(w for *_, w in dw)
    for vn, (bn, f) in VAR.items():
        M.expand = make_expand(f)
        bstart = {"A shared":[b0], "C D,R; W none":[b0,b0]}.get(vn, [b0]*len(bn))
        tp = np.array(cg + bstart)
        I = M.fisher(dw, tp, constrained=True)
        names = ["c_MD","g_MD","c_MW","g_MW","c_DR","g_DR"] + bn
        sp = M.spectrum(I, names)
        Ii = np.linalg.pinv(I); sd = np.sqrt(np.diag(Ii)); R = Ii/np.outer(sd, sd)
        row = dict(cc=cc, N=N, variant=vn, cond=sp["cond"], **{"se_"+k: sp["se"][k] for k in bn})
        for k in bn:   # correlation of each b with its own age slope
            g = {"b":"g_MD","bD":"g_MD","bW":"g_MW","bR":"g_DR","bDiss":"g_MD"}[k]
            row["cor_"+k+"_"+g] = R[names.index(k), names.index(g)]
        out.append(row)
    print(cc, round(N), flush=True)
pd.DataFrame(out).to_csv(sys.argv[3], index=False)
