# MultistageSurv — handover

Last updated 2026-10-07. Covers the sessions of 14–17 Aug 2026.
Project moved out of Dropbox; its home is now `~/may/papers/MultistageSurv` (git repo).

> Status note, 8 October 2026. This handover has been worked through and its decisions are recorded in `DECISIONS.md`, which supersedes it for the current state. The MM claim in §3 was confirmed from the code and data (all M→M episodes had fit = 3) and fixed. Review items 1 to 4, 7, 9, 12 (countries fitted separately), 13, 14 and 16 are addressed in the eight-state rewrite of `prep.Rmd`, `msm.cpp` and `run.r`; items 6, 8, 10 and 11 needed no change; v538 and v511a are deferred. File paths below predate the reorganisation: `paper.qmd` is now in `paper/`.

---

## 1. Where the paper stands

`paper.qmd` → `paper.tex` / `paper.pdf`. Background and Methods are written. Results
contain one figure (fitted vs empirical 1-year transition probabilities) and no numbers.
There is no Discussion, no Limitations, no Conclusion. The abstract is a methods note
rather than an abstract.

Outstanding editorial items, from the first review:

- Write Results: median ages, rates by country and sex, 95% UIs. Use the figures already
  built in `fig/` (`PVXVMXM_*`, `PR`, `PVV`, `PMJ_*`, `vvvxvm`, `xxxm`) — none is cited.
- Delete the stray R output in Results (`List of 1 $ PP: num [1:7,1:7,1:52]`); set that
  chunk to `include: false`. Reconcile the 7×7 array with the six declared states.
- Add Discussion, Limitations, Conclusion; keywords, author contributions, funding,
  ethics, data availability, conflicts. `\journal{medArxiv}` → medRxiv.
- Equation (1): the (6,6) entry is written `\lambda_{66}`; it is `q_{66}`.
- `\tau_{ij,sk} = \exp(\zeta_{ij,s,k} + \zeta_{ij,s,k} \cdot aai)` uses one symbol for
  both intercept and slope.
- The appendix figure documents an EDLL hazard the model no longer uses.
- "EPLD" → ELPD, and report the comparison result.
- **SHASH typo**: the text has γ = τ log(z + (1+z²)^−·⁵) − ν. It should be (1+z²)^+1/2,
  i.e. asinh(z). The code (`logSHASHz`) is correct; only the manuscript is wrong.
- `\bibliography{/Users/knguyen/zotero.bib}` is an absolute path; the `zotero_biblatex.bib`
  symlink will break for co-authors.
- Language pass: "disproportional distributed", "splitted", "formular", "righ-skewed",
  "paricipants", "diversed", "to to access", "The model's parameters was", "dissolute".

---

## 2. The identifiability question

Asked: can the dissolution and remarriage parameters be estimated from DHS at all, what
minimal data would be needed, and what simulations would demonstrate it.

Everything is in `identifiability/` — a 10-page technical note plus runnable code.

### Answer, in order

1. **The model is structurally identified.** Under a correctly specified likelihood all six
   marital intensities are determined by the observable cell probabilities. Proposition 2
   of the note gives a constructive proof by successive differentiation of the occupancy
   curve. The equality constraints (λ34=λ64, λ35=λ65, λ46=λ56) are *not* needed for
   identification — they are variance reduction, and they are testable in principle.

2. **The information is ordered in powers of marital duration.** Dissolution enters the
   observables at O(t), remarriage at O(t²), dissolution of later unions at O(t³),
   remarriage after a later dissolution at O(t⁴). Measured log–log slopes: 0.99, 1.91,
   2.93. At t = 1 the cells are 3.6×10⁴ times less sensitive to α′_D than to α_D.

3. **Practical identifiability.** Pooled over 37 countries (N = 579,616 ever-married women)
   the unconstrained 13-parameter model has Fisher condition number 6.1×10⁶. For the
   median country (Niger, N = 12,834) it is 1.1×10⁷ and the model is unusable: s.e. of
   log ρ_W is 2.30, a 95% interval spanning a factor of ~8,000. Forty-three women in the
   whole Niger sample occupy the widowed-and-remarried cell. **The constrained model is
   comfortably estimable country by country** (s.e. 0.13 / 0.33 / 0.20) — that is the
   positive finding.

4. **Observable cells**, women, all countries, weighted:
   M 424,786 (73.3%) · D₁ 35,591 (6.1%) · W₁ 19,139 (3.3%) · R 84,582 (14.6%) ·
   D₂ 11,416 (2.0%) · W₂ 4,101 (0.7%).

5. **Minimal additional data.** Per country, unconstrained model, s.e. of the log-intensity:

   | | α_W | ρ_W | α′_W |
   |---|---|---|---|
   | DHS as used | 1.41 | 2.30 | 1.25 |
   | + `v538` | 0.26 | 0.56 | 0.77 |
   | + exact dissolution date (all) | 0.18 | 0.39 | 0.81 |

   `v538` (how the previous union ended) is the best value by a wide margin and is already
   collected in roughly 24 DHS samples. It is the only variable that reveals *which*
   dissolution a remarriage followed. Exact dissolution dates for 1–2% of respondents buy
   most of the remaining precision. `v511a` (DHS-8, women, age at current union) alone is
   weak. Nothing DHS can offer identifies widowhood in later unions.

6. **Genuinely unidentified**, regardless of the likelihood: duration-since-*remarriage*
   effects (entry time to R unobserved); third and higher-order unions (`v503` is binary in
   every phase); marital-status-dependent mortality (no death state, so everything is
   conditional on survival); the τ(aai) shape term (within a single-survey country, age at
   interview and birth cohort are exactly collinear — 9 of 37 countries have one survey).

---

## 3. Open and disputed

**The `MM` / `fit = 3` claim.** The note's Lemma 2 says that `prep.Rmd` creates an `MM`
episode only when `ms == 'M'`, that those respondents get `fit = 3`, and that `makeQ`'s
`if (r < size && c < size)` guard then skips `add(2,3,·)` and `add(2,4,·)`, making M
absorbing, so their contribution is `log 1 = 0`. Kinh disputes this. It is not resolved.

Settle it with this test, which needs no agreement about the code:

```r
mm <- with(tdt, A == 2 & Z == 2 & fit == 3)      # the MM episodes
sum(mm); sum(tdt$n[mm])

o1 <- TMB::MakeADFun(modifyList(data, as.list(tdt)),        init, map = map, DLL = "msm")
o2 <- TMB::MakeADFun(modifyList(data, as.list(tdt[!mm, ])), init, map = map, DLL = "msm")

p <- o1$env$last.par.best
o1$fn(p) - o2$fn(p)          # 0  =>  those rows contribute nothing
max(abs(o1$gr(p) - o2$gr(p)))
```

If both are zero, the survival factor exp(−∫(λ_MD+λ_MW)) is not entering the likelihood.
If not, Lemma 2 is wrong and must be withdrawn.

**What depends on it.** Section 2 of the note, Proposition 1, and Figure 1(a) — including
the claim that the dissolution intensities are inflated by e^1.75 ≈ 5.8 and that the
objective is unbounded above. **Lemma 1 does not depend on it**: the ">1 union" D and W
cells are evaluated as P(in D) on a chain that cannot record whether R was visited, so the
"once" and ">once" cells do not partition. That one stands on its own.

---

## 4. Code review — likelihood path (`prep.Rmd` → `run.r` → `msm.cpp`)

### Changes the likelihood value

1. **Missing `v503` is coded as "more than one union".** `prep.Rmd:112` maps ever-married
   respondents with missing `n_union` to `"oneORmore"`; `fit = case_when(...)` has no branch
   for it, so it falls to `otherwise ~ 7`, and `prep.Rmd:141` then classifies them as
   **remarried** (`ms = "R"`). Missing union count inflates the cell that identifies the
   remarriage rate. Check: `cleaned.rds %>% count(is.na(n_union), marriage_age != 0)`.
   A censored respondent should contribute P(M) + P(R), not P(R).
2. **`init$itc = rep(0, 0, 0, 0, 0, 0)` (`run.r:49`) is `numeric(0)`** — `rep`'s second
   argument is `times`. The map expects length 6. Presumably `rep(0, 6)` was meant.
3. **">1 union" dissolution cells are not conditioned on remarriage** (Lemma 1 above).
   Fix by splitting D and W by union order in `makeQ`; costs no parameters.
4. **The same pre-marriage episode gets a different probability depending on `fit`**, because
   `fit` is a person-level attribute applied to all of that person's episodes. Truncate only
   the post-marriage segment.
5. The `MM` term — disputed, see §3.

### Probable bugs, cheap to check

6. `otherwise ~ 7` (`prep.Rmd:151`) is not dplyr syntax — needs `TRUE ~ 7` or `.default`.
7. `if_else(cond, "M", ms)` (l.140–146) on a factor; `dplyr::if_else` is type-strict.
8. `start == 0 → 1` (`run.r:14`) applied before `dur = end - start`: every `VV` episode
   loses its first year of exposure.
9. `fit = if_else(fit != 3, fit - 1, fit)` (`run.r:21`) — this off-by-one patch is a symptom
   that `prep.Rmd` still carries state counts from the old seven-state model in `model.R`.
10. `log(vp(Z) + eps)` with `eps = DBL_EPSILON` floors the penalty at ≈ −36 per unit weight,
    so the likelihood goes flat where a cell probability underflows.
11. `neversex = afs == 0` can be `NA` (if `afs` is NA and `marriage_age == 0`, l.109 never
    repairs it); `group_modify`'s `if (!.y$neversex)` then errors.

### Code does not match the manuscript

12. No hierarchical structure in `msm.cpp`. The manuscript specifies
    μ_{ij,s,k} = μ_{ij,s} + c_{ij,s,k}, c ~ N(0,κ); `run.r` fits one country at a time
    (`CC = 'ST'`) with no random effect and no κ.
13. `tau_aai`, `sigma_aai` and `b_tx` are mapped to `NA` (`run.r:63–68`) — described in the
    paper, not estimated.
14. Kish weights adjust for unequal weighting only; DHS is a two-stage cluster design, so
    the design effect from clustering is unaccounted for and intervals are too narrow.
    If `v001` is retained, inflate `n` by a design effect.
15. 225 respondents with `afs` but no `afm` are dropped (`prep.Rmd:116`). They are
    right-censored, not missing at random.
16. Since exp(c + g·a + b·t_m) = exp(c + (g+b)a − b·afm), the fitted `b_tm` is numerically
    minus the coefficient on age at first marriage. Add `afm` explicitly and report both,
    or the "duration effect" absorbs selection on age at marriage.

`identifiability/PATCH.md` has the concrete code changes for 1–4 and 16, plus the
covariate-centring point: `msm.cpp` puts N(0,1) on `itc` while `age` enters uncentred, so
the intercept is the log-intensity extrapolated to age 0, where that prior is extremely
informative.

---

## 5. Suggested order of work

1. Run the `MM` test in §3. It decides whether §2 of the note survives and whether the
   published dissolution rates need recomputing.
2. Fix review items 1, 2 and 3 — these change the numbers.
3. Refit and compare against the current fit, by country.
4. Decide what to report: occupancy probabilities and all-cause dissolution are estimable
   from DHS without argument; country-specific second-union dissolution and
   widowhood-specific remarriage are not, and should not appear as findings.
5. Add `v538` where available; `v511a` for DHS-8.
6. Then write Results and Discussion.

---

## 6. Files

```
identifiability/
  identifiability.pdf    the note, 10 pp. (compiles from the .tex; needs natbib,
                         cleveref, tcolorbox, mathptmx — not elsarticle)
  identifiability.tex
  fig_ident.pdf          Figure 1
  PATCH.md               concrete code changes
  msm_ident.py           reimplements the post-marriage block of msm.cpp exactly,
                         including the fit-driven submatrix truncation; four
                         observation schemes (dhs / v538 / v511a / exact)
  run_audit.py           produces every number in the note
  plot_figs.py           produces Figure 1
  audit_f.json           all results
```

Reproduce with `cd identifiability && python3 run_audit.py f` (numpy, scipy, pandas;
~4 min). It reads `../data/cc6.csv.bz2` — the path was the cloud session's staging
directory and has been changed to a repo-relative one (2026-10-07).

---

## 7. Caveats on what was done

- The reference parameter θ₀ used throughout the audit is the maximiser of the **corrected**
  likelihood on the real female cell counts, pooled over countries — not the published fit,
  which could not be reused because it comes from the disputed objective.
- The simulation studies S0–S7 in §6 of the note are **proposed designs, not run**. A partial
  run at per-country N was started and abandoned. S0 (likelihood verification) is the gate
  and should be run first.
- The DHS variable facts — no dissolution date in any phase; `v538` in roughly 24 samples;
  `v511a` new in DHS-8 and women-only — were verified against DHS-generated data dictionaries
  (World Bank microdata mirror, IPUMS-DHS) and the DHS-8 questionnaire text, **not** against
  a DHS-8 recode manual, which may not have been published. The ~24-sample list for `v538`
  needs a per-survey check before it goes in print.
- John & Nitsche (2022, *PDR* 48(4):1163–1201) estimate all-cause union dissolution from the
  same DHS data by indirect life-table methods, and state the identical data limitation.
  Their cohort-stationarity assumption is the analogue of the equality constraints here.
  The paper should engage with them.
