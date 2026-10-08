# Minimal patch: make the objective the likelihood of the observed data

Two changes. Neither adds a parameter; both are required before any estimate of
a dissolution or remarriage rate means anything.

## 1. `code/prep.Rmd` — give every ever-married respondent one post-marriage episode

The current code creates the `MM` episode only for `ms == 'M'`, and those
respondents get `fit = 3`, which deletes the dissolution intensities from the
generator. Replace the `MM/MD/MW/MR` block with a single episode whose end state
is the *observed cell*:

```r
# observed post-marriage cell, from v501 (marital_status) and v503 (n_union)
.x$cell <- dplyr::case_when(
  .x$ms == "M" ~ "M",    # in union,  once          -> state M
  .x$ms == "D" & .x$u == "once" ~ "D1",
  .x$ms == "W" & .x$u == "once" ~ "W1",
  .x$ms == "R" ~ "R",    # in union,  > once        -> state R
  .x$ms == "D" ~ "D2",   # div/sep,   > once
  .x$ms == "W" ~ "W2"    # widowed,   > once
)
.x$MPOST <- .x$afm + 1   # one post-marriage episode for EVERY ever-married person
```

Then in `code/run.r` set `fit = 6` for every post-marriage episode: the
post-marriage chain is always the full six-state chain. The `fit`-driven
submatrix truncation was the mechanism of both defects and should be removed for
`A == 'M'` episodes. (Keep it for the pre-marriage episodes, where it correctly
forbids impossible paths.)

Sanity check to add as a unit test: for every `(afm, aai)`,
`sum(P[M, all six cells]) == 1` to machine precision.

## 2. `code/msm.cpp` — separate the first-union and later-union dissolution states

`makeQ` currently folds post-remarriage dissolutions back into the same `D` and
`W` states, so `P[M, D]` cannot distinguish "divorced from a first union" from
"divorced from a second union". Split them:

```cpp
// states: 0 V, 1 X, 2 M, 3 D1, 4 W1, 5 R, 6 D2, 7 W2
template <class Type>
struct makeQ {
  Eigen::SparseMatrix<Type> operator()(const vector<Type> &q, int size) {
    Eigen::SparseMatrix<Type> Q(size, size);
    vector<Type> rs = vector<Type>::Zero(size);
    auto add = [&](int r, int c, Type v) {
      if (r < size && c < size) { Q.coeffRef(r, c) += v; rs(r) += v; }
    };
    add(0, 1, q(0));   // V  -> X
    add(0, 2, q(1));   // V  -> M
    add(1, 2, q(2));   // X  -> M
    add(2, 3, q(3));   // M  -> D1
    add(2, 4, q(4));   // M  -> W1
    add(3, 5, q(5));   // D1 -> R
    add(4, 5, q(6));   // W1 -> R      <-- was constrained equal to q(5)
    add(5, 6, q(7));   // R  -> D2     <-- was constrained equal to q(3)
    add(5, 7, q(8));   // R  -> W2     <-- was constrained equal to q(4)
    add(6, 5, q(5));   // D2 -> R      (3rd+ unions are not observable, so
    add(7, 5, q(6));   // W2 -> R       these must stay tied to q(5), q(6))
    for (int i = 0; i < size; ++i) Q.coeffRef(i, i) = -rs(i);
    Q.makeCompressed();
    return Q;
  }
};
```

To reproduce the manuscript's constrained model, `map` in `run.r` should tie
`q(6) = q(5)`, `q(7) = q(3)`, `q(8) = q(4)` explicitly — so that the constraint
is a visible modelling choice rather than a side effect of a submatrix bound.

## 3. Centre the covariates before applying the N(0,1) priors

`msm.cpp` puts `dnorm(., 0, 1)` on `itc`, `gp_b` and `b_tm`, while `age` enters
uncentred. The intercept is therefore the log-intensity extrapolated to age 0,
where a standard normal prior is extremely informative: the corrected fit needs
`itc` for `M -> W` around `-8.2`, which that prior penalises by ~34 log units.
Centre `age` at, say, 25 and `tm` at 10, or widen the intercept priors.

## 4. Duration since marriage is minus the effect of age at marriage

Because `exp(itc + gp_b * age + b_tm * (age - afm)) = exp(itc + (gp_b + b_tm) *
age - b_tm * afm)`, the fitted `b_tm` is numerically identical to minus the
coefficient on age at first marriage. Add `afm` as an explicit covariate and
report both, or the "duration effect" will absorb any real selection on age at
marriage.

## 5. `tau_aai` / `sigma_aai` / `b_tx`

These are described in the manuscript but mapped to `NA` in `run.r`, so the
fitted model is not the model described. Either drop them from the text or
re-enable them — but see §4 of the technical note: within a single-survey
country, age at interview and birth cohort are exactly collinear, so an
aai-dependent shape parameter is not separable from a cohort trend, and the
right-truncated likelihood already conditions on age at interview.
