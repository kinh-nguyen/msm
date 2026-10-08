# MultistageSurv decisions log

Started 2026-10-08 while working through HANDOVER.md. The aim is a set of transition rates to parameterise a life-course HIV transmission model in the style of EPP-ASM, extended with sexual debut, marriage, divorce and widowhood.

## Decisions

The MM claim (HANDOVER §3) is accepted. Women still in their first union contributed log 1 = 0, because fit = 3 made M absorbing. Their contribution is the survival in M, exp(−∫(λ_MD + λ_MW)).

The model is the middle model. The likelihood uses eight observation states V, X, M, D1, W1, R, D2, W2. Ties are R→D2 = M→D1, R→W2 = M→W1, and W1→R = D2→R = W2→R = D1→R. Remarriage after divorce is identified (per-country s.e. of log ρ about 0.20). Remarriage after widowhood is not separately identified without v538.

The output states for the HIV model are V, X, U, D, W, with U = M + R, D = D1 + D2, W = W1 + W2. The merge is exact under the ties.

The duration effects b_tx (time since debut, X→M) and b_tm (time since first marriage) are kept, together with the age effects. b_tm is to be described as the effect of duration net of selection on age at marriage. In the HIV model, durations are carried either by precomputed duration-averaged hazards by age and cohort (first choice) or by strata of age at debut and age at marriage (check).

Duration coefficients, for reference: three separate coefficients (divorce, widowhood, remarriage) are estimable per country by Fisher information, with median s.e. 0.011, 0.009, 0.013 and worst (ST) 0.037, 0.032, 0.043. See identifiability/fisher_duration.{py,csv}. Fallback is to tie divorce and widowhood. Not yet implemented.

Countries are fitted separately at first, pooled only where it helps identification. Sex is fitted separately or derived.

The cohort term is a linear trend in birth year on the SHASH location μ, for each of the three pre-marriage transitions, in every country. A linear cohort term on the intercept of the post-marriage rates is used only in countries with more than one survey (28 of 37). In the projection the trend is held constant beyond the observed cohorts. First check against the empirical median ages at debut and at marriage by cohort. RW2 later.

Cohort check done (checks/cohort_medians.*, from data/cleaned.rds, weighted Kaplan-Meier medians by five-year cohort). A linear trend in log median age fits well. Median residual is 0.8% of the median age for debut and 1.2% for first marriage, about 0.15 and 0.2 years. A quadratic lowers these to 0.5% and 0.9%. The largest departures are first marriage in GH, ET, SZ and NG and debut in ET and MR, at most about 2.4%. Median age at first marriage rises by 2.2% per decade in the median country, debut by 1.0%, so the interval between debut and marriage widens across cohorts.

tau_aai and sigma_aai are dropped from the model and the paper. Whether cohort should later also act on spread (σ) or skewness (ν) is deferred.

Post-marriage hazards stay log-linear in age (Gompertz form), with age centred at 25 and duration at 10. The intercept prior is N(−4, 3), set in run.r, because centring alone leaves the old N(0,1) prior far from the log rates.

The survey design is handled by a design effect per survey, computed from v001 and v022 for the marital cells and age at marriage, dividing the Kish weights. A cluster bootstrap for one or two countries serves as a check.

v538 (how the previous union ended) is deferred until the base middle model fits. It would then release the tie W1→R = D1→R only in surveys that have it, after a per-survey check of the roughly 24 samples. v511a (DHS-8) is deferred with it.

## Changes applied

First edits, now superseded by the rewrites below: in run.r M→M episodes got fit = 5 and itc = rep(0, 6); in prep.Rmd women with n_union == "oneORmore" (missing v503) are dropped before cleaned.rds is saved (kept).

code/prep.Rmd rewritten for the eight-state model (2026-10-08, untested, backup of the previous version not kept in the repo). It extracts v001/v022 (mv001/mv022), uses completed age floor((doi − dob)/12) instead of round, normalises Kish weights by survey and sex (previously by survey only, pooling the IR and MR files), codes the observed cell M, D1, W1, R, D2, W2, divides the Kish weights by the clustering design effect per survey and sex (written to data/deff.csv), drops the person-level fit, and writes data/cc8/ and data/cc8.csv.bz2. The old data/cc6.csv.bz2 is kept for comparison with the current fit.

code/msm.cpp rewritten for the eight-state model (2026-10-08, not compiled). Generator with explicit ties, tau_aai/sigma_aai and the aai data removed, cohort shift b_coh on the SHASH location of the three pre-marriage rates, b_coh_post on the post-marriage intercepts, post-marriage age and duration centred at age_c and tm_c, intercept prior N(itc_mu, itc_sd) passed as data, SIMULATE reports PP (8×8 annual transition matrices) and QQ (rates by age) for a reference person afs_ref, afm_ref, coh_ref. b_tm has three entries, so run.r chooses shared or separate by the map. Centring alone does not make the old N(0,1) intercept prior harmless, because log rates at age 25 are about −4 (divorce), −6 (widowhood) and −2 (remarriage); hence the prior is now set in run.r.

Folder reorganised on 2026-10-08. code/ holds the current pipeline only (prep.Rmd, run.r, run.sh, msm.cpp). The fitting code on the cluster (msm.cpp, run.r, compile.r, check_fits.r) no longer depends on ktools; ktools.hpp moved to archive/code, where the archived models include it. paper/ holds paper.qmd and its build files. notes/ holds exploratory documents. archive/ holds superseded code and old renders, with empty folders in archive/empty. code/fit/testMW.rds became fit/testMW_code.rds. notes/FPapprox.qmd now compiles notes/dfp2.cpp. archive/code/model_code.cpp and msm_mw_code.ipynb are identical duplicates and can be deleted.

code/run.r rewritten for the eight-state model (2026-10-08, not run). Reads data/cc8.csv.bz2, sets fit per episode (3 before marriage, 5 for episodes ending in M, D1, W1, 8 for R, D2, W2), cohort coh = (yob − 1975)/10, passes age_c = 25, tm_c = 10, itc prior N(−4, 3) and the reference woman (afs 16, afm 18, cohort 1975), estimates b_tx, keeps b_tm shared (separate as a commented alternative), maps b_coh_post off when the country has one survey, saves fit, sdreport and simulation to fit/fit8_<CC><sex>.rds and the comparison plot to fig/fit8_<CC><sex>.pdf. The plotting indices for M→D and M→W were also corrected (the old code took P[3,5], which was W, as MD).

Cluster scripts (2026-10-08). code/run.r takes country and sex as arguments (Rscript code/run.r CC sex) and compiles only when msm.so is missing or older than msm.cpp. code/compile.sh compiles once on a compute node, code/run.sh is a Slurm array with one task per line of code/tasks_f.txt or code/tasks_m.txt (37 countries each). On the Mac, code/hpc_submit.sh syncs code and data to the cluster and submits the compile job and the dependent array; code/hpc_fetch.sh copies fit/fit8_*, fig/fit8_* and the logs back. Defaults HPC_HOST=fuchs and HPC_DIR=/panfs/vdura1/fuchs/fias/knguyen/MultistageSurv, overridable by environment variables.

HPC_INSTRUCTIONS.md (2026-10-08) instructs the cluster to check the setup, compile, fit São Tomé as a test, report, and only after approval run all 37 countries for women, with code/check_fits.r summarising the fits in fit/fit8_summary.csv.

prep.Rmd adam chunk (2026-10-08): the rdhs query had no regional filter and surveyYearStart = 2015, so it pulled every country (India failed to read). It now takes the 100 surveys listed in code/surveys.txt, the surveys present in the existing data/cleaned.rds (99 for women, ZM2002 for men only).

prep.Rmd download and design-effect fixes (2026-10-08, found while rerunning). The DHS login comes from code/rdhs.json, read with rdhs's own reader rdhs:::read_rdhs_config_file() so that data_frame becomes a function, and the client from client_dhs() is passed explicitly to dhs_datasets() and used for client$get_datasets(); the email, project and password-like comment were removed from prep.Rmd (the password is still in earlier commits and should be changed). Flat files whose .DCF dictionary is latin1 (Benin 2006, Cameroon 2004 and others) cannot be parsed by rdhs, which reads them with brio::read_lines as UTF-8 whatever the locale; for those the Stata file of the same recode is used, matched on the upper-case file stem. The design effect uses one 0/1 indicator per cell instead of a factor (a factor fails when a survey has one level) and deff = "replace" (the Kish-scale weights must not be read as a population size); each statistic is wrapped in try().

README.md rewritten for the current project (2026-10-08); the 2022 README is archive/README_2022.md.

## Still to do

- [x] Finish rerunning prep.Rmd (it reached the design effects on 2026-10-08) and check which surveys came from Stata files, count(tmp, cell, u), data/deff.csv and new_d %>% count(A, Z). 
- [ ] Compile and fit one country, then compare with the old fit by country. Decide shared or separate b_tm. 
- [ ] Then RW2 for the cohort term if needed, v538, men, pooling across countries where needed
- [ ] the duration-averaged hazard tables by age and cohort for the HIV model.

## Notes for the Methods section

2,761 ever-married women (0.3%) with missing v503 were excluded, assumed missing at random. Previously they were coded as remarried.
225 women reporting a union but no age at marriage were excluded, assumed missing at random. Their correct contribution would be interval-censored.
All rates are conditional on survival to interview, since there is no death state.
Remarriage is a rate averaged over divorce and widowhood, weighted toward divorce, because of the tie W1→R = D1→R.
Time since dissolution is not observed in any DHS phase, so remarriage depends on age and time since first marriage only.
