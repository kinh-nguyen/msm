# Instructions for running the MultistageSurv fits on the cluster

These instructions are for a Claude session that is already connected to the Slurm cluster. Work only inside the project folder on the cluster, referred to below as `$HPC_DIR` (by default `/scratch/fuchs/fias/knguyen/MultistageSurv`). Read `DECISIONS.md` first; it records why the model has its present form. Your task is to compile and run the model, check the results, and report. It is not to change the model.

## What the code does

`code/msm.cpp` is a TMB model of women's (or men's) sexual and marital careers as an eight-state Markov chain in age: never had sex (V), had sex but never married (X), in first union (M), divorced or widowed after one union (D1, W1), in a later union (R), and divorced or widowed after more than one union (D2, W2). There are six distinct rates, tied as described in the comments of `makeQ`. `code/run.r` fits one country and sex from `data/cc8.csv.bz2` and saves `fit/fit8_<CC><sex>.rds` (the `nlminb` result, the `sdreport`, and simulated annual transition matrices) and, when `data/pd.csv.bz2` and `data/mjs.csv` are present, a comparison plot `fig/fit8_<CC><sex>.pdf`. `code/compile.sh` builds the library once. `code/run.sh` is a Slurm array with one task per line of a task file, each line `<CC> <sex>`. `code/check_fits.r` summarises all fits into `fit/fit8_summary.csv`.

## Rules

Do not edit `code/msm.cpp`, `code/run.r`, `code/prep.Rmd` or anything in `data/`, except as allowed under "Compilation" below. Do not change priors, the `map`, the state order, the ties or the centring constants. Do not delete or overwrite earlier fits; if a fit file already exists, rename it with a date suffix before rerunning. Do not commit or push to git. Ask the user before submitting the full array, before running men, and before requesting more than 8 cores, 16 GB of memory or 8 hours per task. Stop and report rather than improvise when something fails in a way these instructions do not cover.

## Step 1, check the setup

The code comes from the git repository `github.com/kinh-nguyen/msm`. If `$HPC_DIR` is a clone of it, run `git pull` there; if it does not exist, clone the repository into it. The data are not in git: the user copies `data/cc8.csv.bz2` (and, if available, `data/pd.csv.bz2` and `data/mjs.csv`) from their computer. Create an empty file `.here` in `$HPC_DIR` so that `here()` finds the project root.

Confirm that `$HPC_DIR` contains `.here`, `code/msm.cpp`, `code/ktools.hpp`, `code/run.r`, `code/compile.r`, `code/compile.sh`, `code/run.sh`, `code/check_fits.r`, `code/tasks_f.txt` and `data/cc8.csv.bz2`, and create `logs/`, `fit/` and `fig/` if they are missing. Check that the Slurm scripts suit this cluster: the partition (`fuchs`), the `source` line for `.bashrc` and `conda activate kinh`. If the partition or environment name differs, tell the user and correct only those lines in `code/run.sh` and `code/compile.sh`.

Check the R environment on a compute node or in an interactive allocation, not on the login node:

```bash
source ~/.bashrc; conda activate kinh
Rscript -e 'for (p in c("TMB","ktools","here","vroom","dplyr","tidyr","ggplot2","purrr")) cat(p, as.character(packageVersion(p)), "\n")'
```

A missing package is a reason to stop and report, not to install into the shared environment without asking.

## Step 2, compilation

```bash
cd $HPC_DIR
sbatch code/compile.sh
```

When it finishes, `code/msm.so` must exist and `logs/msm8_compile_<jobid>.out` must end without an error. If compilation fails, report the first error message with its line number. You may fix an error yourself only if it is purely a C++ or TMB syntax or type problem whose fix does not change any computed quantity (for example converting a vector type before an assignment). Record every such change, with the before and after lines, in the report.

## Step 3, one test country

```bash
printf "ST f\n" > code/tasks_test.txt
sbatch --array=1-1 code/run.sh code/tasks_test.txt
```

São Tomé is the smallest sample and should finish quickly. In `logs/msm8_<jobid>_1.out` check the following. The table printed by `count(A, Z, fit)` has `fit = 3` only for rows with `A` 0 or 1, `fit = 5` for rows with `A = 2` and `Z` in 2, 3, 4, and `fit = 8` for `Z` in 5, 6, 7. The last line reads `done ST f convergence 0` with a finite objective. Then run `Rscript code/check_fits.r` and check that `max_abs_grad` is below about 1e-3, `pd_hessian` is TRUE, `n_se_missing` is 0, and `max_rowsum_err` is below 1e-8. The last check confirms that each simulated annual transition matrix has rows summing to one.

Look also at the estimates in the log (`fit$par` and `summary(sdr, "fixed")`). The three estimated `itc` values (divorce, widowhood, remarriage) are log rates at age 25 with ten years since first marriage. Values roughly between −7 and −1 are plausible; values pinned near the prior mean of −4 with standard errors near 3 mean the data say almost nothing, and values below −15 or above 2 mean something is wrong. `b_coh` is the shift per decade of birth cohort in log median age; values within about ±0.1 are plausible.

Report the results of Step 3 to the user and wait for approval before Step 4.

## Step 4, all countries, women

After approval:

```bash
jid=$(sbatch --parsable --array=1-$(grep -c . code/tasks_f.txt)%10 code/run.sh code/tasks_f.txt)
echo $jid
```

Monitor with `squeue -u $USER` and `sacct -j $jid --format=JobID,State,Elapsed,MaxRSS`. When all tasks have finished, run `Rscript code/check_fits.r`. A task that hit the time limit or ran out of memory may be resubmitted once alone, with the resources asked for raised after asking the user. A task that fails for any other reason is not resubmitted; collect its error for the report.

## Step 5, report

Write `logs/REPORT_<YYYY-MM-DD>.md` and give the user its content. It should state the Slurm job identifiers; for each country the state, elapsed time, peak memory and the columns of `fit/fit8_summary.csv`; the countries that did not converge, have a non-positive-definite Hessian or missing standard errors, with the relevant log lines; any change made under Step 1 or Step 2; and any estimate outside the plausible ranges in Step 3. Keep interpretation brief and separate from the facts. Men (`code/tasks_m.txt`) are run only if the user asks.
