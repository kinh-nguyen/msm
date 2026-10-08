# MultistageSurv

Rates of sexual debut, first marriage, union dissolution (divorce or separation, and widowhood) and remarriage for 37 sub-Saharan African countries, estimated from 100 DHS surveys. The estimates are meant to parameterise a life-course HIV transmission model in the style of EPP-ASM, extended with the marital states V, X, U, D and W.

`DECISIONS.md` records every modelling decision and change since October 2026 and is the place to start. `HANDOVER.md` is the earlier review (August 2026) that led to them. `HPC_INSTRUCTIONS.md` tells a Claude session on the cluster how to run the fits.

## Model

Each woman (or man) is followed in age through an eight-state Markov chain: V never had sex and never married, X had sex but never married, M in first union, D1 and W1 first union ended by divorce or widowhood without remarriage, R in a later union, and D2 and W2 a later union ended by divorce or widowhood. These are exactly the cells DHS observes through v501 (marital status) and v503 (number of unions). There are six distinct rates. Debut (V→X), marriage at debut (V→M) and marriage after debut (X→M) have sinh-arcsinh hazards in age, with a linear birth-cohort shift of their location and, for X→M, an effect of time since debut. Divorce, widowhood and remarriage are log-linear in age, time since first marriage and, in countries with more than one survey, birth cohort. Later-union rates are tied to the first-union ones (R→D2 = M→D1, R→W2 = M→W1, all returns to union = D1→R), which is what makes the model estimable country by country and is also the structure of the HIV model. Merging M and R into U, D1 and D2 into D, and W1 and W2 into W is exact under these ties.

The likelihood is built from episodes. Before marriage, ages at debut and marriage are treated as reported events; after marriage each woman contributes the probability of her observed cell at interview, computed by matrix exponentials of the generator year by year. Weights are Kish-normalised by survey and sex and divided by a clustering design effect per survey.

## Layout

```
code/             current pipeline
  prep.Rmd        DHS download (rdhs), cleaning, episodes, design effects -> data/cc8.csv.bz2
  msm.cpp         TMB model
  run.r           fit one country and sex: Rscript code/run.r <CC> <sex>
  check_fits.r    summary of all fits -> fit/fit8_summary.csv
  compile.sh, run.sh            Slurm: compile once, then one array task per country
  hpc_submit.sh, hpc_fetch.sh   run on the Mac: sync to the cluster and submit, copy results back
  tasks_f.txt, tasks_m.txt      country and sex per array task
  surveys.txt     the 100 surveys used
  rdhs.json       DHS login and cache settings (git-ignored, holds the password)
identifiability/  technical note on identification, Fisher-information audit, PATCH.md
checks/           empirical median ages at debut and marriage by birth cohort
paper/            paper.qmd and its build files
notes/            exploratory documents
archive/          superseded code and renders, including the 2022 README
data/             DHS extracts and derived data (git-ignored)
fig/, fit/        figures and fitted models
```

## Running

First run `code/prep.Rmd` on a machine with DHS access. The `adam` chunk reads `code/rdhs.json`, downloads the individual and men's recode files of the surveys in `code/surveys.txt` (falling back to the Stata file when rdhs cannot read a flat file), and the later chunks write `data/cleaned.rds`, `data/deff.csv` and `data/cc8.csv.bz2`. Check `count(tmp, cell, u)`, the design effects in `data/deff.csv` (mostly between 1 and 3) and `new_d %>% count(A, Z)`.

To fit one country locally, run `Rscript code/run.r ST f`. To fit all countries on the cluster, copy `data/cc8.csv.bz2` there and either follow `HPC_INSTRUCTIONS.md` or run `bash code/hpc_submit.sh` from the Mac (set `HPC_HOST` and `HPC_DIR` if the defaults do not apply), then `bash code/hpc_fetch.sh` and `Rscript code/check_fits.r`.

## Status (8 October 2026)

The eight-state model is written in `prep.Rmd`, `msm.cpp` and `run.r`. `prep.Rmd` is being rerun; `msm.cpp` has not yet been compiled and no country has been fitted. The open items are listed under "Still to do" in `DECISIONS.md`.
