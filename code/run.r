# Fit the eight-state model (see DECISIONS.md) for one country and sex.
# Interactive: edit the defaults below. Batch: Rscript code/run.r <CC> <sex>
suppressPackageStartupMessages({
  library(TMB)
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# 0 V, 1 X, 2 M, 3 D1, 4 W1, 5 R, 6 D2, 7 W2 -- must match msm.cpp
state_order <- c("V", "X", "M", "D1", "W1", "R", "D2", "W2")

cc <- vroom::vroom(here("data/cc8.csv.bz2"), show_col_types = F)

args <- commandArgs(trailingOnly = TRUE)
CC  <- if (length(args) >= 1) args[1] else 'ST'
SEX <- if (length(args) >= 2) args[2] else 'f'
N_THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
cat("country", CC, "sex", SEX, "threads", N_THREADS, "\n")
COH_REF = 1975 # birth cohort at which the cohort terms are zero

tdt <- cc %>%
  filter(cc == CC) %>% 
  filter(sex == SEX) %>% 
  mutate(start = if_else(start == 0, 1, start)) %>% # SHASH is on log age
  mutate(
    A = match(A, state_order) - 1, 
    Z = match(Z, state_order) - 1, 
    dur = end - start,
    afs = if_else(afs == 0 | is.na(afs), 100, afs), # softplus zero
    afm = if_else(afm == 0 | is.na(afm), 100, afm), # softplus zero
    # chain size per episode: M absorbing before marriage; D1, W1 absorbing
    # for episodes ending in M, D1, W1; full chain for R, D2, W2
    fit = case_when(
      A <= 1 ~ 3,
      Z <= 4 ~ 5,
      TRUE   ~ 8
    ),
    coh = (yob - COH_REF) / 10 # decades
  ) %>%
  select(sv, A, Z, yob, coh, aai, start, end, afs, afm, n, fit, dur) %>%
  filter(dur >= 1)

stopifnot(!anyNA(tdt$A), !anyNA(tdt$Z))
tdt %>% count(A, Z, fit) %>% print(n = Inf)
n_surveys <- n_distinct(tdt$sv)

data <- modifyList(as.list(tdt), list(
  age_c   = 25,  # post-marriage rates: age centred at 25
  tm_c    = 10,  # and duration since first marriage at 10 years
  itc_mu  = -4,  # prior on post-marriage intercepts at the centre: N(-4, 3)
  itc_sd  = 3,
  afs_ref = 16,  # reference woman for SIMULATE
  afm_ref = 18,
  coh_ref = 0
))

init <- list(
  mu         = c(2.6, 2.8, 2.9),
  log_sigma  = c(-2, -2, .1),
  log_nu     = c(.5, -2, -1),
  log_tau    = c(-2, -.5, .8),
  b_coh      = c(0, 0, 0),
  b_tx       = 0,
  itc        = c(0, 0, 0, -4, -6, -2), # V->X, V->M, X->M, divorce, widowhood, remarriage
  gp_b       = c(0, 0, 0),
  b_tm       = c(0, 0, 0),
  b_coh_post = c(0, 0, 0)
)

map <- list(
  itc  = factor(c(NA, NA, NA, 1, 2, 3)), # SHASH carries the level of the first three
  b_tm = factor(rep(1, 3)),              # one duration effect shared by the three post-marriage rates
  # b_tm = factor(1:3),                  # separate: estimable per country (identifiability/fisher_duration.csv)
  b_coh_post = if (n_surveys > 1) factor(1:3) else factor(rep(NA, 3)) # needs more than one survey
)

# compile only if the library is missing or older than the source; on the
# cluster code/compile.sh builds it once before the array starts
dll_file <- TMB::dynlib(here("code/msm"))
if (!file.exists(dll_file) || file.mtime(here("code/msm.cpp")) > file.mtime(dll_file))
  TMB::compile(here("code/msm.cpp"), flags = "-Wno-ignored-attributes", framework = 'TMBad')
base::dyn.load(dll_file)

invisible(TMB::config(tape.parallel = 1, DLL = "msm"))
TMB::openmp(N_THREADS)

obj <- TMB::MakeADFun(data, init, map = map, DLL = "msm", silent = FALSE)

fit <- nlminb(
  obj$par, obj$fn, obj$gr,
  control = list(trace = 0, eval.max = 1000, abs.tol = 1e-12, iter.max = 1000)
)

fit
split(unname(fit$par), names(fit$par))
sdr <- TMB::sdreport(obj)
summary(sdr, "fixed")

rp <- obj$simulate(obj$env$last.par.best)
dir.create(here('fit'), showWarnings = FALSE)
saveRDS(list(fit = fit, sdr = sdr, rp = rp, map = map, data_settings = data[c("age_c", "tm_c", "itc_mu", "itc_sd", "afs_ref", "afm_ref", "coh_ref")]),
        here('fit', paste0('fit8_', CC, SEX, '.rds')))

# fitted one-year transition probabilities for the reference woman
P_yearly = rp$PP # 8 x 8 x age, states in state_order

pe = tibble(
  start = 1:dim(P_yearly)[3],
  VV = P_yearly[1,1,],
  VX = P_yearly[1,2,],
  VM = P_yearly[1,3,],
  XX = P_yearly[2,2,],
  XM = P_yearly[2,3,],
  MM = P_yearly[3,3,],
  MD = P_yearly[3,4,], # M -> D1
  MW = P_yearly[3,5,], # M -> W1
) %>% 
pivot_longer(-start)

if (interactive()) print(
  ggplot(pe, aes(start, value, color = name)) + geom_line() +
  facet_wrap(~name, scales = 'free_y')
)

# compare with the empirical rates (first unions only)
if (file.exists(here('data/pd.csv.bz2')) && file.exists(here('data/mjs.csv'))) {
pd <- vroom::vroom(here('data/pd.csv.bz2'), show_col_types = FALSE)
mj <- vroom::vroom(here('data/mjs.csv'), show_col_types = FALSE)

sex_lab <- if (SEX == 'f') 'female' else 'male'
g <- pd %>% filter(cc == CC, sex == sex_lab) %>% 
  bind_rows(
    mj %>% filter(cc == CC, sex == sex_lab)
  ) %>% 
  mutate(name = paste0(A,Z)) %>%
  left_join(pe, c("start", "name")) %>% 
  mutate(
    name = factor(name, 
    levels = c("VV", "VX", "XX", "VM", "XM", "MM", "MS", "MD", "MW"))
  ) %>% 
  ggplot(aes(start, p_empirical, color = name)) +
  geom_step() +
  geom_line(aes(y = value), linewidth = 1.2) +
  facet_wrap(~name, scale = 'free_y')

ggsave(here('fig', paste0('fit8_', CC, SEX, '.pdf')), g, width = 7, height = 7)
}
cat("done", CC, SEX, "convergence", fit$convergence, "objective", fit$objective, "\n")
