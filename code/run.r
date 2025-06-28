library(ktools)
library("TMB")

state_order <-c("V", "X", "M", "S", "D", "W", "R")

vroom::vroom(here("data/new_dw.csv.bz2")) %>% 
  filter(cc == "CG", sex == 'female', sv == 2011) %>% 
  filter(start > 0) %>%
  mutate(
    A = match(A, state_order) - 1, 
    Z = match(Z, state_order) - 1, 
    fit = case_when(A < 2 ~ 3, otherwise ~ fit), 
    yob = scale(yob), 
    dstart = as.double(start),
    istart = start,
    tx = start - afs, 
    tm = start - afm,
    tx = if_else(tx < 0 | is.na(tx), 0, tx),
    tm = if_else(tm < 0 | is.na(tm), 0, tm),
    aai = scale(aai)
  ) %>%
  select(sv, A, Z, yob, aai, istart, dstart, end, tx, tm, n, fit) %>%
allot(tdt)
# predicting data by age for average duration of tx and tm
tdt %>%
  filter(A == 1 | A == 2) %>%
  group_by(start) %>%
  summarise(
    tx = mean(tx),
    tm = mean(tm),
    .groups = "drop"
  ) %>% 
  bind_rows(
    anti_join(tibble(
      start = 1:(max(tdt$end) + 1),
      tx = 0, tm = 0,
    ), ., by = "start")
  ) %>%
  arrange(start) %>%
allot(sim_data)

N_PAR = 7
data <- list()
data$n_age <- max(tdt$end) + 1 # include zero
data$sim_data = as.matrix(sim_data)

data <- modifyList(data, as.list(tdt))

init <- list(
  mu = c(2.6, 2.8, 2.9, 2.7, 2.7, 2.7, 2.7),
  log_sigma = c(-2, -2, .1, -2, -2, -2, -2),
  log_nu = c(.5, -2, -1, -2, -2, -2, -2),
  log_tau = c(-2, -.5, .8, -2, -2, -2, -2),
  tau_aai = c(-.04, -.05, -.2, -.05, -.05, -.05, -.05),
  b_tx = c(-0.03),
  b_tm = c(-0.03),
  itc = c(3, -.5, 2, 3, 3, 3, 3)
)
TMB::compile("code/msm.cpp", flags = "-Wno-ignored-attributes")
base::dyn.load(TMB::dynlib("code/msm"))
invisible(TMB::config(tape.parallel = 0, DLL = "msm"))
TMB::openmp(20)

obj <- TMB::MakeADFun(
  data, init,
  DLL = "msm", silent = FALSE 
)

fit <- nlminb(
  obj$env$last.par.best, obj$fn, obj$gr,
  control = list(trace = 0, eval.max = 1000, abs.tol = 1e-12, iter.max = 1000)
)
fit
fit$par %>% name2list()
rp <- obj$simulate(obj$env$last.par.best)
plott(rp$PP[1,2,2:51],width = 70)
plott(rp$PP[1,3,2:51],width = 70)
plott(rp$PP[2,3,2:51],width = 70)
dir.create(here('fit'), showWarnings = FALSE)
saveRDS(rp, here('fit/testMW.rds'))

P_yearly = rp$PP

pe = tibble(
  start = 1:dim(P_yearly)[3],
  VV = P_yearly[1,1,],
  VX = P_yearly[1,2,],
  VM = P_yearly[1,3,],
  XX = P_yearly[2,2,],
  XM = P_yearly[2,3,],
) %>% 
pivot_longer(-start)

tdt %>%
  filter(A < 2, Z <= 2) %>%
  group_by(A, start, tx) %>%
  summarise(n_risk = sum(n), .groups = "drop") %>%
  allot(at_risk)

tdt %>%
  filter(A < 2, Z <= 2) %>%
  group_by(A, Z, start, tx) %>%
  summarise(n_move = sum(n), .groups = "drop") %>%
  left_join(at_risk, by = c("A", "start", "tx")) %>%
  mutate(p_empirical = n_move / n_risk) %>% 
  mutate(
    A = case_when(A == 0 ~ 'V', A == 1 ~ 'X', A == 2 ~ "M", A == 3 ~ "S", A == 4 ~ "D", A == 5 ~ "W", A == 6 ~ 'R', otherwise ~ NA_character_),
    Z = case_when(Z == 0 ~ 'V', Z == 1 ~ 'X', Z == 2 ~ "M", Z == 3 ~ "S", Z == 4 ~ "D", Z == 5 ~ "W", Z == 6 ~ 'R', otherwise ~ NA_character_),
    name = paste0(A,Z)) %>% 
  full_join(pe, by = c("start", "name")) %>% 
allot(pd)

g <- pd %>%
  mutate(group = case_when(
    A == "V" ~ 1,
    A == "X" ~ 2,
    A == "M" ~ 3,
  ), 
    p_empirical = ifelse(A == "M", NA_real_, p_empirical),
  ) %>%
  filter(A != "M") %>%
  filter(start > 0) %>% 
  ggplot(aes(start, p_empirical)) +
  geom_point(aes(size = n_risk, color = factor(name)), shape = 21) +
  geom_line(aes(y = value, color = factor(name))) +
  ggh4x::facet_grid2(
    cols = vars(group),
    scales = "free_y",
  ) 

ggsave('g.pdf', g, width = 7, height = 7)
