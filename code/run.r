library(ktools)
library("TMB")

state_order <-c("V", "X", "M", "S", "D", "W", "R")

cc <- vroom::vroom(here("data/cc.csv.bz2"), show_col_types = F)

CC = 'ST'  
SEX = 'f'  

cc %>%
  filter(cc == CC) %>% 
  filter(sex == SEX) %>% 
  mutate(start = if_else(start == 0, 1, start)) %>% 
  mutate(
    A = match(A, state_order) - 1, 
    Z = match(Z, state_order) - 1, 
    afs = if_else(afs == 0 | is.na(afs), 100, afs), # softmax zero
    afm = if_else(afm == 0 | is.na(afm), 100, afm), # softmax zero
  ) %>%
  select(sv, A, Z, yob, aai, start, end, afs, afm, n, fit) %>%
allot(tdt)

tdt <- tdt %>% mutate(across(c(aai, yob), ~ scale(.x)[,1]))

N_PAR = 7
data <- list()

data <- modifyList(data, as.list(tdt))

init <- list(
  # VX VM XM
  mu = c(2.6, 2.8, 2.9),
  log_sigma = c(-2, -2, .1),
  log_nu = c(.5, -2, -1),
  log_tau = c(-2, -.5, .8),
  tau_aai = c(0, 0, 0),
  b_tx = c(-0.03),
  # MS MD MW M...R
  b_tm = c(-0.03),
  itc = c(3, -.5, 2, 3, 3, 3, 3),
  gp_b = c(0, 0, 0, 0)
)

TMB::compile("code/msm.cpp", flags = "-O3 -Wno-ignored-attributes", framework = 'TMBad')
base::dyn.load(TMB::dynlib("code/msm"))
invisible(TMB::config(tape.parallel = 1, DLL = "msm"))
TMB::openmp(20)

obj <- TMB::MakeADFun(
  data, init,
  map = list(
    tau_aai = factor(rep(NA, 3))
  ),
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
  MM = P_yearly[3,3,],
  MS = P_yearly[3,4,],
  MD = P_yearly[3,5,],
  MW = P_yearly[3,6,],
) %>% 
pivot_longer(-start)

# empirical up to M
tdt %>%
  rename(start = dstart) %>% 
  filter(A < 2, Z <= 2) %>%
  group_by(A, start, tx) %>%
  summarise(n_risk = sum(n), .groups = "drop") %>%
  allot(at_risk)

tdt %>%
  rename(start = dstart) %>% 
  filter(A < 2, Z <= 2) %>%
  group_by(A, Z, start, tx) %>%
  summarise(n_move = sum(n), .groups = "drop") %>%
  left_join(at_risk, by = c("A", "start", "tx")) %>%
  mutate(p_empirical = n_move / n_risk) %>% 
allot(pd)

# stepwise empirical for MS MD MW
purrr::map_dfr(3:5, \(x) {
  tdt %>% 
  filter(A == 2, fit != 7) %>% 
    mutate(
      l =  dstart,
      u =  case_when(
        A == Z ~ Inf,
        Z != x ~ Inf,
        Z == x ~ end,
      )
    ) -> tmp
    icenReg::ic_np(tmp[, c('l', 'u')], weights = tmp$n) %>% 
    hz_np()
  }, .id = 'Z') %>% 
  transmute(
    A = 2, Z = as.numeric(Z) + 2, 
    start = interval_start, 
    p_empirical = hazard_rate
  ) %>% 
allot(steps_mj)

pd %>% 
  bind_rows(steps_mj) %>%   
  filter(A != Z, start > 0, start < 45) %>% 
  mutate(
    A = state_order[A + 1],
    Z = state_order[Z + 1],
    name = paste0(A,Z),
  ) %>%
  left_join(pe, char(start, name)) %>% 
  mutate(
    name = factor(name, levels = char(VX, VM, XM, MS, MD, MW))
  ) %>% 
  ggplot(aes(start, p_empirical, color = name)) +
  geom_step() +
  geom_line(aes(y = value), linewidth = 1.2) +
  facet_wrap(~name, scale = 'free_y')

ggsave('g.pdf', g, width = 7, height = 7)
