source("~/Documents/libs.r")
# devtools::install('~/git/adcomp/TMB')

MAX_AGE = 50
N_PAR = 8
N_YOB = 37 # 1965-2001

library(TMB)
TMB::compile("code/model.cpp")
base::dyn.load(TMB::dynlib("code/model"))
invisible(TMB::config(tape.parallel = 0, DLL = "model"))
TMB::openmp(1)
d <- readRDS(here("data/d_agg.rds"))
d

N <- 3000
data <- d %>% 
    mutate(across(everything(), haven::zap_label)) %>%
    mutate(across(everything(), as.integer)) %>%
    mutate(J = J - 1, yob = yob - min(yob)) %>%
    slice_head(n = N) %>% 
    as.list()

# expect_true(max(data$age) <= MAX_AGE)
expect_true(max(data$yob) <= N_YOB)
expect_true(max(data$delta) <= MAX_AGE)
expect_true(all((data$afs + data$n_N) <= data$age))
expect_true(all((data$aam + data$n_N) <= data$age))
expect_true(max(data$J) == 6)

data$sd_q <- c(log(0.001), 0.1) # log normal mean and sd
str(data)

init <- list(
    betas = rnorm(2 * N_PAR, log(0.001), 1),
    yobsm = rnorm(N_YOB * N_PAR, log(0.001), 1),
    pacf = c(0, 1)
    )
str(init)

obj <- TMB::MakeADFun(
    data, 
    init, 
    DLL = "model", silent = TRUE)
fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1, maxit = 500))
rp <- obj$report(obj$env$last.par.best); str(rp)
p_id <- char(VX, XM, MS, MD, MW, SR, DR, WR)
p_lb <- char(
    "Sexual debut", "Marriage", "Separate", "Divorce", "Widow",
    "Separated>>Remarried", "Divorced>>Remarried", "Widowed>>Remarried"
)

rp <- obj$simulate() 
estP <- rp$estP %>% array(c(7,7,51,N_YOB))

tibble(
    age = 1:51,
    VX = estP[1, 2, , 36],
    VM = estP[1, 3, , 36],
    XM = estP[2, 3, , 36],
    MS = estP[3, 4, , 36],
    MD = estP[3, 5, , 36],
    MW = estP[3, 6, , 36],
    SR = estP[4, 7, , 36],
    DR = estP[5, 7, , 36],
    WR = estP[6, 7, , 36],
) %>% 
    pivot_longer(-age) %>%
    mutate(name = factor(name, levels = char(VX, VM, XM, MS, MD, MW, SR, DR, WR))) %>%
    allot(o)

o %>% 
pivot_wider(names_from = name, values_from = value) %>% 
transmute(
    age = age,
    pX = VX,
    pM = (VX * XM) + VM, 
    pNULL = NA,
    pS = pM * MS,
    pD = pM * MD,
    pW = pM * MW,
    pSr = pS * SR,
    pDr = pD * DR,
    pWr = pW * WR,
    ) %>% 
    pivot_longer(-age) %>% 
    mutate(name = factor(name, 
        levels = char(pX, pM, pNULL, pS, pD, pW, pSr, pDr, pWr), 
        labels = char("Debut", "Marriage", "", "Separated|Married", "Divorced|Married", "Widowed|Married", "Remarried|Separated", "Remarried|Divorced", "Remarried|Widowed"), 
        )) %>% 
    ggplot(aes(age, value, color = name)) + geom_line() + facet_wrap(~name, scales = 'free') +
    scale_y_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 7)) +
    theme(
        legend.position = 'none',
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), axis.text.y = element_text(size = 5)) +
    labs(title = "Transition probability in the next year", y = '')

ggsave("fig/prob3000.png", width = 7, height = 4.5)