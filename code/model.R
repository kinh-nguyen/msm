source("~/Documents/libs.r")
# devtools::install('~/git/adcomp/TMB')

d <- vroom::vroom(here("data/new_d.csv"), col_select = -1)
N <- nrow(d)

d %>% 
    mutate(
        cid = as_numeric(cc) - 1, 
        A = case_when(A == 'V' ~ 0, A == "X" ~ 1, A == "M" ~ 2, A == 'S' ~ 3, A == "D" ~ 4, A == "W" ~ 5, A == "R" ~ 6),
        Z = case_when(Z == 'V' ~ 0, Z == "X" ~ 1, Z == "M" ~ 2, Z == 'S' ~ 3, Z == "D" ~ 4, Z == "W" ~ 5, Z == "R" ~ 6),
        across(where(is.numeric), as.integer)
        ) %>% 
    select(-cc, -sex) %T>% 
    print %>% 
    as.list() %>%
    allot(data)

data$cid %>% unique() %>% sort %>% allot(cid)
data$cid %>% max %>% mustbe(36)
(N_CC <- length(cid))

str(data)

data$prior_base <- c(log(0.001), 0.1) # log normal mean and sd
data$prior_t <- c(0, 0.01) # mean and sd
data$prior_cc <- c(0, 0.01) # mean and sd
str(data)

N_PAR = 7
init <- list(
    intercepts = rnorm(N_PAR, data$prior_base[1], data$prior_base[2]), 
    beta_t = rnorm(N_PAR, data$prior_t[1], data$prior_t[2]), 
    cc0 = rep(0, N_CC),
    cc1 = rep(0, N_CC),
    cc2 = rep(0, N_CC),
    cc3 = rep(0, N_CC),
    cc4 = rep(0, N_CC),
    cc5 = rep(0, N_CC),
    cc6 = rep(0, N_CC)
)
str(init)

library(TMB)
TMB::compile("code/model.cpp")
# base::dyn.unload(TMB::dynlib("code/model"))
base::dyn.load(TMB::dynlib("code/model"))
invisible(TMB::config(tape.parallel = 0, DLL = "model"))
TMB::openmp(1)

obj <- TMB::MakeADFun(data, init, DLL = "model", silent = TRUE)
fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1, maxit = 500))
rp <- obj$report(obj$env$last.par.best); str(rp)
p_id <- char(VX, XM, MS, MD, MW, SR, DR, WR)
p_lb <- char(
    "Sexual debut", "Marriage", "Separate", "Divorce", "Widow",
    "Separated>>Remarried", "Divorced>>Remarried", "Widowed>>Remarried"
)

vis <- function(x, rt = F) {
    X <- array(rp[[x]], c(7, 7, 50))
    tibble(
        age = 1:50,
        VX = X[1, 2, ],
        VM = X[1, 3, ],
        XM = X[2, 3, ],
        MS = X[3, 4, ],
        MD = X[3, 5, ],
        MW = X[3, 6, ],
        SR = X[4, 7, ],
        DR = X[5, 7, ],
        WR = X[6, 7, ],
    ) %>%
        pivot_longer(-age) %>%
        mutate(name = factor(name, levels = char(VX, VM, XM, MS, MD, MW, SR, DR, WR))) %>%
        allot(o)
    if (rt) return(o)
    o %>% ggplot(aes(age, value, color = name)) +
        geom_line()
}

vis("KM.masterP") + facet_wrap(~name, scales = "free")

vis("KM.masterP", 1) %>% 
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