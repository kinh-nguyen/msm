source("~/Documents/libs.r")
# devtools::install('~/git/adcomp/TMB')

d <- readRDS(here("data/d_agg.rds"))
N <- nrow(d)

d %>% 
    mutate(
        across(where(is.numeric), as.integer),
        cid = as_numeric(cc),
        J = J - 1) %>% 
    select(-cc) %>% 
    as.list() %>%
    allot(data)

(MAX_AGE = max(data$age))

data$cid %>% unique() %>% sort %>% allot(cid)
(N_CC <- length(cid))

# model matrix including countries random effect
matrix <- crossing(1, 0:MAX_AGE, cid) %>% as.matrix
data$modelmatrix <- cbind(matrix[, 1:2], make_re_matrix(matrix[, 3]))

expect_true(max(data$delta) <= MAX_AGE)
expect_true(all((data$afs + data$n_N) <= data$age))
expect_true(all((data$aam + data$n_N) <= data$age))
expect_true(max(data$J) == 6)

data$prior_base <- c(log(0.001), 0.1) # log normal mean and sd
data$prior_t <- c(0, 0.01) # mean and sd
data$prior_cc <- c(0, 0.01) # mean and sd
str(data)

N_PAR = 7
init <- list(betas = c(
    rnorm(N_PAR, data$prior_base[1], data$prior_base[2]), 
    rnorm(N_PAR, data$prior_t[1], data$prior_t[2]), 
    rnorm(N_PAR*N_CC, data$prior_cc[1], data$prior_cc[2])
))
str(init)

(N_D <- max(data$delta))

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