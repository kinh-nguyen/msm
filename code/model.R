source("~/Documents/libs.r")
# devtools::install('~/git/adcomp/TMB')
vroom::vroom(here("data/new_d.csv"), col_select = -1) %>% 
    mutate(
        cid = as_numeric(cc) - 1, 
        A = case_when(A == 'V' ~ 0, A == "X" ~ 1, A == "M" ~ 2, A == 'S' ~ 3, A == "D" ~ 4, A == "W" ~ 5, A == "R" ~ 6),
        Z = case_when(Z == 'V' ~ 0, Z == "X" ~ 1, Z == "M" ~ 2, Z == 'S' ~ 3, Z == "D" ~ 4, Z == "W" ~ 5, Z == "R" ~ 6),
        across(where(is.numeric), as.integer)
        ) %>% 
    allot(d)

data <- list()
data$prior_base <- c(log(0.001), 0.1) # log normal mean and sd
data$prior_t <- c(0, 0.01) # mean and sd

N_PAR = 7
init <- list(
    intercepts = rnorm(N_PAR, data$prior_base[1], data$prior_base[2]), 
    beta_t = rnorm(N_PAR, data$prior_t[1], data$prior_t[2])
)
str(init)

library(TMB)
TMB::compile("code/model.cpp")
base::dyn.load(TMB::dynlib("code/model"))
invisible(TMB::config(tape.parallel = 0, DLL = "model"))
TMB::openmp(1)

cc <- unique(d$cc)
dir.create(here('fitted'), F)

for (k in cc) {
    for (s in char(male, female)) {
        cat(k, s, "\n")
        tdt <- filter(d, cc == k, sex == s) %>%
            select(A, Z, start, end, n) %>%
            as.list()
        data <- modifyList(data, tdt)
        t0 <- Sys.time()
        obj <- TMB::MakeADFun(data, init, DLL = "model", silent = TRUE)
        fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 0, maxit = 500))
        rp <- obj$report(obj$env$last.par.best)
        save <- modifyList(fit, rp)
        save$runtime <- Sys.time() - t0
        saveRDS(save, here("fitted", paste0(k, s, ".rds")))
    }
}

iso2 <- list.files(here('fitted')) %>% substr(1, 2)
sex <- list.files(here('fitted')) %>% substr(3, 3)
betas <- list.files(here('fitted'), full.names = T) %>% 
    lapply(function(x) {readRDS(x) %$% rbind(intercepts, beta_t)})
modelmatrix <- as.matrix(expand.grid(1, 0:65)) 
PQ <- function(pars) {
    P <- Q <- matrix(0, 7, 7)
    Q[1, 2] = pars[1]
    Q[1, 3] = pars[2]
    Q[2, 3] = pars[3]
    Q[3, 4:6] = pars[4:6]
    Q[4:6, 7] = pars[6]
    Q[7, 4:6] = pars[7]
    for(i in 1:7) Q[i,i] = -sum(Q[i, ])
    P = expm::expm(Q)
    list(P=P, Q=Q)
}

estP <- map_dfr(seq_along(betas), function(z) {
    eta <- exp(modelmatrix %*% betas[[z]])
    PQages <- napply(1:nrow(eta), function(x) PQ(eta[x, ]))
    map_dfr(PQages, function(x) {
        tibble(
            VV = x$P[1,1],
            VX = x$P[1,2],
            VM = x$P[1,3],
            XX = x$P[2,2],
            XM = x$P[2,3],
            MM = x$P[3,3],
            MS = x$P[3,4],
            MD = x$P[3,5],
            MW = x$P[3,6],
            Re = x$P[4,7],
        )}) %>% 
        mutate(age = 0:65, cc = iso2[z], sex = sex[z])
})

estP %>% 
    filter(!(sex == 'f' & age > 49)) %>% 
    mutate(
        region = countrycode::countrycode(cc, 'dhs', 'un.regionsub.name'),
        cc = countrycode::countrycode(cc, 'dhs', 'country.name'),
        debut = VX, 
        married = (VX * XM) + VM,
        divorce = married * MD, 
        widowed = married * MW, 
        separated = married * MS, 
        remarried = married * Re, 
    ) %>% 
    allot(estPP)

ktheme <- theme(
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank())

estPP %>% 
    ggplot(aes(age, VV, color = sex)) + 
    geom_line() + 
    facet_wrap(~cc) +
    ktheme +
    guides(color = guide_legend(,direction = 'horizontal', ncol = 2)) +
    theme(legend.position = c(.9, 1.07)) +
    labs(title = "Probability of staying virgin", y = '', linetype = "Transition")
ggsave(here('fig/PVV.pdf'), width = 7, height = 7)

estPP %>% 
filter(sex == 'm') %>% 
    select(age, sex, cc, `Debut|Virgin` = VX, `Married|Virgin` = VM, `Married|Debut` = XM) %>% 
    pivot_longer(4:6) %>% 
    ggplot(aes(age, value, color = name)) + 
    geom_line() + 
    facet_wrap(~cc) +
    ktheme +
    guides(color = guide_legend(,direction = 'horizontal', ncol = 3)) +
    theme(legend.position = "bottom") +
    labs(title = "Probability of sexual debut or married - Male", y = '', color = "Transition") -> PVXVMXM_male
ggsave(here('fig/PVXVMXM_male.pdf'), PVXVMXM_male, width = 9, height = 7.5)
  
estPP %>% 
filter(sex == 'f') %>% 
    select(age, sex, cc, `Debut|Virgin` = VX, `Married|Virgin` = VM, `Married|Debut` = XM) %>% 
    pivot_longer(4:6) %>% 
    ggplot(aes(age, value, color = name)) + 
    geom_line() + 
    facet_wrap(~cc) +
    ktheme +
    guides(color = guide_legend(,direction = 'horizontal', ncol = 3)) +
    theme(legend.position = "bottom") +
    labs(title = "Probability of sexual debut or married - female", y = '', color = "Transition") -> PVXVMXM_female
ggsave(here('fig/PVXVMXM_female.pdf'), PVXVMXM_female, width = 9, height = 7.5)
  
estPP %>% 
filter(sex == 'm') %>% 
    select(age, sex, cc, `Separate|Married` = MS, `Divorce|Married` = MD, `Widowed|Married` = MW) %>% 
    pivot_longer(4:6) %>% 
    ggplot(aes(age, value, color = name)) + 
    geom_line() + 
    facet_wrap(~cc) +
    ktheme +
    guides(color = guide_legend(,direction = 'horizontal', ncol = 3)) +
    theme(legend.position = "bottom") +
    scale_y_continuous(trans = 'log') +
    labs(title = "Probability of marital dissolution - male", y = '', color = "Transition") -> PMJ_male
ggsave(here('fig/PMJ_male.pdf'), PMJ_male, width = 9, height = 7.5)
    
estPP %>% 
filter(sex == 'f') %>% 
    select(age, sex, cc, `Separate|Married` = MS, `Divorce|Married` = MD, `Widowed|Married` = MW) %>% 
    pivot_longer(4:6) %>% 
    ggplot(aes(age, value, color = name)) + 
    geom_line() + 
    facet_wrap(~cc) +
    ktheme +
    guides(color = guide_legend(,direction = 'horizontal', ncol = 3)) +
    theme(legend.position = "bottom") +
    scale_y_continuous(trans = 'log') +
    labs(title = "Probability of marital dissolution - female", y = '', color = "Transition") -> PMJ_female
ggsave(here('fig/PMJ_female.pdf'), PMJ_female, width = 9, height = 7.5)
    
estPP %>% 
    select(age, sex, cc, `Remarried` = Re) %>% 
    ggplot(aes(age, Remarried, color = sex)) + 
    geom_line() + 
    facet_wrap(~cc) +
    ktheme +
    guides(color = guide_legend(,direction = 'horizontal', ncol = 2)) +
    theme(legend.position = "bottom") +
    scale_y_continuous(trans = 'log') +
    labs(title = "Probability of remarried from either of the dissolution states", y = '', color = "Sex") -> PR
ggsave(here('fig/PR.pdf'), PR, width = 9, height = 7.5)
    
pdftools::pdf_combine(list.files(here('fig'), "^P.*\\.pdf", full.names = T))
