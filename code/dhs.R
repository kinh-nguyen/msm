#' ---
#' title: "INLA fit to Poisson vs GAM"
#' bibliography: "../zotero_biblatex.bib"
#' output:
#'   pdf_document:
#'     toc: true
#'     number_sections: true
#'     keep_tex: true
#'     include:
#'       in_header: "~/templates/floatHforRmarkdown.tex"
#' documentclass: article
#' linkcolor: blue
#' urlcolor: blue
#' citecolor: red
#' ---
#'
# -----------------------------------------------------------------------------

#+ packages and config, include=FALSE
library(rdhs)
library(tidyverse)
library(ggplot2)
library(janitor)
library(survival)
library(Epi)
library(popEpi)
library(mgcv)
tabyl <- function(dat, ...) janitor::tabyl(dat, ..., show_missing_level = FALSE)
library(knitr)
opts_chunk$set(echo = FALSE, cache = FALSE, out.extra = "")
set.seed(1)

library(INLA)
remotes::install_github("kklot/ktools")
library(ktools)
remotes::install_github("kklot/inlar")
library(inlar)

# Pull with `rdhs`
# set_rdhs_config(email = "ath19@ic.ac.uk", project = "Statistics and Machine Learning for HIV")
# datasets <- dhs_datasets("MW", fileFormat = "flat", fileType = c("IR")) %>%
#     filter(SurveyYear == 2015)  %>%
#     mutate(sex = if_else(FileType == "Individual Recode", "female", "male"))
# datasets
datasets <- readRDS("~/dhs_datasets.rds")
# downloads <- get_datasets(datasets$FileName)
# names(downloads)
downloads <- readRDS("~/dhs_downloads.rds")

x <- 1 # adapt from multiple datasets
o <- readRDS(downloads[[x]]) %>%
    as_tibble() %>%
    mutate(lab = datasets$SurveyId[x])

# Extract and converting date time 
dta <- o %>%
    select(
        psu = v021, strata = v023, weights = v005, dob = v011, doi = v008, afs = v531, lab,
        marriage_age = v511,
        n_union = v503,
        month_1st_union = v507,
        year_1st_union = v508,
        cmc_uninon = v509,
        marital_status = v501
    ) %>%
    mutate(
        marriage_age =  if_else(is.na(marriage_age), 0L, marriage_age),
        marriage_age2 = (cmc_uninon - dob) / 12,
        iso = substr(lab, 1, 2),
        yob = cmc_to_year(dob), svy = cmc_to_year(doi),
        age = svy - yob, sex = datasets$sex[x]
    )

#' Data cleaning
#'
#' - afs or marriage age greater than age
#' - married but no afs
#' 
dta %<>% filter(!(afs > age | marriage_age > age))
dta %<>% filter(!(marital_status != 0 & afs == 0))
# dta %<>% filter(afs <= marriage_age)
dta <- dta |>
    mutate(across(c(marital_status, n_union), as_factor)) %>%
    mutate(marital_status = fct_recode(
        marital_status,
        separated = "no longer living together/separated",
        union = "living with partner"
    ))

#' # Age at marriage (AAM) model
#'
#' Five possible scenarios outlined in \cref{fig:aamModel} where time *at risk of
#' getting married* is highlighted.
#'
#+ aamModel, fig.cap = "AAM risk set model", fig.height=3
bind_rows(
    tibble(start = char(birth, afs, aam), end = char(afs, aam, aai), a = 0:2, z = 1:3, risk = c(0, 1, 0)),
    tibble(start = char(birth, afs), end = char(afs, aai), a = c(0, 1), z = c(1, 3), risk = c(0, 1)),
    tibble(start = char(birth, aam, afs), end = char(aam, afs, aai), a = c(0, .7, 1), z = c(.7, 1, 3), risk = c(1, 0, 0)),
    tibble(start = char(birth, aam), end = char(aam, aai), a = c(0, .7), z = c(.7, 3), risk = c(1, 0)),
    #' no AFS
    tibble(start = char(birth), end = char(aai), a = c(0), z = c(3), risk = c(1)),
    .id = 'set'
) %>%
    mutate(risk = factor(risk, labels = char(No, Yes)), 
    across(where(is.character), str_to_upper)) %>%
    ggplot() +
        geom_segment(aes(a, 1, xend = z, yend = 1, color=risk), arrow = arrow(type = "closed")) +
        geom_text(aes(a, 1.1, label = start)) +
        geom_text(aes(z, 1.1, label = end)) +
        theme_void() +
        scale_color_manual(values = c("grey", "#FD7C02")) +
        coord_cartesian(ylim = c(.5, 1.5)) +
        facet_grid(vars(set)) +
        guides(color = 'none')
#' The first case is the typical states transition with observed marriage event;
#' the second case is right censored at the time of interview. The third and
#' fourth case are special in that the time "at risk of getting marriage" are
#' counted from birth. To reflect the difference between this exposure time and
#' the two standard cases (1 and 2), respondent's age at the exposure time are
#' taken into account.
#' 
lexis <- dta %>%
    select(age, marriage_age, afs) %>%
    mutate(
        across(where(is.integer), as.numeric),
        event = as.numeric(marriage_age != 0),
        age_0 = case_when(
            event == 1 & afs != 0 & afs <= marriage_age ~ afs, # age count from afs
            event == 0 & afs != 0 ~ afs, # age count from afs
            event == 1 & afs != 0 & afs > marriage_age ~ 0, # age count from birth
            event == 1 & afs == 0 ~ 0, # age count from birth
            event == 0 & afs == 0 ~ 0 # age count from a
        ),
        time = case_when(
            event == 1 & afs != 0 & afs <= marriage_age ~ marriage_age - afs, 
            event == 0 & afs != 0 ~ age - afs, 
            event == 1 & afs != 0 & afs > marriage_age ~ marriage_age, 
            event == 1 & afs == 0 ~ marriage_age,
            event == 0 & afs == 0 ~ age
        )
    )

#' ## Time-split data format
#'
#' The exposure time is divided by equidistant intervals of length $\tau = 1$,
#' from zero year to maximum 46 in this data. Status of the respondents at all
#' periods prior to the end of the study period (\cref{fig:aamModel}) was set to
#' single.
#'
#' Respondent's age at each of the intervals is calculated accordingly and
#' modelled as a covariate. To avoid dropping of samples with immediate
#' transition (time at risk is zero), a month was added to the age at marriage
#' if age at marriage equal age at first sex (or adding 1/12 to exposure time
#' where exposure time is zero). The example is shown in \cref{tab:data_lexis}.
#'
#+ message=F, warnings=F
#' lexis %<>% mutate(time = if_else(time == 0, time + 1 / 12, time)) # add a month
Lx <- Epi::Lexis(
    exit = list(tar = time),
    exit.status = factor(event, labels = c("Single", "Married")),
    data = lexis
)
sL <- splitMulti(Lx, tar = seq(0, 46, 1))
sL %<>% mutate(aar = age_0 + tar)

#+ data_lexis, results = "asis"
sL %>%
    head() %>%
    kable(caption = "\\label{tab:data_lexis}Lexis data format: 
    lex.id: respondent's id,
    tar: time at risk,
    lex.dur: at risk duration,
    aar: age at risk,
    age: age at interview,")

summary(sL)

#' ## Poisson model - counts data format
#'
#' Counting the events and person-year by time at risk and age at risk
#' (\cref{tab:dataCountFormat}). The empirical rates of marriage are shown in
#' \cref{fig:empirical_rate} where time since at risk is grouped for
#' visualization.
#'
psdata <- sL %>%
    group_by(tar, aar) %>%
    summarise(
        count = sum(lex.Xst == "Married"),
        person_year = n()
    ) %>%
    ungroup() %>% as.data.frame()

#+ results = "asis"
psdata %>%
    head() %>%
    kable(caption = "\\label{tab:dataCountFormat}Count data format")

caption_long <-
    "Rate of marriage per person-year by exposure time and age at risk based on
    count format data. Thick line is median across exposure time."

psdata %>%
    # grouping for plotting only
    mutate(tar2 = findInterval2(tar, c(0:10, 15, 20))) %>%
    group_by(tar2, aar) %>%
    mutate(count = sum(count), person_year = sum(person_year)) %>%
    ungroup() %>%
    mutate(rate = count / person_year) %>%
    group_by(aar) %>%
    mutate(med = median(rate)) -> empirical_rate

#+ empirical_rate, fig.cap = caption_long
empirical_rate %>%
    ggplot() +
    geom_line(aes(aar, rate, color = factor(tar2))) +
    geom_line(aes(aar, med), size = 2) +
    scale_color_viridis_d() +
    theme(legend.position = 'bottom') +
    labs(color = "Time since at risk", x = "Age at risk")

#' ## Smoothed age effect and tar - GAM model
#' 
#' To reflect the observed nonlinear age effect in \cref{fig:empirical_rate}.
#' Age effect is smoothed with `gam` as follows:
#' 
mp2 <- gam(count ~ s(aar) + s(tar), poisson(link = "log"), psdata, offset = log(person_year))
#' The model with smoothed age effect shows significant improvement, as seen in
#' in \Cref{fig:poisson_rate_separate}. Rate of marriage increases with age and taping
#' off from age 25. Consider adolescent age, marriage are more likely to occur
#' within a year exposed to risk (sexual debut in this case). 

# prediction data frame
nd <- crossing(aar = 10:30, tar = c(0,1,2,3,5,7,10))

#+ poisson_rate_merge, fig.cap = 'Estimate rate of marriage by age'
rom <- ci.pred(mp2, nd) %>%
    as_tibble() %>%
    rename_with(~ char(est, lo, up)) %>%
    bind_cols(nd) %>%
    left_join(psdata, char(aar, tar))  %>% 
    ggplot() +
    geom_point(aes(aar, count/person_year, color = factor(tar))) +
    geom_line(aes(aar, est, color = factor(tar))) +
    geom_ribbon(aes(aar, est, ymin=lo, ymax=up, fill = factor(tar)), alpha = .5) +
    labs(
        title = "Poisson model with fixed effect of TAR",
        x = "Age", y = "Rate", color = "Time since at risk"
    ) + guides(fill = 'none')
rom

#+ poisson_rate_separate, fig.cap = 'Estimate rate of marriage by age'
rom + facet_wrap(~tar, scales = 'fixed')

#' ## Poisson vs Zero-inflated poisson with INLA
#'
#' The likelihood is defined as
#'
#' $$\text{Prob}(y \mid \ldots ) = p \times 1_{[y=0]} + (1-p)\times \text{Poisson}(y)$$
#'
#' where $p$ is a hyperparameter where
#'
#' $$p = \frac{\exp(\theta)}{1+\exp(\theta)}$$
#'
#' and the Gaussian prior is defined on $\theta$.
result0 = inla(
    count ~ f(aar, model = "rw1") + f(tar, model = "rw1") + offset(ln_py),
    family = "zeroinflatedpoisson1",
    data = psdata %>% bind_rows(nd %>% mutate(count = NA)) %>% mutate(ln_py = log(person_year)),
    # data = psdata %>% filter(person_year != 1) %>% bind_rows(nd %>% mutate(count = NA)) %>% mutate(ln_py = log(person_year)),
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(config = TRUE, waic = T)
)

result1 = inla(
    count ~ f(aar, model = "rw1") + f(tar, model = "rw1") + offset(ln_py),
    family = "poisson",
    data = psdata %>% bind_rows(nd %>% mutate(count = NA)) %>% mutate(ln_py = log(person_year)),
    # data = psdata %>% filter(person_year != 1) %>% bind_rows(nd %>% mutate(count = NA)) %>% mutate(ln_py = log(person_year)),
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(config = TRUE, waic = T)
)

#+ WAIC, results = 'asis'
bind_rows(
    result1$waic[1:2] %>% bind_cols() %>% mutate(model = "Poisson"),
    result0$waic[1:2] %>% bind_cols() %>% mutate(model = "Zero-inflated")
) %>%
    relocate(model, .before = waic) %>%
    knitr::kable(caption = "\\label{tab:waic}WAIC Poisson vs Zero-inflated Poisson")

#' Estimate the probability of having no case of marriage is 7%. This is not
#' what we want. The exceed zeros in is the zero in time, not zero in number of
#' marriages in a period.
#' 
#' Let work on this to get INLA into the story. Mainly how to deal with
#' prediction without the offset: fitted values will include `offset` but does
#' not include `E`
#' 
#+ hyperpar_poisson1, results = 'asis'
result0$summary.hyperpar %>%
    select(1:2)  |> 
    knitr::kable(caption='\\label{tab:hyperpar_poisson}Poisson zero-inflated hyperparameters')
result1$summary.hyperpar %>%
    select(1:2)  |> 
    knitr::kable(caption='\\label{tab:hyperpar_poisson}Poisson hyperparameters')

#' Effect of AAR and TAR are shown in \cref{fig:inla_poisson_re}.
#'
#+ inla_poisson_re, fig.cap = 'RW2 models for age and time at risk'
autoplot(result0, "summary", "random") + facet_wrap(~var, scales = "free") +
    geom_hline(yintercept = 0)

#+ inla_poisson_fitted, fig.cap = 'Fitted values vs median - INLA'
autoplot(result0, "summary", "fitted.values") +
    geom_abline(slope = 1, intercept = 0)

#' Sample from posteriors and get the prediction for the combinations of TAR and AAR.
#' 
r.samples <- inla.posterior.sample(n = 100, result = result0, num.threads = 4)
post_latent <- r.samples %>%
    as_tibble(var = "latent") %>%
    filter(term == "Predictor") %>%
    slice(nrow(psdata) + 1:nrow(nd)) %>%
    group_by(term_id) %>%
    summarise(value = mean(value), up = quantile(value, probs = .975), lo = quantile(value, probs = 0.025)) %>%
    bind_cols(nd)

r.samples <- inla.posterior.sample(n = 100, result = result1, num.threads = 4)
post_latent1 <- r.samples %>%
    as_tibble(var = "latent") %>%
    filter(term == "Predictor") %>%
    slice(nrow(psdata) + 1:nrow(nd)) %>%
    group_by(term_id) %>%
    summarise(value = mean(value), up = quantile(value, probs = .975), lo = quantile(value, probs = 0.025)) %>%
    bind_cols(nd)

psdata %>%
    group_by(tar, aar) %>%
    mutate(count = sum(count), person_year = sum(person_year)) %>%
    ungroup() %>%
    mutate(rate = count / person_year) %>%
    group_by(aar) %>%
    mutate(med = median(rate)) -> empirical_rate

post_latent %>%
    bind_rows(post_latent1, .id = "model") %>%
    mutate(model = if_else(model == 1, "Poisson", "Zero-inflated")) %>%
    left_join(empirical_rate, char(tar, aar))  %>%
    ggplot() +
    geom_line(aes(aar, exp(value), color = factor(tar), linetype = model)) +
    geom_point(aes(aar, rate, color = factor(tar))) +
    facet_wrap(~tar)

# -----------------------------------------------------------------------------
#' 
#' # References {-}
#' 
#' <div id="refs"></div>
#' 
#' \setcounter{section}{0}
#' \renewcommand{\thesection}{\Alph{section}}
#' \setcounter{table}{0}
#' \renewcommand{\thetable}{A\arabic{table}}
#' \setcounter{figure}{0}
#' \renewcommand{\thefigure}{A\arabic{figure}}
#
#' # Appendix
#' 
# -----------------------------------------------------------------------------
#' ## Splitting time - collecting...
#'
#' The common survival analysis views time as the response variable. If survival
#' model is viewed from the demographic aspect, the basic observation is not one
#' time to event (or censoring) for each individual, but rather many periods of
#' follow up from each individual. This means modeling of rates rather than time
#' to response, where response is a dichotomous outcome in each interval.
#'
#' In this set up, time is a covariate and exposure time is response. Stratified
#' analysis becomes an interaction between time and a categorical covariate, and
#' time-varying coefficients becomes interactions between time and a continuous
#' covariate. Finally, the modeling can be done with Poisson regression
#' framework.
#'
#' The Cox-model is a special case of a Poisson model where the time covariate
#' is modelled with one parameter per failure time. Poisson modeling of disease
#' rates and follow-up studies are usually restricted to constant rate over
#' time-spans, If follow-up time are in small intervals, smoothing of rates can
#' be done by splitting the follow-up in small intervals, and modeling datasets
#' in which each record represents a period of follow-up time for a person.
#'
#' Takes the survival $X$ as response variable and censoring by time $Z$. In
#' survival data we only observe the time $\min(X,Z)$ and the event indicator
#' $\delta = 1{X < Z}$.  In a life-table, differences on the time scale are
#' accumulated as risk time whereas the position on the time scale (age) for
#' these are used as a covariate.
#'
#' Consider a study where the follow-up time for each individual is divided into
#' small intervals of equal length $y$ and each with an exit status recorded (0
#' for most of the intervals and 1 for the last interval for individuals
#' experiencing an event). Each interval contributes an observation of an
#' empirical rate, $(d, y)$, where $d$ is the number of events in the interval,
#' and $y$ is the length of the interval, i.e. the risk time. This definition is
#' slightly different from the traditional as $d/y$ (or $\sum d/\sum y$) in that
#' it keep the entire information content in the observation even if the number
#' of events is 0. This makes it usable as a response variable in all
#' situations.
#'
#' The rate of event occurrence is defined as a function of some timescale, $t$:
#'
#' $$\lambda(t) = \lim_{h\rightarrow 0} \frac{P{\text{event in} (t, t + h]|
#' \text{at risk at time} t}}{h}$$
#'
#' The rate may depend on covariates. In this formulation time $t$ is a
#' covariate and $h$ is risk time, namely the difference between two points on
#' the time ($t + h$ and $t$).
#'
#' The likelihood contribution from an observed $(d, y)$ for an interval with
#' constant rate $\lambda$ is a Bernoulli likelihood with probability $\lambda
#' y$
#'
#' $$L(\lambda | (d, y)) = (\lambda y)^d (1 - \lambda y)^{1 - d} =
#' ( \frac{\lambda y}{1 - \lambda y} )^d (1 - \lambda y)$$
#'
#' $$\log(L) = l(\lambda | (d, y)) = d \log( \frac{\lambda y}{1 - \lambda y} ) +
#' \log(1 - \lambda y) $$
#'
#' The contributions to the likelihood from one individual are conditionally
#' independent, i.e, the total likelihood from one individual across intervals
#' $t_0, t_1, t_2, t_3, t_4$ is

#' \begin{flalign}
#' P{\text{event at} t_4| alive at t0} &= P(\text{event at} t_4| alive at t_3)\\
#' &\times P(survive (t_2, t_3)| alive at t_2) \\
#' &\times P(survive (t_1, t_2)| alive at t_1) \\
#' &\times P(survive (t_0, t_1)| alive at t_0)
#' \end{flalign}

#' which is similar to independent Poisson observations likelihood. The amount
#' and spacing of events limits how detailed the rates can be modelled. The
#' model depends on how large intervals of constant rate is acceptable to the
#' study context.
#'
#' The Cox model specifies the intensity $\lambda$ as a function of time ($t$)
#' and the covariates $(x_1,...,x_p)$ through the linear predictor $\eta_i =
#' \beta_1 x_{1i} + \beta_p x_{pi}$ as:
#'
#' $$\lambda(t, x_i) = \lambda_0(t) \exp(\eta_i)$$
#'
#' with the partial log-likelihood
#'
#' $$l(\beta) = \sum_{\text{death times}} \log(\frac{e^{\eta_{death}}{\sum_{i\in
#' R_t}e^{\eta_i}}})$$
#'
#' where $R_t$ is the risk set at time $t$ including all individuals at risk at
#' time $t$
#'
#' The assumption behind the Poisson approach is essentially only the assumption
#' that a model with constant rates in each small interval gives an adequate
#' description of data. So in practice we would split data in small equidistant
#' intervals