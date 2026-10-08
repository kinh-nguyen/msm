#' ---
#' title: "Married, divorce, widowed in women - Malawi 2015"
#' bibliography: "../zotero_biblatex.bib"
#' output:
#'   pdf_document:
#'     toc: true
#'     number_sections: true
#'     keep_tex: false
#'     include:
#'       in_header: "~/templates/floatHforRmarkdown.tex"
#' documentclass: article
#' link-citations: yes
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
library(ggfortify)
library(janitor)
tabyl <- function(dat, ...) janitor::tabyl(dat, ..., show_missing_level = FALSE)
devtools::load_all("~/Code/R/ktools/")
library(knitr)
opts_chunk$set(echo = FALSE, cache = FALSE, out.extra = "")
set.seed(1)

library(INLA)
remotes::install_github("kklot/ktools")
library(ktools)
remotes::install_github("kklot/inlar")
library(inlar)

library(mgcv, help, pos = 2, lib.loc = NULL)

here::i_am("code/dhs.R")
library(here)

#' # Review of mixture survival model
#'
#' The most well-known survival mixture model is the 'cure' model
#' (@farewellMixtureModelsSurvival1986,@berksonSurvivalCurveCancer1952), where
#' a fraction of the sample (censored sample in practice) is assumed to never
#' experience the event (thus a point mass at $T = \infty$).
#'
#' Given the incidence rate $\mu_i$, Weibull's shape parameters $\nu$ and a cure
#' fraction $\rho \in (0, 1)$, the event times $T_i$ is modelled as
#'
#' $$\begin{cases} T_i \sim \text{Weibull}(\mu_i, \nu) \qquad \text{with
#'  probability} 1-\rho\\
#'  T_i = \infty \qquad \text{with probability} \rho \end{cases}$$
#'
#' The model can be fitted with `smcure` package using EM algorithm. The same
#' model can be fitted with INLA using `weibullcure` family. These versions do
#' not allow to model covarirates on the probability of being cure.
#'
#' With INLA, the rate can be modelled adding covariate and spatial correllation
#' as @jiangGeostatisticalSurvivalModels2014
#'
#' $$\mu_i = \exp(\alpha + \beta X_i + W(s_i) + U(s_i))$$
#'
#' where $W$ the factors at individual $i$'s location of residence $s_i$ and
#' $U(s_i)$ a smooth spatial surface centred at zero reflecting possible
#' additional factors that are not measured and spatial correlation between
#' subjects.
#'
#' The likelihood in the cure rate model with measurements $Y_i = (L, R, T)$
#' (left-truncated, right-censored, event time) and event indicator $Z_i$ is
#' defined as
#'
#' $$P(Y|\cdot) = \prod_{i, Z_i = 1} P(T = T_i|T>L_i,\mu_i) \prod_{i, Z_i = 0}
#' P(T > R_i|T>L_i,\mu_i)
#' $$
#'
#' with
#'
#' $$P(T = T|T>L,\mu) = \frac{(1-\rho)f(.)}{\rho + (1 - \rho)
#' \int_{L}^{\infty} f(.)du }$$
#'
#' and
#'
#' $$P(T > R|T>L,\mu) = \frac{\rho + (1 - \rho)
#' \int_{R}^{\infty} f(.)du }{\rho + (1 - \rho) \int_{L}^{\infty} f(.)du }$$
#'
#' which have closed-form solutions.
#'
#' The AAM model is the opposite and simpler than the cure model. It has exceed
#' number of 'deaths' at time zero and all those were observed in the data. The
#' censored observations did not experience marriage at time zero. So instead to
#' having a cure rate, we have a zero-inflated rate.
#'
#' But since the probability density will generally support zeros as well, it
#' will need to be modified to truncate zero value (?) if keeping the same
#' approach as cure rate model. @louzadaZeroinflatedNonDefault2018 did
#' this to model loan (straight to default group - loan time is zero), in
#' addition to "cure". In particular,
#'
#' \begin{equation}
#' S(t) = p_1 + (1 - p_0 - p_1) S_0(t) \\
#' F(t) = p_0 + (1 - p_0 - p_1) F_0(t) \\
#' F(0) = p_0, \qquad \lim_{t\rightarrow \infty} F(t) = 1-p_1
#' \end{equation}
#'
#' $$f(t) = \begin{cases}
#' p_0, \qquad t=0\\
#' (1 - p_0 - p_1) f_0(t), \qquad t > 0
#' \end{cases}
#' $$
#'
#' where $p_0$ the proportion of zero-inflated times, $p_1$ the proportion of
#' cure rate. The censored observations in this model have likelihood
#' contributions
#'
#' $$p_1 + (1 - p_0 - p_1) S_0(t)$$
#'
#' Covariates effect on the proportions was modelled with a multinomial
#' logistic regression. There was no mention of truncated distribution in the
#' paper.
#'
#' In AAM case, the cure part can be removed and the model becomes similar to
#' cure rate model with exception that *the distribution needs to be truncated*,
#' i.e. $T | T > 0 \sim \text{Weibull}()$ (do we?),
#' @defreitascostaZeroinflatedcensoredWeibullGamma2021 did exactly this but used
#' the normal Weibull distribution in their likelihood formula.
#'  
#' $$\begin{aligned}
#' \begin{aligned}
#' \log
#' \{L( \pmb {\theta };D_\text {obs})\} &= 
#' \sum_{\begin{array}{c} i:t_i = 0 \end{array}}\log (p_{0i}) + 
#' \sum_{\begin{array}{c} i:t_i > 0 \end{array}}\delta_i
#' \Bigg [
#' \log \bigg (
#' \frac{\alpha _{\mathrm{w}}}{\lambda _{\mathrm{w}}}
#' \bigg ) + (\alpha _{\mathrm{w}}-1)
#' \log \bigg (
#' \frac{t_i}{\lambda _{\mathrm{w}}}
#' \bigg )\Bigg ] \\& 
#' \quad + 
#' \sum_{\begin{array}{c} i:t_i>0 \end{array}}
#' \log (1-p_{0i}) - 
#' \sum_{\begin{array}{c} i:t_i>0 \end{array}}
#' \left(\frac{t_i}{\lambda _{\mathrm{w}}}\right)^{\alpha _{\mathrm{w}}}.
#' \end{aligned} 
#' \end{aligned}$$
#' 
#' > In summary, in the case of observed AFS and AFS is equal or greater than
#' > AAM, this approach can be used. 
#' 
#' 
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
        psu = v021, strata = v023, weights = v005, dob = v011, doi = v008, lab,
        afs = v531, 
        marriage_age = v511,
        n_union = v503,
        month_1st_union = v507,
        year_1st_union = v508,
        cmc_uninon = v509,
        marital_status = v501,
        union_married = v502,
        ever_in = v535
    ) %>%
    mutate(
        marriage_age =  if_else(is.na(marriage_age), 0L, marriage_age),
        marriage_age2 = (cmc_uninon - dob) / 12,
        iso = substr(lab, 1, 2),
        yob = cmc_to_year(dob), svy = cmc_to_year(doi),
        age = svy - yob, sex = datasets$sex[x]
    )

#' # Descriptive statistics
#' 
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
#'
#' ## Key variables
#'
#' All married individuals responded to number of union question. Among those
#' with more than once unions, more than 19% had more than one union in any of
#' the marital statuses. The following differences are needed to be considered
#' when modelling:
#'
#' - The number of currently married with more than one unions reflects
#'   *partially* the number of remarried. Because remarried can also be
#'   currently in other states, including living with partner, widowed,
#'   divorced, and separated. Using this as an estimate of remarried rate would
#'   underestimate it.
#' - Among more than once union respondents, how was the 1st union ended is
#'   unknown in any of the states.

#+ echo=F, results="asis"
dta %>% tabyl(marital_status) |> adorn_pct_formatting() |> kable(caption = "Marital distribution")

dta %>% tabyl(n_union) |> adorn_pct_formatting() |> kable(caption = "Number of union")

dta %>%
    tabyl(marital_status, n_union) %>%
    adorn_percentages() %>%
    adorn_pct_formatting() %>%
    adorn_ns("front") %>%
    kable(caption = "Percent rowwise")

dta %>%
    tabyl(marital_status, n_union) %>%
    adorn_percentages("col") %>%
    adorn_pct_formatting() %>%
    adorn_ns("front") %>%
    kable(caption = "Percent colwise")

#' ## Union - time since sexual debut to first union
#'
#' ### Delay sexual intercourse in child marriage
#' 
dta %<>%
    mutate(
        time_since_debut_c = marriage_age2 - afs,
        marriage_age_d = findInterval2(marriage_age2, seq(15, 30, 5)),
        afs_d = findInterval2(afs, seq(15, 30, 5)),
        time_since_debut_d = findInterval2(time_since_debut_c, c(-2, -1, 0, 1, 2, 5, 10)), 
        time_since_debut_d = fct_recode(time_since_debut_d, '<-3'='-12--3')
    )
#' Very young age of marriage has larger proportion of having *sexual debut
#' after marriage* (\@ref(fig:time_since_db_fig)), up to 50% in those married
#' under 14 (\@ref(tab:time_since_db_tb)). This period is mostly in a year or
#' two; three or more years delay of sexual intercourse after marriage is
#' visible only in those married at age under 14. We might assume the delay is
#' correctly reported.
#'
#' This implies risk of STDs transmission in very young age should not be
#' inferred based on married status but AFS. Also, this group is already in a
#' stable relationship. These suggest excluding this group from the analyses of
#' union formation is reasonable and standard survival model can be used.
#'
#+ time_since_db_fig, fig.cap="Time since debut to married by age of marriage"
dta %>%
    group_by(time_since_debut_d, marriage_age_d) %>%
    count() %>%
    filter(marriage_age_d != "NA-NA") %>%
    ggplot() +
    geom_col(aes(marriage_age_d, n, fill = time_since_debut_d), position = position_fill()) +
    theme(axis.text.x = element_text(angle = 30)) +
    scale_fill_manual(values = ktools::gen_colors(okabe, 8)) +
    labs(x = "Age at marriage", title = "Time since AFS to first union")

#+ fig.cap="Time since debut to married by age of marriage - filtered"
dta %>%
    filter(time_since_debut_c >= 0) %>%
    group_by(time_since_debut_d, marriage_age_d) %>%
    count() %>%
    filter(marriage_age_d != "NA-NA") %>%
    ggplot() +
    geom_col(aes(marriage_age_d, n, fill = time_since_debut_d), position = position_fill()) +
    theme(axis.text.x = element_text(angle = 30)) +
    scale_fill_manual(values = ktools::gen_colors(okabe, 8)) +
    labs(x = "Age at marriage", title="Time since AFS to first union - exclude delay sex group")

#' ### Late sexual debut - short time to marriage
#' 
#' Plot time since debut by AFS shows late sexual debut shorten the time to
#' marriage. Might need to include this in the model as covariate.
#'
#+ fig.cap="Time since debut to married by age at first sex"
dta %>%
    filter(time_since_debut_c >= 0) %>%
    filter(afs != 0) %>%
    group_by(time_since_debut_d, afs_d) %>%
    count() %>%
    filter(time_since_debut_d != "NA-NA") %>%
    ggplot() +
    geom_col(aes(afs_d, n, fill = time_since_debut_d), position = position_fill()) +
    theme(axis.text.x = element_text(angle = 30)) +
    scale_fill_manual(values = ktools::gen_colors(okabe, 9)) +
    labs(x = "Age at first sex", title = "Time since AFS to first union")

dta %<>%
    rownames_to_column("id") %>%
    mutate(afs = as.numeric(haven::zap_labels(afs)), marriage_age = as.numeric(marriage_age)) %>%
    mutate(
        t_virgin = if_else(afs == 0, age, afs),
        t_debut = if_else(afs != 0, afs, NA_real_),
        t_union = if_else(marriage_age != 0, marriage_age, NA_real_),
        t_union = if_else(marital_status == 2, age, t_union), 
        t_widowed = if_else(marital_status == 3, age, NA_real_), 
        t_separate = if_else(marital_status %in% 4:5, age, NA_real_),
        censored = marital_status %in% 2:5
    )

#' 
#' ## Dissolution - time since married
#' 
dta %<>%
    mutate(
        time_since_married_c = (doi - cmc_uninon) / 12,
        time_since_married_d = findInterval2(time_since_married_c, seq(0, 20, 5))
    )
#' ### Once union - competing risk model with right censored data
#'
#' **Among those with only once union**, a majority of the population have
#' remained in the same marital status for more than 5 years. ~~Nearly 95% of
#' the widowed respondents has not remarried for 5 to 20+ year after the first
#' union~~. But since we do not know the time of the event, e.g., time of
#' partner's death, the above interpretation is only correct for the married and
#' partnered group. The others groups, widowed, divorced, and separated, the
#' numbers would more likely reflect the varying in the time of the events.
#'
#' Limited to this set of data, a survival model using the time since married
#' would treat separated, divorce, and widowed as competing events while those
#' in married or partnered states will be right-censored.

#+ fig.cap="Time since married"
dta %>%
    group_by(marital_status, time_since_married_d, n_union) %>%
    count() %>%
    drop_na() %>%
    ggplot() +
    geom_col(aes(marital_status, n, fill = time_since_married_d), position = position_fill()) +
    facet_wrap(~n_union) +
    theme(axis.text.x = element_text(angle = 30))

#' All those in union have records of time of union
#+ union_by_marriage_age, results = "asis"
dta %>%
    filter(n_union == "once") %>%
    tabyl(marriage_age_d, marital_status) %>%
    kable(caption = "Time since union by marital statue")

#' 
#' ## More than once union - left-censoring of both event and time
#' 
#' **Among those with more than once unions**, it is unknown what states a
#' respondent has passed through before reaching the current state. Inexact
#' number of unions also does not allow to narrow down the possible pathways. If
#' we are interested in hazard of union dissolution from the first union, where
#' union dissolution is defined as either widowed, divorce, or separate, we can
#' treat all the states as left-censored at the current time since married.
#'
#' However, if we are interested in in specific events of widowed or divorce,
#' how to specify which state a respondent is censored? Could a probability of
#' the first union ended in widowed or divorce be used as a weight in specifying
#' one of the censored events (likelihood contribution). Event if this works,
#' there are still cases where both states widowed and divorce has been passed
#' through. So we could only limit the analyses at best to hazard of an event
#' *after the first union*.
#' 
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
    tibble(start = char(birth), end = char(aai), a = c(0), z = c(3), risk = c(1)),
    .id = 'set'
) %>%
    mutate(risk = factor(risk, labels = char(No, Yes)), 
    across(where(is.character), str_to_upper)) %>%
    ggplot() +
        geom_segment(aes(a, 1, xend = z, yend = 1, color = risk), arrow = arrow(type = "closed")) +
        geom_text(aes(a, 1.1, label = start)) +
        geom_text(aes(z, 1.1, label = end)) +
        theme_void() +
        scale_color_manual(values = c("grey", "#FD7C02")) +
        coord_cartesian(ylim = c(.5, 1.5)) +
        facet_grid(vars(set)) +
        guides(color = "none")

#' The first case is the typical states transition with observed marriage event;
#' the second case is right censored at the time of interview. The third and
#' fourth case are special in that the time "at risk of getting marriage" are
#' counted from birth. The difference between this exposure time and the two
#' standard cases (1 and 2) lies in the AFS and respondent's age at the time.
#' There are also many cases in scenario 1 which has AFS equals to AAM.
#'
#' Using cases 1 and 2 data, the modeled's output is marriage rate after sexual
#' debut -- this can be used in EPPASM for example. But cases 3, 4, 5 can not be
#' used as they either do not have information.
#'
#' To use all data, we can formulate the AFS as time varying covariate -- need
#' to use time splitting format.
#'
#' ## AAM - since sexual debut
#'
#' \cref{tab:aam_equal_afs} shows among sexually debuted respondents, nearly 40%
#' had also married at the same age. This proportion needs to be modelled
#' similar to @defreitascostaZeroinflatedcensoredWeibullGamma2021.
dta |>
    filter(afs != 0, marriage_age >= afs) %>%
    mutate(aam_equal_afs = marriage_age == afs) %>%
    tabyl(aam_equal_afs) %>%
    knitr::kable(caption = "\\label{tab:aam_equal_afs}Proportion of sexually debuted sample with AFS equals AAM.")

#' However, if computing AFS as $\text{date of birth} + \text{AFS} \times 12$,
#' proportion reduces drastically, from 40% to 1% and 4% in female and male respectively.
#'
here("data", "marriage_epis.rds") %>%
    readRDS() %>%
    group_by(sex) %>%
    summarise(n = sum(duration == 0), pc = n / n()) %>%
    knitr::kable(caption = "\\label{tab:aam_afs_reduce}Percent with AFS equals AAM when doing AFS by CMC")

#' Compare with "cure model", mixture model with a point mass at zero is not as
#' needed: the cure model estimates unknown cure rate from censored observations
#' which change the hazard estimates, i.e.,
#'
#' $$h(t) = \dfrac{(1 - \pi)f(t)}{\pi + (1 - \pi) S(t)}$$
#'
#' with $\pi$ the cure rate. The appearance of $\pi$ in the denominator is
#' because it represents the theoretical probability of an event happened at
#' infinite time. In zero-mass model, however, there a no hidden process as the
#' zeroed observations have been observed. In addition, the zero-mass proportion
#' does not contribute to the hazard.
#'
#' $$h(t) = \dfrac{(1 - \pi_0)f(t)}{(1 - \pi_0) S(t)}, \; t > 0$$
#'
#' where $\pi_0$ the proportion that the event occurred at time zero, thus no
#' contributing towards survival probability in the denominator.
#'
#' In addition, from sexual debut's perspective, individuals who have AAM the
#' same AFS and who have AAM after AFS could be considered as two independent
#' processes. The former is less proactive in making sexual debut than the
#' latter.
#'
#' We can fit independent likelihood for the whole data as
#'
#' $$L(\theta, D) = \prod_{t_i = 0} \pi_0 \prod_{t_i \neq 0} (1 - \pi_0)
#' \prod_{t_i \neq 0} f(t_i)^{\delta_i} S(t_i)^{1 - \delta_i}$$
#'
#' with $Z = 1$ if $t_i = 0$ and $Z =0$ otherwise.
#'
#' $$Z \sim Binomial(n, \pi_0)$$
#' 
#' $$logit(\pi_0) =  \beta_0 + fn(age) $$
#'
#' $fn(age)$ here is a smoothed function of age, the model is fitted separately
#' for each sex.
#'
#' Consider the AAM sample that are observed from the AFS to the occurrence
#' marriage. Observation can be right censored if the last follow-up time is
#' before the occurrence of the event. By partitioning the time axis into small
#' intervals, the number of events in each interval is approximately Poisson
#' with mean $\mu(t) = h(t) P(t)$, where $P(t)$ is the total person-time in the
#' interval. The time at risk of each subject can thus be split into $n$ bins or
#' time intervals common to all subjects (using `ktools::surv_split` to keep
#' zero time subjects). The bin of six months was chosen so that hazard can be
#' considered approximately constant within each interval.
#'
#' Let $t$ denotes the vector of the bins (`episode` variable), $t_i$
#' representing the time for bin $i = 1,..., n$. The vectors $y(t)$ and $P(t)$
#' represent the total number of events observed and the total person-year in
#' each interval. Using the Poisson likelihood to approximate the general
#' likelihood for survival data, the hazard can be estimated by modeling the
#' expected number of events, $\mu(t)$, in each interval as a Poisson variable
#' by using $P(t)$ as an offset term:
#'
#' $$\log[\mu(t)] =  fn(t) + fn(age) + \log[P(t)]$$,
#'
#' where $fn(t)$ denotes the logarithm of the hazard and $fn(age)$ denotes the
#' smoothed effect of age (cubic splines). Similarly, this model is also fitted
#' separately for each sex.
here('data', 'marriage_epis_durspl.rds') %>%
    readRDS() %>%
    mutate(
        py =  if_else(t_end == t_start, log(1/365), log((t_end - t_start)/12)),
        episode = as.numeric(episode),
        age_start = age_start / 12, 
        zero = as.integer(t_end == t_start)
    ) -> sdata

#+ fit_zero, cache=TRUE
zero <- gam(
    zero ~ s(age_start) + s(age_start, by = sex),
    family = binomial,
    data = sdata %>% mutate(sex = if_else(sex=='male', 0, 1))
)

newd <- crossing(sex = 0:1, episode = 1:10, age_start = 10:40)
zero_d <- tibble(newd,
    zero = predict(zero, newd, type = "response"), 
    hz = zero / (1 - zero)
)

#+ fit_gam, cache = TRUE
tst <- gam(
    event ~ 1 + s(episode, bs = 'cr') + s(age_start, bs = 'cr') + 
                s(episode, by = sex, bs = 'cr') + s(age_start, by = sex, bs = 'cr'),
    family = "poisson",
    data = sdata %>% mutate(sex = if_else(sex == 'male', 0, 1)), 
    subset = t_end != t_start,
    offset = py
)

pst <- bind_cols(hz = predict(tst, newd, type = "response"), newd)

#' \cref{fig:Rate_of_marriage} (top rows) shows that sexual debuted when
#' marriage was more prevalence in female than male, happened in younger age in
#' female than male (peaked at 15 vs 20 for female and male, respectively).
#' 
#' Rate of marriage estimates by age shows that, given an individual had ever
#' had sex, rate of marriage is highest in age 20 and 25 in female and male,
#' respectively. 
#' 
#' Rate of marriage is the highest within one year of sexual debut and there
#' were little differences when year of sexual debut was 5+ years, stabilizing
#' at approximately 0.2/year in both sex.
#' 
#' For EPPASM: use the probability estimate to move immediately the newly sexual
#' debuts and the rate to move previously sexually debutted to married compartment.
#' 
#+ Rate_of_marriage, fig.cap = 'Rate of marriage and probability of AFS == AAM'
zero_d %>%
    mutate(episode = 0, type = "Probability of AAM = AFS") %>%
    bind_rows(pst %>% mutate(type = "Rate of marriage by year after sexual debut")) %>%
    mutate(sex = if_else(sex == 0, "Male", "Female")) %>%
    ggplot(aes(age_start, hz, color = factor(episode))) +
    facet_wrap(~ type + sex, scales = "free") +
    geom_line() +
    labs(
        x = "Age", y = "Marriage rate", color = "Years after sexual debut",
        title = "Marriage rate among sexually debutted"
    ) +
    theme(legend.position = 'bottom') + guides(color = guide_legend(nrow=1))

#' 
#' # Methods
#'
#' ## Estimate ever divorce by time since married weighted
#'
#' @clarkDivorceSubSaharanAfrica2015 estimated the probability of the first
#' union ending in divorce or widowed using those numbers among *respondents with
#' only one union*.
#'
#' $$p_{\text{divorce}, t} = \frac{\text{currently divorce}}{\text{currently divorce +
#'   currently widowed}} $$
#' 
#' The proportion ever divorce at time $T$ after married was estimated as
#'
#' $$p_{\text{ever divorce, T}} = p_{\text{currently divorce}, T} +
#' p_{\text{remarried}, T} \times
#' \sum_{t}^T p_{\text{divorce}, t} p_{\text{union dissolution}, t}$$
#'
#' where
#'
#' $$\sum_{t}^T p_{\text{divorce}, t} p_{\text{union dissolution}, t}$$
#'
#' is the cumulative probability of having a prior divorce up to time $T$, with
#'
#' $$p_{\text{union dissolution}, t} = p_{\text{ever divorce},t} - p_{\text{ever
#' union dissolution}, t-1}$$
#'
#' The probability of widowhood was estimated as
#'
#' $$p_{\text{widowhood}} = p_{\text{union dissolution}} - p_{\text{divorce}}$$
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
