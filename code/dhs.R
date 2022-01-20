#' ---
#' title: "Poisson format - splitting time"
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
library(ggfortify)
library(janitor)
library(survival)
library(Epi)
library(popEpi)
library(mgcv)
tabyl <- function(dat, ...) janitor::tabyl(dat, ..., show_missing_level = FALSE)
devtools::load_all("~/Code/R/ktools/")
library(knitr)
opts_chunk$set(echo = FALSE, cache = FALSE, out.extra = "")
set.seed(1)

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
    tibble(start = char(birth, afs, aam), end = char(afs, aam, aai), a = 0:2, z = 1:3, set = 1, risk = c(0, 1, 0)),
    tibble(start = char(birth, afs), end = char(afs, aai), a = c(0, 1), z = c(1, 3), set = 2, risk = c(0, 1)),
    tibble(start = char(birth, aam, afs), end = char(aam, afs, aai), a = c(0, .7, 1), z = c(.7, 1, 3), set = 3, risk = c(1, 0, 0)),
    tibble(start = char(birth, aam), end = char(aam, aai), a = c(0, .7), z = c(.7, 3), set = 4, risk = c(1, 0)),
    #' no AFS
    tibble(start = char(birth), end = char(aai), a = c(0), z = c(3), set = 5, risk = c(1))
) %>%
    mutate(risk = factor(risk, labels = char(No, Yes)), 
    across(where(is.character), str_to_upper)) %>%
    ggplot() +
        geom_segment(aes(a, 1, xend = z, yend = 1, color=risk), 
        arrow = arrow(type = "closed")) +
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
#' ## Cox model - original data
#' 
#' Age at the start of exposure time is used as covariate, not age at interview.
#' In particular, the data is recoded as
#' 
#' ```r
#' age_0 = case_when(
#'     event == 1 & afs != 0 & afs <= marriage_age ~ afs, # age count from afs
#'     event == 0 & afs != 0 ~ afs, # age count from afs
#'     event == 1 & afs != 0 & afs > marriage_age ~ 0, # age count from birth
#'     event == 1 & afs == 0 ~ 0, # age count from birth
#'     event == 0 & afs == 0 ~ 0 # age count from a
#' )
#' ```
#' 
#' which can be see in \ref{tab:data_cox_classic}.
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
#+ data_cox_classic, results = "asis"
lexis %>%
    head() %>%
    kable(caption = "\\label{tab:data_cox_classic}Data format for Cox model with age at risk recoded")

#+ echo=TRUE
cox_mod_ori <- coxph(Surv(time, event) ~ age_0, lexis)

#' Fitted Cox model show age increases rate by 14%.
ci.exp(cox_mod_ori)
#'
#' ## Cox model - time-split data format
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
lexis %<>% mutate(time = if_else(time == 0, time + 1 / 12, time)) # add a month
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

#' A Cox model for this data is fitted with interval format as
#+ echo=T
cox_mod_int <- coxph(Surv(tar, tar + lex.dur, lex.Xst=="Married") ~ aar, data = sL)

#+ echo=F
ci.exp(cox_mod_int)

#' which gives the same estimate of the age at risk coefficient.
#'
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

#+ empirical_rate, fig.cap = caption_long
psdata %>%
    # grouping for plotting only
    mutate(tar2 = findInterval2(tar, c(0:10, 15, 20))) %>%
    group_by(tar2, aar) %>%
    mutate(count = sum(count), person_year = sum(person_year)) %>%
    ungroup() %>%
    mutate(rate = count / person_year) %>%
    group_by(aar) %>%
    mutate(med = median(rate)) %>%
    ggplot() +
    geom_line(aes(aar, rate, color = factor(tar2))) +
    geom_line(aes(aar, med), size = 2) +
    scale_color_viridis_d() +
    theme(legend.position = 'bottom') +
    labs(color = "Time since at risk", x = "Age at risk")

#' Fitting Poisson model with time at risk (TAR) as a covariate and
#' $\log(\text{person-year}))$ as offset gave a ~~slightly larger effect of age at
#' risk, 16% vs. 14%~~ (due to removal of 30+ time at risk) the same estimate of
#' age at risk effect (\cref{tab:poisson_model}).
#' 
#' Now the model 
#' 
mp1 <- gam(count ~ aar + tar, poisson(link = "log"), psdata, offset = log(person_year))

#+ poisson_model, results = 'asis'
ci.exp(mp1) |> knitr::kable(caption='\\label{tab:poisson_model}Poisson model with fixed TAR effect')

#' ## Smoothed age effect
#' 
#' To reflect the observed nonlinear age effect in \cref{fig:empirical_rate}.
#' Age effect is smoothed with `gam` as follows:
#' 
mp2 <- gam(count ~ s(aar) + s(tar), poisson(link = "log"), psdata, offset = log(person_year))
#' The model with smoothed age effect shows significant improvement, as seen in
#' in \Cref{fig:poisson_rate_separate}. Rate of marriage increases with age and taping
#' off from age 25. Consider adolescent age, marriage are more likely to occur
#' within a year exposed to risk (sexual debut in this case). 

#+ echo=F
anova(mp1, mp2, test = "Chisq")

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
    ) + guides(fill = FALSE)
rom

#+ poisson_rate_separate, fig.cap = 'Estimate rate of marriage by age'
rom + facet_wrap(~tar, scales = 'fixed')

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