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