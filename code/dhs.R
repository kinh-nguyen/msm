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