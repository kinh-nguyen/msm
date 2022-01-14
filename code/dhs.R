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
#' linkcolor: blue
#' urlcolor: blue
#' citecolor: red
#' ---
#'

library(rdhs)
library(tidyverse)
library(ktools)
set_rdhs_config(email = "ath19@ic.ac.uk", project = "Statistics and Machine Learning for HIV")
#' get one sex one country first
datasets <- dhs_datasets("MW", fileFormat = "flat", fileType = c("IR")) %>%
    filter(SurveyYear == 2015)  %>%
    mutate(sex = if_else(FileType == "Individual Recode", "female", "male"))
datasets
downloads <- get_datasets(datasets$FileName)
names(downloads)

#' Extract 
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

#' cleaning
dta %<>% filter(!(afs > age | marriage_age > age))

dta$marriage_age %>% table(useNA = "a")
dta$afs %>% table(useNA = "a")

#' AFS and AAm
dta %>%
    filter(afs %in% 10:55 & marriage_age != 0) %>%
    summarise(mean(afs > marriage_age))

dta %>%
    filter(afs %in% 10:55 & marriage_age != 0) %>%
    filter(afs > marriage_age) %>%
    summarise(mean(afs - marriage_age))

#' married but no afs
dta %<>% filter(!(marital_status != 0 & afs == 0))

#' cleaning
#' todo let assume aam >= afs for now
# dta %<>% filter(afs <= marriage_age)

#' states indicators
dta %>%
    pull(marital_status) %>%
    table(useNA = "a")

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
dta

library("flexsurv")
n_trans <- max(tmat, na.rm = TRUE)
#' plot raw data
library(ggpubr, help, pos = 2, lib.loc = NULL)

four_state_raw <- ggarrange(
    dta %>% ggplot() +
        geom_density(aes(t_debut)) +
        labs(title = "t_debut", x = "Time"),
    dta %>%
        mutate(diff = if_else(censored, t_union, t_union - t_debut)) %>%
        filter(diff >= 0) %>% 
        ggplot() +
        geom_density(aes(diff, fill = factor(censored)), alpha=.7) +
        labs(title = "Debut before Married", fill = "Censored", x = "Time"),
    dta %>%
        mutate(diff = if_else(censored, t_union, t_union - t_debut)) %>%
        filter(diff < 0) %>% 
        ggplot() +
        geom_density(aes(diff, fill = factor(censored)), alpha=.7) +
        labs(title = "Debut after Married", fill = "Censored", x = "Time"),
    
    dta %>%
        mutate(
            diff = if_else(marriage_age != 0, t_separate - t_union, t_separate),
            know_aam = marriage_age != 0
        ) %>%
        ggplot() +
            geom_density(aes(diff, fill = factor(know_aam)), alpha = .7) +
            labs(title = "Separated", fill = "Known age at married?", x = "Time"),
    
    dta %>%
        mutate(
            diff = if_else(marriage_age != 0, t_widowed - t_union, t_widowed),
            know_aam = marriage_age != 0
        ) %>%
        ggplot() +
            geom_density(aes(diff, fill = factor(know_aam)), alpha = .7) +
            labs(title = "Widowed", fill = "Known age at married?", x = "Time"),
   
    ncol = 3, nrow = 2
)
ggsave("img/four_state_raw.png")

#' Prepare msm format
msdta <-
dta %>% 
    select(age, marital_status, marriage_age, starts_with("t_")) %>%
    filter(t_union >= t_debut) %>%
    mutate(
        pid = 1:n(),
        t_debut = if_else(t_debut == t_virgin, t_virgin + 1 / 12, t_debut),
        t_union = if_else(t_union <= t_debut, t_debut + 1 / 12, t_union),
        t_separate = if_else(t_separate <= t_union, t_union + 1 / 12, t_separate),
        t_widowed = if_else(t_widowed <= t_union, t_union + 1 / 12, t_widowed),
    ) %>%
    pivot_longer(starts_with("t_"), names_to = "state", values_to = "time", names_prefix = "t_") %>%
    mutate(
        state = factor(state, levels = char(virgin, debut, union, separate, widowed)),
        state_numeric = as.numeric(state)
    ) %>% drop_na()
    
msdta

saveRDS(msdta, "data/mw2015.rds")
