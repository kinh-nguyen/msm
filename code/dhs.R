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
#' # Data cleaning
#'
#' Remove those
#'
#' - afs or marriage age greater than age
#' - married but no afs

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

#' # Descriptive statistics
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

#' ### Surviving "single" curve by age of respondent
#'
#' Excluding delaying sexual debut after married respondents, standard survival
#' methods can be applied.
#'
#' This is not good. Both the time points are mostly in the past, current age
#' should not have any role. Replacing by age of first sex would be more
#' appropriate.
#'
#+ cache=TRUE
library(survival)
library(ggfortify) # for autoplot

svvdt <- dta %>%
    filter(afs > 0) %>%
    mutate(
        time_since_debut_c = if_else(marriage_age == 0, age - afs, time_since_debut_c),
        event = as.numeric(marital_status != "never in union"), 
        agegr = findInterval2(age, seq(15, 50, 5))
    ) %>%
    filter(time_since_debut_c >= 0) %>%
    select(time_since_debut_c, agegr, event)

svvm1_fit <- survfit(Surv(time_since_debut_c, event) ~ agegr, data = svvdt)
#+ surviving_single_age, fig.cap = "Estimate survival curves by age groups"
autoplot(svvm1_fit) + labs(title = "Survival single by age", x = "Year since AFS") 

#' 
#' ### Surviving "single" curve by age at first sex
#' 
svvdt <- dta %>%
    filter(afs > 0) %>%
    mutate(
        time_since_debut_c = if_else(marriage_age == 0, age - afs, time_since_debut_c),
        event = as.numeric(marital_status != "never in union"), 
        agegr = findInterval2(afs, seq(15, 30, 5))
    ) %>%
    filter(time_since_debut_c >= 0) %>%
    select(time_since_debut_c, agegr, event) %>%
    mutate(agegr = fct_recode(agegr, '25-29' = '30-42'))

library(survival)
svvm1_fit <- survfit(Surv(time_since_debut_c, event) ~ agegr, data = svvdt)
#+ surviving_single_afs, fig.cap = "Estimate survival curves by age at first sex"
library(ggfortify) # for autoplot
autoplot(svvm1_fit) + labs(title = "Survival single by AFS", x = "Year since AFS") 
    mutate(
        state = factor(state, levels = char(virgin, debut, union, separate, widowed)),
        state_numeric = as.numeric(state)
    ) %>% drop_na()
    
msdta

saveRDS(msdta, "data/mw2015.rds")
