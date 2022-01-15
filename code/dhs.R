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
#' ### Competing risk model for one union data subset
#'
#' Package `cmprsk` do the @fineProportionalHazardsModel1999 model.
#' 
#' Lines in the figure is
#'
#' > exp(-B(t)), where B(t) is the estimated cumulative sub-distribution
#' hazard obtained for the specified covariate values, obtained from the
#' Breslow-type estimate of the underlying hazard and the estimated regression
#' coefficients
#' 
competing_risk_data <- dta %>%
    filter(n_union == "once") %>%
    select(age, time = time_since_married_c, marital_status) %>%
    drop_na() %>%
    mutate(
        agegr = findInterval2(age, seq(15, 50, 5)),
        event = case_when(
            marital_status == "married" ~ 0, # censored
            marital_status == "union" ~ 0, # censored
            marital_status == "divorce" ~ 1,
            marital_status == "separated" ~ 1,
            marital_status == "widowed" ~ 2
        )
    )

competing_risk_data %<>%
    mutate(agegr = fct_recode(agegr, "45-49" = "50-51"), agegr_c = as.numeric(agegr))
agegr_lab <- competing_risk_data  %>% tabyl(agegr) %>% pull(agegr)


competing_risk_model_divorce <- cmprsk::crr(
    ftime = competing_risk_data$time,
    fstatus = competing_risk_data$event,
    cov1 = competing_risk_data$agegr_c, failcode = 1
)
summary(competing_risk_model_divorce)

competing_risk_model_widowed <- cmprsk::crr(
    ftime = competing_risk_data$time, 
    fstatus = competing_risk_data$event,
    cov1 = competing_risk_data$agegr_c, failcode = 2
)
summary(competing_risk_model_widowed)

predict_agr <- matrix(1:7, nrow=7)
competing_risk_predict_divorce <- cmprsk::predict.crr(competing_risk_model_divorce, predict_agr)
competing_risk_predict_widowed <- cmprsk::predict.crr(competing_risk_model_widowed, predict_agr)

#+ fig.cap = "Competing risk model with `cmprsk` package on one union data subset"
as_tibble(unclass(competing_risk_predict_divorce)) %>%
    bind_rows(as_tibble(unclass(competing_risk_predict_widowed)), .id = "risk") %>%
    mutate(risk = if_else(risk == "1", "Divorce", "Widowed")) %>%
    pivot_longer(3:9, names_prefix = "V") %>%
    mutate(agegr = agegr_lab[as.numeric(name) - 1]) %>%
    ggplot(aes(V1, 1 - value, color = agegr)) +
    geom_line() +
    facet_wrap(~risk, scales = "fixed") +
    labs(
        title = "Competing risks model on one union data subset",
        x = "Time since married", y = "Hazard"
    )

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
#' ## Descriptive statistics
#' 
#+ results="asis"
dta %>% tabyl(afs) |> adorn_pct_formatting() |> kable(caption = "AFS distribution")

dta %>% tabyl(marriage_age) |> adorn_pct_formatting() |> kable(caption = "Marriaged age distribution")

#'
#' ## Union - time since sexual debut to first union
#'
#+ time_since_db_tb, results="asis"
dta %>%
    filter(marriage_age_d != "NA-NA") %>%
    tabyl(time_since_debut_d, marriage_age_d) %>%
    adorn_totals() %>%
    adorn_totals("col") %>%
    adorn_percentages("col") %>%
    adorn_pct_formatting() %>%
    adorn_ns("front") %>%
    kable(caption = "Time since debut to first union")

#+ results="asis"
dta %>%
    filter(time_since_debut_c >= 0) %>%
    filter(afs != 0) %>%
    tabyl(time_since_debut_d, afs_d) %>%
    adorn_totals() %>%
    adorn_totals("col") %>%
    adorn_percentages("col") %>%
    adorn_pct_formatting() %>%
    adorn_ns("front") %>%
    kable(caption = "Time since debut to first union by AFS")

#' 
#' ## Dissolution - time since married
#' 
# Three way cross-table
three_ways <- dta %>%
    tabyl(marital_status, time_since_married_d, n_union) %>%
    adorn_totals() %>%
    adorn_totals("col") %>%
    adorn_percentages() %>%
    adorn_pct_formatting() %>%
    adorn_ns("front")

#+ results="asis"
three_ways$once %>% kable(caption = "Time since married - one union")
three_ways$`more than once` %>% kable(caption = "Time since married - more than once union")