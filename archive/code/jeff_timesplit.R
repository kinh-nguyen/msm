library(rdhs)
library(survival)
library(dplyr)
library(haven)
library(tidyr)
library(ktools)

#' Malawi 2015 DHS

ir <- readRDS(downloads[[1]])
mr <- readRDS(downloads[[2]])

df <- ir %>%
  transmute(
    caseid, 
    psu = v021,
    strata = v023,
    weights = v005,
    dob_cmc = v011,
    doi_cmc = v008,
    afs = v525,
    afs_recode = v531,
    afs_flag = as_factor(v532),
    first_birth_cmc = v211,
    marital_status = as_factor(v501),
    marriage_age = v511,
    first_union_cmc = v509,
    n_union = v503, 
    sex = 'female'
  )

df <- mr %>%
  transmute(
    caseid = mcaseid,
    psu = mv021,
    strata = mv023,
    weights = mv005,
    dob_cmc = mv011,
    doi_cmc = mv008,
    afs = mv525,
    afs_recode = mv531,
    afs_flag = as_factor(mv532),
    first_birth_cmc = NA,
    marital_status = as_factor(mv501),
    marriage_age = mv511,
    first_union_cmc = mv509,
    n_union = mv503,
    sex = "male"
  ) %>%
  bind_rows(df)

#' Minimal dataset; no faffing with weights, etc.
df <- df %>%
  mutate(
    eversex = afs > 0,
    afs = if_else(afs == 0, NA_real_, as.double(afs)),
    first_sex_cmc = if_else(eversex, dob_cmc + afs * 12, NA_real_)
  ) %>%
  select(caseid, sex, dob_cmc, doi_cmc, eversex, first_sex_cmc, first_union_cmc)

#' Define start and end of calendar periods, duration periods and age periods
marriage_epis <- df %>%
  filter(eversex) %>%
  mutate(
    time_start = first_sex_cmc,
    time_end = pmin(first_union_cmc, doi_cmc, na.rm = TRUE),
    event = as.integer(!is.na(first_union_cmc) & first_union_cmc == time_end),
    duration = time_end - time_start,
  )

#' survSplit() can't handle cases where tstart == tstop (first sex coincides with first union)
#' ~~Drop~~ keep for zero-inflated.

# remove marriage before AFS
marriage_epis <- filter(marriage_epis, time_end >= time_start)


#' Split episodes by single_year
cuts <- seq(6, max(marriage_epis$duration), 6)

marriage_epis_durspl <- ktools::surv_split(
  marriage_epis, duration = "duration", event = "event", cuts = cuts
)

#' Need to update start/stop of time and age episodes
marriage_epis_durspl %<>%
  mutate(
    duration =  t_end - t_start,
    time_start = first_sex_cmc + t_start,
    time_end = time_start + duration,
    age_start = time_start - dob_cmc,
    age_end = age_start + duration
  )

saveRDS(marriage_epis, "data/marriage_epis.rds")
saveRDS(marriage_epis_durspl, "data/marriage_epis_durspl.rds")
