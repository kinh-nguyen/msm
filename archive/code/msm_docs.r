library(tidyverse)
library(msm)
cav %>% head()
statetable.msm(state = state, subject = PTNUM, data = cav)

Q <- rbind(
    c(0, 0.25, 0, 0.25),
    c(0.166, 0, 0.166, 0.166),
    c(0, 0.25, 0, 0.25),
    c(0, 0, 0, 0)
)
Q

# get crude Q
Q.crude <- crudeinits.msm(state ~ years, PTNUM, data = cav, qmatrix = Q)
Q.crude

cav.msm <- msm(state ~ years,
    subject = PTNUM, data = cav, qmatrix = Q, deathexact = 4,
    # exacttimes = TRUE, 
    # obstype = 2
)

plot(cav.msm)

cav.msm <- msm(state ~ years,
    subject = PTNUM, data = cav, qmatrix = Q, deathexact = 4,
    covariates = ~ sex
    # exacttimes = TRUE, 
    # obstype = 2
)

plot(cav.msm)

hazard.msm(cav.msm)
pmatrix.msm(cav.msm, t=1)
