# Interval-censored single-visit msm example
# - hidden truth simulated at exact times
# - one observed visit per subject (age at report)
# - msm data in interval-censored format (obstype = 2)

library(dplyr)
library(msm)

set.seed(42)
N0 <- 2500

# True hazard parameters (log-logistic scale parameter = typical time scale)
true_params <- list(
  "1_to_2" = c(shape = 1.5, scale = 4),
  "1_to_3" = c(shape = 1.2, scale = 15),
  "2_to_4" = c(shape = 2.0, scale = 6),
  "3_to_4" = c(shape = 1.8, scale = 11),
  "4_to_2" = c(shape = 1.5, scale = 4),
  "4_to_3" = c(shape = 1.2, scale = 15)
)

# Correct log-logistic random generator (inverse CDF)
rllogis_custom <- function(n, shape, scale) {
  u <- runif(n)
  scale * (u / (1 - u))^(1 / shape)
}

# Baseline sample and ages
base <- data.frame(id = 1:N0)
base$current_age <- runif(N0, 15, 49)            # age at interview/report
base$afm <- rgamma(N0, shape = 50, rate = 2.5)  # age first married (example)
# Keep only those with afm < current_age
base <- base[base$afm < base$current_age, ]
N <- nrow(base)
base$id <- 1:N

# Simulate exact hidden transition histories and store final observed state at report
final_states <- integer(N)
# We'll also keep the exact hidden event list in case you want it
hidden_histories <- vector("list", N)

for (i in seq_len(N)) {
  current_state <- 1L
  time_at_entry_age <- base$afm[i]   # absolute age of last entry
  subj_events <- data.frame(age = time_at_entry_age, state = current_state) # record start
  
  repeat {
    if (current_state == 1L) {
      t_to_2 <- rllogis_custom(1, true_params$"1_to_2"["shape"], true_params$"1_to_2"["scale"])
      t_to_3 <- rllogis_custom(1, true_params$"1_to_3"["shape"], true_params$"1_to_3"["scale"])
      dt <- min(t_to_2, t_to_3); next_state <- c(2L, 3L)[which.min(c(t_to_2, t_to_3))]
    } else if (current_state == 2L) {
      dt <- rllogis_custom(1, true_params$"2_to_4"["shape"], true_params$"2_to_4"["scale"]); next_state <- 4L
    } else if (current_state == 3L) {
      dt <- rllogis_custom(1, true_params$"3_to_4"["shape"], true_params$"3_to_4"["scale"]); next_state <- 4L
    } else { # state 4
      t_to_2 <- rllogis_custom(1, true_params$"4_to_2"["shape"], true_params$"4_to_2"["scale"])
      t_to_3 <- rllogis_custom(1, true_params$"4_to_3"["shape"], true_params$"4_to_3"["scale"])
      dt <- min(t_to_2, t_to_3); next_state <- c(2L, 3L)[which.min(c(t_to_2, t_to_3))]
    }
    age_next <- time_at_entry_age + dt
    # If next event occurs before reported age, record it and continue
    if (age_next < base$current_age[i]) {
      current_state <- next_state
      time_at_entry_age <- age_next
      subj_events <- rbind(subj_events, data.frame(age = age_next, state = current_state))
    } else {
      # Event would occur after report age => final observed state is current_state at report
      final_states[i] <- current_state
      hidden_histories[[i]] <- subj_events
      break
    }
  }
}

# Observed time since process start (time from afm to interview)
obs_time <- base$current_age - base$afm

# Build msm interval-censored dataset:
# For each subject create two rows:
#  - baseline known state at time 0 (exact)
#  - one follow-up observation at time2 = obs_time, with interval (time1 = 0, time2 = obs_time)
#    and reported state = final_states[i]
msm_rows <- vector("list", N)
for (i in seq_len(N)) {
  msm_rows[[i]] <- data.frame(
    subject = base$id[i],
    # First row: exact known start at time 0
    time = 0,
    time1 = 0,
    time2 = 0,
    state = 1L
  )
  # Second row: single visit at obs_time; interval is (0, obs_time]
  msm_rows[[i]] <- rbind(
    msm_rows[[i]],
    data.frame(
      subject = base$id[i],
      time = obs_time[i],    # time variable (msm will use time1/time2 because obstype=2)
      time1 = 0,
      time2 = obs_time[i],
      state = final_states[i]
    )
  )
}
msm_data_interval <- bind_rows(msm_rows) %>% arrange(subject, time)

# Map numeric states to labels if you want readable output (msm uses numeric internally)
state_labels <- c("Married", "Divorced", "Widowed", "Remarried")
msm_data_interval$state_label <- factor(msm_data_interval$state,
                                        levels = 1:4,
                                        labels = state_labels)

# Check the resulting table (one start + one visit per subject)
print(head(msm_data_interval, 6))

# Define Q matrix (numbers indicate shared parameters)
Q_constrained <- rbind(
  c(0, 1, 2, 1),  # From Married -> Divorced (1), Widowed (2)
  c(0, 0, 0, 3),  # From Divorced -> Remarried (3)
  c(0, 0, 0, 4),  # From Widowed -> Remarried (4)
  c(0, 1, 2, 0)   # From Remarried -> Divorced (1), Widowed (2) same as Married
)
dimnames(Q_constrained) <- list(state_labels, state_labels)

# Fit msm with interval censoring (obstype = 2 uses time1 & time2)
# No qinit/init to avoid version conflicts; let msm generate initial values
fit_msm_interval <- msm(
  state ~ time,
  subject = subject,
  data = msm_data_interval,
  qmatrix = Q_constrained,
  obstype = 2    # interval censored using time1 & time2
)

summary(fit_msm_interval)

# msm_data_interval as in your dataset: one row time0 and one time2 per subject
library(dplyr)
obs <- msm_data_interval %>%
  group_by(subject) %>%
  summarize(start = first(state), end = last(state)) %>%
  ungroup()
table(obs$start, obs$end)

# Jacobian & Krylov identifiability check for single-visit interval data
# Assumes you already have `msm_data_interval` with columns:
#   subject, time, time1, time2, state  (as produced earlier)
# and state codes 1..4 where start state = 1 at time1=0 and time2 > 0 for the visit.
#
# Outputs:
#  - m: number free parameters
#  - S: number distinct observation times (excluding time 0)
#  - k: number states
#  - counting bound check: m <= S*(k-1)
#  - Jacobian rank (finite diff) and whether full column rank (=m)
#  - Krylov rank (span dimension of e1' Q^n)

library(expm)    # for matrix exponential
library(Matrix)  # for rank estimation if needed

# --- Model-to-Q function (adjust if you use different parametrization) ---
# Parameter theta = c(q1, q2, q3, q4)
# q1 = M/R -> D, q2 = M/R -> W, q3 = D -> R, q4 = W -> R
Q_from_theta <- function(theta) {
  q1 <- theta[1]; q2 <- theta[2]; q3 <- theta[3]; q4 <- theta[4]
  Q <- matrix(0, 4, 4)
  # from Married (1)
  Q[1,2] <- q1; Q[1,3] <- q2; Q[1,1] <- -(q1 + q2)
  # from Divorced (2)
  Q[2,4] <- q3; Q[2,2] <- -q3
  # from Widowed (3)
  Q[3,4] <- q4; Q[3,3] <- -q4
  # from Remarried (4) shares q1,q2 with Married
  Q[4,2] <- q1; Q[4,3] <- q2; Q[4,4] <- -(q1 + q2)
  return(Q)
}

# --- Prepare observed distinct times t_s (exclude time 0) ---
# user must have msm_data_interval in environment
if (!exists("msm_data_interval")) stop("msm_data_interval not found in workspace.")
times2 <- unique(msm_data_interval$time2)
times2 <- sort(times2[!is.na(times2) & times2 > 0])
S <- length(times2)
k <- 4
m <- 4  # number of free params in this parametrization

cat("k =", k, ", m =", m, ", distinct visit times S =", S, "\n")
cat("Counting bound: m <=", S*(k-1), " -> ", m <= S*(k-1), "\n\n")

# --- Observation map F(theta): stack p_j(t_s) for j=1..k-1 and s=1..S ---
F_of_theta <- function(theta, times = times2) {
  Q <- Q_from_theta(theta)
  out <- numeric(length(times)*(k-1))
  idx <- 1
  for (t in times) {
    P <- expm(Q * t)
    # take first row probabilities p1..pk, but drop last since they sum to 1
    row1 <- as.numeric(P[1, ])
    out[idx:(idx + k - 2)] <- row1[1:(k-1)]
    idx <- idx + (k-1)
  }
  return(out)
}

# --- Numeric Jacobian via central finite differences ---
numeric_jacobian <- function(theta, eps = 1e-6) {
  f0 <- F_of_theta(theta)
  p <- length(theta)
  nout <- length(f0)
  J <- matrix(0, nout, p)
  for (j in seq_len(p)) {
    e <- numeric(p); e[j] <- 1
    th_p <- theta + eps*e
    th_m <- theta - eps*e
    f_p <- F_of_theta(th_p)
    f_m <- F_of_theta(th_m)
    J[, j] <- (f_p - f_m) / (2*eps)
  }
  return(J)
}

# --- Krylov span dimension from e1' Q^n  (n=0..k-1) ---
krylov_rank <- function(Q) {
  rows <- matrix(0, nrow = k, ncol = k)
  v <- matrix(0, nrow = 1, ncol = k)
  v[1] <- 1  # e1'
  for (n in 0:(k-1)) {
    rows[n+1, ] <- as.numeric(v %*% (Q %^% n))  # Q^n
  }
  # compute rank via svd tolerance
  s <- svd(rows)$d
  tol <- max(dim(rows)) * max(s) * .Machine$double.eps
  r <- sum(s > tol)
  return(list(rank = r, singular.values = s))
}

# --- Run checks at a supplied theta (default: equal small rates) ---
theta0 <- c(0.1, 0.1, 0.1, 0.1)
cat("Using theta0 =", paste(round(theta0,4), collapse = ", "), "\n\n")

J <- numeric_jacobian(theta0, eps = 1e-6)
sv <- svd(J)$d
tolJ <- max(dim(J)) * max(sv) * .Machine$double.eps
rankJ <- sum(sv > tolJ)

cat("Jacobian: output dim =", nrow(J), ", params =", ncol(J), "\n")
cat("Singular values (Jacobian):\n"); print(sv)
cat("Jacobian numeric rank =", rankJ, " (full rank required = ", m, ")\n\n")

Q0 <- Q_from_theta(theta0)
Kres <- krylov_rank(Q0)
cat("Krylov singular values:\n"); print(Kres$singular.values)
cat("Krylov rank =", Kres$rank, " (max =", k, ")\n\n")

# --- Summary interpretation ---
cat("Interpretation:\n")
if (rankJ < m) {
  cat("- Jacobian is rank-deficient: local non-identifiability (some parameter directions invisible).\n")
} else {
  cat("- Jacobian is full column-rank: local identifiability NOT ruled out by linearized test.\n")
}
if (Kres$rank < k) {
  cat("- Krylov span < k: the first-row functional information cannot distinguish all directions in Q.\n")
} else {
  cat("- Krylov span = k: first-row moments span full state-space (necessary but not sufficient for identifiability).\n")
}
