loglikelihood <- function(p) {
  beta <- p[1:3]
  xi <- p[4:6] # Contribution per unit:
  loglik <- 0
  for(i in 1:N){
    # Data for subject i:
    data.i <- dta.split[[i]]
    O <- data.i$state; t <- data.i$age
    # Loop over observations for subject i:
    for (j in 2:length(O)) {
      # Q and P matrix:
      Q <- matrix(0, 3, 3)
      t1 <- t[j - 1]
      t2 <- t[j]
      Q[1, 2] <- exp(beta[1] + xi[1] * t1)
      Q[1, 3] <- exp(beta[2] + xi[2] * t1)
      Q[1, 1] <- -sum(Q[1, ])
      Q[2, 3] <- exp(beta[3] + xi[3] * t1)
      Q[2, 2] <- -sum(Q[2, ])
      P <- MatrixExp(mat = Q, t = t2 - t1)
      # Likelihood contribution:
      death <- as.numeric(O[j] == D)
      loglik <- loglik + log(
        (1 - death) * P[O[j - 1], O[j]] +
          death * (
            P[O[j - 1], 1] * Q[1, D] +
            P[O[j - 1], 2] * Q[2, D]
          )
      )
    }
  }
  return(-loglik)
}

