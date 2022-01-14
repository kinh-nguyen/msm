# simulated data to test
set.seed(10)
ftime <- rexp(200)
fstatus <- sample(0:2, 200, replace = TRUE)
cov <- matrix(runif(200), nrow = 200)
dimnames(cov)[[2]] <- c("x")

library("cmprsk")

z <- crr(ftime, fstatus, cov, failcode = 1, cencode = 0)
summary(z)
plot(z)
z.p <- predict(z, matrix(runif(1), nrow = 1))
plot(z.p, lty = 1, color = 2:3)
