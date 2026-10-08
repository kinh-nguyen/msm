# Summarise every fit/fit8_*.rds: convergence, gradient, Hessian, standard
# errors and the row sums of the simulated annual transition matrices.
# Usage: Rscript code/check_fits.r            (writes fit/fit8_summary.csv)
files <- Sys.glob(here::here("fit", "fit8_*.rds"))
rows <- lapply(files, function(f) {
  x <- readRDS(f)
  se <- summary(x$sdr, "fixed")[, "Std. Error"]
  PP <- x$rp$PP
  rs <- sapply(2:dim(PP)[3], function(a) max(abs(rowSums(PP[, , a]) - 1)))
  data.frame(
    file          = basename(f),
    convergence   = x$fit$convergence,
    message       = x$fit$message,
    objective     = x$fit$objective,
    iterations    = x$fit$iterations,
    max_abs_grad  = max(abs(x$sdr$gradient.fixed)),
    pd_hessian    = x$sdr$pdHess,
    n_par         = length(se),
    n_se_missing  = sum(!is.finite(se)),
    max_se        = suppressWarnings(max(se[is.finite(se)])),
    max_rowsum_err = max(rs)
  )
})
out <- do.call(rbind, rows)
print(out, row.names = FALSE)
write.csv(out, here::here("fit", "fit8_summary.csv"), row.names = FALSE)
