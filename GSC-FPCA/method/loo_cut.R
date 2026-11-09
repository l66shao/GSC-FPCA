# Leave-one-out estimate under the cut posterior. ll is a units-by-draws matrix from the sampler with PSLOO = 1.
# Controls inform the components, so they use PSIS on the marginal log likelihood. Treated units inform only nu, and
# their rows hold the exact leave-one-out predictive log density given the components, so the estimate for a treated
# unit is the log of its average over the draws.
elpd_cut <- function(ll, is_control) {
  lme <- function(v) { m <- max(v); m + log(mean(exp(v - m))) }
  psis <- suppressWarnings(loo::loo(t(ll[is_control, , drop = FALSE])))
  e <- numeric(nrow(ll))
  e[is_control] <- psis$pointwise[, "elpd_loo"]
  if (any(!is_control)) e[!is_control] <- apply(ll[!is_control, , drop = FALSE], 1, lme)
  n <- length(e); et <- e[!is_control]
  list(elpd_i = e, is_control = is_control, elpd = sum(e), se = sqrt(n) * sd(e), looic = -2 * sum(e), looic_se = 2 * sqrt(n) * sd(e),
       elpd_ctl = sum(e[is_control]), se_ctl = sqrt(sum(is_control)) * sd(e[is_control]),
       elpd_trt = sum(et), se_trt = if (length(et) > 1) sqrt(length(et)) * sd(et) else NA_real_,
       p_loo_ctl = psis$estimates["p_loo", 1], p_loo_ctl_se = psis$estimates["p_loo", 2], pareto_k = psis$diagnostics$pareto_k, psis = psis)
}
# Difference in elpd between two fits of the same data: a minus b, with its standard error.
elpd_diff <- function(a, b) { d <- a$elpd_i - b$elpd_i; c(diff = sum(d), se = sqrt(length(d)) * sd(d)) }
