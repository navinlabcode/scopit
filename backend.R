invert.pmultinom.stoppable <- function (lower = -Inf, upper = Inf, probs, target.prob, method, maxval)
{
  if (isTRUE(method != "exact"))
    stop("Method must be specified. Only available method right now is \"exact\"")
  stopifnot(is.numeric(lower) | all(is.na(lower)))
  stopifnot(is.numeric(upper) | all(is.na(upper)))
  stopifnot(is.numeric(probs) | all(is.na(probs)))
  stopifnot(is.numeric(target.prob) | all(is.na(target.prob)))
  stopifnot(all(sapply(0 <= probs & probs <= 1, isTRUE) | is.na(probs)))
  stopifnot(all(sapply(0 <= target.prob & target.prob <= 1,
                       isTRUE) | is.na(target.prob)))
  tryCatch(mapply(function(x, y, z) NULL, lower, upper, probs),
           warning = function(w) stop(sprintf("Error generated from warning in %s: %s",
                                              format(conditionCall(w)), conditionMessage(w))))
  if (all(is.na(target.prob)) | any(is.na(c(lower, upper, probs)))) {
    return(NA)
  }
  if (!isTRUE((all(is.infinite(lower) & lower < 0) | all(is.infinite(upper) &
                                                         upper > 0)))) {
    stop("Lower and upper both given. Inverse is not guaranteed to exist under these conditions, and correct behavior has not yet been implemented.")
  }
  if (isTRUE(all(is.infinite(lower) & lower < 0) & all(is.infinite(upper) &
                                                       upper > 0))) {
    return(0)
  }
  if (all(is.infinite(lower) & lower < 0))
    increasing <- FALSE
  if (all(is.infinite(upper) & upper > 0))
    increasing <- TRUE
  r.lower <- mapply(function(x, y, z) x, lower, upper, probs)
  r.upper <- mapply(function(x, y, z) y, lower, upper, probs)
  r.probs <- mapply(function(x, y, z) z, lower, upper, probs)
  if (!isTRUE(all.equal(sum(r.probs), 1)))
    stop(sprintf("Probabilities do not sum to 1 after recycling: %s",
                 paste(r.probs, collapse = ", ")))
  is.zero.prob <- sapply(r.probs, function(p) isTRUE(all.equal(p,
                                                               0)))
  is.lower.bound <- r.lower >= 0
  causes.infinite.loop <- sapply(is.zero.prob & is.lower.bound,
                                 isTRUE)
  if (any(causes.infinite.loop)) {
    stop(sprintf("Requiring more than %.3f in a category with probability %.3f will cause an infinite loop",
                 r.lower[causes.infinite.loop][1], r.probs[causes.infinite.loop][1]))
  }
  a <- 3
  v <- a
  n <- v^2
  distribution <- pmultinom:::sum.tpois.pmf(r.lower, r.upper, n * r.probs,
                                n)
  path <- pmultinom:::prob.between(r.lower, r.upper, r.probs, 0:n, distribution,
                       n)
  while ((increasing & path[length(path)] < max(target.prob,
                                                na.rm = TRUE)) | !increasing & path[length(path)] > min(target.prob,
                                                                                                        na.rm = TRUE)) {
    if (n > maxval) return(NA)
    old.n <- n
    v <- v + a
    n <- v^2 - 1
    distribution <- pmultinom:::sum.tpois.pmf(r.lower, r.upper, n * r.probs,
                                  n)
    path <- c(path, pmultinom:::prob.between(r.lower, r.upper, r.probs,
                                 (old.n + 1):n, distribution, n))
  }
  index.condition.satisfied <- sapply(target.prob, function(tp) which(increasing &
                                                                        path >= tp | !increasing & path <= tp)[1])
  sample.size.condition.satisfied <- index.condition.satisfied -
    1
  return(sample.size.condition.satisfied)
}
