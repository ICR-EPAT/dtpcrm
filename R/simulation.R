
#' @export
applied_crm_sim <- function(true_tox, prior, target,
                            max_sample_size, first_dose,
                            num_sims, cohort_size = 1,
                            dose_func = applied_crm,
                            ...) {

  iterations <- list()
  for(i in 1:num_sims) {
    # Start afresh.
    # TODO insert previously observed outcomes here for
    # simulations of partially-observed trials
    tox <- c()
    level <- c()
    dose <- first_dose
    stop <- FALSE
    stop_reason <- NULL
    while(!stop & length(tox) < max_sample_size) {
      # Simulate outcomes for a cohort
      cohort_tox = stats::rbinom(n = cohort_size, size = 1, prob = true_tox[dose])
      cohort_level = rep(dose, cohort_size)
      # Accumulate data
      tox <- c(tox, cohort_tox)
      level <- c(level, cohort_level)
      # Update the model
      x <- dose_func(prior = prior, target = target, tox = tox, level = level,
                     ...)
      dose <- x$mtd
      stop <- ifelse(is.null(x$stop), FALSE, x$stop)
      stop_reason <- x$stop_reason
    }

    print(i)
    iterations[[i]] <- list(tox = tox, level = level, mtd = dose,
                            stop = stop, stop_reason = stop_reason)
    # TODO: further reporting
  }

  # Summarise
  dose_selections = sapply(iterations, function(x) x$mtd)
  doses_given = unlist(sapply(iterations, function(x) x$level))
  summary = list(
    # Echo what you were given
    true_tox = true_tox, prior = prior, target = target,
    max_sample_size = max_sample_size, first_dose = first_dose,
    num_sims = num_sims, cohort_size = cohort_size,
    # Summarise trial outcomes
    prob_stop = mean(sapply(iterations, function(x) x$stop)),
    mtd = sapply(1:length(prior), function(d)
      sum(dose_selections == d, na.rm = TRUE) / num_sims),
    # Summarise doses givne to patients
    doses_given = sapply(1:length(prior), function(d)
      sum(doses_given == d, na.rm = TRUE) / num_sims),
    prob_dose_given = sapply(1:length(prior), function(d)
      sum(doses_given == d, na.rm = TRUE) / length(doses_given))
    # TODO: further reporting
  )


  return(list(summary = summary, iterations = iterations))
}

