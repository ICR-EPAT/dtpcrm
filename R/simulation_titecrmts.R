# obswin - observation period in units of time (days)
# minfu - minimum follow-up required for each patient in cohort to conduct model update
# recrate - # of pts recruited per obswin e.g. 5
# initdes - vector specifying the doses to be assisned to pts as per the initial design e.g. c(1,1,2,2,3,3,4,4,4,4)

#' @export
applied_titecrmts_sim <- function(true_tox, prior, target,
                                max_sample_size,
                                num_sims, cohort_size = 1,
                                obswin, minfu,recrate,
                                initdes,
                                dose_func = applied_titecrm,
                                ...) {

  iterations <- list()
  for(i in 1:num_sims) {
    # Start afresh.
    # TODO insert previously observed outcomes here for
    # simulations of partially-observed trials
    tox <- c()
    level <- c()
    fu <- c()
    stop <- FALSE
    stop_reason <- NULL
    rectime <- obswin / recrate

    # simulate data for initial design
    tox <- stats::rbinom(n = length(initdes), size = 1, prob = true_tox[initdes])
    level <- initdes
    if(any(tox == 1)){
      pos <- cohort_size * ceiling(which(tox == 1)[1] / cohort_size)
      tox <- tox[1:pos]
      level <- level[1:pos]

      fu = (rectime * (length(tox) - 1)) - (rectime * c(0:(length(tox) - 1))) + minfu
      fu[tox == 1] <- obswin
      fu <- pmin(fu, obswin)

      x <- dose_func(prior = prior, target = target, tox = tox, level = level,
                     followup = fu, obswin = obswin, ...)
      dose <- x$mtd
      stop <- ifelse(is.null(x$stop), FALSE, x$stop)
      stop_reason <- x$stop_reason

    }

    if(all(tox == 0)){
      dose = utils::tail(level, 1)
    } else{

    while(!stop & length(tox) < max_sample_size) {


      # Simulate outcomes for a cohort (after initial des)
      cohort_tox = stats::rbinom(n = cohort_size, size = 1, prob = true_tox[dose])
      cohort_level = rep(dose, cohort_size)
      cohort_fu = (rectime * (cohort_size - 1)) - (rectime * c(0:(cohort_size-1))) + minfu # follow-up based on fixed accrual from recrate TODO allow for non fixed accrual
      cohort_fu[cohort_tox == 1] <- obswin # weight of 1 for dlt patients


      # Accumulate data
      tox <- c(tox, cohort_tox)
      level <- c(level, cohort_level)
      fu <- fu + (cohort_size * rectime) + minfu # add on additional follow-up for previous patients
      fu <- c(fu, cohort_fu)
      fu <- pmin(fu, obswin) # fix follow-up to maximum of observational period


      # Update the model
      x <- dose_func(prior = prior, target = target, tox = tox, level = level,
                     followup = fu, obswin = obswin, ...)
      dose <- x$mtd
      stop <- ifelse(is.null(x$stop), FALSE, x$stop)
      stop_reason <- x$stop_reason
    }}


    if(!stop){
      x <- dose_func(prior = prior, target = target, tox = tox, level = level,   # run model final time with full follow-up
                     followup = rep(obswin, length(tox)), obswin = obswin, ...)
      dose <- x$mtd
    }

    print(i)
    iterations[[i]] <- list(tox = tox, level = level, mtd = dose,
                            stop = stop, stop_reason = stop_reason)
    #print(iterations[[i]])
    # TODO: further reporting
  }

  # Summarise
  dose_selections = sapply(iterations, function(x) x$mtd)
  doses_given = unlist(sapply(iterations, function(x) x$level))
  summary = list(
    # Echo what you were given
    true_tox = true_tox, prior = prior, target = target,
    max_sample_size = max_sample_size, initdes = initdes,
    num_sims = num_sims, cohort_size = cohort_size,
    # Summarise trial outcomes
    prob_stop = table(substr(unlist(sapply(iterations, function(x) x$stop_reason)), 1, 15)) / num_sims,
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
