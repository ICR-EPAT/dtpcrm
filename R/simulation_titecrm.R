# obswin - observation period in units of time (days)
# minfu - minimum follow-up required for each patient in cohort to conduct model update
# recrate - # of pts recruited per obswin e.g. 5

#' @title Simulate TITE-CRM trials using specified design options
#'
#' @description applied_titecrm_sim is used to simulate trials using the
#'     time-to-event continual reassessment method with specified design options
#'     to determine the operating characteristics.
#'
#' @usage applied_titecrm_sim(true_tox, prior, target, max_sample_size,
#'     first_dose, num_sims, cohort_size = 1, obswin, minfu, recrate, dose_func
#'     = applied_titecrm, ...)
#'
#' @param true_tox A vector of 'true' underlying rates of toxicity for each of
#'     the dose levels.
#' @param prior A vector of prior estimates of toxicity probabilties for the
#'     dose levels.
#' @param target The target DLT rate.
#' @param max_sample_size The maximum number of patients to be recruited in any
#'     simulation.
#' @param first_dose The first dose level to tested.
#' @param num_sims The total number of simulations to be run.
#' @param cohort_size The size of the patient cohorts. Default is 1.
#' @param obswin The observation period for total patient follow up.
#' @param minfu The minimum amount of follow-up required for each patient.
#' @param recrate The number of patients recruited per obswin.
#' @param dose_func The function to be employed in executing the CRM. Default is
#'     applied_titecrm.
#' @param ... Any other arguements detailed in dtp::applied_titecrm.
#'
#' @return A list containg two further lists. The first of these lists contains
#'     the operating charateristics of the design, the second contains the
#'     underlying data for each of the simulation iterations.
#'
#' @references O'Quigley, J. O., Pepe, M., and Fisher, L. (1990). Continual
#'     reassessment method: A practical design for phase I clinical trials in
#'     cancer. Biometrics 46:33-48.
#'
#' Cheung, Y. K. (2011). Dose Finding by the Continual Reassessment Method. New
#' York: Chapman & Hall/CRC Press.
#'
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#'
#' @examples
#' prior  <- c(0.1, 0.3, 0.5)
#' target <- 0.2
#' true_tox <- c(0.05, 0.2, 0.35)
#' first_dose <- 1
#' num_sims <- 1000
#' obswin = 80
#'
#' applied_titecrm_sim(true_tox = true_tox, prior = prior, target = target,
#'                               max_sample_size = 21, first_dose = first_dose, num_sims = num_sims, cohort_size = 3,
#'                               obswin = obswin, minfu = 20, recrate = 3, dose_func = applied_titecrm)
#'
#' @keywords CRM Simulations TITE dtpcrm
#'
#' @export
applied_titecrm_sim <- function(true_tox, prior, target,
                            max_sample_size, first_dose,
                            num_sims, cohort_size = 1,
                            obswin, minfu,recrate,
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
    dose <- first_dose
    stop <- FALSE
    stop_reason <- NULL
    rectime <- obswin / recrate # time to recruit one patient



    while(!stop & length(tox) < max_sample_size) {

      # Simulate outcomes for a cohort
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
    }

    if(!stop){
    x <- dose_func(prior = prior, target = target, tox = tox, level = level,   # run model final time with full follow-up
                   followup = rep(obswin, length(tox)), obswin = obswin, ...)
    dose <- x$mtd
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
