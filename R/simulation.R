#' @title Simulate CRM trials using specified design options
#'
#' @description applied_crm_sim is used to simulate trials using the continual
#'     reassessment method with specified design options to determine the
#'     operating characteristics.
#'
#' @usage applied_crm_sim(true_tox, prior, target, max_sample_size, first_dose,
#'     num_sims, cohort_size = 1, dose_func = applied_crm, ...)
#'
#' @param true_tox A vector of 'true' underlying rates of toxicity for each of
#'     the dose levels.
#' @param prior A vector of prior estimates of toxicity probabilties for the
#'     dose levels.
#' @param target The target DLT rate.
#' @param max_sample_size The maximum number of subjects to be recruited in any
#'     simulation.
#' @param first_dose The first dose level to tested.
#' @param num_sims The total number of simulations to be run.
#' @param cohort_size The size of the cohorts. Default is 1.
#' @param dose_func The function to be employed in executing the CRM. Default is
#'     applied_crm.
#' @param ... Any other arguements detailed in dtpcrm::applied_crm.
#'
#' @return A list containing two further lists. The first of these lists contains
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
#' @examples
#' # It may take quite long for large num_sims
#' prior  <- c(0.1, 0.3, 0.5)
#' target <- 0.2
#' true_tox <- c(0.15, 0.25, 0.45)
#' first_dose <- 1
#' num_sims <- 5  # recommend doing 5000 simulations for the final design
#'
#' applied_crm_sim(true_tox, prior, target, max_sample_size = 30, first_dose,
#'                 num_sims, cohort_size = 1, dose_func = applied_crm)
#'
#' @keywords CRM Simulations dtpcrm
#'
#' @export
applied_crm_sim <- function(true_tox, prior, target,
                            max_sample_size, first_dose,
                            num_sims, cohort_size = 1,
                            dose_func = applied_crm,
                            ...) {

  iterations <- list()
  for(i in 1:num_sims) {
    # Start afresh.
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
    prob_stop = table(substr(unlist(sapply(iterations, function(x) x$stop_reason)), 1, 15)) / num_sims,
    mtd = sapply(1:length(prior), function(d)
      sum(dose_selections == d, na.rm = TRUE) / num_sims),
    # Summarise doses given to subjects
    doses_given = sapply(1:length(prior), function(d)
      sum(doses_given == d, na.rm = TRUE) / num_sims),
    prob_dose_given = sapply(1:length(prior), function(d)
      sum(doses_given == d, na.rm = TRUE) / length(doses_given))
  )


  return(list(summary = summary, iterations = iterations))
}
