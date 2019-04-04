

.conduct_dose_finding_cohorts <- function(next_dose, tox_counts, cohort_sizes,
                                          prev_tox = c(), prev_dose = c(),
                                          dose_func = applied_crm, ...) {
  # This is a helper function to dtps
  # It calculates the doses recommended to a set of future cohorts that have
  # already observed prev_tox outcomes at doses prev_dose.
  # Dose decisions are made using dose_func, taking args tox, level & ...

  # next_dose, the dose that will be given to the very next cohort
  # tox_counts, vector of toxicity counts in cohorts.
  # cohort_sizes, vector of future cohort sizes
  # prev_tox, vector of (bool) toxicity events already observed
  # prev_dose, vector of dose-levels already given
  # dose_func, function that will perform the dose-finding calculation
  #            func should take tox and level as args, plus ...
  #            func should return the next dose (int) or an object with name mtd
  if (length(prev_dose) != length(prev_tox)) {
    stop ('prev_doses and prev_tox should be same length.')
  }

  toxes       = prev_tox
  doses       = prev_dose
  dose_recs   = integer(length(tox_counts))
  dose        = next_dose
  num_cohorts = length(tox_counts)

  for (i in 1:num_cohorts) {
    these_toxes = c(rep(1, tox_counts[i]), rep(0, cohort_sizes[i] - tox_counts[i]))
    toxes       = c(toxes, these_toxes)
    these_doses = rep(dose, cohort_sizes[i])
    doses       = c(doses, these_doses)
    x           = dose_func(tox=toxes, level=doses, ...)

    if ("mtd" %in% names(x)) {
      dose = x$mtd
    } else {
      dose = x
    }

    dose_recs[i] = dose

    if ("stop" %in% names(x)) {
      if (x[['stop']]) {
        dose_recs[i:num_cohorts] = NA
        break
      }
    }
  }

  return(dose_recs)
}

#' @title Produce the Dose Transition Pathways
#'
#' @description calculate_dtps is used to produce the dose transition pathways
#'     for the continual reassessment method with specified design
#'     options. These pathways present the possible model recommendations based
#'     on all permumations of trial outcomes.
#'
#' @usage calculate_dtps(next_dose, cohort_sizes, prev_tox = c(), prev_dose =
#'     c(), dose_func = applied_crm, ...)
#'
#' @param next_dose An integer value representing the dose to be assigned to the
#'     first cohort of subjects in the pathways.
#' @param cohort_sizes A vector of cohort sizes representing the size of the
#'     cohorts to be treated with the recommended dose at each decision point.
#' @param prev_tox A vector of previous subject outcomes; 1 indicates toxicity,
#'     0 otherwise.
#' @param prev_dose A vector of previous subject doses; The length of prev_dose
#'     must be equal to that of prev_tox.
#' @param dose_func A function such as applied_crm which produces an object of
#'     class 'mtd'. To be used for calculation of the next recommended dose for
#'     each pathway permutation.
#' @param ... Any other arguments to be passed to dose_func; for specific
#'     arguments related to applied_crm see.
#'
#' @return Produces a dataframe containing all possible permutations of outcomes
#'     for each cohort based on cohort_sizes and the recommended
#'     doses for such permutations.
#'
#' @examples
#' prior  <- c(0.1, 0.2, 0.5)
#' target <- 0.15
#' prev_tox <- c(0, 0, 0)
#' prev_dose <- c(2, 2, 2)
#' cohort_sizes <- c(2, 3)
#'
#' next_dose = applied_crm(prior = prior, target = target,
#'                         tox = prev_tox, level = prev_dose)$mtd
#'
#' dose_func <- applied_crm
#'
#' DTP = calculate_dtps(next_dose, cohort_sizes, prev_tox = prev_tox,
#'                       prev_dose = prev_dose, dose_func = applied_crm,
#'                       prior = prior, target = target)
#'
#' @keywords CRM DTP dtpcrm
#'
#' @export
calculate_dtps = function(next_dose, cohort_sizes, prev_tox = c(),
                          prev_dose = c(), dose_func = applied_crm, ...) {

  # Helper functions
  # 1) This function produces a row for the dtp data.frame
  # next_dose is the immediately recommended next dose, d0 say
  # tox_counts is a vector of the counts of DLTs observed, (t1, t2, t3) say
  # dose_recs is a congruent vector of the doses recommended, (d1, d2, d3) say
  # This yields c(d0, t1, d1, t2, d2, t3, d3)
  .make_dtp_row = function(next_dose, tox_counts, dose_recs) {
    return (c(next_dose, as.vector(rbind(tox_counts, dose_recs))))
  }

  num_cohorts <- length(cohort_sizes)
  feasible_tox_counts <- lapply(cohort_sizes, function(x) 0:x)
  paths <- expand.grid(feasible_tox_counts)
  # Order by first col, then second col, etc
  paths <- paths[do.call(order, as.data.frame(paths)),]
  # Reset row names
  row.names(paths) <- 1:nrow(paths)

  # Invoke DTP calculation on each permutation of the toxicity counts
  dtps <- apply(paths, 1, function(x) .make_dtp_row(
    next_dose, x, .conduct_dose_finding_cohorts(
      next_dose, x, cohort_sizes, prev_tox = prev_tox, prev_dose = prev_dose,
      dose_func = dose_func, ...)))
  dtps <- t(dtps)
  dtps <- data.frame(dtps)
  colnames(dtps) <- c('D0', as.vector(rbind(paste0('T', 1:num_cohorts),
                                            paste0('D', 1:num_cohorts))))

  dtps[t(apply(is.na(dtps), 1, cumsum)) > 0 ] <- NA # change to NA for all columns after first NA
  return(dtps)
}
