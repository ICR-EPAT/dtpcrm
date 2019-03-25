

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

#' @title Calculate Dose-Transition Pathways (DTPs)
#'
#' @description Calculate Dose-Transition Pathways (DTPs) for future cohorts in a CRM-like trial.
#'
#' @param next_dose the dose that will be given to the very next cohort
#' @param cohort_sizes vector of cohort sizes for future paths, e.g. c(2, 3) considers a cohort of 2 followed by a cohort of 3
#' @param prev_tox vector of (bool) toxicity events already observed
#' @param prev_dose, vector of dose-levels already given
#' @param dose_func function that will perform the dose-finding calculation. It should take tox and level as args. Other args are passed through ... ,
#            It should return the next dose (int) or object like dfcrm:: mtd. Default function is applied_crm
#' @param ... other arg
#' @export
calculate_dtps = function(next_dose, cohort_sizes, prev_tox = c(),
                          prev_dose = c(), dose_func = applied_crm, ...) {
  # Calculate Dose-Transition Pathways (DTPs) for future cohorts in a CRM-like
  #   trial. The first cohort will be receive next_dose, conditional
  #   on having already observed prev_tox outcomes at prev_dose doses (optional)
  # Dose decisions are made using dose_func, taking args tox and level, & ...
  # dose_func should return either a dose-selection or an dfcrm::mtd-like object
  # When using dose_func = applied_crm or dfcrm::crm, prior & target should
  # be provided via ...
  #
  # Params:
  # next_dose, the dose that will be given to the very next cohort
  # cohort_sizes, vector of cohort sizes for future paths,
  #   e.g. c(2, 3) considers a cohort of 2 followed by a cohort of 3
  # prev_tox, vector of (bool) toxicity events already observed
  # prev_dose, vector of dose-levels already given
  # dose_func, function that will perform the dose-finding calculation
  #            It should take tox and level as args.
  #            Other args are passed through ...
  #            It should return the next dose (int) or object like dfcrm:: mtd
  #            Default function is applied_crm

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

# order(paths[,1], paths[,2], paths[,3])
# paths[order(paths[,1], paths[,2], paths[,3]), ]
#
# do.call(order, as.data.frame(paths))
# paths[do.call(order, as.data.frame(paths)),]
