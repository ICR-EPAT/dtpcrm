# Start: Remove (eventually)
# install.packages('testthat')
library(testthat)
library(dfcrm)

source('R/stopping_delegates.R')
source('R/applied_crm.R')
source('R/dtps.R')
# End: Remove


test_that("dtps are correct with previous outcomes and no stopping, v1", {

  source('R/dtps.R')

  prior  <- c(0.1, 0.2, 0.5)
  target <- 0.15
  prev_tox <- c(0, 0, 0)
  prev_dose <- c(2, 2, 2)
  cohort_sizes <- c(2, 3)
  next_dose = applied_crm(prior = prior, target = target,
                          tox = prev_tox, level = prev_dose)$mtd
  dose_func <- applied_crm
  dtps1 = calculate_dtps(next_dose, cohort_sizes, prev_tox = prev_tox,
                         prev_dose = prev_dose, dose_func = applied_crm,
                         prior = prior, target = target)

  # Reproduce from first principles
  dtps2 = matrix(nrow = prod(cohort_sizes + 1),
                 ncol = 1 + 2 * length(cohort_sizes))
  cohort1_dose <- rep(next_dose, cohort_sizes[1])

  # Row 1
  cohort1_tox <- c(0, 0)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 0)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[1, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 2
  cohort1_tox <- c(0, 0)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[2, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 3
  cohort1_tox <- c(0, 0)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[3, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 4
  cohort1_tox <- c(0, 0)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(1, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[4, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 5
  cohort1_tox <- c(0, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 0)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[5, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 6
  cohort1_tox <- c(0, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[6, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 7
  cohort1_tox <- c(0, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[7, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 8
  cohort1_tox <- c(0, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(1, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[8, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 9
  cohort1_tox <- c(1, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 0)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[9, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 10
  cohort1_tox <- c(1, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[10, ] = c(next_dose,
                  sum(cohort1_tox), cohort2_mtd,
                  sum(cohort2_tox), cohort3_mtd)

  # Row 11
  cohort1_tox <- c(1, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[11, ] = c(next_dose,
                  sum(cohort1_tox), cohort2_mtd,
                  sum(cohort2_tox), cohort3_mtd)

  # Row 12
  cohort1_tox <- c(1, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose))$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(1, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose))$mtd
  dtps2[12, ] = c(next_dose,
                  sum(cohort1_tox), cohort2_mtd,
                  sum(cohort2_tox), cohort3_mtd)

  # Compare
  expect_equal(all(dtps1 == dtps2), TRUE)
})

test_that("dtps are correct with previous outcomes and stopping, v1", {

  source('R/dtps.R')
  source('R/stopping_delegates.R')

  prior  <- c(0.1, 0.2, 0.5)
  target <- 0.15
  prev_tox <- c(0, 0, 0)
  prev_dose <- c(2, 2, 2)
  cohort_sizes <- c(2, 3)
  next_dose = applied_crm(prior = prior, target = target,
                          tox = prev_tox, level = prev_dose)$mtd

  stop_func <- function(x) {
    x = stop_for_excess_toxicity_empiric(x, tox_lim = 0.15, prob_cert = 0.8, dose = 1, nsamps = 10000)
  }

  dose_func <- applied_crm
  dtps1 = calculate_dtps(next_dose, cohort_sizes, prev_tox = prev_tox,
                         prev_dose = prev_dose, dose_func = applied_crm,
                         prior = prior, target = target, stop_func = stop_func)

  # Reproduce from first principles
  dtps2 = matrix(nrow = prod(cohort_sizes + 1),
                 ncol = 1 + 2 * length(cohort_sizes))
  cohort1_dose <- rep(next_dose, cohort_sizes[1])

  # Row 1
  cohort1_tox <- c(0, 0)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 0)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[1, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 2
  cohort1_tox <- c(0, 0)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[2, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 3
  cohort1_tox <- c(0, 0)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[3, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 4
  cohort1_tox <- c(0, 0)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(1, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[4, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 5
  cohort1_tox <- c(0, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 0)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[5, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 6
  cohort1_tox <- c(0, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[6, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 7
  cohort1_tox <- c(0, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[7, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 8
  cohort1_tox <- c(0, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(1, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[8, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 9
  cohort1_tox <- c(1, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 0)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[9, ] = c(next_dose,
                 sum(cohort1_tox), cohort2_mtd,
                 sum(cohort2_tox), cohort3_mtd)

  # Row 10
  cohort1_tox <- c(1, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 0, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[10, ] = c(next_dose,
                  sum(cohort1_tox), cohort2_mtd,
                  sum(cohort2_tox), cohort3_mtd)

  # Row 11
  cohort1_tox <- c(1, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(0, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[11, ] = c(next_dose,
                  sum(cohort1_tox), cohort2_mtd,
                  sum(cohort2_tox), cohort3_mtd)

  # Row 12
  cohort1_tox <- c(1, 1)
  cohort2_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox),
                             level = c(prev_dose, cohort1_dose), stop_func = stop_func)$mtd
  cohort2_dose <- rep(cohort2_mtd, cohort_sizes[2])
  cohort2_tox <- c(1, 1, 1)
  cohort3_mtd <- applied_crm(prior = prior, target = target,
                             tox = c(prev_tox, cohort1_tox, cohort2_tox),
                             level = c(prev_dose, cohort1_dose, cohort2_dose), stop_func = stop_func)$mtd
  dtps2[12, ] = c(next_dose,
                  sum(cohort1_tox), cohort2_mtd,
                  sum(cohort2_tox), cohort3_mtd)

  # Compare
  expect_equal(all(dtps1 == dtps2, na.rm = T), TRUE)
})
