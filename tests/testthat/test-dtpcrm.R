context("test-dtpcrm")

######################
######## DTPS ########
######################

test_that("dtps are correct with previous outcomes and no stopping, v1", {

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

#############################
######## APPLIED_CRM ########
#############################

test_that("applied_crm matches dfcrm::crm under neutral conditions, v1", {

  # By neutral conditions, I mean no logic to avoid skipping, enforce
  # coherence, or opt to stop

  prior  <- c(0.1, 0.2, 0.3)
  target <- 0.2
  tox    <- c(0, 0, 0, 0, 0, 0, 0, 1, 1)
  level  <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  obj1 <- dfcrm::crm(prior = prior, target = target, tox = tox, level = level)
  obj2 <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                      no_skip_esc = FALSE, no_skip_deesc = FALSE,
                      global_coherent_esc = FALSE, stop_func = NULL)
  expect_equal(obj1$mtd, obj2$mtd)
  expect_equal(obj1$ptox, obj2$ptox, tolerance = 0.001)
})

test_that("applied_crm matches dfcrm::crm under neutral conditions, v2", {

  # By neutral conditions, I mean no logic to avoid skipping, enforce
  # coherence, or opt to stop

  prior  <- c(0.224, 0.4203, 0.6478)
  target <- 0.378
  tox    <- c(0, 0, 0, 0, 0, 1, 0, 0, 0, 1)
  level  <- c(1, 2, 1, 1, 3, 3, 1, 2, 3, 3)
  obj1 <- dfcrm::crm(prior = prior, target = target, tox = tox, level = level)
  obj2 <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                      no_skip_esc = FALSE, no_skip_deesc = FALSE,
                      global_coherent_esc = FALSE, stop_func = NULL)
  expect_equal(obj1$mtd, obj2$mtd)
  expect_equal(obj1$ptox, obj2$ptox, tolerance = 0.001)
})

test_that("applied_crm omits stop data when stop_func omitted", {

  prior  <- c(0.1, 0.2, 0.3)
  target <- 0.2
  tox    <- c(0, 0, 1)
  level  <- c(1, 2, 3)
  obj <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                     no_skip_esc = FALSE, no_skip_deesc = FALSE,
                     global_coherent_esc = FALSE, stop_func = NULL)
  expect_null(obj$stop)
  expect_null(obj$stop_reason)
})

test_that("applied_crm does not skip dose in escalation, v1", {

  prior  <- c(0.05, 0.15, 0.3)
  target <- 0.3
  tox    <- c(0, 0, 0)
  level  <- c(1, 1, 1)
  obj <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                     no_skip_esc = TRUE, no_skip_deesc = FALSE,
                     global_coherent_esc = FALSE, stop_func = NULL)

  expect_lte(obj$mtd, max(level) + 1)
})

test_that("applied_crm does not skip dose in de-escalation, v1", {

  prior  <- c(0.1, 0.25, 0.4)
  target <- 0.05
  tox    <- c(0, 0, 1, 1)
  level  <- c(3,3, 3,3)
  obj <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                     no_skip_esc = FALSE, no_skip_deesc = TRUE,
                     global_coherent_esc = FALSE, stop_func = NULL)

  expect_gte(obj$mtd, min(level) - 1)
})


# test_that("applied_crm does not violate global coherency in escalation", {
#   # DS TODO
#   # For a valid test, the scenario chosen should want incoherent escalation.
#   # That could be tricky to find...
# })

test_that("applied_crm stops when max sample size is reached", {

  prior  <- c(0.1, 0.2, 0.3)
  target <- 0.2
  tox    <- c(0,0,0, 0,0,0, 0,1,1)
  level  <- c(1,1,1, 2,2,2, 3,3,3)

  stop_func = function(mtd) stop_for_sample_size(mtd, max_sample_size = 9)
  obj <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                     stop_func = stop_func)
  expect_equal(obj$stop, TRUE)
  expect_gt(nchar(obj$stop_reason), 1)
})

test_that("applied_crm does not stop when max sample size is not reached", {

  prior  <- c(0.1, 0.2, 0.7)
  target <- 0.333
  tox    <- c(0,0,0, 0,0,0)
  level  <- c(1,1,1, 2,2,2)

  stop_func = function(mtd) stop_for_sample_size(mtd, max_sample_size = 12)
  obj <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                     stop_func = stop_func)
  expect_equal(obj$stop, FALSE)
  expect_null(obj$stop_reason)
})

test_that("applied_crm stops when excess toxicity detected by empiric method", {

  prior  <- c(0.1, 0.2, 0.7)
  target <- 0.15
  tox    <- c(1,0, 1,1)
  level  <- c(1,1, 1,1)
  stop_func <- function(y) stop_for_excess_toxicity_empiric(y,
                                                            tox_lim = 0.25,
                                                            prob_cert = 0.9,
                                                            dose = 1)
  obj <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                     stop_func = stop_func)

  expect_equal(obj$stop, TRUE)
  expect_gt(nchar(obj$stop_reason), 1)
})

test_that("applied_crm does not stop when excess toxicity not detected by
          empiric method", {

            prior  <- c(0.1, 0.2, 0.7)
            target <- 0.15
            tox    <- c(0,0, 1,1)
            level  <- c(1,1, 2,2)
            stop_func <- function(y) stop_for_excess_toxicity_empiric(y, tox_lim = 0.25,
                                                                      prob_cert = 0.9,
                                                                      dose = 1)
            obj <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                               stop_func = stop_func)
            expect_equal(obj$stop, FALSE)
            expect_null(obj$stop_reason)
            })

test_that("applied_crm stops when excess toxicity detected by logit method", {

  prior  <- c(0.1, 0.2, 0.7)
  target <- 0.15
  tox    <- c(1,0, 1,1)
  level  <- c(1,1, 1,1)
  stop_func <- function(y) stop_for_excess_toxicity_logistic(y, tox_lim = 0.25,
                                                             prob_cert = 0.9,
                                                             dose = 1)
  obj <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                     stop_func = stop_func, model = 'logistic')

  expect_equal(obj$stop, TRUE)
  expect_gt(nchar(obj$stop_reason), 1)
})

test_that("applied_crm does not stop when excess toxicity not detected by
          logit method", {

            prior  <- c(0.1, 0.2, 0.7)
            target <- 0.15
            tox    <- c(0,0, 1,1)
            level  <- c(1,1, 2,2)
            stop_func <- function(y) stop_for_excess_toxicity_empiric(y, tox_lim = 0.25,
                                                                      prob_cert = 0.9,
                                                                      dose = 1)
            obj <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                               stop_func = stop_func, model = 'logistic')
            expect_equal(obj$stop, FALSE)
            expect_null(obj$stop_reason)
            })

test_that("applied_crm stops when consensus reached", {

  prior  <- c(0.1, 0.2, 0.5)
  target <- 0.15
  tox    <- c(0,0,0, 1,0,1, 0,0,1)
  level  <- c(1,1,1, 2,2,2, 1,1,1)
  stop_func <- function(y) stop_for_consensus_reached(y, req_at_mtd = 6)
  obj <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                     stop_func = stop_func)
  expect_equal(obj$stop, TRUE)
  expect_gt(nchar(obj$stop_reason), 1)
})

test_that("applied_crm does not stop when consensus not reached", {

  prior  <- c(0.1, 0.2, 0.5)
  target <- 0.15
  tox    <- c(0,0,0, 1,0,0, 0,0,0)
  level  <- c(1,1,1, 2,2,2, 1,1,1)
  stop_func <- function(y) stop_for_consensus_reached(y, req_at_mtd = 6)
  obj <- applied_crm(prior = prior, target = target, tox = tox, level = level,
                     stop_func = stop_func)
  expect_equal(obj$stop, FALSE)
  expect_null(obj$stop_reason)
})

