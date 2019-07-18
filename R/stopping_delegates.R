#' @title Stopping for sample size reached
#'
#' @description This is a function for use with applied_crm for the stop_func
#'     argument. The rule will suggest stopping in the scenario that a maximum
#'     number of subjects has been recruited.
#'
#' @usage stop_for_sample_size(x, max_sample_size)
#'
#' @param x An object of class 'mtd'.
#' @param max_sample_size An integer; specifying the maxmium number of subjects
#'     to be recruited.
#'
#' @details This function is an example of a possible stopping function to be
#'     used with applied_crm, it will modifiy the 'mtd' class object produced by
#'     applied_crm to include a logical value under the name 'stop' indicting
#'     whether or not the trial should stop. The package dtpcrm contains a few
#'     of these functions for possible use with applied_crm.
#'
#' @examples
#' prior  <- c(0.1, 0.3, 0.5)
#' target <- 0.2
#' tox    <- c(0, 0, 1, 0, 1, 1)
#' level  <- c(1, 1, 1, 2, 2, 2)
#'
#' stop_rule <- function(x){
#'   x <- stop_for_sample_size(x, max_sample_size = 20)
#' }
#'
#' crm <- applied_crm(prior, target, tox, level, no_skip_esc = TRUE, no_skip_deesc = TRUE,
#'                    global_coherent_esc = TRUE, stop_func = stop_rule)
#'
#' @keywords CRM stop dtpcrm
#'
#' @export
stop_for_sample_size <- function(x, max_sample_size) {

  # x is an object isomorphic to that returned by dfcrm:crm

  stop_decision = length(x$level) >= max_sample_size
  x$stop = stop_decision
  if(stop_decision) {
    x$stop_reason = paste("Maximum sample size of", max_sample_size, "reached.")
  }
  return(x)
}

#' @title Stopping for excess toxicity - Empiric method
#'
#' @description This is a function for use with applied_crm for the stop_func
#'     argument. The rule will suggest stopping in the scenario that the
#'     probability of toxicity being greater than a specifed value at a defined
#'     dose is greater than some further specified certainty value.
#'
#' @usage stop_for_excess_toxicity_empiric(x, tox_lim, prob_cert, dose = 1,
#'     nsamps = 10^6, suppress_dose = TRUE)
#'
#' @param x An object of class 'mtd'.
#' @param tox_lim A numeric; specifying the value for which the estimated
#'     toxicity at the selcted dose is not to exceed.
#' @param prob_cert A numeric; specifying the probability value to be used when
#'     assessing the certainty required that toxicty at the specificed dose
#'     exceeds tox_lim.
#' @param dose An integer; the dose to be assessed.
#' @param nsamps number of samples used for beta in the underlying normal
#'     sampling of beta.
#' @param suppress_dose A logical value indicating if the MTD should be set to
#'     NA if trial should stop.
#'
#' @details This function is an example of a possible stopping function to be
#'     used with applied_crm, it will modifiy the 'mtd' class object produced by
#'     applied_crm to include a logical value under the name 'stop' indicating
#'     whether or not the trial should stop. The package dtpcrm contains a few
#'     of these functions for possible use with applied_crm.
#'
#' @examples
#' prior  <- c(0.1, 0.3, 0.5)
#' target <- 0.2
#' tox    <- c(0, 0, 1, 0, 1, 1)
#' level  <- c(1, 1, 1, 2, 2, 2)
#'
#' stop_rule <- function(x){
#'   x <- stop_for_excess_toxicity_empiric(x, tox_lim = 0.25, prob_cert = 0.85)
#' }
#'
#' crm <- applied_crm(prior, target, tox, level, no_skip_esc = TRUE, no_skip_deesc = TRUE,
#'                    global_coherent_esc = TRUE, stop_func = stop_rule)
#'
#' @keywords CRM stop dtpcrm
#'
#' @export
stop_for_excess_toxicity_empiric <- function(x, tox_lim, prob_cert, dose = 1,
                                             nsamps=10^6,
                                             suppress_dose = TRUE) {

  # If x was estimated with est.var=F this will fail cos x$post.var will be NULL
  post_beta_mean = x$estimate
  post_beta_var  = x$post.var
  post_beta_samp = stats::rnorm(n = nsamps, mean = post_beta_mean, sd = sqrt(post_beta_var))
  post_prob_tox_samp = x$prior[dose] ^ exp(post_beta_samp)
  prob_too_toxic = mean(post_prob_tox_samp > tox_lim)
  stop_decision = prob_too_toxic > prob_cert
  x$stop = stop_decision
  if(stop_decision) {
    # stop_reason = paste0("Probability of toxicity being > ", tox_lim,
    #                      " at dose ", dose, " is > ", prob_cert)
    x$stop_reason = paste0("Prob(Prob(Tox[", dose,"]) > ", tox_lim, ") = ",
                           round(prob_too_toxic, 3), " > ", prob_cert)
    if(suppress_dose)
      x$mtd = NA
  }
  return(x)
}

#' @title Stopping for excess toxicity - Logistic method
#'
#' @description This is a function for use with applied_crm for the stop_func
#'     argument. The rule will suggest stopping in the scenario that the
#'     probability of toxicity being greater than a specifed value at a defined
#'     dose is greater than some further specified certainty value.
#'
#' @usage stop_for_excess_toxicity_logistic(x, tox_lim, prob_cert, dose = 1,
#'     nsamps = 10^6, suppress_dose = TRUE)
#'
#' @param x An object of class 'mtd'.
#' @param tox_lim A numeric; specifying the value for which the estimated
#'     toxicity at the selcted dose is not to exceed.
#' @param prob_cert A numeric; specifying the probability value to be used when
#'     assessing the certainty required that toxicty at the specificed dose
#'     exceeds tox_lim.
#' @param dose An integer; the dose to be assessed.
#' @param nsamps Number of samples used for beta in the underlying normal
#'     sampling of beta.
#' @param suppress_dose A logical value indicating if the MTD should be set to
#'     NA if trial should stop.
#'
#' @details This function is an example of a possible stopping function to be
#'     used with applied_crm, it will modifiy the 'mtd' class object produced by
#'     applied_crm to include a logical value under the name 'stop' indicating
#'     whether or not the trial should stop. The package dtpcrm contains a few
#'     of these functions for possible use with applied_crm.
#'
#' @examples
#' prior  <- c(0.1, 0.3, 0.5)
#' target <- 0.2
#' tox    <- c(0, 0, 1, 0, 1, 1)
#' level  <- c(1, 1, 1, 2, 2, 2)
#'
#' stop_rule <- function(x){
#'   x <- stop_for_excess_toxicity_logistic(x, tox_lim = 0.25, prob_cert = 0.85)
#' }
#'
#' crm <- applied_crm(prior, target, tox, level, no_skip_esc = TRUE, no_skip_deesc = TRUE,
#'                    global_coherent_esc = TRUE, stop_func = stop_rule)
#'
#' @keywords CRM stop dtpcrm
#'
#' @export
stop_for_excess_toxicity_logistic <- function(x, tox_lim, prob_cert, dose = 1,
                                              nsamps = 10^6,
                                              suppress_dose = TRUE) {
  # If x was estimated with est.var=F this will fail cos x$post.var will be NULL
  post_beta_mean = x$estimate
  post_beta_var  = x$post.var
  post_beta_samp = stats::rnorm(nsamps, post_beta_mean, sqrt(post_beta_var))
  post_prob_tox_samp = exp(x$intcpt + exp(post_beta_samp) * x$dosescaled[dose]) /
    (1 + exp(x$intcpt + exp(post_beta_samp) * x$dosescaled[dose]))
  prob_too_toxic = mean(post_prob_tox_samp > tox_lim)
  stop_decision = prob_too_toxic > prob_cert
  x$stop = stop_decision
  if(stop_decision) {
    x$stop_reason = paste0("Probability of toxicity being > ", tox_lim,
                           " at dose ", dose, " is > ", prob_cert)
    x$mtd = NA
  }
  return(x)
}

#' @title Stopping for consensus
#'
#' @description This is a function for use with applied_crm for the stop_func
#'     argument. The rule will suggest stopping in the scenario that a
#'     particular number of patients has already been treated at the current
#'     recommended MTD.
#'
#' @usage stop_for_consensus_reached(x, req_at_mtd)
#'
#' @param x An object of class 'mtd'.
#' @param req_at_mtd An integer; the number of patients required at current
#'     estimate of MTD to suggest stopping for consensus.
#'
#' @details This function is an example of a possible stopping function to be
#'     used with applied_crm, it will modifiy the 'mtd' class object produced by
#'     applied_crm to include a logical value under the name 'stop' indicting
#'     whether or not the trial should stop. The package dtpcrm contains a few
#'     of these functions for possible use with applied_crm.
#'
#' @examples
#' prior  <- c(0.1, 0.3, 0.5)
#' target <- 0.2
#' tox    <- c(0, 0, 1, 0, 1, 1)
#' level  <- c(1, 1, 1, 2, 2, 2)
#'
#' stop_rule <- function(x){
#'   x <- stop_for_consensus_reached(x, req_at_mtd = 6)
#' }
#'
#' crm <- applied_crm(prior, target, tox, level, no_skip_esc = TRUE, no_skip_deesc = TRUE,
#'                    global_coherent_esc = TRUE, stop_func = stop_rule)
#'
#' @keywords CRM stop dtpcrm
#' @export
stop_for_consensus_reached <- function(x, req_at_mtd) {

  cur_est   = x$mtd
  num_treat = sum(x$level == cur_est)

  stop_decision = num_treat >= req_at_mtd
  x$stop = stop_decision
  if(stop_decision) {
    x$stop_reason = paste0("Consensus Reached - ", req_at_mtd, " treated at dose ", cur_est)
  }
  return(x)
}
