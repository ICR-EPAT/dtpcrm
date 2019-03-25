

# stop_for_sample_size <- function(x, max_sample_size) {
#
#   # x is an object isomorphic to that returned by dfcrm:crm
#
#   stop_decision = length(x$level) >= max_sample_size
#   if(stop_decision) {
#     stop_reason = paste("Maximum sample size of", max_sample_size, "reached.")
#     return(list(stop_decision, stop_reason))
#   } else {
#     return(stop_decision)
#   }
# }
#
# stop_for_excess_toxicity_empiric <- function(x, tox_lim, prob_cert, dose = 1,
#                                              nsamps=10^5) {
#
#   # If x was estimated with est.var=F this will fail cos x$post.var will be NULL
#   post_beta_mean = x$estimate
#   post_beta_var  = x$post.var
#   post_beta_samp = rnorm(n = nsamps, mean = post_beta_mean, sd = post_beta_var)
#   post_prob_tox_samp = x$prior[dose] ^ exp(post_beta_samp)
#   prob_too_toxic = mean(post_prob_tox_samp > tox_lim)
#   stop_decision = prob_too_toxic > prob_cert
#   if(stop_decision) {
#     # stop_reason = paste0("Probability of toxicity being > ", tox_lim,
#     #                      " at dose ", dose, " is > ", prob_cert)
#     stop_reason = paste0("Prob(Prob(Tox[", dose,"]) > ", tox_lim, ") = ",
#                          round(prob_too_toxic, 3), " > ", prob_cert)
#     return(list(stop_decision, stop_reason))
#   } else {
#     return(stop_decision)
#   }
# }
#
#
# stop_for_excess_toxicity_logistic <- function(x, tox_lim, prob_cert, dose = 1,
#                                               nsamps = 10^5) {
#   # If x was estimated with est.var=F this will fail cos x$post.var will be NULL
#   post_beta_mean = x$estimate
#   post_beta_var  = x$post.var
#   post_beta_samp = rnorm(nsamps, post_beta_mean, post_beta_var)
#   post_prob_tox_samp = exp(x$intcpt + exp(post_beta_samp) * x$prior[dose]) /
#     (1 + exp(x$intcpt + exp(post_beta_samp) * x$prior[dose]))
#   prob_too_toxic = mean(post_prob_tox_samp > tox_lim)
#   stop_decision = prob_too_toxic > prob_cert
#   if(stop_decision) {
#     stop_reason = paste0("Probability of toxicity being > ", tox_lim,
#                          " at dose ", dose, " is > ", prob_cert)
#     return(list(stop_decision, stop_reason))
#   } else {
#     return(stop_decision)
#   }
# }
#
# stop_for_consensus_reached <- function(x, req_at_mtd) {
#
#   cur_est   = x$mtd
#   num_treat = sum(level == cur_est)
#
#   stop_decision = num_treat >= req_at_mtd
#   if(stop_decision) {
#     stop_reason = paste0("Consensus Reached - ", req_at_mtd, " treated at dose ", cur_est)
#     return(list(stop_decision, stop_reason))
#   } else {
#     return(stop_decision)
#   }
# }


# I had to change this pattern because it did not work:
# it did not differentiate good and bad stopping.
# In bad stopping, you want to select no dose.
# Thus, these delegates now decorate x, rather than returning
# just stopping information. I am not surte it is right or wrong
# but it allows differentiation between good & bad stopping.
# I checked tests; they run OK.
# The dtp calls may need to change.
# Apologies for disruption. KB

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

#' @export
stop_for_excess_toxicity_empiric <- function(x, tox_lim, prob_cert, dose = 1,
                                             nsamps=10^5,
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

#' @export
stop_for_excess_toxicity_logistic <- function(x, tox_lim, prob_cert, dose = 1,
                                              nsamps = 10^5,
                                              suppress_dose = TRUE) {
  # If x was estimated with est.var=F this will fail cos x$post.var will be NULL
  post_beta_mean = x$estimate
  post_beta_var  = x$post.var
  post_beta_samp = stats::rnorm(nsamps, post_beta_mean, sqrt(post_beta_var))
  post_prob_tox_samp = exp(x$intcpt + exp(post_beta_samp) * x$prior[dose]) /
    (1 + exp(x$intcpt + exp(post_beta_samp) * x$prior[dose]))
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
