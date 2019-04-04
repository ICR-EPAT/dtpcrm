#' @title Execute the TITE-CRM
#'
#' @description applied_titecrm is used to execute the time-to-event continual
#'     reassessment method with specified design options to determine the dose
#'     for the next subject.
#'
#' @usage applied_titecrm(prior, target, tox, level, followup, obswin,
#'     no_skip_esc = TRUE, no_skip_deesc = TRUE, global_coherent_esc = TRUE,
#'     stop_func = NULL, ...)
#'
#' @param prior A vector of prior estimates of toxicity probabilties for the
#'     dose levels.
#' @param target The target DLT rate.
#' @param tox A vector of subject outcomes; 1 indicates toxicity, 0 otherwise.
#' @param level A vector of dose levels assigned to subjects. The length of
#'     level must be equal to that of tox.
#' @param followup A vector of follow up times of subjects. The length must be equal
#'     to that of tox.
#' @param obswin The observation period with respect to which DLT is assessed. 
#' @param no_skip_esc If FALSE, the method will not enforce no skipping of doses
#'     in escalation. Default is TRUE.
#' @param no_skip_deesc If FALSE, the method will not enforce no skipping of
#'     doses in de-escalation. Default is TRUE.
#' @param global_coherent_esc If FALSE, the method will not enforce global
#'     coherent escalation, that is, escalation if the overall rate of toxicity
#'     seen at the current dose level is above the target rate. Default is TRUE.
#' @param stop_func An optional argument to provide a function which will
#'     utilised alongside the TITE-CRM to determine if the trial should be stopped.
#' @param ... Any other arguments detailed in dfcrm::titecrm.
#'
#' @details The adaptive weighting scheme is given in Cheung and Chappell (2000)
#'     given in the reference list.
#'
#' @return An object of class "mtd" is returned as per package "dfcrm",
#'     additional information is provided if a stopping function is used.
#'
#' \item{prior}{Initial guesses of toxicity rates.}  \item{target}{The target
#' probability of toxicity at the MTD.}  \item{ptox}{Updated estimates of
#' toxicity rates.}  \item{ptoxL}{Lower confidence/probability limits of
#' toxicity rates.}  \item{ptoxU}{Upper confidence/probability limits of
#' toxicity rates.}  \item{mtd}{The updated estimate of the MTD.}
#' \item{prior.var}{The variance of the normal prior.}  \item{post.var}{The
#' posterior variance of the model parameter.}  \item{estimate}{Estimate of the
#' model parameter.}  \item{method}{The method of estimation.}  \item{model}{The
#' working model.}  \item{dosescaled}{The scaled doses obtained via backward
#' substitution.}  \item{tox}{subjects' toxicity indications.}
#' \item{level}{Dose levels assigned to subjects.}  \item{followup}{Follow-up
#' times of subjects.}  \item{obswin}{Observation window with respect to which
#' DLT is assessed.}  \item{weights}{Weights assigned to subjects.}
#' \item{entry}{Entry times of subjects.}  \item{exit}{Exit times of subjects.}
#' \item{scheme}{Weighting scheme.}  \item{stop}{A logical variable detailing if
#' the trial should be stopped; TRUE to stop, FALSE otherwise}
#' \item{stop_reason}{A detailed reason for why the trial should be
#' stopped. Only provided if stop is TRUE}
#'
#' @references O'Quigley, J. O., Pepe, M., and Fisher, L. (1990). Continual
#'     reassessment method: A practical design for phase I clinical trials in
#'     cancer. Biometrics 46:33-48.
#'
#' Cheung, Y. K. (2011). Dose Finding by the Continual Reassessment Method. New
#' York: Chapman & Hall/CRC Press.
#'
#' Cheung, Y. K. and Chappell, R. (2000). Sequential designs for phase I
#' clinical trials with late-onset toxicities. Biometrics 56:1177-1182.
#'
#' @examples
#' prior    <- c(0.1, 0.3, 0.5)
#' target   <- 0.2
#' tox      <- c(0, 0, 1, 0, 1, 1)
#' level    <- c(1, 1, 1, 2, 2, 2)
#' followup <- c(96, 82, 77, 60, 51, 44)
#' obswin   <- 80
#'
#' applied_titecrm(prior = prior, target = target, tox = tox,level = level,
#'                 followup = followup, obswin = obswin)
#'
#' @keywords CRM TITE
#'
#' @export
applied_titecrm <- function(prior, target, tox, level, followup, obswin,
                        no_skip_esc = TRUE, no_skip_deesc = TRUE,
                        global_coherent_esc = TRUE,
                        stop_func = NULL, ...) {

  # Start with the dfcrm decision
  x <- dfcrm::titecrm(prior = prior, target = target, tox = tox, level = level,
                      followup = followup, obswin = obswin,
                      var.est = TRUE, ...)

  # Avoid skipping doses in escalation, if required
  if (no_skip_esc & x$mtd > (max(level) + 1)) {
    x$mtd <- max(level) + 1
  }

  # Avoid skipping doses in de-escalation, if required
  if (no_skip_deesc & x$mtd < (min(level) - 1)) {
    x$mtd <- min(level) - 1
  }

  # Global coherence in escalation means not escalating from a dose with an
  # observed toxicity rate exceeding the target
  if (global_coherent_esc) {
    last_dose <- utils::tail(level, 1)
    tox_rate_last_dose <- sum(tox[level == last_dose]) / sum(level == last_dose)
    if(tox_rate_last_dose > target) {
      x$mtd <- min(x$mtd, last_dose)
    }
  }

  # Invoke decision to determine whether trial should stop if stop_func is given
  if(!is.null(stop_func)) {
    x = stop_func(x)  # Let stopping delegate decorate x
  }

  return(x)
}
