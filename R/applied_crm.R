
#' @title Execute the CRM
#'
#' @description applied_crm is used to execute the continual reassessment method
#'     with specified design options to determine the dose for the next subject.
#'
#' @param prior A vector of prior estimates of toxicity probabilties for the
#'     dose levels.
#' @param target The target DLT rate.
#' @param tox A vector of subject outcomes; 1 indicates toxicity, 0 otherwise.
#' @param level A vector of dose levels assigned to subjects. The length of
#'     level must be equal to that of tox.
#' @param no_skip_esc If FALSE, the method will not enforce no skipping of doses
#'     in escalation. Default is TRUE.
#' @param no_skip_deesc If FALSE, the method will not enforce no skipping of
#'     doses in de-escalation. Default is TRUE.
#' @param global_coherent_esc If FALSE, the method will not enforce global
#'     coherent escalation, that is, escalation if the overall rate of toxicity
#'     seen at the current dose level is above the target rate. Default is TRUE.
#' @param stop_func An optional argument to provide a function which will
#'     utilised alongside the CRM to determine if the trial should be stopped.
#' @param ... Any other arguments detailed in dfcrm::crm.
#'
#' @details For maximum likelihood estimation, the variance of the estimate of beta
#'     (post.var) is approximated by the posterior variance of beta with a
#'     dispersed normal prior.
#'
#' The empiric model is specified as F(d, beta) = d^{exp(beta)}. The logistic model is
#' specified as logit (F(d,beta)) = intcpt + exp(beta) * d. For method="bayes", the
#' prior on beta is normal with mean 0. Exponentiation of beta ensures an increasing
#' dose-toxicity function.
#'
#' This function is largely a wrapper for the dfcrm function crm.  It provides
#' functionality for additional design choices for the CRM including global coherency
#' and stopping for excess toxicity and stopping when sufficient number of subjects are dosed at MTD.
#'
#'
#' @return An object of class "mtd" is returned as per package "dfcrm",
#'     additional information is provided if a stopping function is used.
#'     \item{prior}{Initial guesses of toxicity rates.}  \item{target}{The
#'     target probability of toxicity at the MTD.}  \item{ptox}{Updated
#'     estimates of toxicity rates.}  \item{ptoxL}{Lower confidence/probability
#'     limits of toxicity rates.}  \item{ptoxU}{Upper confidence/probability
#'     limits of toxicity rates.}  \item{mtd}{The updated estimate of the MTD.}
#'     \item{prior.var}{The variance of the normal prior.}  \item{post.var}{The
#'     posterior variance of the model parameter.}  \item{estimate}{Estimate of
#'     the model parameter.}  \item{method}{The method of estimation.}
#'     \item{model}{The working model.}  \item{dosescaled}{The scaled doses
#'     obtained via backward substitution.}  \item{tox}{Patients' toxicity
#'     indications.}  \item{level}{Dose levels assigned to patients.}
#'     \item{stop}{A logical variable detailing if the trial should be stopped;
#'     TRUE to stop, FALSE otherwise} \item{stop_reason}{A detailed reason for
#'     why the trial should be stopped. Only provided if stop is TRUE}
#'
#' @references O'Quigley, J. O., Pepe, M., and Fisher, L. (1990). Continual
#'     reassessment method: A practical design for phase I clinical trials in
#'     cancer. Biometrics 46:33-48.
#'
#' Cheung, Y. K. (2011). Dose Finding by the Continual Reassessment Method. New
#' York: Chapman & Hall/CRC Press.
#'
#' @examples
#' prior  <- c(0.1, 0.3, 0.5)
#' target <- 0.2
#' tox    <- c(0, 0, 1, 0, 1, 1)
#' level  <- c(1, 1, 1, 2, 2, 2)
#' applied_crm(prior, target, tox, level, no_skip_esc = TRUE, no_skip_deesc = TRUE,
#'             global_coherent_esc = TRUE, stop_func = NULL)
#'
#' @keywords CRM
#'
#' @export
applied_crm <- function(prior, target, tox, level,
                        no_skip_esc = TRUE, no_skip_deesc = TRUE,
                        global_coherent_esc = TRUE,
                        stop_func = NULL, ...) {

  # Start with the dfcrm decision
  x <- dfcrm::crm(prior = prior, target = target, tox = tox, level = level,
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

#' @title Provide a summary of applied_crm output
#'
#' @description summary_crm is used to return a dataframe of the summary of the
#'     output from applied_crm.
#'
#' @usage summary_crm(x)
#'
#' @param x An object assigned to be the output from applied_crm.
#'
#' @details This function takes an object of class "mtd" and produces a
#'     dataframe containing a summary of information within the
#'     object. Specifically it shows the dose levels, prior probabilities, number
#'     of evaluable patients, number of DLTs and the posterior probability
#'     estimates along with confidence/probability intervals if estimated in the
#'     underlying object.
#'
#' @return Dataframe of the summary of the output from applied_crm.
#'
#' @references O'Quigley, J. O., Pepe, M., and Fisher, L. (1990). Continual
#'     reassessment method: A practical design for phase I clinical trials in
#'     cancer. Biometrics 46:33-48.
#'
#' Cheung, Y. K. (2011). Dose Finding by the Continual Reassessment Method. New
#' York: Chapman & Hall/CRC Press.
#'
#' @examples
#' prior  <- c(0.1, 0.3, 0.5)
#' target <- 0.2
#' tox    <- c(0, 0, 1, 0, 1, 1)
#' level  <- c(1, 1, 1, 2, 2, 2)
#'
#' crm_obj <- applied_crm(prior, target, tox, level, no_skip_esc = TRUE, no_skip_deesc = TRUE,
#'                        global_coherent_esc = TRUE, stop_func = NULL)
#'
#' summary_crm(crm_obj)
#'
#' @keywords CRM dtpcrm
#'
#' @export
summary_crm <- function(x) {
  summary <- data.frame('Dose.level' = c(1:length(x$prior)), 'Prior.Prob(DLT)' = x$prior,
                        'Number.of.Evaluable.Patients' = rep(NA, length(x$prior)),
                        'Number.of.DLTs' = rep(NA, length(x$prior)),
                        'Posterior.Prob' = round(x$ptox, digits = 3))

  for(i in 1:length(x$prior)) {
    summary$Number.of.Evaluable.Patients[i] = sum(x$level == i)
    summary$Number.of.DLTs[i] = sum(x$tox[x$level == i])
  }

  summary$`Posterior.Prob` <- as.character(summary$`Posterior.Prob`)
  for(i in 1:length(x$ptox)) {
    summary$`Posterior.Prob`[i] = paste0(summary$`Posterior.Prob`[i], ' (',
                                                 round(x$ptoxL[i], digits = 3), ', ',
                                                 round(x$ptoxU[i], digits = 3), ')')
  }

  names(summary) <- c('Dose level', 'Prior Prob(DLT)', 'Number of Evaluable Patients',
                        'Number of DLTs', 'Posterior Prob(DLT)')
  return(summary)
}


#' @title Plot of posterior estimates from the CRM
#'
#' @description Provides functionality for plotting the posterior estimates of
#'     probabilities of toxicity at each dose level for both the most recent
#'     update and for past cohort updates if specified.
#'
#' @usage plot_crm(crm, dose_labels, cohort_sizes = NULL, file = NULL,
#'                 height = 600, width = 750, dose_func = NULL, ...,
#'                 ylim = c(0, 1), lwd = 1, cex.axis = 1, cex.lab = 1,
#'                 cex = 1, cohort.last = F)
#'
#' @param crm An object of class 'mtd' produced by applied_crm to be plotted.
#' @param dose_labels A vector of character strings detailing the labels to be
#'     used for each dose level in the plot.
#' @param cohort_sizes An optional vector of cohort sizes; if provided the
#'     previous estimates for each cohort will be plotted in addition.
#' @param file An optional string for the file name; if provided the plot will
#'     be saved as a .PNG to the current working directory under the provided
#'     file name.
#' @param height A numeric value specifying the vertical pixel count of the
#'     plot. Default is 600.
#' @param width A numeric value specifiying the horizontal pixel count of the
#'     plot. Default is 750.
#' @param dose_func Must be provided if cohort_sizes is provided. The function
#'     to be used to when implementing the CRM for previous cohorts.
#' @param ... Arguments to be provided to dose_func detailing CRM specification.
#'     See applied_crm.
#' @param ylim The y-axis range. Default is c(0, 1)
#' @param lwd line width relative to the default (default=1). 2 is twice as wide. Default is 1.
#' @param cex.axis The magnification to be used for axis annotation relative to the current setting of cex. Default is 1.
#' @param cex.lab The magnification to be used for x and y labels relative to the current setting of cex. Default is 1.
#' @param cex A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. Default is 1.
#' @param cohort.last If TRUE, the last cohort will have lwd = 6 for emphasis. Default is FALSE.
#'
#' @details Produces a plot of current dose-toxicity estimates including the
#'     priors and outputs a .png of plot to current directory if 'file' is
#'     provided. Potential for history of estimates by cohort if cohort.sizes is
#'     provided; dose_func is required to do this.
#'
#' @examples
#' prior  <- c(0.1, 0.3, 0.5)
#' target <- 0.2
#' tox    <- c(0, 0, 1, 0, 1, 1)
#' level  <- c(1, 1, 1, 2, 2, 2)
#'
#' crm <- applied_crm(prior, target, tox, level, no_skip_esc = TRUE, no_skip_deesc = TRUE,
#'                    global_coherent_esc = TRUE, stop_func = NULL)
#'
#' plot_crm(crm, dose_labels = c("1", "2", "3"))
#'
#' @keywords CRM plot dtpcrm
#'
#' @export
plot_crm <- function(crm, dose_labels, cohort_sizes = NULL, file = NULL,
                     height = 600, width = 750, dose_func = NULL, ..., ylim = c(0, 1),
                     lwd = 1, cex.axis = 1, cex.lab = 1, cex = 1,
                     cohort.last = F) {

  if(is.null(cohort_sizes)) {

    if(!is.null(file)) {
      grDevices::png( file = paste0(file, '.png'), height = height, width = width)
      graphics::par(cex.axis = cex.axis, cex.lab = cex.lab)
    }

    doses     <- c(1:length(crm$prior))
    postprob  <- crm$ptox
    priorprob <- crm$prior

    graphics::plot(x = doses, y = postprob, col = 2, type = 'b', xlab = 'Dose', xaxt = 'n',
         ylab = 'Probability of Dose Limiting Toxicity', ylim = ylim, lwd = lwd)
    graphics::points(x= doses, y =priorprob, type = 'b', col = 1, lwd = lwd)
    graphics::abline(h = crm$target, lty = 2)
    graphics::axis(1, at = doses, labels = dose_labels)
    graphics::legend(x = max(postprob, priorprob) + 0.15, col = c(1, 2), lty = 1, lwd = lwd,
           legend = c('Prior Curve', 'Posterior Curve'), cex = cex)

    if(!is.null(file)) {
      grDevices::dev.off()
    }

  } else {

    if(is.null(dose_func)) {
      stop('dose_func required for cohort history plot')
    }

    if(!is.null(file)) {
      grDevices::png( file = paste0(file, '.png'), height = height, width = width )
      graphics::par(cex.axis = cex.axis, cex.lab = cex.lab)
    }

    colours <- c(1, 2, 'chartreuse4', 'darkgoldenrod2', 'hotpink1',
                 'royalblue2', 'chocolate1', 'mediumorchid4', 'brown', 'aquamarine1',
                 'darkmagenta', 'darkolivegreen', 'deepskyblue4', 'dimgray', 'darksalmon',
                 'darkseagreen', 'darkslateblue', 'darkslategray1','cyan4', 'coral2')
    doses     <- c(1:length(crm$prior))
    priorprob <- crm$prior

    graphics::plot(x = doses, y = priorprob, col = 1, type = 'b', xlab = 'Dose', xaxt = 'n',
         ylab = 'Probability of Dose Limiting Toxicity', ylim = ylim, lwd = lwd)
    graphics::abline(h = crm$target, lty = 2)
    graphics::axis(1, at = doses, labels = dose_labels)

    j = 1
    k = 1
    legend_label = c('Prior Curve')
    legend_pch = c(NA)

    for(i in cohort_sizes) {

      legend_label = c(legend_label, paste0('Cohort ', j))
      legend_pch   = c(legend_pch, paste0(j))

      loopcohort_crm <- dose_func(prior = crm$prior, target = crm$target, tox = crm$tox[1:(k+i-1)],
                                  level = crm$level[1:(k+i-1)], ...)
      graphics::points(x = doses , y = loopcohort_crm$ptox, type = 'b', col = colours[1+j], pch = paste0(j),
             lwd = lwd)

      k = k + i
      j = j + 1
    }

    graphics::legend(x = max(priorprob) + 0.15, col = colours[1:(length(cohort_sizes)+1)],
           lty = rep(1, length(cohort_sizes)+1), lwd = lwd, pch = legend_pch,
           legend = legend_label, cex = cex)

    if(cohort.last){
      graphics::points(x = doses , y = crm$ptox, type = 'b', col = colours[j], pch = paste0(j - 1),
             lwd = 6)
    }

    if(!is.null(file)) {
      grDevices::dev.off()
    }
  }
}
