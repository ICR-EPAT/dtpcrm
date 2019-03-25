
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' @export
applied_crm <- function(prior, target, tox, level,
                        # The above signature should be mimicked by specialisations
                        no_skip_esc = TRUE, no_skip_deesc = TRUE,
                        # coherent_esc = FALSE, coherent_deesc = FALSE,
                        global_coherent_esc = TRUE,
                        # TODO: I don't imagine anyone wants global_coherent_deesc
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

  # # Coherence, or local coherence, in escalation means not escalating
  # # immediately after a toxicity event
  # if (coherent_esc & tail(tox, 1) == 1) {
  #   x$mtd <- min(x$mtd, tail(level, 1))
  # }

  # # Coherence, or local coherence, in de-escalation means not de-escalating
  # # immediately after a non-toxicity event
  # if (coherent_deesc & tail(tox, 1) == 0) {
  #   x$mtd <- max(x$mtd, tail(level, 1))
  # }

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
  # if(!is.null(stop_func)) {
  #   stop_decision = stop_func(x)
  #   if(is.list(stop_decision)) {
  #     x$stop <- stop_decision[[1]]
  #     x$stop_reason <- stop_decision[[2]]
  #   } else {
  #     x$stop <- stop_decision
  #   }
  # } else {
  #   x$stop <- FALSE
  #   x$stop_reason <- NULL
  # }
  if(!is.null(stop_func)) {
    x = stop_func(x)  # Let stopping delegate decorate x
  }

  return(x)
}

# provides summary dataframe of applied_crm dose level information - useful for reports ----

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


# Produces a plot of current dose-toxicity curve against prior for applied_crm obj
# and outputs png of plot to current directory if 'file' is given. Potential for
# histroy of curves by cohort if cohort.sizes is provided - dose_func is required to do this.
# height and width are for pixels of output plot.
# current colour limit for number of cohorts is 10 couldn't be bothered to add more....
# dose_labels is required as people have a certain habit of not naming doses 1:n

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

    graphics::plot(x = doses, y = postprob, col = 2, type = 'b', xlab = 'Dose Level', xaxt = 'n',
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

    colours <- c(1, 2, 'chartreuse4', 'darkgoldenrod2', 'hotpink1', 'royalblue2', 'chocolate1',
                 'mediumorchid4', 'brown', 'aquamarine1')
    doses     <- c(1:length(crm$prior))
    priorprob <- crm$prior

    graphics::plot(x = doses, y = priorprob, col = 1, type = 'b', xlab = 'Dose Level', xaxt = 'n',
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


