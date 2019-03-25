#' @export
applied_titecrm <- function(prior, target, tox, level, followup, obswin,
                        # The above signature should be mimicked by specialisations
                        no_skip_esc = TRUE, no_skip_deesc = TRUE,
                        # coherent_esc = FALSE, coherent_deesc = FALSE,
                        global_coherent_esc = TRUE,
                        # TODO: I don't imagine anyone wants global_coherent_deesc
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
