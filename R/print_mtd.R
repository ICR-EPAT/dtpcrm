#' @export
print.mtd <- function (x, dgt = 3, model.detail = x$model.detail, patient.detail = x$patient.detail,
          ...){
  cat("Today: ", date(), "\n")
  n <- length(x$pid)
  used <- rep(0, n)
  used[x$include] <- 1
  ptox <- round(x$ptox, digits = dgt)
  ptoxL <- round(x$ptoxL, digits = dgt)
  ptoxU <- round(x$ptoxU, digits = dgt)
  if (x$tite) {
    wts <- round(x$weights, digits = dgt)
    if (is.null(x$followup)) {
      followup <- rep("N/A", n)
    }
    else {
      followup <- round(x$followup, digits = dgt)
    }
    if (patient.detail) {
      cat("DATA SUMMARY (TITE-CRM) \n")
      cat("PID", "\t", "Level", "\t", "Toxicity", "\t",
          "f/u", "\t", "Weight", "\t", "Included", "\n")
      for (i in 1:n) cat(x$pid[i], "\t", x$level[i], "\t",
                         x$tox[i], "\t\t", followup[i], "\t", wts[i],
                         "\t\t", used[i], "\n")
    }
  }
  else {
    if (patient.detail) {
      cat("DATA SUMMARY (CRM) \n")
      cat("PID", "\t", "Level", "\t", "Toxicity", "\t",
          "Included", "\n")
      for (i in 1:n) {
        cat(x$pid[i], "\t", x$level[i], "\t", x$tox[i],
            "\t\t", used[i], "\n")
      }
    }
  }
  cat("\nToxicity probability update (with", x$conf.level *
        100, "percent probability interval):", "\n")
  if (is.null(x$dosename)) {
    cat("Level", "\t", "Prior", "\t", "n", "\t", "total.wts",
        "\t", "total.tox", "\t", "Ptox", "\t", "LoLmt", "\t",
        "UpLmt", "\n")
    K <- length(x$prior)
    for (k in 1:K) {
      expt <- which(x$level == k & used == 1)
      if (!x$tite) {
        totwts <- length(expt)
      }
      else {
        totwts <- round(sum(x$weights[expt]), digits = dgt)
      }
      cat(k, "\t", x$prior[k], "\t", length(expt), "\t",
          totwts, "\t\t", sum(x$tox[expt]), "\t\t", ptox[k],
          "\t", ptoxL[k], "\t", ptoxU[k], "\n")
    }
    if("stop" %in% x){         
              if(x$stop){
                        cat("It is recommended to stop the trial\n")
              } else{cat("Next recommended dose level:", x$mtd, "\n")}
    } else{cat("Next recommended dose level:", x$mtd, "\n")}
  }
  else {
    cat("Dose", "\t\t", "Level", "\t", "Prior", "\t", "n",
        "\t", "total.wts", "\t", "total.tox", "\t", "Ptox",
        "\t", "LoLmt", "\t", "UpLmt", "\n")
    K <- length(x$prior)
    for (k in 1:K) {
      expt <- which(x$level == k & used == 1)
      if (!x$tite) {
        totwts <- length(expt)
      }
      else {
        totwts <- round(sum(x$weights[expt]), digits = dgt)
      }
      cat(x$dosename[k], "\t", k, "\t", x$prior[k], "\t",
          length(expt), "\t", totwts, "\t\t", sum(x$tox[expt]),
          "\t\t", ptox[k], "\t", ptoxL[k], "\t", ptoxU[k],
          "\n")
    }
    if("stop" %in% x){         
              if(x$stop){
                        cat("It is recommended to stop the trial\n")
              } else{cat("Next recommended dose level:", x$mtd, "\n")}
    } else{cat("Next recommended dose level:", x$mtd, "\n")}
  }

  cat("Recommendation is based on a target toxicity probability of",
      x$target, "\n")
  if (model.detail) {
    cat("\nEstimation details:\n")
    if (x$model == "empiric") {
      cat("Empiric dose-toxicity model: p = dose^{exp(beta)}\n")
    }
    else {
      cat("Logistic dose-toxicity model: p = {1 + exp(-a-exp(beta)*dose)}^{-1} with a =",
          x$intcpt, "\n")
    }
    cat("dose =", round(x$dosescaled, digits = dgt), "\n")
    if (x$method == "bayes") {
      cat("Normal prior on beta with mean 0 and variance",
          x$prior.var, "\n")
      cat("Posterior mean of beta:", round(x$estimate,
                                           digits = dgt), "\n")
      cat("Posterior variance of beta:", round(x$post.var,
                                               digits = dgt), "\n")
    }
    else if (x$method == "mle") {
      cat("MLE of beta:", round(x$estimate, digits = dgt),
          "\n")
      cat("Approximate variance of beta:", round(x$post.var,
                                                 digits = dgt), "\n")
    }
  }
}
