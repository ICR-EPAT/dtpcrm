
.M.mat <- function(numrow, c, M, c.sz){
  for(i in 1:(c.sz + 1)){
    M[numrow + i - 1, c] <- as.character(i - 1)
  }
  return(M)
}

# .col.linearrow.func <- function(newdose, start.dose, arr.lcolour, ind, ind.col){
#   arr.lcolour[ind,ind.col][newdose > start.dose] <- "green"
#   arr.lcolour[ind,ind.col][newdose == start.dose] <- "orange"
#   arr.lcolour[ind,ind.col][newdose < start.dose | is.na(newdose)] <- "red"
#   return(arr.lcolour)
# }


#' @title Produce DTP flow diagram
#'
#' @description dtpflow will produce a flow diagram of the possible paths for
#'     the next three cohorts of subjects.
#'
#' @usage dtpflow(dtptable, cohort.labels = c('C1', 'C2', 'C3'))
#'
#' @param dtptable a dataframe produced by calculate_dtps where cohort_sizes was of length 3.
#' @param cohort.labels A vector of length 3, containing character strings for
#'     the cohort labels.
#'
#' @details The function will produce a visual flow diagram for the first three
#'     cohorts of the provided dataframe.
#'
#' @examples
#'
#' prior  <- c(0.1, 0.2, 0.5)
#' target <- 0.15
#' prev_tox <- c(0, 0, 0)
#' prev_dose <- c(2, 2, 2)
#' cohort_sizes <- c(2, 3, 3)
#'
#' next_dose = applied_crm(prior = prior, target = target,
#'                         tox = prev_tox, level = prev_dose)$mtd
#'
#' dose_func <- applied_crm
#'
#' DTP = calculate_dtps(next_dose, cohort_sizes, prev_tox = prev_tox,
#'                       prev_dose = prev_dose, dose_func = applied_crm,
#'                       prior = prior, target = target)
#'
#' dtpflow(dtptable = DTP, cohort.labels = c('C1', 'C2', 'C3'))
#'
#' @keywords CRM DTP dtpcrm
#'
#' @export
dtpflow <- function(dtptable, cohort.labels = c('C1', 'C2', 'C3')){

  ### Check dtptable for correct number of cohorts
  if(ncol(dtptable) != 7) {stop('dtptable is required to project doses for 3 cohorts, where length(cohort_sizes) = 3')}

  ### m: vector of cohort size for each cohort, e.g. m<-c(3,3,3) for cohort sizes of 3  ####

    m <- c(max(dtptable[ , 2], na.rm = T),
           max(dtptable[ , 4], na.rm = T),
           max(dtptable[ , 6], na.rm = T))

  ### Obtain unique doses for Cohort 2 and 3    #####

  start.dose <- dtptable[1, 1]

  ### Cohort 2 ####

  C2.dose <- unique(dtptable[ , 3][dtptable[ , 2] == 0])

  for(i in 1:m[1]){
    C2.dose <- c(C2.dose, unique(dtptable[ , 3][dtptable[ , 2] == i]))
    if(length(C2.dose) > (m[1] + 1))
      {stop('There are inconsistent recommendations for the same pathways, unable to produce flow chart')}
  }


  ### Cohort 3 ####

  ind <- 0
  j <- 1
  C3.alldose <- rep(NA, length = length(C2.dose) * (m[2] + 1))

  for(y in 1:(m[1] + 1)){
    ind <- (ind + 1):((m[2] + 1) * (m[3] + 1) * y)
    C3.dose <- unique(dtptable[ind, 5][dtptable[ind, 4] == 0])
    for (i in 1:m[2]){
      C3.dose <- c(C3.dose, unique(dtptable[ind, 5][dtptable[ind, 4] == i]))
      if(length(C3.dose) > (m[2] + 1))
      {stop('There are inconsistent recommendations for the same pathways, unable to produce flow chart')}
    }
    C3.alldose[j:((m[2] + 1) * y)] <- C3.dose
    j <- j + (m[2] + 1)
    ind <- max(ind)

  }

  final.dose <- c(start.dose, C2.dose, C3.alldose)

  names <- paste("d",  "(", c(start.dose, C2.dose, C3.alldose), ")", sep = "")
  length(names)

  M <- matrix(nrow = length(final.dose),
              ncol = length(final.dose),
              byrow = TRUE, data = 99)
  col.line <- matrix(nrow = length(final.dose),
                     ncol = length(final.dose),
                     byrow = TRUE, data = "black")


  ### C1 to C2  ------------------
  ind.row <- 2
  ind.col = 1

  M <- .M.mat(ind.row, c = ind.col, M = M, c.sz = m[1])

  # col.line <- .col.linearrow.func(newdose = C2.dose, start.dose = start.dose,
  #                                 arr.lcolour = 'black',
  #                                 ind = (ind.row):(ind.row + m[1]),
  #                                 ind.col = ind.col)


  ### C2 to C3 ( Number of different sets = length(C2.dose) )  -----------------------
  rn <- m[1]
  for(ind.col in 2:(m[1] + 2)){
    M <- .M.mat(ind.row + (rn + 1), c = ind.col, M = M, c.sz = m[2])

    # col.line <- .col.linearrow.func(newdose = C3.alldose[((ind.col - 2) * ((m[2] + 1) + 1)):((ind.col - 1) * (m[2] + 1))],
    #                                 start.dose = C2.dose[ind.col - 1],
    #                                 arr.lcolour = 'black',
    #                                 ind = (ind.row + (rn + 1)):(ind.row + (rn + 1) + m[2]),
    #                                 ind.col = ind.col)

    ind.row <- ind.row + (rn + 1)
    rn <- m[2]

  }

  box.colour <- rep("white", length(final.dose))
  box.colour[is.na(final.dose)] <- "red"
  names[which(names == "d(NA)")] <- "STOP"

  ## PLOT DTP FLOW DIAGRAM --------------

  graphics::par(mfrow=c(1, 1))

  diagram::plotmat(M, pos = c(1, length(C2.dose), length(C3.alldose)),
          curve = 0, name = names, lwd = 2,absent=99,
          cex=1,dtext=0.08,box.lwd = 1.5, cex.txt = 0.9,
          box.type = "square", box.prop = 0.49,
          box.cex=c(1, rep(1,length(C2.dose)),
                    rep(0.75, length(C3.alldose))),
          box.size=c(0.08, rep(0.06,length(C2.dose)),
                     rep(0.029, length(C3.alldose))),
          box.col=box.colour,
          #arr.lcol=col.line,
          #arr.col=col.line,
          shadow.size=0)

  ## Adding Cohort Labels --------------
  graphics::text(x = 0.02, y = 0.85, cohort.labels[1], cex = 0.8)
  graphics::text(x = 0.02, y = 0.5, cohort.labels[2], cex = 0.8)
  graphics::text(x = 0.02, y = 0.25, cohort.labels[3], cex = 0.8)

}




