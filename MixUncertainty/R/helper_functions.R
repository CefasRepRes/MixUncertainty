# ---
# title: 'Define helper functions'
# author: 'Various'
# date: 'August 2023'
# ---

#' Crude classification of 'true' and 'false' zeros in historic data
#'
#' Function takes a time-series of metier effort-share data as a matrix and
#' crudely classified any zero values as 'true' or 'false' zeros. A true zero
#' implies that an observation of zero effort-share is likely given historic
#' non-zero trends. A false zero implies that an observation of zero effort-share
#' is unlikely and may be due to missing data in that year. False zeros are
#' converted to \code{NA}.
#'
#' @param dat A matrix of historic effortshare where rows = years and cols = metiers

findTrueZeros <- function(dat) {

  ## Check for presence of structural and false zeros
  if (any(dat==0)) {

    ## does zero element fall outsize 95% CI of non zero elements?
    datMissing <- sapply(1:nrow(t(dat)), function(y) {

      datSd <- sd(t(dat)[y,t(dat)[y,] > 0])
      datMu <- mean(t(dat)[y,t(dat)[y,] > 0])
      datn0 <- sum(t(dat)[y,] > 0)

      sapply(1:ncol(t(dat)), function(x){

        res <- FALSE
        if(t(dat)[y,x] == 0 & datn0 >=5) {
          if((t(dat)[y,x] > (datMu + 2*datSd)) | (t(dat)[y,x] < (datMu - 2*datSd)))
            res <- TRUE
        }
        return(res)
      }, simplify = "array")
    }, simplify = "array")

    ## Remove zeros corresponding to missing data
    dat[datMissing] <- NA
  }
  return(dat)
}

#' Simple function to draw random Dirichlet distributed draws
#'
#' Function takes a vector of alpha parameters and samples $n$ random draws
#' from a Dirichlet distribution.
#'
#' @param n Number of draws
#' @param alpha a numeric vector of alpha parameters.

## draw random samples from a Dirichlet distribution
rdirichlet <- function(n ,alpha) {
  x <- replicate(n, sapply(alpha, function(a) rgamma(n = 1, a, 1)))
  p <- sweep(x, 2, colSums(x), "/")
  return(p)
}

#' Multinomial inverse logit transformation
#'
#' Multinomial inverse logit transformation of a vector.
#'
#' @param x a numeric vector on a natural scale

## multinomial inverse logit function
invlogit <- function(x) {
  xsize <- length(x)
  expx  <- exp(x)
  sumx  <- sum(expx)

  p <- expx / (1 + sumx)
  q <- 1 - sum(p)
  p[xsize+1] <- q
  return(p)
}

#' Check optimisation of a fitted model
#'
#' Check optimisation. Return code.
#'
#' @param x a list

checkOpt <- function(x, verbose, makeLog) {

  ## check optimisation
  if(x$opt$convergence > 0) {
    if (verbose)
      cat(" no Dirichlet convergence |")
    if (makeLog)
      logs <- "Dir - no convergence"
    else
      logs <- NULL
  } else {
    if (verbose)
      cat(" Dirichlet success |")
    if (makeLog)
      logs <- "Dir - success"
    else
      logs <- NULL
  }
  return(logs)
}

#' Check for model fitting failure
#'
#' Check for failure. Return code.
#'
#' @param x a list
#' @param verbose (logical) should progress/summaries be printed? Defaults to
#'                \code{TRUE}

checkFail <- function(x, verbose) {

  if (is.null(x$opt$convergence)) {
    if (verbose)
      cat(" Dirichlet failed |")
    rerun  <- TRUE
    fail   <- TRUE
    pdHess <- FALSE
    nans   <- FALSE
    conv   <- FALSE

  } else if (!x$sdr$pdHess) {
    if (verbose)
      cat(" Hessian not positive-definite |")
    rerun  <- TRUE
    fail   <- FALSE
    pdHess <- FALSE

    xplsd   <- as.list(x$sdr,"Std")
    nans   <- any(is.nan(xplsd$rw))
    conv   <- x$opt$convergence == 0

  } else {

    rerun  <- FALSE
    fail   <- FALSE
    pdHess <- TRUE

    xplsd  <- as.list(x$sdr,"Std")
    nans   <- any(is.nan(xplsd$rw))
    conv   <- x$opt$convergence == 0
  }

  return(list(rerun  = rerun,
              fail   = fail,
              pdHess = pdHess,
              nans   = nans,
              conv   = conv))
}

#' Calculate AIC from an optimisation summary
#'
#' ...
#'
#' @param opt a list

simpleAIC <- function(opt) {

  npar <- length(opt[["par"]])
  nll  <- opt[["objective"]]
  return(npar * 2 + 2 * nll)
}

#' Small function to impute values to cases of zero and one in proportional data
#'
#' Function identifies years where a metier has zero effort.
#' In these cases, impute 1e-6 to the proportional
#' effort-share of each zero effort metier. Then deduct the
#' added effort-share from non-zero effort metiers.
#'
#' @param dat A data frame containing numeric data.

impute_cases <- function(dat) {

  dat2 <- sapply(1:nrow(dat),
                 function(x) {

                   ## Check if 1 occurs in row
                   if(any(dat[x,] == 0, na.rm = TRUE)) {

                     ## find column where zero occurs
                     datcol    <- which(dat[x,] == 0)
                     nnonzero  <- ncol(dat) - length(datcol)

                     ## impute small values
                     dat[x, datcol] <- 1e-6

                     ## distribute the cost evenly across remaining data
                     # dat[x, -datcol]  <- dat[x, -datcol] - (sum(dat[x, datcol])/nnonzero)

                     ## distribute the cost proportionally according the value of remaining data
                     dat[x, -datcol]  <- dat[x, -datcol] - (sum(dat[x, datcol]) * (dat[x, -datcol]/sum(dat[x, -datcol])))

                   }

                   return(dat[x,])
                 })

  ## edit output to match effort-share matrix orientation
  colnames(dat2) <- rownames(dat)
  return(t(dat2))
}

