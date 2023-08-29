# ---
# title: 'Forecast effortshare from a Dirichlet distributed noise and multivariate log-normal random walk'
# author: 'Various'
# date: 'November 2022'
# ---

#' Forecast catchability from a state-space time-series model with Dirichlet distributed observation noise
#'
#' Function takes a time-series of metier effort-share data as a matrix and
#' fits a state-space model with Dirichlet distributed observation error. The
#' specific model structure will depend on the observation data. The final data
#' year is simulated to generate a distribution based on observation uncertainty.
#' A forecast is carried out for each simulated catchability value that propagates
#' uncertainty estimated from the random walk.
#'
#' @param dat A matrix of historic effortshare where rows = years and cols = metiers
#' @param nit (Integer) the number of replicate parameter trajectories to sample
#' @param fillyear Vector of intermediate and forecast years to be sampled for.
#' @param verbose (Logical) Should the function print progress?
#' @param makeLog (Logical) Should the function generate a log of the methods
#'                applied to each variable (and their success)?
#' @param makePlots (Logical) Should the function generate plots showing model
#'                  fit to data for each stock?
#'
#' @return An matrix of sampled parameter values
#'
#' @export

TMB_DirMissingAR1Hurdle <- function(dat,
                                    nit,
                                    fillyear,
                                    verbose,
                                    makeLog,
                                    makePlots) {

  require(TMB)

  # -------------------#
  # Filter excess data
  # -------------------#

  ## remove rows (years) before first data point
  dat <- dat[min(which(rowSums(dat) > 0)):nrow(dat),,drop = FALSE]

  ## Prepare random walk matrix here
  rw <- t(dat)*0
  rw <- rw[-nrow(rw),,drop = FALSE]

  # ---------------------------------#
  # Identify missing values
  # ---------------------------------#

  ## Check for presence of structural and false zeros
  dat <- findTrueZeros(dat)

  # --------------------------------------------#
  # Build TMB objective functions and fit models
  # --------------------------------------------#

  hurdle <- any(dat==0, na.rm = TRUE)

  ## Fit hurdle if true zeros are present
  if (hurdle) {
    out <- fitMVNhurdle(dat, rw)
  } else {
    out <- fitMVNDirrw(dat, rw)
  }

  # --------------------------------------------#
  # Check optimisation
  # --------------------------------------------#
  #
  # Re-run optimisation with a Gaussian distributed intervals if hessian is not
  # positive-definite

  ## Check identifiability of parameters
  checkout <- checkFail(out, verbose)
  rerun  <- checkout$rerun
  fail   <- checkout$fail
  pdHess <- checkout$pdHess
  nans   <- checkout$nans
  conv   <- checkout$conv

  if(rerun) {
    if (verbose)
      cat(" RERUN |")

    if (hurdle) {
      out <- fitNhurdle(dat, rw)
    } else {
      out <- fitNDirrw(dat, rw)
    }

    ## Check identifiability of parameters
    checkout <- checkFail(out, verbose)
    rerun  <- checkout$rerun

  }

  ## check optimisation
  if(!rerun & (verbose | makeLog))
    logs <- checkOpt(out, verbose, makeLog)

  # If random walk hessian is not positive definite with Gaussian distributed
  # intervals on the latent process, check if key parameters are identifiable in
  # either MVN or univariate normal models.

  if(rerun & !hurdle) {

    ## if both multivariate and univariate methods cannot generate forecasts
    ## then fail!
    if ((fail | nans) & (checkout$fail | checkout$nans)) {
      if (makeLog) {
        logs <- "Dir - failed"
      } else {
        logs <- NULL
      }
      return(list(res  = NULL,
                  logs = logs,
                  plots = NULL))
    }

    ## otherwise, choose the better of the two methods
    if(which.max(c(sum(pdHess, conv), sum(checkout$pdHess, checkout$conv))) == 1) {
      out <- fitMVNDirrw(dat, rw)
    }

    ## check optimisation
    if(verbose | makeLog)
      logs <- checkOpt(out, verbose, makeLog)
  }

  # If hurdle hessian is not positive definite with Gaussian distributed intervals
  # on the latent process, use fixed values for AR1 process OR Bernoulli
  # observation. Select model with lowest AIC

  ## Check identifiability of parameters
  if(rerun & hurdle) {
    if (verbose)
      cat(" RERUN |")

    out <- list(out1 = fitNhurdle_fixARrho(dat, rw, logitARrho = 1.75),
                out2 = fitNhurdle_fixb0b1(dat, rw, logb0 = 2.3, logb1 = 2))

    checklist <- c(checkFail(out$out1, FALSE)$rerun, checkFail(out$out2, FALSE)$rerun)

    ## If both fail, stop function here
    if (sum(checklist) == 2) {
      if (verbose)
        cat(" Dirichlet failed |")

      if (makeLog) {
        logs <- "Dir - failed"
      } else {
        logs <- NULL
      }
      return(list(res  = NULL,
                  logs = logs,
                  plots = NULL))
    }

    ## If only one succeeds, select this one
    if (sum(checklist) == 1) {
      out <- out[[which(!checklist)]]

      ## check optimisation
      if(verbose | makeLog)
        logs <- checkOpt(out, verbose, makeLog)
    }

    # If both positive-definite hessian, choose model with lowest AIC
    if (sum(checklist) == 0) {

      outIdx<- which.min(c(simpleAIC(out$out1$opt),
                           simpleAIC(out$out2$opt)))
      out <- out[[outIdx]]

      ## check optimisation
      if(verbose | makeLog)
        logs <- checkOpt(out, verbose, makeLog)
    }
  }

  # ---------------------------------#
  # Extract fitted parameters
  # ---------------------------------#

  pl     <- as.list(out$sdr,"Est")
  plsd   <- as.list(out$sdr,"Std")
  MVNcov <- out$obj$report()$MVNcov

  ## Add information to logs if NaN occur in parameter estimate variances
  if(any(is.nan(unlist(plsd)))) {
    logs <- sapply(logs, function(x) paste0(x, " - NaN"))
  }

  # ---------------------------------#
  # Forecast from model
  # ---------------------------------#

  ## Define number of forecast years
  nyear_forecast <- length(fillyear)

  if (hurdle) {
    resvariates <- forecast_from_hurdle(dat, pl, plsd, MVNcov, fillyear, nit)
  } else {
    resvariates <- forecast_from_rw(dat, pl, plsd, MVNcov, fillyear, nit)
  }

  ## Check if forecast has failed
  if(any(is.nan(resvariates))) {
    if (verbose)
      cat(" Forecast failed |")

    if (makeLog)
      logs <- "Dir - failed"

    return(list(res  = NULL,
                logs = logs,
                plots = NULL))
  }

  ## make diagnostic plots to explore quality of fit and forecast
  if (makePlots) {

    pred_quantiles <- array(NA, dim = c(nyear_forecast, ncol(dat), 3))
    for(i in 1:dim(pred_quantiles)[1]) {
      for(j in 1:dim(pred_quantiles)[2]) {
        pred_quantiles[i,j,] <- quantile(resvariates[i,j,], na.rm = TRUE,
                                         prob = c(0.025,0.50,0.975))
      }
    }

    p1 <- tryCatch(plot_forecast_Dir(list(dat = t(dat)),
                                     pl,
                                     plsd,
                                     pred_quantiles,
                                     fillyear,
                                     invlogitfun = invlogit),
                   error = function(e)e)

  } else {
    p1 <- NULL
  }

  return(list(res  = resvariates,
              logs = logs,
              plots = p1))
}

