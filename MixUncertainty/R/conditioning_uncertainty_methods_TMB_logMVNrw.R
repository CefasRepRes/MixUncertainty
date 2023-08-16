# ---
# title: 'Forecast catchability from a multivariate log-normal random walk'
# author: 'Various'
# date: 'November 2022'
# ---

#' Forecast catchability from fitted multivariate log-normal random walk
#'
#' Function takes a time-series of stock catchability data as a matrix and
#' fits a multivariate normal random walk to the time-series. The final data
#' year is simulated to generate a distribution based on observation uncertainty.
#' A forecast is carried out for each simulated catchability value that propagates
#' uncertainty estimated from the random walk.
#'
#' A useful feature of this method is that missing values are easily accommodated.
#'
#' @param qs A matrix of historic catchability where rows = years and cols = stocks
#' @param meter_mt An object of class \code{FLMetier}, \code{FLMetier},
#'                 \code{FLMetierExt} or a \code{FLMetierExt}.
#' @param fillyear (Integer) Vector of years to be filled with forecasted catchability
#' @param deterministic (Logical) Should the results for first iteration be
#'                      a simple mean of the historical period? Defaults
#'                      to \code{TRUE}
#' @param verbose (Logical) Should the function print progress? Defaults
#'                to \code{TRUE}
#' @param makeLog (Logical) Should the function generate a log of the methods
#'                applied to each stock (and their success)? Defaults
#'                to \code{TRUE}
#' @param makePlots (Logical) Should the function generate a figure of model fit?
#'                  Defaults to \code{TRUE}
#'
#' @return An updated FLR metier object.

TMB_logMVNrw <- function(qs,
                         metier_mt,
                         fillyear,
                         deterministic,
                         detMethod,
                         verbose = TRUE,
                         makeLog = TRUE,
                         makePlots = TRUE) {

  # ==========================================================#
  # Prepare objects and arguments
  # ==========================================================#

  ## Find number of iterations in object
  nit <- ifelse(deterministic, dims(metier_mt)$iter-1, dims(metier_mt)$iter)

  ## prepare log file
  if(makeLog) {
    logs <- matrix(NA, nrow = 1, ncol = ncol(qs))
    colnames(logs) <- colnames(qs)
  }
  if(makePlots) {
    plots <- vector(mode = "list", ncol(qs))
    names(plots) <- colnames(qs)
  }

  # ==========================================================#
  # Prepare Data
  # ==========================================================#

  # ----------------------------------------------------------#
  # Filter excess data
  # ----------------------------------------------------------#

  ## extract vector of years
  qs_years <- rownames(qs)

  ## convert NaN, Inf and NA into zeros
  qs[is.na(qs)|is.nan(qs)] <- 0
  qs[is.infinite(qs)]      <- 0

  ## exit function if all data are zero
  if (all(qs == 0)) {
    if (verbose) {
      cat(" no data available |")
    }
    if (makeLog) {
      logs[] <- "no data"
      return(list(res   = metier_mt,
                  logs  = logs,
                  plots = NULL))
    } else {
      return(list(res   = metier_mt,
                  logs  = NULL,
                  plots = NULL))
    }
  }

  ## remove rows (years) before first data point
  qs <- qs[min(which(rowSums(qs) > 0)):nrow(qs),,drop = FALSE]

  ## remove cols (stocks) with zero catchability
  qs <- qs[,colSums(qs) > 0, drop = FALSE]

  # ==========================================================#
  # Fit models
  # ==========================================================#

  # ---------------------------------------------------------#
  # Apply appropriate analysis
  # ---------------------------------------------------------#

  # In this section, we apply a univariate normal random walk if only one stock
  # is caught. If more than one stock is caught, a multivariate normal random
  # walk is applied to all stocks with >50% coverage of the time-series. A
  # univariate normal random walk is applied to each stock with <50% coverage
  # of the time-series. Stocks with < 3 data points are skipped.

  ## univariate normal random walk if only 1 stock caught
  if (ncol(qs) < 2) {
    out <- fit_N_AR1(qs,
                     metier_mt,
                     qs_years,
                     fillyear,
                     deterministic,
                     detMethod,
                     verbose,
                     makeLog,
                     makePlots)

    ## Check if model fitted well
    if (!is.null(out$logs)) {
      if (verbose) {
        cat(paste0(" ",out$logs," |"))
      }
      fail <- TRUE
    } else {

      ## check for model fitting issues
      checkout <- checkFail(out, verbose = FALSE)
      fail <- checkout$fail

      if(verbose) {
        if (fail) {
          cat(" Model failed |")
        }
      }

      ## make log
      if (makeLog) {
        logs <- checkOpt(out, FALSE, makeLog)
      }
    }

    ## proceed with forecast if model has not failed
    if (!fail) {
      ## extract fitted parameters
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

      ## number of years to forecast
      nyear_forecast <- length(fillyear)

      ## generate predictions
      resvariates <- tryCatch(forecast_from_AR1_logN(qs,
                                                     pl,
                                                     plsd,
                                                     MVNcov,
                                                     fillyear,
                                                     nit), error = function(e) e)

      # ---------------------------------#
      # Insert Forecast into FLR objects
      # ---------------------------------#

      ## if successful, insert forecast into FLR object
      if (class(resvariates)[1] == "array") {

        if (verbose){
          cat(" N model success |")
        }

        metier_mt <- insert_forecast(resvariates,
                                     qs,
                                     qs_years,
                                     metier_mt,
                                     deterministic,
                                     detMethod,
                                     nit)

        # ---------------------------------#
        # (Optional) make diagnostic plots
        # ---------------------------------#

        ## make diagnostic plots to explore quality of fit and forecast
        if (makePlots & !any(is.nan(plsd$rw))) {

          pred_quantiles <- array(NA, dim = c(nrow(pl$rw), nyear_forecast, 3))
          for(x in 1:dim(pred_quantiles)[1]) {
            for(y in 1:dim(pred_quantiles)[2]) {
              pred_quantiles[x,y,] <- quantile(resvariates[x,y,], na.rm = TRUE,
                                               prob = c(0.025,0.50,0.975))
            }
          }

          plots[[colnames(qs)]] <- plot_forecast_MVN(list(dat = t(qs)),
                                                     pl,
                                                     plsd,
                                                     pred_quantiles,
                                                     fillyear)

        } # END if makePlots
      } else {
        if(verbose) {
          cat(" Model failed |")
        }
        if (makeLog) {
          logs[,i] <- "Model failed"
        }
      } # END if forecast successful
    } # END if model successful
  }

  ## more nuanced approach if more than 1 stock caught
  if (ncol(qs) >=2) {

    ## proportion of observation data per stock
    pObs <- colSums(qs > 0)/nrow(qs) >= .7

    ## vector of stocks for multivariate normal
    stocks_mvn <- colnames(qs)[pObs]

    ## vector of stocks for univariate normal
    stocks_norm <- colnames(qs)[!pObs]

    ## if <2 stocks available for multivariate normal, just use univariate for all
    if (length(stocks_mvn) < 2) {
      stocks_norm <- c(stocks_norm, stocks_mvn)
      stocks_mvn  <- vector("numeric", length = 0)
    }

    ## fit and forecast multivariate normal random walk
    if (length(stocks_mvn) > 0) {
      out <- fitMVN_AR1(qs[,stocks_mvn, drop = FALSE],
                        metier_mt,
                        qs_years,
                        fillyear,
                        deterministic,
                        detMethod,
                        verbose,
                        makeLog,
                        makePlots)

      ## Check if model fitted well
      if (!is.null(out$logs)) {

        ## insufficient data for MVN model
        rerun <- TRUE
      } else {

        ## check for model fitting issues
        checkout <- checkFail(out, verbose = FALSE)
        rerun <- checkout$rerun

        if(verbose) {
          if (checkout$fail) {
            cat(" Model failed |")
          }
        }

        ## make log
        if (makeLog) {
          logs[,stocks_mvn] <- sapply(stocks_mvn, function(x) checkOpt(out, FALSE, makeLog), USE.NAMES = TRUE)
        }
      }

      ## If MVN model fitting fails, use univariate normal model instead
      ## Otherwise, forecast from the model
      if (rerun) {

        if (verbose)
          cat(" RERUN |")

        ## re-make stock name vector for univariate model
        stocks_norm <- c(stocks_norm, stocks_mvn)
        stocks_mvn  <- vector("numeric", length = 0)

      } else {

        ## extract fitted parameters
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

        ## number of years to forecast
        nyear_forecast <- length(fillyear)

        ## generate predictions
        resvariates <- tryCatch({forecast_from_AR1_logN(qs[,stocks_mvn, drop = FALSE],
                                                        pl,
                                                        plsd,
                                                        MVNcov,
                                                        fillyear,
                                                        nit)},
                                error = function(e)e)

        # ---------------------------------#
        # Insert Forecast into FLR objects
        # ---------------------------------#

        ## if successful, insert forecast into FLR object
        if (class(resvariates)[1] == "array") {

          if (verbose) {
            cat(" MVN model success |")
          }

          ## if successful, insert forecast into FLR object
          metier_mt <- insert_forecast(resvariates,
                                       qs[,stocks_mvn, drop = FALSE],
                                       qs_years,
                                       metier_mt,
                                       deterministic,
                                       detMethod,
                                       nit)

          # ---------------------------------#
          # (Optional) make diagnostic plots
          # ---------------------------------#

          ## make diagnostic plots to explore quality of fit and forecast
          if (makePlots & !any(is.nan(plsd$rw))) {

            pred_quantiles <- array(NA, dim = c(nrow(pl$rw), nyear_forecast, 3))
            for(x in 1:dim(pred_quantiles)[1]) {
              for(y in 1:dim(pred_quantiles)[2]) {
                pred_quantiles[x,y,] <- quantile(resvariates[x,y,], na.rm = TRUE,
                                                 prob = c(0.025,0.50,0.975))
              }
            }

            ## Generate plot
            plots[[stocks_mvn[1]]] <- plot_forecast_MVN(list(dat = t(qs[,stocks_mvn, drop = FALSE])),
                                                             pl,
                                                             plsd,
                                                             pred_quantiles,
                                                             fillyear)
          } # END makePlots
        } else {
          if(verbose) {
            cat(" Model failed |")
          }
          if (makeLog) {
            logs[,stocks_mvn] <- "Model failed"
          }
        } # END if else forecast success
      } # END if else rerun
    } # END if fit MVN stocks

    ## fit and forecast univariate normal random walk
    if (length(stocks_norm) > 0) {

      ## loop over each stock
      for(i in stocks_norm) {
        out <- fit_N_AR1(qs[,i,drop = FALSE],
                         metier_mt,
                         qs_years,
                         fillyear,
                         deterministic,
                         detMethod,
                         verbose,
                         makeLog,
                         makePlots)

        ## Check if model fitted well
        if (!is.null(out$logs)) {
          if (verbose) {
            cat(paste0(" ",out$logs," |"))
          }
          fail <- TRUE
        } else {

          ## check for model fitting issues
          checkout <- checkFail(out, verbose = FALSE)
          fail <- checkout$fail

          if(verbose) {
            if (fail) {
              cat(" Model failed |")
            }
          }

          ## make log
          if (makeLog) {
            logs[,i] <- checkOpt(out, FALSE, makeLog)
          }
        }

        ## proceed with forecast if model has not failed
        if (!fail) {
          ## extract fitted parameters
          pl     <- as.list(out$sdr,"Est")
          plsd   <- as.list(out$sdr,"Std")
          MVNcov <- out$obj$report()$MVNcov

          ## Add information to logs if NaN occur in parameter estimate variances
          if(any(is.nan(unlist(plsd)))) {
            logs[,i] <- sapply(logs[,i], function(x) paste0(x, " - NaN"))
          }

          # ---------------------------------#
          # Forecast from model
          # ---------------------------------#

          ## number of years to forecast
          nyear_forecast <- length(fillyear)

          ## generate predictions
          resvariates <- tryCatch(forecast_from_AR1_logN(qs[,i, drop = FALSE],
                                                         pl,
                                                         plsd,
                                                         MVNcov,
                                                         fillyear,
                                                         nit), error = function(e) e)

          # ---------------------------------#
          # Insert Forecast into FLR objects
          # ---------------------------------#

          ## Forecast fails if contains NaNs
          forecastsuccess <- ifelse(class(resvariates)[1] == "array",
                                    ifelse(any(is.nan(resvariates)),
                                           FALSE,
                                           TRUE),
                                    FALSE)

          ## if successful, insert forecast into FLR object
          if (forecastsuccess) {

            if (verbose) {
              cat(" N model success |")
            }

            metier_mt <- insert_forecast(resvariates,
                                         qs[,i, drop = FALSE],
                                         qs_years,
                                         metier_mt,
                                         deterministic,
                                         detMethod,
                                         nit)

            # ---------------------------------#
            # (Optional) make diagnostic plots
            # ---------------------------------#

            ## make diagnostic plots to explore quality of fit and forecast
            if (makePlots & !any(is.nan(plsd$rw))) {

              pred_quantiles <- array(NA, dim = c(nrow(pl$rw), nyear_forecast, 3))
              for(x in 1:dim(pred_quantiles)[1]) {
                for(y in 1:dim(pred_quantiles)[2]) {
                  pred_quantiles[x,y,] <- quantile(resvariates[x,y,], na.rm = TRUE,
                                                   prob = c(0.025,0.50,0.975))
                }
              }

              plots[[i]] <- plot_forecast_MVN(list(dat = t(qs[,i, drop = FALSE])),
                                              pl,
                                              plsd,
                                              pred_quantiles,
                                              fillyear)

            } # END if makePlots
          } else {
            if(verbose) {
              cat(" Model failed |")
            }
            if (makeLog) {
              logs[,i] <- "Model failed"
            }
          } # END if else forecast successful
        } # END if model successful
      } # END loop over stocks
    } # END if >0 valid stock
  } # END if >1 stock caught

  return(list(res   = metier_mt,
              logs  = logs,
              plots = plots))
}
