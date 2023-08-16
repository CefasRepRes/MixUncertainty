# ---
# title: 'Define wrapper functions to fit TMB log-normal models'
# author: 'Various'
# date: 'August 2023'
# ---

#' Wrapper function to fit latent AR1 with with multivariate normal intervals
#'
#' A wrapper function to construct data and parameter inputs to build a TMB
#' objective function.
#'
#' @param qs A matrix of time-series observations for catchability data where
#'            cols = years and rows = stocks
#' @param metier_mt description
#' @param qs_years description
#' @param fillyear description
#' @param deterministic description
#' @param detMethod description
#' @param verbose description
#' @param makeLog description
#' @param makePlots description
#'
#' @return A list containing:
#'
#'         - 'res' = ...
#'         - 'logs' = ...
#'         - 'plots' = ...
#'
#' @export

fitMVN_AR1 <- function(qs,
                       metier_mt,
                       qs_years,
                       fillyear,
                       deterministic,
                       detMethod,
                       verbose,
                       makeLog,
                       makePlots) {
  
  ## remove rows (years) before first data point
  qs <- qs[min(which(rowSums(qs) > 0)):nrow(qs),,drop = FALSE]
  
  # ---------------------------------#
  # Build data and parameter objects
  # ---------------------------------#
  
  ## Prepare random walk matrix here
  rw <- t(qs)*0
  
  ## convert zeros back into NA for log-transformation
  qs[qs == 0] <- NA
  
  ## translate catchability to a log-scale
  logqs <- log(qs)
  
  ## Check availability of catchability data
  checkout <- checkCatchability(logqs, qs_years, verbose, makeLog)
  
  if(checkout$run == FALSE) {
    return(list(logs = checkout$logs))
  }
  
  ## prepare objects for TMB
  data <- list (code = "A",
                dat = t(logqs))
  
  ## build parameters object
  parameters <- list(rw = rw,
                     logSdMVN = rep(0,nrow(data$dat)),
                     logitRho = 0,
                     mu         = rep(0,nrow(data$dat)),
                     logitARrho = rep(0,nrow(data$dat)),
                     logSdObs = rep(0,nrow(data$dat)),
                     logTau = 0,
                     logb0  = 0,
                     logb1  = 0,
                     dat_missing = rep(0, sum(is.na(data$dat))))
  
  map1 <- list(logTau = as.factor(NA), 
               logb0  = as.factor(NA), 
               logb1  = as.factor(NA))
  
  # ---------------------------------#
  # Build model and optimise
  # ---------------------------------#
  
  obj <- tryCatch({MakeADFun(data, 
                             parameters, 
                             random = "rw",
                             DLL = "MixUncertainty", 
                             map = map1,
                             silent = TRUE)},
                  error = function(e)e)
  
  opt <- tryCatch({nlminb(obj$par, obj$fn, obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)
  
  sdr <- tryCatch({sdreport(obj)},
                  error = function(e)e)
  
  return(list(sdr = sdr,
              obj = obj,
              opt = opt))
}

#' Wrapper function to fit latent AR1 with univariate normal intervals
#'
#' A wrapper function to construct data and parameter inputs to build a TMB
#' objective function.
#'
#' @param qs A matrix of time-series observations for catchability data where
#'            cols = years and rows = stocks
#' @param metier_mt description
#' @param qs_years description
#' @param fillyear description
#' @param deterministic description
#' @param detMethod description
#' @param verbose description
#' @param makeLog description
#' @param makePlots description
#'
#' @return A list containing:
#'
#'         - 'sdr': a TMB report summary for the optimised model parameters,
#'         - 'opt': a summary returned by the optimisation function,
#'         - 'obj': the optimised TMB objective function
#'
#' @export

fit_N_AR1 <- Norm_option <- function(qs,
                                     metier_mt,
                                     qs_years,
                                     fillyear,
                                     deterministic,
                                     detMethod,
                                     verbose,
                                     makeLog,
                                     makePlots) {
  
  ## remove rows (years) before first data point
  qs <- qs[min(which(rowSums(qs) > 0)):nrow(qs),,drop = FALSE]
  
  # ---------------------------------#
  # Build data and parameter objects
  # ---------------------------------#
  
  ## Prepare random walk matrix here
  rw <- t(qs)*0
  
  ## convert zeros back into NA for log-transformation
  qs[qs == 0] <- NA
  
  ## translate catchability to a log-scale
  logqs <- log(qs)
  
  ## Check availability of catchability data
  checkout <- checkCatchability(logqs, qs_years, verbose, makeLog)
  
  if(checkout$run == FALSE) {
    return(list(logs = checkout$logs))
  }
  
  ## prepare objects for TMB
  data <- list (code = "B",
                dat = t(logqs))
  
  ## build parameters object
  parameters <- list(rw = rw,
                     logSdMVN = rep(0,nrow(data$dat)),
                     logitRho = 0,
                     mu         = rep(0,nrow(data$dat)),
                     logitARrho = rep(0,nrow(data$dat)),
                     logSdObs = rep(0,nrow(data$dat)),
                     logTau = 0,
                     logb0  = 0,
                     logb1  = 0,
                     dat_missing = rep(0, sum(is.na(data$dat))))
  
  map1 <- list(logitRho = as.factor(NA),
               logTau = as.factor(NA), 
               logb0  = as.factor(NA), 
               logb1  = as.factor(NA))
  
  # ---------------------------------#
  # Build model and optimise
  # ---------------------------------#
  
  obj <- tryCatch({MakeADFun(data, 
                             parameters, 
                             random = "rw",
                             DLL = "MixUncertainty", 
                             map = map1,
                             silent = TRUE)},
                  error = function(e)e)
  
  opt <- tryCatch({nlminb(obj$par, obj$fn, obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)
  
  sdr <- tryCatch({sdreport(obj)},
                  error = function(e)e)
  
  return(list(sdr = sdr,
              obj = obj,
              opt = opt))
}
