# ---
# title: 'Define wrapper functions to fit TMB Dirichlet models'
# author: 'Various'
# date: 'August 2023'
# ---

#' Wrapper function to fit latent AR1 MVN with hurdle Dirichlet and Bernoulli observation error
#'
#' A wrapper function to construct data and parameter inputs to build a TMB
#' objective function.
#'
#' @param dat A matrix of time-series observations for compositional data where
#'            cols = years and rows = compositional elements
#' @param rw A matrix of latent-state parameters to be optimised
#'
#' @return A list containing:
#'
#'         - 'sdr': a TMB report summary for the optimised model parameters,
#'         - 'opt': a summary returned by the optimisation function,
#'         - 'obj': the optimised TMB objective function
#'
#' @export

fitMVNhurdle <- function(dat, rw) {

  ## build data object
  data <- list (dat = t(dat))

  ## build parameters object
  parameters <- list(logSdMVN = rep(0,nrow(rw)),
                     logitRho = 0,
                     logTau   = 0,
                     mu         = rep(0,nrow(rw)),
                     logitARrho = rep(0,nrow(rw)),
                     logb0      = 0,
                     logb1      = 0,
                     rw = rw)

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "fit_Dirichlet_MVN_AR1_missing_hurdle_nminus1traj",
                             silent = TRUE)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- sdreport(obj)

  return(list(sdr = sdr,
              opt = opt,
              obj = obj))

}

#' Wrapper function to fit latent AR1 Gaussian with hurdle Dirichlet and Bernoulli observation error
#'
#' A wrapper function to construct data and parameter inputs to build a TMB
#' objective function.
#'
#' @param dat A matrix of time-series observations for compositional data where
#'            cols = years and rows = compositional elements
#' @param rw A matrix of latent-state parameters to be optimised
#'
#' @return A list containing:
#'
#'         - 'sdr': a TMB report summary for the optimised model parameters,
#'         - 'opt': a summary returned by the optimisation function,
#'         - 'obj': the optimised TMB objective function
#'
#' @export

fitNhurdle <- function(dat, rw) {

  ## build data object
  data <- list (dat = t(dat))

  ## build parameters object
  parameters <- list(logSdMVN = rep(0,nrow(rw)),
                     logTau   = 0,
                     mu         = rep(0,nrow(rw)),
                     logitARrho = rep(0,nrow(rw)),
                     logb0      = 0,
                     logb1      = 0,
                     rw = rw)

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "fit_Dirichlet_N_AR1_missing_hurdle_nminus1traj",
                             silent = TRUE)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- sdreport(obj)

  return(list(sdr = sdr,
              opt = opt,
              obj = obj))

}

#' Wrapper function to fit latent AR1 Gaussian with fixed AR1 correlation and hurdle Dirichlet and Bernoulli observation error
#'
#' A wrapper function to construct data and parameter inputs to build a TMB
#' objective function.
#'
#' @param dat A matrix of time-series observations for compositional data where
#'            cols = years and rows = compositional elements
#' @param rw A matrix of latent-state parameters to be optimised
#'
#' @return A list containing:
#'
#'         - 'sdr': a TMB report summary for the optimised model parameters,
#'         - 'opt': a summary returned by the optimisation function,
#'         - 'obj': the optimised TMB objective function
#'
#' @export

fitNhurdle_fixARrho <- function(dat, rw, logitARrho) {

  ## build data object
  data <- list (dat = t(dat))

  ## build parameters object
  parameters <- list(logSdMVN = rep(0,nrow(rw)),
                     logTau   = 0,
                     mu         = rep(0,nrow(rw)),
                     logitARrho = rep(logitARrho,nrow(rw)),
                     logb0      = 0,
                     logb1      = 0,
                     rw = rw)

  map1 <- list(logitARrho = rep(as.factor(NA), nrow(rw)))

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "fit_Dirichlet_N_AR1_missing_hurdle_nminus1traj",
                             silent = TRUE,
                             map = map1)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- sdreport(obj)

  return(list(sdr = sdr,
              opt = opt,
              obj = obj))

}

#' Wrapper function to fit latent AR1 Gaussian with hurdle Dirichlet observation error and fixed Bernoulli parameters
#'
#' A wrapper function to construct data and parameter inputs to build a TMB
#' objective function.
#'
#' @param dat A matrix of time-series observations for compositional data where
#'            cols = years and rows = compositional elements
#' @param rw A matrix of latent-state parameters to be optimised
#'
#' @return A list containing:
#'
#'         - 'sdr': a TMB report summary for the optimised model parameters,
#'         - 'opt': a summary returned by the optimisation function,
#'         - 'obj': the optimised TMB objective function
#'
#' @export

fitNhurdle_fixb0b1 <- function(dat, rw, logb0, logb1) {

  ## build data object
  data <- list (dat = t(dat))

  ## build parameters object
  parameters <- list(logSdMVN = rep(0,nrow(rw)),
                     logTau   = 0,
                     mu         = rep(0,nrow(rw)),
                     logitARrho = rep(0,nrow(rw)),
                     logb0      = logb0,
                     logb1      = logb1,
                     rw = rw)

  map1 <- list(logb0 = as.factor(NA), logb1 = as.factor(NA))

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "fit_Dirichlet_N_AR1_missing_hurdle_nminus1traj",
                             silent = TRUE,
                             map = map1)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- sdreport(obj)

  return(list(sdr = sdr,
              opt = opt,
              obj = obj))

}

#' Wrapper function to fit latent Random Walk with MVN distributed intervals and Dirichlet observation error
#'
#' A wrapper function to construct data and parameter inputs to build a TMB
#' objective function.
#'
#' @param dat A matrix of time-series observations for compositional data where
#'            cols = years and rows = compositional elements
#' @param rw A matrix of latent-state parameters to be optimised
#'
#' @return A list containing:
#'
#'         - 'sdr': a TMB report summary for the optimised model parameters,
#'         - 'opt': a summary returned by the optimisation function,
#'         - 'obj': the optimised TMB objective function
#'
#' @export

fitMVNDirrw <- function(dat, rw) {

  ## build data object
  data <- list (dat = t(dat))

  ## build parameters object
  parameters <- list(logSdMVN = rep(0,nrow(rw)),
                     logitRho = 0,
                     logTau   = 0,
                     rw = rw)

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "fit_Dirichlet_rw_missing",
                             silent = TRUE)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- sdreport(obj)

  return(list(sdr = sdr,
              opt = opt,
              obj = obj))

}

#' Wrapper function to fit latent Random Walk with Gaussian distributed intervals and Dirichlet observation error
#'
#' A wrapper function to construct data and parameter inputs to build a TMB
#' objective function.
#'
#' @param dat A matrix of time-series observations for compositional data where
#'            cols = years and rows = compositional elements
#' @param rw A matrix of latent-state parameters to be optimised
#'
#' @return A list containing:
#'
#'         - 'sdr': a TMB report summary for the optimised model parameters,
#'         - 'opt': a summary returned by the optimisation function,
#'         - 'obj': the optimised TMB objective function
#'
#' @export

fitNDirrw <- function(dat, rw) {

  ## build data object
  data <- list (dat = t(dat))

  ## build parameters object
  parameters <- list(logSdMVN = rep(0,nrow(rw)),
                     logTau   = 0,
                     rw = rw)

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "fit_Dirichlet_rw_missing_nocorr",
                             silent = TRUE)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- sdreport(obj)

  return(list(sdr = sdr,
              opt = opt,
              obj = obj))

}
