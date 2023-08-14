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
  data <- list (code = "E",
                dat = t(dat))

  ## build parameters object
  parameters <- list(rw = rw,
                     logSdMVN = rep(0,nrow(rw)),
                     logitRho = 0,
                     mu         = rep(0,nrow(rw)),
                     logitARrho = rep(0,nrow(rw)),
                     logSdObs = 0,
                     logTau   = 0,
                     logb0      = 0,
                     logb1      = 0,
                     dat_missing = vector(mode = "numeric"))

  map1 <- list(logSdObs = as.factor(NA))

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "MixUncertainty",
                             map = map1,
                             silent = TRUE)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- tryCatch({sdreport(obj)},
                  error = function(e)e)

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
  data <- list (code = "F",
                dat = t(dat))

  ## build parameters object
  parameters <- list(rw = rw,
                     logSdMVN = rep(0,nrow(rw)),
                     logitRho = 0,
                     mu         = rep(0,nrow(rw)),
                     logitARrho = rep(0,nrow(rw)),
                     logSdObs = 0,
                     logTau = 0,
                     logb0  = 0,
                     logb1  = 0,
                     dat_missing = vector(mode = "numeric"))

  map1 <- list(logitRho = as.factor(NA),
               logSdObs = as.factor(NA))

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "MixUncertainty",
                             map = map1,
                             silent = TRUE)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- tryCatch({sdreport(obj)},
                  error = function(e)e)

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
  data <- list (code = "F",
                dat = t(dat))

  ## build parameters object
  parameters <- list(rw = rw,
                     logSdMVN = rep(0,nrow(rw)),
                     logitRho = 0,
                     mu         = rep(0,nrow(rw)),
                     logitARrho = rep(logitARrho,nrow(rw)),
                     logSdObs = 0,
                     logTau   = 0,
                     logb0      = 0,
                     logb1      = 0,
                     dat_missing = vector(mode = "numeric"))

  map1 <- list(logitRho = as.factor(NA),
               logitARrho = rep(as.factor(NA), nrow(rw)),
               logSdObs = as.factor(NA))

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "MixUncertainty",
                             silent = TRUE,
                             map = map1)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- tryCatch({sdreport(obj)},
                  error = function(e)e)

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
  data <- list (code = "F",
                dat = t(dat))

  ## build parameters object
  parameters <- list(rw = rw,
                     logSdMVN = rep(0,nrow(rw)),
                     logitRho = 0,
                     mu         = rep(0,nrow(rw)),
                     logitARrho = rep(0,nrow(rw)),
                     logSdObs = 0,
                     logTau   = 0,
                     logb0      = logb0,
                     logb1      = logb1,
                     dat_missing = vector(mode = "numeric"))

  map1 <- list(logitRho = as.factor(NA),
               logSdObs = as.factor(NA),
               logb0 = as.factor(NA),
               logb1 = as.factor(NA))

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "MixUncertainty",
                             silent = TRUE,
                             map = map1)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- tryCatch({sdreport(obj)},
                  error = function(e)e)

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
  data <- list (code = "C",
                dat = t(dat))

  ## build parameters object
  parameters <- list(rw = rw,
                     logSdMVN = rep(0,nrow(rw)),
                     logitRho = 0,
                     mu = 0,
                     logitARrho = 0,
                     logSdObs = 0,
                     logTau   = 0,
                     logb0 = 0,
                     logb1 = 0,
                     dat_missing = vector(mode = "numeric"))

  map1 <- list(mu = as.factor(NA),
               logitARrho = as.factor(NA),
               logSdObs = as.factor(NA),
               logb0 = as.factor(NA),
               logb1 = as.factor(NA))

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "MixUncertainty",
                             map = map1,
                             silent = TRUE)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- tryCatch({sdreport(obj)},
                  error = function(e)e)

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
  data <- list (code = "D",
                dat = t(dat))

  ## build parameters object
  parameters <- list(rw = rw,
                     logSdMVN = rep(0,nrow(rw)),
                     logitRho = 0,
                     mu = 0,
                     logitARrho = 0,
                     logSdObs = 0,
                     logTau   = 0,
                     logb0 = 0,
                     logb1 = 0,
                     dat_missing = vector(mode = "numeric"))

  map1 <- list(logitRho = as.factor(NA),
               mu = as.factor(NA),
               logitARrho = as.factor(NA),
               logSdObs = as.factor(NA),
               logb0 = as.factor(NA),
               logb1 = as.factor(NA))

  obj <- tryCatch({MakeADFun(data,
                             parameters,
                             random = "rw",
                             DLL = "MixUncertainty",
                             map = map1,
                             silent = TRUE)},
                  error = function(e)e)

  opt <- tryCatch({nlminb(obj$par,
                          obj$fn,
                          obj$gr,
                          control = list(eval.max = 1000, iter.max = 1000))},
                  error = function(e)e)

  sdr <- tryCatch({sdreport(obj)},
                  error = function(e)e)

  return(list(sdr = sdr,
              opt = opt,
              obj = obj))
}
