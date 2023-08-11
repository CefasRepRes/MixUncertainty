# ---
# title: 'Define helper functions to carry out stochastic forecasts'
# author: 'Various'
# date: 'August 2023'
# ---

#' Stochastic forecast from an AR1 latent process
#'
#' ...
#'
#' @param dat
#' @param pl
#' @param plsd
#' @param MVNcov
#' @param fillyear
#' @param nit

## forecast from hurdle
forecast_from_hurdle <- function(dat, pl, plsd, MVNcov, fillyear, nit) {

  ## Define number of forecast years
  nyear_forecast <- length(fillyear)

  ## convert AR1 correlation parameter
  ARrho <- 2 * exp(pl$logitARrho)/(exp(pl$logitARrho) + 1) - 1

  ## calculate parameter sigma
  ARsigma <- sqrt(1 - (ARrho^2))

  ## Generate matrix of predictions
  resvariates <- replicate(nit, {

    ## define a random walk matrix and a prediction matrix for the forecast period
    predMatrix  <- matrix(NA, nrow = nrow(pl$rw)+1, ncol = nyear_forecast+1)
    alphaMatrix <- matrix(NA, nrow = nrow(pl$rw)+1, ncol = nyear_forecast+1)
    rwMatrix    <- matrix(NA, nrow = nrow(pl$rw), ncol = nyear_forecast+1)

    ## extract starting random walk
    rwMatrix[,1] <- sapply(1:nrow(pl$rw), function(x) rnorm(1,
                                                            pl$rw[x,ncol(pl$rw)],
                                                            plsd$rw[x,ncol(pl$rw)]))

    ## calculate the starting value for each variable
    predMatrix[,1] <- invlogit(rwMatrix[,1])

    for(i in 2:(nyear_forecast+1)){

      ## sample random error interval
      randomInterval <- MASS::mvrnorm(1, rep(0, nrow(MVNcov)), MVNcov)

      ## adjust data given mean and AR1 correlation parameters
      rwMatrix[,i] <- pl$mu + ARrho* (rwMatrix[,i-1] - pl$mu) + (ARsigma * randomInterval)

      # ## calculate probability that any observation is zero
      # l <- (-exp(pl$logb0)) + exp(pl$logb1) * invlogit(rwMatrix[,i])
      # p <- 1 / (1 + exp(-l))
      #
      # ## sample zero observations
      # zeros <- rbinom(n = length(predMatrix[,i]), size = rep(1,length(predMatrix[,i])), prob = p)

      ## Sample a value for each component
      predMatrix[,i] <- invlogit(rwMatrix[,i]) #* zeros

    }
    return(t(predMatrix[,-1]))
  })

  ## Reinsert year and stock dimension names
  dimnames(resvariates) <- list(year  = fillyear,
                                stock = dimnames(dat)[[2]],
                                iter  = 1:nit)

  return(resvariates)

}

#' Stochastic forecast from a random walk model
#'
#' ...
#'
#' @param dat
#' @param pl
#' @param plsd
#' @param MVNcov
#' @param fillyear
#' @param nit

## forecast from random walk
forecast_from_rw <- function(dat, pl, plsd, MVNcov, fillyear, nit) {

  ## Define number of forecast years
  nyear_forecast <- length(fillyear)

  ## Generate matrix of predictions
  resvariates <- replicate(nit, {

    ## define a random walk matrix and a prediction matrix for the forecast period
    predMatrix  <- matrix(NA, nrow = nrow(pl$rw)+1, ncol = nyear_forecast+1)
    alphaMatrix <- matrix(NA, nrow = nrow(pl$rw)+1, ncol = nyear_forecast+1)
    rwMatrix    <- matrix(NA, nrow = nrow(pl$rw), ncol = nyear_forecast+1)

    ## extract starting random walk
    rwMatrix[,1] <- sapply(1:nrow(pl$rw), function(x) rnorm(1,
                                                            pl$rw[x,ncol(pl$rw)],
                                                            plsd$rw[x,ncol(pl$rw)]))

    ## calculate the starting value for each variable
    predMatrix[,1] <- invlogit(rwMatrix[,1])

    for(i in 2:(nyear_forecast+1)){

      ## sample random interval
      rwMatrix[,i] <- MASS::mvrnorm(1, rwMatrix[,i-1], MVNcov)

      ## Sample a value for each component
      predMatrix[,i] <- invlogit(rwMatrix[,i])
    }
    return(t(predMatrix[,-1]))
  })

  ## Reinsert year and stock dimension names
  dimnames(resvariates) <- list(year  = fillyear,
                                stock = dimnames(dat)[[2]],
                                iter  = 1:nit)

  return(resvariates)

}
