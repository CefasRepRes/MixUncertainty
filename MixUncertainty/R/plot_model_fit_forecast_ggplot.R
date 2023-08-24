# ---
# title: 'Functions to generate illustrive plots for model fits'
# author: 'Various'
# date: 'August 2023'
# ---

#' Plot latent-process model fit and forecast for a random walk and noise state-space model
#'
#' Plot the marginal maximum likelihood estimate and median stochastic forecast
#' for the latent process of a random walk and observation noise state-space model.
#' The 95% confidence intervals around the estimate and forecast are shown as
#' dashed lines. Observation data are shown as open points.
#'
#' @param data named list input to TMB containing observation data
#' @param pl   named list output from TMB containing parameter estimates
#' @param plsd named list output from TMB containing parameter standard error
#' @param pred_quantiles array of forecast value summaries at the 2.5, 50, 97.5
#'                       centiles.
#' @param years vector of years associated with observations and forecasts
#'
#' @importFrom tidyr gather
#' @import ggplot2

plot_forecast_MVN <- function(data, pl, plsd, pred_quantiles, years) {

  ## Check that dplyr, tidyr and ggplot2 are available
  if (!all(requireNamespace("dplyr", quietly = TRUE),
           requireNamespace("tidyr", quietly = TRUE),
           requireNamespace("ggplot2", quietly = TRUE))) {
    stop("packages 'dplyr', 'tidyr' and 'ggplot2' are needed for plotting")
  }

  makelong <- function(data, dref) {
    dd <- as.data.frame(t(data))
    names(dd) <- rownames(dref)
    dd$years <- as.numeric(rownames(dd))
    dd <- gather(dd, metier, effortshare, -years)
  }

  makelongfwd <- function(data) {
    dd <- as.data.frame(data)
    dd$years <- as.numeric(rownames(dd))
    dd <- gather(dd, metier, effortshare, -years)
  }

  ## compute historic confidence intervals -- best to sample these because we need
  ## to inverse-logit BEFORE calculating the quantiles
  rw_sampled <- sapply(1:nrow(pl$rw), function(x)
    sapply(1:ncol(pl$rw), function(y)
      rnorm(50000, pl$rw[x,y], plsd$rw[x,y]),simplify = "array"), simplify = "array")

  est   <- pl$rw
  lower <- t(apply(rw_sampled, MARGIN = c(2,3), quantile, prob = c(0.025)))
  upper <- t(apply(rw_sampled, MARGIN = c(2,3), quantile, prob = c(0.975)))

  ## get y-range
  ylim <- c(min(min(lower,na.rm = TRUE),
                min(data$dat, na.rm = TRUE),
                min(pred_quantiles, na.rm = TRUE)),
            max(max(upper,na.rm = TRUE),
                max(data$dat, na.rm = TRUE),
                max(pred_quantiles, na.rm = TRUE)))

  ## rearrange data
  data$dat[data$dat == 0] <- NA
  d <- as.data.frame(t(log(data$dat)))
  d$years <- as.numeric(rownames(d))
  d <- tidyr::gather(d, metier, effortshare, -years)

  colnames(lower) <- colnames(pl$rw)
  colnames(upper) <- colnames(pl$rw)

  dest <- makelong(data = est, dref = data$dat)
  dlow <- makelong(data = lower, dref = data$dat)
  dupp <- makelong(data = upper, dref = data$dat)

  ## rearrange forecast
  if (dim(lower)[1] > 1) {
    fwd_low <- cbind(t(lower)[ncol(pl$rw),], pred_quantiles[,,1])
    fwd_est <- cbind(t(est)[ncol(pl$rw),], pred_quantiles[,,2])
    fwd_upp <- cbind(t(upper)[ncol(pl$rw),], pred_quantiles[,,3])

    colnames(fwd_low) <- colnames(fwd_est) <- colnames(fwd_upp) <-
      (as.numeric(tail(colnames(data$dat),1))) : (as.numeric(tail(colnames(data$dat),1)) + dim(pred_quantiles)[1])
    rownames(fwd_low) <- rownames(fwd_est) <- rownames(fwd_upp) <- rownames(data$dat)

    fwd_est <- makelongfwd(data = t(fwd_est))
    fwd_low <- makelongfwd(data = t(fwd_low))
    fwd_upp <- makelongfwd(data = t(fwd_upp))

  } else {
    fwd_low <- t(as.matrix(c(t(lower)[ncol(pl$rw),], pred_quantiles[,,1])))
    fwd_est <- t(as.matrix(c(t(est)[ncol(pl$rw),], pred_quantiles[,,2])))
    fwd_upp <- t(as.matrix(c(t(upper)[ncol(pl$rw),], pred_quantiles[,,3])))

    colnames(fwd_low) <- colnames(fwd_est) <- colnames(fwd_upp) <- c(colnames(lower)[ncol(pl$rw)], years)
    rownames(fwd_low) <- rownames(fwd_est) <- rownames(fwd_upp) <- rownames(data$dat)

    fwd_est <- makelongfwd(data = t(fwd_est))
    fwd_low <- makelongfwd(data = t(fwd_low))
    fwd_upp <- makelongfwd(data = t(fwd_upp))
  }

  ## plot fit
  p1 <- ggplot(mapping = aes(x = years, y = effortshare, colour = metier)) +
    geom_point(data = d, shape = 21) +
    geom_line(data = dest) +
    geom_line(data = dupp, linetype = 2) +
    geom_line(data = dlow, linetype = 2) +
    geom_line(data = fwd_est, linewidth = 0.9) +
    geom_line(data = fwd_low, linewidth = 0.9, linetype = 2) +
    geom_line(data = fwd_upp, linewidth = 0.9, linetype = 2) +
    scale_y_continuous("catchability") +
    scale_x_continuous("year") +
    scale_colour_viridis_d(option = "D") +
    theme_classic()

  return(p1)
}

#' Plot latent-process model fit for a random walk/AR1 and Dirichlet noise state-space model
#'
#' Plot the marginal maximum likelihood estimate for the latent process of a
#' random walk (or AR1) state-space model with Dirichlet observation noise.
#' The 95% confidence intervals around the estimate are shown as dashed lines.
#' Observation data are shown as open points.
#'
#' @param data named list input to TMB containing observation data
#' @param pl   named list output from TMB containing parameter estimates
#' @param plsd named list output from TMB containing parameter standard error
#' @param years vector of years associated with observations and forecasts
#' @param invlogitfun function to perform multinomial inverse logit transformation
#'                    of latent process values.
#'
#' @importFrom tidyr gather
#' @import ggplot2

plot_fit_Dir <- function(data, pl, plsd, years, invlogitfun) {

  ## Check that dplyr, tidyr and ggplot2 are available
  if (!all(requireNamespace("dplyr", quietly = TRUE),
           requireNamespace("tidyr", quietly = TRUE),
           requireNamespace("ggplot2", quietly = TRUE))) {
    stop("packages 'dplyr', 'tidyr' and 'ggplot2' are needed for plotting")
  }

  makelong <- function(data, dref) {
    dd <- as.data.frame(t(data))
    names(dd) <- rownames(dref)
    dd$years <- as.numeric(rownames(dd))
    dd <- gather(dd, metier, effortshare, -years)
  }

  ## compute historic confidence intervals -- best to sample these because we need
  ## to inverse-logit BEFORE calculating the quantiles
  rw_sampled <- sapply(1:nrow(pl$rw), function(x)
    sapply(1:ncol(pl$rw), function(y)
      rnorm(50000, pl$rw[x,y], plsd$rw[x,y]),simplify = "array"), simplify = "array")

  rw_sampled <- apply(rw_sampled, MARGIN = c(1,2), invlogitfun)

  est   <- apply(t(pl$rw), MARGIN = c(1), invlogitfun)
  lower <- apply(rw_sampled, MARGIN = c(1,3), quantile, prob = c(0.025))
  upper <- apply(rw_sampled, MARGIN = c(1,3), quantile, prob = c(0.975))

  ## get y-range
  ylim <- c(min(min(lower,na.rm = TRUE), min(data$dat, na.rm = TRUE)),
            max(max(upper,na.rm = TRUE), max(data$dat, na.rm = TRUE)))

  ## rearrange data
  d <- as.data.frame(t(data$dat))
  d$years <- as.numeric(rownames(d))
  d <- gather(d, metier, effortshare, -years)

  colnames(lower) <- colnames(data$dat)
  colnames(upper) <- colnames(data$dat)

  dest <- makelong(data = est, dref = data$dat)
  dlow <- makelong(data = lower, dref = data$dat)
  dupp <- makelong(data = upper, dref = data$dat)

  ## plot fit
  p1 <- ggplot(mapping = aes(x = years, y = effortshare, colour = metier)) +
    geom_point(data = d, shape = 21) +
    geom_line(data = dest) +
    geom_line(data = dupp, linetype = 2) +
    geom_line(data = dlow, linetype = 2) +
    scale_y_continuous("proportional effort-share") +
    scale_x_continuous("year") +
    scale_colour_viridis_d(option = "D") +
    theme_classic()

  return(p1)

}

#' Plot latent-process model fit and forecast for a random walk/AR1 and Dirichlet noise state-space model
#'
#' Plot the marginal maximum likelihood estimate and median stochastic forecast
#' for the latent process of a random walk (or AR1) state-space model with Dirichlet
#' observation noise. The 95% confidence intervals around the estimate and forecast
#' are shown as dashed lines. Observation data are shown as open points.
#'
#' @param data named list input to TMB containing observation data
#' @param pl   named list output from TMB containing parameter estimates
#' @param plsd named list output from TMB containing parameter standard error
#' @param pred_quantiles array of forecast value summaries at the 2.5, 50, 97.5
#'                       centiles.
#' @param years vector of years associated with observations and forecasts
#' @param invlogitfun function to perform multinomial inverse logit transformation
#'                    of latent process values.
#'
#' @importFrom tidyr gather
#' @import ggplot2

plot_forecast_Dir <- function(data, pl, plsd, pred_quantiles, years, invlogitfun) {

  ## Check that dplyr, tidyr and ggplot2 are available
  if (!all(requireNamespace("dplyr", quietly = TRUE),
           requireNamespace("tidyr", quietly = TRUE),
           requireNamespace("ggplot2", quietly = TRUE))) {
    stop("packages 'dplyr', 'tidyr' and 'ggplot2' are needed for plotting")
  }

  makelong <- function(data, dref) {
    dd <- as.data.frame(t(data))
    names(dd) <- rownames(dref)
    dd$years <- as.numeric(rownames(dd))
    dd <- gather(dd, metier, effortshare, -years)
  }

  makelongfwd <- function(data) {
    dd <- as.data.frame(data)
    dd$years <- as.numeric(rownames(dd))
    dd <- gather(dd, metier, effortshare, -years)
  }

  ## number of forecast years
  nyears <- dim(pred_quantiles)[2] - 1

  ## compute historic confidence intervals -- best to sample these because we need
  ## to inverse-logit BEFORE calculating the quantiles
  rw_sampled <- sapply(1:nrow(pl$rw), function(x)
    sapply(1:ncol(pl$rw), function(y)
      rnorm(50000, pl$rw[x,y], plsd$rw[x,y]),simplify = "array"), simplify = "array")

  rw_sampled <- apply(rw_sampled, MARGIN = c(1,2), invlogitfun)

  est   <- apply(t(pl$rw), MARGIN = c(1), invlogitfun)
  lower <- apply(rw_sampled, MARGIN = c(1,3), quantile, prob = c(0.025))
  upper <- apply(rw_sampled, MARGIN = c(1,3), quantile, prob = c(0.975))

  ## get y-range
  ylim <- c(min(min(lower,na.rm = TRUE),
                min(data$dat, na.rm = TRUE),
                min(pred_quantiles, na.rm = TRUE)),
            max(max(upper,na.rm = TRUE),
                max(data$dat, na.rm = TRUE),
                max(pred_quantiles, na.rm = TRUE)))

  ## rearrange data
  d <- as.data.frame(t(data$dat))
  d$years <- as.numeric(rownames(d))
  d <- gather(d, metier, effortshare, -years)

  colnames(lower) <- colnames(data$dat)
  colnames(upper) <- colnames(data$dat)

  dest <- makelong(data = est, dref = data$dat)
  dlow <- makelong(data = lower, dref = data$dat)
  dupp <- makelong(data = upper, dref = data$dat)

  ## rearrange forecast
  fwd_low <- rbind(t(lower)[ncol(pl$rw),], pred_quantiles[,,1])
  fwd_est <- rbind(t(est)[ncol(pl$rw),], pred_quantiles[,,2])
  fwd_upp <- rbind(t(upper)[ncol(pl$rw),], pred_quantiles[,,3])

  rownames(fwd_low) <- rownames(fwd_est) <- rownames(fwd_upp) <-
    (as.numeric(tail(colnames(data$dat),1))) : (as.numeric(tail(colnames(data$dat),1)) + dim(pred_quantiles)[1])
  colnames(fwd_low) <- colnames(fwd_est) <- colnames(fwd_upp) <- rownames(data$dat)

  fwd_est <- makelongfwd(data = fwd_est)
  fwd_low <- makelongfwd(data = fwd_low)
  fwd_upp <- makelongfwd(data = fwd_upp)

  ## plot fit
  p1 <- ggplot(mapping = aes(x = years, y = effortshare, colour = metier)) +
    geom_point(data = d, shape = 21) +
    geom_line(data = dest) +
    geom_line(data = dupp, linetype = 2) +
    geom_line(data = dlow, linetype = 2) +
    geom_line(data = fwd_est, linewidth = 0.9) +
    geom_line(data = fwd_low, linewidth = 0.9, linetype = 2) +
    geom_line(data = fwd_upp, linewidth = 0.9, linetype = 2) +
    scale_y_continuous("proportional effort-share") +
    scale_x_continuous("year") +
    scale_colour_viridis_d(option = "D") +
    theme_classic()

  return(p1)
}
