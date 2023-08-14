# ---
# title: 'Functions to generate illustrive plots for model fits'
# author: 'Various'
# date: 'August 2023'
# ---

#' Plot latent-process model fit for a random walk and noise state-space model
#'
#' Plot the marginal maximum likelihood estimate for the latent process of a
#' random walk and observation noise state-space model together with 95%
#' confidence intervals around the estimate and observation points.
#'
#' @param data named list input to TMB containing observation data
#' @param pl   named list output from TMB containing parameter estimates
#' @param plsd named list output from TMB containing parameter standard error
#' @param years vector of years associated with observation data

plot_fit_MVN <- function(data, pl, plsd, years) {

  require(ggplot)

  ## compute confidence intervals
  lower <- t(pl$rw)+2*t(plsd$rw)
  upper <- t(pl$rw)-2*t(plsd$rw)

  ## get y-range
  ylim <- c(min(min(lower,na.rm = TRUE), min(data$dat, na.rm = TRUE)),
            max(max(upper,na.rm = TRUE), max(data$dat, na.rm = TRUE)))

  layout(matrix(c(1,2),nrow=1), width=c(4,2))
  par(mar=c(5,4,4,0)) #No margin on the right side
  matplot(t(data$dat), pch = 1, ylim = ylim, xlab = "year", ylab = "log-catchability",axes=F)
  axis(2)
  axis(side = 1, at = 1:ncol(data$dat), labels = years)
  matplot(t(pl$rw),   type = "l", lty = 1, add = TRUE)
  matplot(lower, type = "l", lty = 2, add = TRUE)
  matplot(upper, type = "l", lty = 2, add = TRUE)
  par(mar=c(5,0,4,2)) #No margin on the left side
  plot(c(0,1), type = "n", axes = F, xlab = "", ylab = "")
  legend("center", colnames(t(data$dat)), col = seq_len(ncol(t(data$dat))),
         cex = 0.8, fill = seq_len(ncol(t(data$dat))))

}

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

plot_forecast_MVN <- function(data, pl, plsd, pred_quantiles, years) {

  require(ggplot)

  ## number of forecast years
  nyears <- dim(pred_quantiles)[2] - 1

  ## compute historic confidence intervals
  lower <- t(pl$rw)+2*t(plsd$rw)
  upper <- t(pl$rw)-2*t(plsd$rw)

  ## get y-range
  ylim <- c(min(min(lower,na.rm = TRUE),
                min(data$dat, na.rm = TRUE),
                min(pred_quantiles, na.rm = TRUE)),
            max(max(upper,na.rm = TRUE),
                max(data$dat, na.rm = TRUE),
                max(pred_quantiles, na.rm = TRUE)))

  ## Generate fit plot
  layout(matrix(c(1,2),nrow=1), width=c(4,2))
  par(mar=c(5,4,4,0)) #No margin on the right side
  matplot(t(data$dat), pch = 1,
          xlim = c(0, ncol(pl$rw)+nyears), ylim = ylim,
          xlab = "year", ylab = "log-catchability",axes=F)
  axis(2)
  axis(side = 1, at = 1:(ncol(pl$rw)+nyears), labels = years)
  matplot(t(pl$rw),   type = "l", lty = 1, add = TRUE)
  matplot(upper, type = "l", lty = 2, add = TRUE)
  matplot(lower, type = "l", lty = 2, add = TRUE)

  ## Add forecast
  matplot(rbind((t(pl$rw)*NA)[1:(nyrs-1),], t(pred_quantiles[,,1])),
          add = TRUE, type = "l", lwd = 2, lty = 2)
  matplot(rbind((t(pl$rw)*NA)[1:(nyrs-1),], t(pred_quantiles[,,2])),
          add = TRUE, type = "l", lwd = 2, lty = 1)
  matplot(rbind((t(pl$rw)*NA)[1:(nyrs-1),], t(pred_quantiles[,,3])),
          add = TRUE, type = "l", lwd = 2, lty = 2)
  par(mar=c(5,0,4,2)) #No margin on the left side
  plot(c(0,1), type = "n", axes = F, xlab = "", ylab = "")
  legend("center", colnames(t(data$dat)), col = seq_len(ncol(t(data$dat))),
         cex = 0.8, fill = seq_len(ncol(t(data$dat))))

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

plot_fit_Dir <- function(data, pl, plsd, years, invlogitfun) {

  require(ggplot)

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

plot_forecast_Dir <- function(data, pl, plsd, pred_quantiles, years, invlogitfun) {

  require(ggplot)

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
