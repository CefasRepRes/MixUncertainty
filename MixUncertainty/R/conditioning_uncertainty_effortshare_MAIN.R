# ---
# title: 'Functions to stochastically forecast fleet-metier effort-share'
# author: 'Various'
# date: 'November 2022'
# ---

#' Stochastically forecast future metier effort-share from historic estimates
#'
#' Function takes a \code{FLFleet}, \code{FLFleets}, \code{FLFleetExt} or a
#' \code{FLFleetsExt} object and fits a statistical model to historic
#' observations of metier effort-share for each fleet. The metier effort-share
#' for future years is a stochastic forecast from fitted model parameters.
#'
#' @param fleets An object of class \code{FLFleet}, \code{FLFleets},
#'               \code{FLFleetExt} or a \code{FLFleetsExt}.
#' @param method (Character) The algorithm used to fit the statistical model and
#'               generate stochastic forecasts. Using the default is highly
#'               recommended.
#' @param datayears A vector of historic years used to estimate future fleet
#'                  effort-share
#' @param TACyear (Integer) The projection year in the fleets object
#' @param deterministic (Logical) Should the results for first iteration be
#'                      a simple mean of the historical period? Defaults
#'                      to \code{TRUE}
#' @param deterministic_yrs (Integer) The number of recent data years to use in
#'                          the deterministic calculation. Defaults to 3.
#' @param detMethod (Character) The method to calculate the deterministic average.
#'                  If "mean" then a simple mean of the most recent data years
#'                  is used. Otherwise, the median of the stochastic forecast is
#'                  used.
#' @param verbose (Logical) Should the function print progress? Defaults
#'                to \code{TRUE}
#' @param makeLog (Logical) Should the function record model fitting and forecasting
#'                performance for each fleet? Defaults to \code{TRUE}
#' @param makePlots (Logical) Should the function generate a figure showing
#'                  observations and model fit and forecasts with uncertainty
#'                  intervals for each fleet? Defaults to \code{TRUE}
#' @param removeNonRecent (Logical) Should the function remove fleets that have not
#'                        registered landings in recent years? Defaults to \code{FALSE}
#'
#' @return a named list containing:
#'         - **fleets** - An updated \code{FLFleet}, \code{FLFleets},
#'               \code{FLFleetExt} or \code{FLFleetsExt} object
#'         - **logs** - An optional summary of model fitting and forecasting
#'                      performance for each fleets
#'         - **plots** - An optional figure showing the model maximum likelihood
#'                       fit (+ 95% CI) and median (+ 95% CI) forecasts against
#'                       observations for each fleet.
#'
#' @export

setGeneric("uncertainty_effortshare", function(fleets,
                                               method = "TMB_DirMissingAR1Hurdle",
                                               datayears     = NULL,
                                               TACyear       = NULL,
                                               deterministic = TRUE,
                                               deterministic_yrs = 3,
                                               detMethod     = "mean",
                                               verbose       = TRUE,
                                               makeLog       = TRUE,
                                               makePlots     = TRUE,
                                               removeNonRecent = FALSE) {
  standardGeneric("uncertainty_effortshare")
})

## if fleets == "FLFleetExt"
#' @rdname uncertainty_effortshare
setMethod(f = "uncertainty_effortshare",
          signature = signature(fleets = "FLFleetExt"),
          definition = function(fleets,
                                method = "TMB_DirMissingAR1Hurdle",
                                datayears     = NULL,
                                TACyear       = NULL,
                                deterministic = TRUE,
                                deterministic_yrs = 3,
                                detMethod     = "mean",
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE,
                                removeNonRecent = FALSE) {

            ## Function expects that multiple iterations are present
            if(dims(fleets[[f]])$iter < 2)
              stop("input should have > 1 iterations to store sampled uncertainty")

            stop("Methods not yet implemented for 'FLFleetExt'")

          })

## if fleets == "FLFleetsExt"
#' @rdname uncertainty_effortshare
setMethod(f = "uncertainty_effortshare",
          signature = signature(fleets = "FLFleetsExt"),
          definition = function(fleets,
                                method = "TMB_DirMissingAR1Hurdle",
                                datayears     = NULL,
                                TACyear       = NULL,
                                deterministic = TRUE,
                                deterministic_yrs = 3,
                                detMethod     = "mean",
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE,
                                removeNonRecent = FALSE) {

            ## make a log file
            if (makeLog == TRUE) {
              logs <- vector(mode = "list", length(fleets))
              names(logs) <- names(fleets)
            }
            if (makePlots) {
              plots <- vector(mode = "list", length(fleets))
              names(plots) <- names(fleets)
            }

            ## loop over fleets
            for(f in names(fleets)) {
              if(verbose)
                cat(f)

              out_f <- uncertainty_effortshare(fleets[[f]],
                                               datayears     = datayears,
                                               TACyear       = TACyear,
                                               deterministic = deterministic,
                                               deterministic_yrs = deterministic_yrs,
                                               detMethod     = detMethod,
                                               verbose       = verbose,
                                               makeLog       = makeLog,
                                               makePlots     = makePlots,
                                               removeNonRecent = removeNonRecent)

              fleets[[f]] <- out_f$fleets
              logs[[f]]   <- out_f$logs
              plots[[f]]  <- out_f$plots

            }

            ## return updated fleets
            return(list(fleets = fleets,
                        logs   = logs,
                        plots  = plots))
          })

## if fleets == "FLFleet"
#' @rdname uncertainty_effortshare
setMethod(f = "uncertainty_effortshare",
          signature = signature(fleets = "FLFleet"),
          definition = function(fleets,
                                method = "TMB_DirMissingAR1Hurdle",
                                datayears     = NULL,
                                TACyear       = NULL,
                                deterministic = TRUE,
                                deterministic_yrs = 3,
                                detMethod     = "mean",
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE,
                                removeNonRecent = FALSE) {

            ## Find number of iterations in object
            nit <- min(dims(fleets)$iter)

            ## Function expects that multiple iterations are present
            if(nit < 2)
              stop("input should have > 1 iterations to store sampled uncertainty")

            if(is.null(datayears))
              stop("argument 'datayears' cannot be NULL")

            if(is.null(TACyear))
              TACyear <- as.integer(tail(datayears,1)) + 1

            if(TACyear > dims(fleets)$maxyear)
              stop("argument 'TACyear' exceeds available years")

            ## If TACyear > datayear + 1, then fill intermediate years too
            fillyear <- (as.integer(tail(datayears,1))+1):TACyear

            # ----------------#
            # Extract data
            # ----------------#

            ## Extract metier effortshare for fleet f
            effshare_f <- sapply(fleets@metiers@names, function(mt) {
              fleets@metiers[[mt]]@effshare[,ac(datayears),,,,1]
            }, simplify = "array")[,,,,,,,drop = TRUE]

            ## if only one metier, coerce to matrix
            if(!is.matrix(effshare_f)) {
              effshare_f <- as.matrix(effshare_f)
              colnames(effshare_f) <- fleets@metiers@names
            }

            # -----------------------------------------------------------------#
            # Handle cases where historical effort-share is consistently zero
            # -----------------------------------------------------------------#

            ## consider consistant historical zero as zero in all iterations
            if (removeNonRecent) {
              effidx <- colSums(tail(effshare_f, max(deterministic_yrs,3)), na.rm = TRUE) > 0
            } else {
              effidx <- colMeans(effshare_f) > 0
            }
            effdata <- effshare_f[, effidx, drop = FALSE]
            effdata <- sweep(effdata, 1, rowSums(effdata, na.rm = TRUE), "/")
            effdata_mt <- colnames(effdata)

            ## Only run if multiple metiers are present
            if (length(effdata_mt) > 1) {

              # --------------------#
              # Handle 0/1 cases
              # --------------------#
              #
              ## The Dirichlet distribution does not like perfect 1/0 proportions.
              ## If a metier == 0: fill 0 metiers with 1e-6 and subtract proportionally
              ## from non-zero metiers.

              if (method == "TMB_Dirrw") {
                if (any(effdata == 0)) {
                  if (verbose)
                    cat(" Metier effort-share = 0: imputing 1e-6 to each zero effort metier.")

                  effdata <- impute_cases(effdata)
                }
              }

              ## convert into dataframe for model
              effdata <- as.data.frame(effdata)

              # ---------------------------#
              # Fit Dirichlet model
              # ---------------------------#

              out <- do.call(method, list(dat      = as.matrix(effdata),
                                          nit      = ifelse(deterministic, nit-1, nit),
                                          fillyear = fillyear,
                                          verbose  = verbose,
                                          makeLog  = makeLog,
                                          makePlots= makePlots))

              if(makeLog) {
                logs  <- out$logs
                plots <- out$plots
              } else {
                logs  <- NULL
                plots <- out$plots
              }

              resvariates <- out$res

              # ----------------------------------#
              # Initial deterministic calculation
              # ----------------------------------#
              #
              # To handle the possibility of a failed forecast, we first fill
              # all slots with a simple mean of historic years.

              for (mt in colnames(effdata)) {
                fleets@metiers[[mt]]@effshare[,ac(fillyear)] <- mean(effdata[tail(1:nrow(effdata), deterministic_yrs),mt],na.rm = TRUE)
              }

              if (is.null(resvariates)) {
                return(list(fleets = fleets,
                            logs   = logs,
                            plots  = plots))
              }

              # ----------------------------------#
              # Insert stochastic forecast
              # ----------------------------------#

              ## Insert data into correct slots
              if (deterministic) {

                ## Loop over each metier
                for(mt in colnames(effdata)) {

                  ## Fill all iterations with deterministic calculations
                  if (detMethod == "mean") {
                    fleets@metiers[[mt]]@effshare[,ac(fillyear),,,,1] <- mean(effdata[tail(1:nrow(effdata), deterministic_yrs),mt],na.rm = TRUE)
                  } else {
                    fleets@metiers[[mt]]@effshare[,ac(fillyear),,,,1] <- apply(exp(resvariates[,mt,]), 1, median, na.rm = TRUE)
                  }

                  ## Overwrite 2+ iterations with stochastically sampled values
                  fleets@metiers[[mt]]@effshare[,ac(fillyear),,,,2:nit] <- resvariates[,mt == colnames(resvariates),]
                } ## END loop over metiers

              } else {

                ## Overwrite all iterations with stochastically sampled values
                for(mt in colnames(resvariates)) {
                  fleets@metiers[[mt]]@effshare[,ac(fillyear)] <- resvariates[,mt == colnames(resvariates),]
                }
              }
            } else { ### If only 1 metier

              ## loop because we might have removed metiers with consistent zero allocation
              for(ii in fleets@metiers@names) {
                fleets@metiers[[ii]]@effshare[,ac(fillyear)] <- 1
              }

              if (verbose)
                cat(" < 2 metier |")

              if(makeLog) {
                logs  <- "< 2 metier"
                plots <- NULL
              } else {
                logs  <- NULL
                plots <- NULL
              }

            } ## End metier if/else statement

            ## return updated fleets
            return(list(fleets = fleets,
                        logs   = logs,
                        plots  = plots))
          })

## if fleets == "FLFleets"
#' @rdname uncertainty_effortshare
setMethod(f = "uncertainty_effortshare",
          signature = signature(fleets = "FLFleets"),
          definition = function(fleets,
                                method = "TMB_DirMissingAR1Hurdle",
                                datayears     = NULL,
                                TACyear       = NULL,
                                deterministic = TRUE,
                                deterministic_yrs = 3,
                                detMethod     = "mean",
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE,
                                removeNonRecent = FALSE) {

            ## Find number of iterations in object
            nit <- min(dims(fleets)$iter)

            ## Function expects that multiple iterations are present
            if(nit < 2)
              stop("input should have > 1 iterations to store sampled uncertainty")

            if(is.null(datayears))
              stop("argument 'datayears' cannot be NULL")

            ## Define TAC year if not provided
            if(is.null(TACyear))
              TACyear <- as.integer(tail(datayears,1)) + 1

            if(TACyear > dims(fleets)$maxyear)
              stop("argument 'TACyear' exceeds available years")

            ## make a log file
            if (makeLog) {
              logs <- vector(mode = "list", length(fleets))
              names(logs) <- names(fleets)
            }
            if (makePlots) {
              plots <- vector(mode = "list", length(fleets))
              names(plots) <- names(fleets)
            }


            ## loop over fleets
            for(f in names(fleets)) {
              if (verbose)
                cat("\n", f)

              out_f <- uncertainty_effortshare(fleets[[f]],
                                               method = method,
                                               datayears     = datayears,
                                               TACyear       = TACyear,
                                               deterministic = deterministic,
                                               deterministic_yrs = deterministic_yrs,
                                               detMethod     = detMethod,
                                               verbose       = verbose,
                                               makeLog       = makeLog,
                                               makePlots     = makePlots,
                                               removeNonRecent = removeNonRecent)

              fleets[[f]] <- out_f$fleets
              logs[[f]]   <- out_f$logs
              plots[[f]]  <- out_f$plots
            }

            ## return updated fleets
            return(list(fleets = fleets,
                        logs   = logs,
                        plots  = plots))
          })

