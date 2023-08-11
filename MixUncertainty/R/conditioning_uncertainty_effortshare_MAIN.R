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
#' @param datayear (Integer) The final data year in the fleets object
#' @param TACyear (Integer) The projection year in the fleets object
#' @param nyrs (Integer) The number of historic years used to estimate
#'             future effort-share
#' @param deterministic (Logical) Should the results for first iteration be
#'                      a simple mean of the historical period? Defaults
#'                      to \code{TRUE}
#' @param verbose (Logical) Should the function print progress? Defaults
#'                to \code{TRUE}
#' @param makeLog (Logical) Should the function record model fitting and forecasting
#'                performance for each fleet? Defaults to \code{TRUE}
#' @param makePlots (Logical) Should the function generate a figure showing
#'                  observations and model fit and forecasts with uncertainty
#'                  intervals for each fleet? Defaults to \code{TRUE}
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
                                               method = "TMB_Dirrw",
                                               datayear      = NULL,
                                               TACyear       = NULL,
                                               nyrs          = 4,
                                               deterministic = TRUE,
                                               verbose       = TRUE,
                                               makeLog       = TRUE,
                                               makePlots     = TRUE) {
  standardGeneric("uncertainty_effortshare")
})

## if fleets == "FLFleetExt"
#' @rdname uncertainty_effortshare
setMethod(f = "uncertainty_effortshare",
          signature = signature(fleets = "FLFleetExt"),
          definition = function(fleets,
                                method = "TMB_Dirrw",
                                datayear      = NULL,
                                TACyear       = NULL,
                                nyrs          = 4,
                                deterministic = TRUE,
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE) {

            ## Function expects that multiple iterations are present
            if(dims(fleets[[f]])$iter < 2)
              stop("input should have > 1 iterations to store sampled uncertainty")



          })

## if fleets == "FLFleetsExt"
#' @rdname uncertainty_effortshare
setMethod(f = "uncertainty_effortshare",
          signature = signature(fleets = "FLFleetsExt"),
          definition = function(fleets,
                                method = "TMB_Dirrw",
                                datayear      = NULL,
                                TACyear       = NULL,
                                nyrs          = 4,
                                deterministic = TRUE,
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE) {

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
                                               datayear      = datayear,
                                               TACyear       = TACyear,
                                               nyrs          = nyrs,
                                               deterministic = deterministic,
                                               verbose       = verbose,
                                               makeLog       = makeLog,
                                               makePlots     = makePlots)

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
                                method = "TMB_Dirrw",
                                datayear      = NULL,
                                TACyear       = NULL,
                                nyrs          = 4,
                                deterministic = TRUE,
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE) {

            ## Find number of iterations in object
            nit <- min(dims(fleets)$iter)

            ## Function expects that multiple iterations are present
            if(nit < 2)
              stop("input should have > 1 iterations to store sampled uncertainty")

            if(is.null(datayear))
              stop("argument 'datayear' cannot be NULL")

            if(is.null(TACyear))
              TACyear <- datayear + 1

            if(TACyear > dims(fleets)$maxyear)
              stop("argument 'TACyear' exceeds available years")

            ## If TACyear > datayear + 1, then fill intermediate years too
            fillyear <- (datayear+1):TACyear

            ## Only run if multiple metiers are present
            if (length(names(fleets@metiers)) > 1) {

              # ----------------#
              # Extract data
              # ----------------#

              ## Extract metier effortshare for fleet f
              effshare_f <- sapply(fleets@metiers@names, function(mt) {
                fleets@metiers[[mt]]@effshare[,ac((datayear-(nyrs-1)):datayear),,,,1]
              }, simplify = "array")[,,,,,,,drop = TRUE]

              # -----------------------------------------------------------------#
              # Handle cases where historical effort-share is consistently zero
              # -----------------------------------------------------------------#

              ## consider consistant historical zero as zero in all iterations
              effdata <- effshare_f[, colMeans(effshare_f) > 0]
              effdata_mt <- colnames(effdata)

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

              ## Insert data into correct slots
              ## Remember! 1st iteration should be untouched for compatability
              ## with deterministic conditioning
              for(mt in colnames(resvariates)) {
                fleets@metiers[[mt]]@effshare[,ac(fillyear),,,,2:nit] <- resvariates[,mt == colnames(resvariates),]
              }

            } else { ### If only 1 metier

              fleets@metiers[[1]]@effshare[,ac(fillyear)] <- 1

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
                                method = "TMB_Dirrw",
                                datayear      = NULL,
                                TACyear       = NULL,
                                nyrs          = 4,
                                deterministic = TRUE,
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE) {

            ## Find number of iterations in object
            nit <- min(dims(fleets)$iter)

            ## Function expects that multiple iterations are present
            if(nit < 2)
              stop("input should have > 1 iterations to store sampled uncertainty")

            if(is.null(datayear))
              stop("argument 'datayear' cannot be NULL")

            ## Define TAC year if not provided
            if(is.null(TACyear))
              TACyear <- datayear + 1

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
                                               datayear      = datayear,
                                               TACyear       = TACyear,
                                               nyrs          = nyrs,
                                               deterministic = deterministic,
                                               verbose       = verbose,
                                               makeLog       = makeLog,
                                               makePlots     = makePlots)

              fleets[[f]] <- out_f$fleets
              logs[[f]]   <- out_f$logs
              plots[[f]]  <- out_f$plots
            }

            ## return updated fleets
            return(list(fleets = fleets,
                        logs   = logs,
                        plots  = plots))
          })

