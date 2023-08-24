# ---
# title: 'Functions to stochastically sample metier-stock catchability'
# author: 'Various'
# date: 'November 2022'
# ---

#' Stochastically sample future catchability from historic estimates
#'
#' Function takes a \code{FLFleet}, \code{FLFleets}, \code{FLFleetExt} or a
#' \code{FLFleetsExt} object and fits a statistical model to historic
#' observations of metier-stock catchability for each fleet. The catchability
#' for future years is a stochastic forecast from fitted model parameters.
#'
#' @param fleets An object of class \code{FLFleet}, \code{FLFleets},
#'               \code{FLFleetExt} or a \code{FLFleetsExt}.
#' @param method (Character) The algorithm used to fit the statistical model and
#'               generate stochastic forecasts. Using the default is highly
#'               recommended.
#' @param datayear (Integer) The final data year in the fleets object
#' @param TACyear (Integer) The projection year in the fleets object
#' @param nyrs (Integer) The number of historic years used to estimate
#'             future catchability. Defaults to 3. Recommend > 10 if
#'             using a random walk model.
#' @param deterministic (Logical) Should the results for first iteration be
#'                      a deterministic average of the historical period? Defaults
#'                      to \code{TRUE}
#' @param detMethod (Character) The method to calculate the deterministic average.
#'                  If "mean" then a simple mean of the most recent three data years
#'                  is used. Otherwise, the median of the stochastic forecast is
#'                  used.
#' @param verbose (Logical) Should the function print progress? Defaults
#'                to \code{TRUE}
#' @param makeLog (Logical) Should the function generate a log of the methods
#'                applied to each stock (and their success)? Defaults
#'                to \code{TRUE}
#' @param makePlots (Logical) Should the function generate plots showing model
#'                  fit to data for each stock? Defaults to \code{TRUE}
#'
#' @return named list
#'
#' @export

setGeneric("uncertainty_catchability", function(fleets,
                                                method        = "TMB_logMVNrw",
                                                datayear      = NULL,
                                                TACyear       = NULL,
                                                nyrs          = 3,
                                                deterministic = TRUE,
                                                detMethod     = "mean",
                                                verbose       = TRUE,
                                                makeLog       = TRUE,
                                                makePlots     = TRUE,
                                                parallel      = FALSE) {
  standardGeneric("uncertainty_catchability")
})

# NOTE:
#  CURRENTLY DOES NOT WORK WITH HISTORICAL UNCERTAINTY. HOW WOULD WE ESTIMATE
#  FUTURE CATCHABILITY UNCERTAINTY IF THERE IS HISTORICAL CATHABILITY
#  UNCERTAINTY? ONE OPTION MIGHT BE TO SIMPLY RESAMPLE FROM THE HISTORICAL
#  PERIOD RATHER THAN FITTING A DISTRIBUTION. ALTERNATIVELY, WE MIGHT SIMPLY
#  FIT A DISTRIBUTION TO EACH ITERATION AND THEN SAMPLE A SINGLE VALUE FROM THE
#  DISTRIBUTION FOR EACH FUTURE YEAR.

#  CURRENTLY DOES NOT WORK WITH MULTIPLE UNITS, SEASONS OR AREAS.

## if fleets == "FLFleetExt"
#' @rdname uncertainty_catchability
setMethod(f = "uncertainty_catchability",
          signature = signature(fleets = "FLFleetExt"),
          definition = function(fleets,
                                method        = "TMB_logMVNrw",
                                datayear      = NULL,
                                TACyear       = NULL,
                                nyrs          = 3,
                                deterministic = TRUE,
                                detMethod     = "mean",
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE) {

            ## Function expects that multiple iterations are present
            if(dims(fleets[[f]])$iter < 2)
              stop("input should have > 1 iterations to store sampled uncertainty")

            stop("Methods not yet implemented for 'FLFleetExt'")
          })

## if fleets == "FLFleetsExt"
#' @rdname uncertainty_catchability
setMethod(f = "uncertainty_catchability",
          signature = signature(fleets = "FLFleetsExt"),
          definition = function(fleets,
                                method        = "TMB_logMVNrw",
                                datayear      = NULL,
                                TACyear       = NULL,
                                nyrs          = 3,
                                deterministic = TRUE,
                                detMethod     = "mean",
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE,
                                parallel      = FALSE) {

            ## make a log file
            if (makeLog) {
              logs <- vector(mode = "list", length(fleets))
              names(logs) <- names(fleets)
            }
            ## make a plots list
            if (makePlots) {
              plots <- vector(mode = "list", length(fleets))
              names(plots) <- names(fleets)
            }

            if(parallel) {

            } else {

              ## loop over fleets
              for(f in names(fleets)) {
                if(verbose)
                  print(f)

                out_f <- uncertainty_catchability(fleets[[f]],
                                                  method        = method,
                                                  datayear      = datayear,
                                                  TACyear       = TACyear,
                                                  nyrs          = nyrs,
                                                  deterministic = deterministic,
                                                  detMethod     = detMethod,
                                                  verbose       = verbose,
                                                  makeLog       = makeLog,
                                                  makePlots     = makePlots)
                fleets[[f]] <- out_f$fleets
                logs[[f]]   <- out_f$logs
                plots[[f]]  <- out_f$plots
              }
            }

            ## return updated fleets
            return(list(fleets = fleets,
                        logs   = logs,
                        plots  = plots))
          })

## if fleets == "FLFleet"
#' @rdname uncertainty_catchability
setMethod(f = "uncertainty_catchability",
          signature = signature(fleets = "FLFleet"),
          definition = function(fleets,
                                method        = "TMB_logMVNrw",
                                datayear      = NULL,
                                TACyear       = NULL,
                                nyrs          = 3,
                                deterministic = TRUE,
                                detMethod     = "mean",
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE) {

            ## Find number of iterations in object
            nit <- dims(fleets)$iter

            ## Function expects that multiple iterations are present
            if(nit < 2)
              stop("input should have > 1 iterations to store sampled uncertainty")

            if(is.null(datayear))
              stop("argument 'datayear' cannot be NULL")

            if(is.null(TACyear))
              TACyear <- datayear + 1

            if(TACyear > dims(fleets)$maxyear)
              stop("argument 'TACyear' exceeds available years")

            ## make a log file
            if (makeLog) {
              logs <- vector(mode = "list", length(fleets@metiers))
              names(logs) <- names(fleets@metiers)
            }
            ## make a plots list
            if (makePlots) {
              plots <- vector(mode = "list", length(fleets@metiers))
              names(plots) <- names(fleets@metiers)
            }

            ## If TACyear > datayear + 1, then fill intermediate years too
            fillyear <- (datayear+1):TACyear

            ## loop over metiers
            for(mt in fleets@metiers@names) {
              if (verbose)
                cat("\n","-- ",mt)

              ## Extract data for the mt'th metier
              metier_mt <- fleets@metiers[[mt]]

              ## Extract vector of data years
              years <- dimnames(metier_mt@catches[[1]])$year

              ## Define conditioned year
              if(!is.null(nyrs))
                catchabilityyears <- which(sapply(years, function(x) x %in%
                                                    ((datayear-(nyrs-1)):datayear)))

              ## extract catchability for recent years
              qs <- sapply(names(metier_mt@catches),
                           function(ct) metier_mt@catches[[ct]]@catch.q[,,,,,1],
                           simplify = "array")[,,,,,,,drop = TRUE]

              ## if only one stock caught, coerce into matrix
              if(length(names(metier_mt@catches)) == 1) {
                qs <- matrix(qs, ncol = 1, nrow = length(qs))
                colnames(qs) <- names(metier_mt@catches)
                rownames(qs) <- years
              }

              ## extract recent historical period to condition catchability from
              if(!is.null(nyrs))
                qs <- qs[catchabilityyears,, drop=FALSE]  ## conditioned years

              ## Apply method to fit a distribution to historic catchability and generate a forecast
              out_mt <- do.call(method, list(qs,
                                             metier_mt,
                                             fillyear,
                                             deterministic,
                                             detMethod,
                                             verbose,
                                             makeLog,
                                             makePlots))

              ## re-insert updated metier data
              fleets@metiers[[mt]] <- out_mt$res

              ## extract logfile
              logs[[mt]]  <- out_mt$logs
              plots[[mt]] <- out_mt$plots

            } # END loop over metiers
            return(list(fleets = fleets,
                        logs   = logs,
                        plots  = plots))
          })

## if fleets == "FLFleets"
#' @rdname uncertainty_catchability
setMethod(f = "uncertainty_catchability",
          signature = signature(fleets = "FLFleets"),
          definition = function(fleets,
                                method        = "TMB_logMVNrw",
                                datayear      = NULL,
                                TACyear       = NULL,
                                nyrs          = 3,
                                deterministic = TRUE,
                                detMethod     = "mean",
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE,
                                parallel      = FALSE) {

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

            ## make a log file
            if (makeLog) {
              logs <- vector(mode = "list", length(fleets))
              names(logs) <- names(fleets)
            }
            if (makePlots) {
              plots <- vector(mode = "list", length(fleets))
              names(plots) <- names(fleets)
            }

            if(parallel == TRUE | is.numeric(parallel)) {

              ## set up parallel environment
              cl <- beginParallel(parallel)

              ## stop cluster on exit
              on.exit({
                parallel::stopCluster(cl)
                foreach::registerDoSEQ()
              }, add = TRUE)

              out <- foreach(fleet_f = fleets,
                             .errorhandling = "pass") %dopar% {

                uncertainty_catchability(fleet_f,
                                         method        = method,
                                         datayear      = datayear,
                                         TACyear       = TACyear,
                                         nyrs          = nyrs,
                                         deterministic = deterministic,
                                         detMethod     = detMethod,
                                         verbose       = verbose,
                                         makeLog       = makeLog,
                                         makePlots     = makePlots)
              }


              fleets <- lapply(out, "[[", "fleets")
              logs   <- lapply(out, "[[", "logs")
              plots  <- lapply(out, "[[", "plots")


            } else {

              ## loop over fleets
              for(f in names(fleets)) {
                if (verbose)
                  cat("\n", f)

                out_f <- uncertainty_catchability(fleets[[f]],
                                                  method        = method,
                                                  datayear      = datayear,
                                                  TACyear       = TACyear,
                                                  nyrs          = nyrs,
                                                  deterministic = deterministic,
                                                  detMethod     = detMethod,
                                                  verbose       = verbose,
                                                  makeLog       = makeLog,
                                                  makePlots     = makePlots)

                fleets[[f]] <- out_f$fleets
                logs[[f]]   <- out_f$logs
                plots[[f]]  <- out_f$plots
              }
            }

            ## return updated fleets
            return(list(fleets = fleets,
                        logs   = logs,
                        plots  = plots))
          })
