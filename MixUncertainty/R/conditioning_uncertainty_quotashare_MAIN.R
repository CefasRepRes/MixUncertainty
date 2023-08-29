# ---
# title: 'Functions to stochastically sample fleet quota-share'
# author: 'Various'
# date: 'November 2022'
# ---

#' Stochastically sample future fleet quota-share from historic estimates
#'
#' Function takes a \code{FLFleet}, \code{FLFleets}, \code{FLFleetExt} or a
#' \code{FLFleetsExt} object and fits a statistical model to historic
#' observations of stock quota-share for each fleet. The stock quota-share
#' for future years is a stochastic forecast from fitted model parameters.
#'
#' @param fleets An object of class \code{FLFleet}, \code{FLFleets},
#'               \code{FLFleetExt} or a \code{FLFleetsExt}.
#' @param method (Character) The algorithm used to fit the statistical model and
#'               generate stochastic forecasts. Using the default is highly
#'               recommended.
#' @param quotashare (Optional) A data frame containing: year, stock, fleet
#'                   and landingshare.
#' @param datayears A vector of historic years used to estimate future stock
#'                  quota-share
#' @param TACyear (Integer) The projection year in the fleets object
#' @param deterministic (Logical) Should the results for first iteration be
#'                      a simple mean of the historical period? Defaults
#'                      to \code{TRUE}
#' @param deterministic_yrs (Integer) The number of recent data years to use in
#'                          the deterministic calculation. Defaults to 3.
#' @param verbose (Logical) Should the function print progress? Defaults
#'                to \code{TRUE}
#' @param makeLog (Logical) Should the function record model fitting and forecasting
#'                performance for each fleet? Defaults to \code{TRUE}
#' @param makePlots (Logical) Should the function generate a figure showing
#'                  observations and model fit and forecasts with uncertainty
#'                  intervals for each fleet? Defaults to \code{TRUE}
#'
#' @return
#'
#' @export

setGeneric("uncertainty_quotashare", function(fleets,
                                              method = "TMB_DirMissingAR1Hurdle",
                                              quotashare    = NULL,
                                              datayears     = NULL,
                                              TACyear       = NULL,
                                              deterministic = TRUE,
                                              deterministic_yrs = 3,
                                              verbose       = TRUE,
                                              makeLog       = TRUE,
                                              makePlots     = TRUE) {
  standardGeneric("uncertainty_quotashare")
})

## if fleets == "FLFleetsExt"
#' @rdname uncertainty_quotashare
setMethod(f = "uncertainty_quotashare",
          signature = signature(fleets = "FLFleetsExt"),
          definition = function(fleets,
                                method = "TMB_DirMissingAR1Hurdle",
                                quotashare    = NULL,
                                datayears     = NULL,
                                TACyear       = NULL,
                                deterministic = TRUE,
                                deterministic_yrs = 3,
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE) {

            stop("Methods not yet implemented for 'FLFleetExt'")

          })

## if fleets == "FLFleets"
#' @rdname uncertainty_quotashare
setMethod(f = "uncertainty_quotashare",
          signature = signature(fleets = "FLFleets"),
          definition = function(fleets,
                                method = "TMB_DirMissingAR1Hurdle",
                                quotashare    = NULL,
                                datayears     = NULL,
                                TACyear       = NULL,
                                deterministic = TRUE,
                                deterministic_yrs = 3,
                                verbose       = TRUE,
                                makeLog       = TRUE,
                                makePlots     = TRUE) {


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

            # ----------------------------#
            # Extract Fleet Landings-share
            # ----------------------------#

            if (is.null(quotashare)) {
              fleetcatches <- calculate_landingshare(fleets, datayears)
            } else {
              fleetcatches <- quotashare
            }

            if(any(fleetcatches$landingshare < 0))
              stop("< 0 values detected in observation data")

            ## make a log file
            if (makeLog) {
              logs <- vector(mode = "list", length(unique(fleetcatches$stock)))
              names(logs) <- unique(fleetcatches$stock)
            }

            ## make a plots list
            if (makePlots) {
              plots <- vector(mode = "list", length(unique(fleetcatches$stock)))
              names(plots) <- unique(fleetcatches$stock)
            }

            # -----------------------------------------#
            # Model future landing-share & uncertainty
            # -----------------------------------------#
            #
            # In this section, we loop over each stock and predict stochastically the
            # future landings share of each fleet.

            ## extract fleet landings as done in FCube
            sl. <- eval(parse("",text="landings"))

            fleets_land <-lapply(fleets, function (x) {
              mt. <- lapply(x@metiers, function(x1) {
                res. <- do.call(rbind,lapply(names(sl.(x1)), function(n) {
                  m <- as.data.frame.table(sl.(x1)[[n]])
                  m$qname <- n
                  return(m)
                }))
                names(res.)[which(names(res.)=="Freq")]   <- "landings"
                res.$fleet  <- x@name
                res.$metier <- x1@name
                return(res.)})
              mt. <- eval(parse(text=paste('rbind(mt.[[',paste(seq(length(mt.)),collapse=']],mt.[['),']])',sep='')))
            })
            fleets_land <- eval(parse(text=paste('rbind(fleets_land[[',
                                                 paste(seq(length(fleets_land)),
                                                       collapse=']],fleets_land[['),
                                                 ']])',
                                                 sep='')))

            ## filter for data years
            fleets_land <- fleets_land[fleets_land$year %in% tail(datayears, deterministic_yrs),]
            quotashare  <-  tapply(fleets_land$landings,
                                 list(fleets_land$qname,
                                      fleets_land$fleet,
                                      fleets_land$iter),
                                 sum,
                                 na.rm=T)
            if (nit>1 && all(is.na(quotashare[,,2:nit])))  quotashare[,,2:nit] <- quotashare[,,1]

            quotashare[is.na(quotashare)] <- 0
            quotashare <- sweep(quotashare,c(1,3),apply(quotashare,c(1,3),sum),FUN="/")

            quotashare <- quotashare[rownames(quotashare[sort(c(unique(fleetcatches$stock))),,]),,]

            ## loop over stocks
            for(stk in unique(fleetcatches$stock)) {
              if (verbose)
                cat("\n", stk)

              # At this point, we need to define two methods. One method where we
              # condition future quota based on the historic share of landings,
              # and another where the distribution among countries is defined
              # by FIDES and the distribution among country-level fleets is
              # based on historic landings share.

              # ---------------------------#
              # Extract stock data
              # ---------------------------#

              ## select stock
              stkcatches <- fleetcatches[fleetcatches$stock %in% stk,
                                         c("year","fleet","landingshare")]

              ## spread fleet across columns
              stkcatches <- tidyr::spread(data  = stkcatches,
                                          key   = fleet,
                                          value = landingshare)

              ## remove fleets that catch zero stock across timeseries
              stkcatches      <- stkcatches[,colnames(stkcatches) != "year", drop = FALSE]
              landsharefleets <- which(colSums(stkcatches, na.rm = TRUE) > 0)
              landsharedata   <- stkcatches[,landsharefleets, drop = FALSE]

              ## Consider NAs to be zeros
              landsharedata[is.na(landsharedata)] <- 0

              ## If more than one fleet catches this stock...
              if (ncol(landsharedata) > 1) {

                ## Check if there are any cases of zero landings share
                if (any(landsharedata == 0) & method == "TMB_Dirrw") {
                  cat(" Fleet landings-share = 0: imputing 1e-6 to each zero effort fleet")
                  landsharedata <- as.data.frame(impute_cases(as.matrix(landsharedata)))
                }

                # ---------------------------#
                # Fit Dirichlet model
                # ---------------------------#

                ## use datayears as rownames
                rownames(landsharedata) <- datayears

                ## fit model
                out <- do.call(method, list(dat      = as.matrix(landsharedata),
                                            nit      = ifelse(deterministic, nit-1, nit),
                                            fillyear = fillyear,
                                            verbose  = verbose,
                                            makeLog  = makeLog,
                                            makePlots= makePlots))

                if(makeLog) {
                  logs[[stk]]  <- out$logs
                  plots[[stk]] <- out$plots
                } else {
                  logs[[stk]]  <- NULL
                  plots[[stk]] <- out$plots
                }

                resvariates <- out$res


                ## Insert data into correct slots
                ## Remember! 1st iteration should be untouched for compatability
                ## with deterministic conditioning
                for(f in colnames(resvariates)) {
                  if(deterministic) {
                    quotashare[stk,f,-1] <- resvariates[1,colnames(resvariates) == f,]
                  } else {
                    quotashare[stk,f,] <- resvariates[1,colnames(resvariates) == f,]
                  }
                }

              } else { # if only one fleet catches this stock...

                quotashare[stk,landsharefleets, ] <- 1

              } # END if/else > 1 fleet catches stock
            } # END loop over stocks

            ## return updated object
            return(list(quotashare = quotashare,
                        logs       = logs,
                        plots      = plots))
          })


#' Calculate proportional fleet share of stock landings
#'
#' Function takes an \code{FLFleets} object and calculates the annual share of
#' overall stock landings for each fleet. It is assumed that there is no
#' stochasticity in landings estimates, and values are drawn from the first
#' iteration of each landings \code{FLQuant}.
#'
#' @param fleets An \code{FLFleets} object
#' @param datayears A vector of historic data years to include in calculations
#'
#' @return A data frame with variables: 'year', 'stock', 'fleet' and 'landingshare'

calculate_landingshare <- function(fleets, datayears) {

  # ------------------------------#
  # Calculate fleet stock landings
  # ------------------------------#

  fleetcatches <- lapply(names(fleets), function(ff) {

    ## extract metiers
    fl <- fleets[[ff]]@metiers
    fl_catches <- lapply(names(fl), function(mm) {

      ## extract catches
      mt <- fl[[mm]]

      ## extract landings for each stock as an array
      mt_catches <- sapply(names(mt@catches), function(cc) {
        mt@catches[[cc]]@landings[,,,,,1]
      }, simplify = "array", USE.NAMES = TRUE)

      ## convert array to dataframe
      mt_catches <- as.data.frame.table(mt_catches)[,c("year","Var7","Freq")]
      colnames(mt_catches) <- c("year","stock","data")
      mt_catches$year <- as.integer(as.character(mt_catches$year))
      mt_catches$stock <- as.character(mt_catches$stock)

      return(mt_catches)
    })

    fl_catches <- do.call(rbind, fl_catches)
    fl_catches$fleet <- ff
    return(fl_catches)
  })

  # -------------------------------------#
  # Calculate proportional stock landings
  # -------------------------------------#

  ## calculate stock landings by fleet
  fleetcatches <- do.call(rbind, fleetcatches)
  fleetcatches <- aggregate(data ~ year + stock + fleet,
                            data = fleetcatches,
                            FUN = sum)

  ## calculate proportional stock landings by fleet
  totalcatches <- aggregate(data ~ year + stock,
                            data = fleetcatches,
                            FUN = sum)
  fleetcatches <- merge(fleetcatches, totalcatches, by = c("year", "stock"))
  fleetcatches$landingshare <- fleetcatches$data.x/fleetcatches$data.y

  ## select recent years
  fleetcatches <- fleetcatches[fleetcatches$year %in% ac(datayears),]

  ## remove unnecessary columns
  fleetcatches <- fleetcatches[,!(colnames(fleetcatches) %in% c("data.x", "data.y"))]

  return(fleetcatches)
}
