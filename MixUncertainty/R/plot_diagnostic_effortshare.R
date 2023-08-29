# ---
# title: 'Functions to generate diagnostic plots for sampled metier-stock effort-share'
# author: 'Various'
# date: 'April 2023'
# ---

setGeneric("diagnostic_effortshare", function(fleets,
                                              datayears = NULL,
                                              TACyear   = NULL,
                                              fl        = NULL) {
  standardGeneric("diagnostic_effortshare")
})

## if fleets == "FLFleetExt"
#' @rdname diagnostic_effortshare
setMethod(f = "diagnostic_effortshare",
          signature = signature(fleets = "FLFleetExt"),
          definition = function(fleets,
                                datayears = NULL,
                                TACyear   = NULL) {

          })


## if fleets == "FLFleetsExt"
#' @rdname diagnostic_effortshare
setMethod(f = "diagnostic_effortshare",
          signature = signature(fleets = "FLFleetsExt"),
          definition = function(fleets,
                                datayears = NULL,
                                TACyear   = NULL,
                                fl        = NULL) {

          })

## if fleets == "FLFleet"
#' @rdname diagnostic_effortshare
setMethod(f = "diagnostic_effortshare",
          signature = signature(fleets = "FLFleet"),
          definition = function(fleets,
                                datayears = NULL,
                                TACyear   = NULL) {

            require(dplyr)

            ## Extract metier effortshare for fleet f
            effshare_f <- sapply(fleets@metiers@names, function(mt) {
              fleets@metiers[[mt]]@effshare[,ac(datayears)]
            }, simplify = "array")[,,,,,,,drop = TRUE]

            ## Convert into dataframe
            effdata <- as.data.frame.table(effshare_f)
            colnames(effdata) <- c("year", "iter", "metier", "data")

            effdata %>%
              group_by(year, metier) %>%
              summarise(p50 = median(data),
                        p05 = quantile(data, p = 0.05),
                        p25 = quantile(data, p = 0.25),
                        p75 = quantile(data, p = 0.75),
                        p95 = quantile(data, p = 0.95)) %>%
              mutate(year = as.numeric(as.character(year))) %>%
              ggplot2::ggplot() +
              geom_line(aes(x = year, y = p50)) +
              geom_ribbon(aes(x = year, ymin = p25, ymax = p75), alpha = 0.20) +
              geom_ribbon(aes(x = year, ymin = p05, ymax = p95), alpha = 0.10) +
              geom_point(aes(x = year, y = data),
                         data = data.frame(data = sapply(fleets@metiers@names,
                                                         function(mt) {
                                                           c(fleets@metiers[[mt]]@effshare[,ac(TACyear),,,,1])

                                                         }),
                                           metier = fleets@metiers@names, year = TACyear)) +
              geom_vline(aes(xintercept = as.integer(tail(datayears,1))), linetype = 2) +
              scale_y_continuous("Metier effort-share") +
              facet_wrap(~metier, scales = "free_y") +
              ggtitle(fleets@name) +
              theme_bw()

          })

## if fleets == "FLFleetsExt"
#' @rdname diagnostic_effortshare
setMethod(f = "diagnostic_effortshare",
          signature = signature(fleets = "FLFleets"),
          definition = function(fleets,
                                datayears = NULL,
                                TACyear   = NULL,
                                fl        = NULL) {

            if (is.null(fl))
              stop("Fleet name 'fl' must be specified")

            if (!(fl %in% names(fleets)))
              stop("'fl' must be the name of a fleet in 'fleets'")

            ## subset for specific fleet
            fleets_fl <- fleets[[fl]]

            ## run plotting function
            diagnostic_effortshare(fleets    = fleets_fl,
                                   datayears = datayears,
                                   TACyear   = TACyear)

          })
