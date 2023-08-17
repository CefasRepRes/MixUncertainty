# ---
# title: 'Functions to generate diagnostic plots for sampled metier-stock catchability'
# author: 'Various'
# date: 'April 2023'
# ---

setGeneric("diagnostic_catchability", function(fleets,
                                               nyrs     = 3,
                                               datayear = NULL,
                                               TACyear  = NULL,
                                               fl       = NULL,
                                               mt       = NULL,
                                               c1       = NULL,
                                               c2       = NULL) {
  standardGeneric("diagnostic_catchability")
})

## if fleets == "FLFleetExt"
#' @rdname diagnostic_catchability
setMethod(f = "diagnostic_catchability",
          signature = signature(fleets = "FLFleetExt"),
          definition = function(fleets,
                                nyrs     = 3,
                                datayear = NULL,
                                TACyear  = NULL,
                                mt       = NULL,
                                c1       = NULL,
                                c2       = NULL) {

          })


## if fleets == "FLFleetsExt"
#' @rdname diagnostic_catchability
setMethod(f = "diagnostic_catchability",
          signature = signature(fleets = "FLFleetsExt"),
          definition = function(fleets,
                                nyrs     = 3,
                                datayear = NULL,
                                TACyear  = NULL,
                                fl   = NULL,
                                mt   = NULL,
                                c1   = NULL,
                                c2   = NULL) {

          })

## if fleets == "FLFleet"
#' @rdname diagnostic_catchability
setMethod(f = "diagnostic_catchability",
          signature = signature(fleets = "FLFleet"),
          definition = function(fleets,
                                nyrs     = 3,
                                datayear = NULL,
                                TACyear  = NULL,
                                mt       = NULL,
                                c1       = NULL,
                                c2       = NULL) {

            ## Extract vector of data years
            years <- dimnames(fleets@metiers[[mt]]@catches[[1]])$year

            ## Define conditioned year
            catchabilityyears <- which(sapply(years, function(x) x %in% ((datayear-(nyrs-1)):datayear)))

            ## extract sampled catchability
            plot_df <- data.frame(c1 = c(fleets@metiers[[mt]]@catches[[c1]]@catch.q[,ac(TACyear)]),
                                  c2 = c(fleets@metiers[[mt]]@catches[[c2]]@catch.q[,ac(TACyear)]))
            plot_df <- setNames(plot_df, c(c1,c2))

            ## extract historic catchability
            qs_df <- data.frame(c1 = c(fleets@metiers[[mt]]@catches[[c1]]@catch.q[,catchabilityyears]),
                                c2 = c(fleets@metiers[[mt]]@catches[[c2]]@catch.q[,catchabilityyears]))
            qs_df <- setNames(qs_df, c(c1,c2))

            qs_mean_df <- data.frame(c1 = mean(fleets@metiers[[mt]]@catches[[c1]]@catch.q[,catchabilityyears]),
                                     c2 = mean(fleets@metiers[[mt]]@catches[[c2]]@catch.q[,catchabilityyears]))
            qs_mean_df <- setNames(qs_mean_df, c(c1,c2))

            p1 <- ggplot2::ggplot(plot_df, aes(x = get(c1), y = get(c2))) +
              geom_point(aes(colour = "random draw"),alpha = .5) +
              geom_density_2d() +
              geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
              geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
              geom_point(data = qs_df, aes(x = get(c1), y = get(c2), colour = "data"), size = 3) +
              geom_point(data = qs_mean_df, aes(x = get(c1), y = get(c2), colour = "data mean"), size = 3) +
              xlab(c1) + ylab(c2) + ggtitle("multivariate random draws") +
              theme_classic() +
              theme(legend.position = "none")

            p2 <- ggplotFL::plot(fleets@metiers[[mt]]@catches[[c1]]@catch.q,
                                 probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) + ggtitle(c1)

            p3 <- ggplotFL::plot(fleets@metiers[[mt]]@catches[[c2]]@catch.q,
                                 probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) + ggtitle(c2)

            p4 <- cowplot::plot_grid(p2,p3, ncol = 2)

            legend1 <- cowplot::get_legend(p1 + theme(legend.position = "right"))
            legend2 <- NULL

            biv    <- cowplot::plot_grid(p1,p4, nrow = 2, rel_widths = c(2,1,1))
            legend <- cowplot::plot_grid(legend1, legend2, nrow = 2, rel_widths = c(1,1))

            combined <- cowplot::plot_grid(biv, legend, ncol = 2, rel_widths = c(3,1))
            return(combined)

          })


## if fleets == "FLFleets"
#' @rdname diagnostic_catchability
setMethod(f = "diagnostic_catchability",
          signature = signature(fleets = "FLFleets"),
          definition = function(fleets,
                                nyrs     = 3,
                                datayear = NULL,
                                TACyear  = NULL,
                                fl   = NULL,
                                mt   = NULL,
                                c1   = NULL,
                                c2   = NULL) {

            if (is.null(fl))
              stop("Fleet name 'fl' must be specified")

            if (c1 == c1)
              stop("catch names 'c1' and 'c2' must be different")

            if (!(fl %in% names(fleets)))
              stop("'fl' must be the name of a fleet in 'fleets'")

            ## subset for specific fleet
            fleets_fl <- fleets[[fl]]

            ## run plotting function
            diagnostic_catchability(fleets = fleets_fl,
                                    datayear = datayear,
                                    TACyear = TACyear,
                                    nyrs = nyrs,
                                    mt = mt,
                                    c1 = c1,
                                    c2 = c2)
          })
