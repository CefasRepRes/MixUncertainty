# ---
# title: 'Functions to generate diagnostic plots for sampled fleet stock quota-share'
# author: 'Various'
# date: 'April 2023'
# ---

setGeneric("diagnostic_quotashare", function(fleets,
                                             quotashare,
                                             nyrs     = 4,
                                             datayear = NULL,
                                             TACyear  = NULL,
                                             stk      = NULL) {
  standardGeneric("diagnostic_quotashare")
})

## if fleets == "FLFleetsExt"
#' @rdname diagnostic_quotashare
setMethod(f = "diagnostic_quotashare",
          signature = signature(fleets = "FLFleetsExt"),
          definition = function(fleets,
                                quotashare,
                                nyrs     = 4,
                                datayear = NULL,
                                TACyear  = NULL,
                                stk      = NULL) {

          })

## if fleets == "FLFleets"
#' @rdname diagnostic_quotashare
setMethod(f = "diagnostic_quotashare",
          signature = signature(fleets = "FLFleets"),
          definition = function(fleets,
                                quotashare,
                                nyrs     = 4,
                                datayear = NULL,
                                TACyear  = NULL,
                                stk      = NULL) {

            require(dplyr)

            # ----------------------------------#
            # Process estimated quota-share data
            # ----------------------------------#

            ## Convert array to dataframe
            quotashare <- as.data.frame.table(quotashare)

            ## add variable names
            colnames(quotashare) <- c("stock", "fleet", "iter", "quotashare")

            ## filter for stock
            quotashare <- quotashare[quotashare$stock == stk,]

            ## Add TAC year
            quotashare$year <- TACyear

            # ----------------------------#
            # Extract Fleet Landings-share
            # ----------------------------#

            fleetcatches <- calculate_landingshare(fleets = fleets,
                                                   datayear = datayear,
                                                   nyrs = nyrs)

            ## reformat fleetcatches to match quotashare
            fleetcatches$iter <- 1
            fleetcatches <- fleetcatches[fleetcatches$stock == stk,]
            colnames(fleetcatches)[colnames(fleetcatches) == "landingshare"] <- "quotashare"

            ## consider missing historical quotashare as zero
            fleetcatches <- fleetcatches %>%
              tidyr::spread(fleet, quotashare) %>%
              tidyr::gather(fleet, quotashare, -year, -stock, -iter) %>%
              dplyr::mutate(quotashare = ifelse(is.na(quotashare),0,quotashare))

            # -----#
            # Plot
            # -----#

            ## summarise data and plot
            rbind(fleetcatches, quotashare) %>%
              dplyr::group_by(year, fleet) %>%
              dplyr::summarise(p50 = median(quotashare, na.rm = TRUE),
                        p05 = quantile(quotashare, p = 0.05, na.rm = TRUE),
                        p25 = quantile(quotashare, p = 0.25, na.rm = TRUE),
                        p75 = quantile(quotashare, p = 0.75, na.rm = TRUE),
                        p95 = quantile(quotashare, p = 0.95, na.rm = TRUE)) %>%
              dplyr::mutate(year = as.numeric(year)) %>%
              ggplot() +
              geom_line(aes(x = year, y = p50)) +
              geom_ribbon(aes(x = year, ymin = p25, ymax = p75), alpha = 0.20) +
              geom_ribbon(aes(x = year, ymin = p05, ymax = p95), alpha = 0.10) +
              geom_point(aes(x = year, y = quotashare),
                         data = quotashare[quotashare$iter == 1,]) +
              geom_vline(aes(xintercept = datayear), linetype = 2) +
              scale_y_continuous("Fleet landings-share") +
              facet_wrap(~fleet, scales = "free_y") +
              ggtitle(stk) +
              theme_bw()
          })
