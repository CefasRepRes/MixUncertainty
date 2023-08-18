# MixUncertainty
An R package to stochastically forecast parameter values for mixed fisheries models from historic observations. Simple state-space time-series models are fitted to historic fleet parameter observations and future values are stochastically projected from the fitted latent process. MixUncertainty allows for the easy conditioning of fleet structures with uncertainty around key parameters:

- metier-stock **catchability**
- metier proportional **effort-share**
- fleet proportional stock **quota-share**

MixUncertainty uses FLR libraries. Fleet data are expected to follow the `FLFleet`(s), `FLMetier`(s), `FLCatch`(es) hierarchy defined by the `FLFleet` package. Uncertainty conditioning methods depend on the `TMB` package. Optional figures are generated using `tidyr`, `plyr`, `ggplot2` and `cowplot`.

# Installation
MixUncertainty can be installed using the 'devtools' package

```{r}
install.packages("devtools")
devtools::install_github("pacematt/MixUncertainty/MixUncertainty")
```

# Examples
The example below uses a simple synthetic data set containing two fleets, three metiers and two stocks. Historic observations of catchability, effort-share and landings numbers are available for years 1 to 10, and estimates of future catchability and effort-share in years 11 and 12 are a simple mean of values for years 8 to 10. 

``` {r}
## load package
library(MixUncertainty)

## load example dataset
data("example_fleets")

## check fleet data dimensions and structure
dim(example_fleets)
```

### Catchability
Fit simple state-space time-series model to 10 years of historic observations. Specify that year 10 is the final data year and forecasts should extend to year 12.

``` {r}
out1 <- uncertainty_catchability(example_fleets, datayear = 10, TACyear = 12, nyrs = 10)

## check model fitting logs for each catch
out1$logs

## check time-series of model fit (maximum likelihood estimate + 95% CI)
out1$plots

# generate figure showing covariance in observations and random draws for two stocks harvested by a give metier
diagnostic_catchability(out1$fleets, nyrs = 10, datayear = 10, TACyear = 12,
                        fl = "I", mt = "1ab", c1 = "A", c2 =  "B")
```
### Metier effort-share
To combine uncertainty around catchability and effort-share conditioning, use the updated fleet object as input.

```{r}
out2 <- uncertainty_effortshare(out1$fleets, datayear = 10, TACyear = 12, nyrs = 10)

## check model fitting logs for each fleet
out2$logs

## check time-series of model fit (maximum likelihood estimate + 95% CI)
out2$plots

## generate figure showing times-series of historic data and stochastic forecast
diagnostic_effortshare(out2$fleets$I, nyrs = 10, datayear = 10, TACyear = 12)
```
### Fleet quota-share
The historic proportional share of stock quota allocated to each fleet is not typically known and mixed fisheries models typically use the share of stock landings as a proxy for quota-share. If quota-share is `NULL`, a time-series model is fitted to the historic proportional share of stock landings.

```{r}
out3 <- uncertainty_quotashare(out2$fleets, quotashare = NULL, datayear = 10, TACyear = 12, nyrs = 10)

## a three-dimensional array (stock, fleet, iteration)
dim(out3$quotashare)

## check model fitting logs for each fleet
out3$logs

## generate figure showing times-series of historic data and stochastic forecast
out3$plots

## generate figure showing times-series of historic data and stochastic forecast
diagnostic_quotashare(fleets = out2$fleets, quotashare = out3$quotashare, stk = "B",
                      nyrs = 10, datayear = 10, TACyear = 12)

```

# Methods
## Catchability

## Metier effort-share

## Fleet quota-share
