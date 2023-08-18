# MixUncertainty
An R package to stochastically forecast parameter values for mixed fisheries models from historic observations. Simple state-space time-series models are fitted to historic fishing fleet parameter observations and future values are stochastically projected from the fitted latent process. MixUncertainty allows for the easy conditioning of fleet structures with uncertainty around three key parameters:

- catchability of each stock for each metier (**catchability**)
- proportional share of fleet effort across metiers (**effort-share**)
- proportional share of stock quota across fleets (**quota-share**)

MixUncertainty uses FLR libraries. Fleet data are expected to follow the `FLFleet`(s), `FLMetier`(s), `FLCatch`(es) hierarchy defined by the `FLFleet` package. Uncertainty conditioning methods depend on the `TMB` package. Optional figures are generated using `tidyr`, `plyr`, `ggplot2` and `cowplot`.

## Installation
MixUncertainty can be installed using the 'devtools' package

```{r}
install.packages("devtools")
devtools::install_github("pacematt/MixUncertainty/MixUncertainty")
```

## Application
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

## generate figure showing covariance in observations and random draws for two stocks
## harvested by a given metier
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

## Methods
### Catchability
Within a given metier, catchability might be expected to co-vary among stocks, reflecting how changes in gear efficiency might have correlated impacts on stocks with similar morphologies and life-histories. The underlying vector of "true" catchabilities $\mathbf{q}^{\ast}\_t$ at time $t$, where $\mathbf{q_t} = (q_{1,t}, ..., q_{s,t})$ for $s$ exploited stocks, follows an AR1 process on a log-scale with multivariate normal distributed noise $\eta$:

$$\log \mathbf{q}^{\ast}\_t = \mathbf{\mu} + \phi \left(\log \mathbf{q}^{\ast}\_{t-1} - \mathbf{\mu} \right) + \rho \mathbf{\eta_t}, \text{where } \mathbf{\eta_t} \sim \mathbf{N}(0, \Sigma),$$

$\Sigma$ is the variance-covariance matrix for the multivariate normal distribution, $\mu$ is the stationary mean for each catchability element, $\phi$ is the AR1 correlation parameter, and $\rho = \sqrt(1 - \phi^2)$. The observation error $\epsilon_{i,t}$ is assumed to take a univariate normal distribution:

$$\log q_t = \log q^{\ast}\_{i,t} + \epsilon_{i,t}, \text{where } \epsilon_{i,t} \sim N(0, \sigma_i)$$

In cases of model convergence issues, the model is simplified by fixing covariances to zero and therefore relaxing the assumption of multivariate dependency in process noise.

### Metier effort-share and fleet quota-share
Metier effort-share and fleet quota-share are both compositional time-series. The data at each time increment are a vector $\mathbf{y}\_t = (y_{1,t}, \dots, y_{K,t})$, where $y_{1,t} + \dots + y_{K,t} = 1$ and $y_{i,t} \in \[0,1\]$. Two different models are adopted, depending on whether the data contain zero observations. 

If all observations are non-zero, then the multinomial logit transformed expectation $E(\mathbf{y)\_t})$ is assumed to follow a random walk with multivariate normal distributed increments $\nu$:

$$\text{logit~}\text{E}(\mathbf{y}\_t) = \text{logit~}\text{E}(\mathbf{y}\_{t-1}) + \mathbf{\nu}\_t, \text{where } \mathbf{\nu}\_t \sim \mathbf{N}(0, \Sigma),$$

$\Sigma$ is the variance-covariance matrix for the multivariate normal distribution. Observations are modelled using the Dirichlet distribution

$$\mathbf{y}\_t = \text{Dir}(\mathbf{\alpha}_t), \text{where } \mathbf{\alpha}_t = E(\mathbf{y}_t \cdot \tau)$$

$\tau$ is the time-invariant concentration parameter.

One short-coming is that the Dirichlet distribution cannot accommodate zero observations. Hence, when data contain zero observations, a hurdle approach is used to model the observation process and an AR1 process is used to model the expectation.

Similar to models for catchability, in cases of convergence issues, the models are simplified by fixing covariances to zero and therefore relaxing the assumption of multivariate dependency in process noise.
 
