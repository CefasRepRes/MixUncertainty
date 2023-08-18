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

# Overview
## Catchability

## Metier effort-share

## Fleet quota-share
