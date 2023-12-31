% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/conditioning_uncertainty_methods_TMB_DirMissingAR1HurdleRevised.R
\name{TMB_DirMissingAR1Hurdle}
\alias{TMB_DirMissingAR1Hurdle}
\title{Forecast catchability from a state-space time-series model with Dirichlet distributed observation noise}
\usage{
TMB_DirMissingAR1Hurdle(dat, nit, fillyear, verbose, makeLog, makePlots)
}
\arguments{
\item{dat}{A matrix of historic effortshare where rows = years and cols = metiers}

\item{nit}{(Integer) the number of replicate parameter trajectories to sample}

\item{fillyear}{Vector of intermediate and forecast years to be sampled for.}

\item{verbose}{(Logical) Should the function print progress?}

\item{makeLog}{(Logical) Should the function generate a log of the methods
applied to each variable (and their success)?}

\item{makePlots}{(Logical) Should the function generate plots showing model
fit to data for each stock?}
}
\value{
An matrix of sampled parameter values
}
\description{
Function takes a time-series of metier effort-share data as a matrix and
fits a state-space model with Dirichlet distributed observation error. The
specific model structure will depend on the observation data. The final data
year is simulated to generate a distribution based on observation uncertainty.
A forecast is carried out for each simulated catchability value that propagates
uncertainty estimated from the random walk.
}
