% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DETMfull.R
\name{DETMfull}
\alias{DETMfull}
\title{DETM fitting with full parameters}
\usage{
DETMfull(
  posdata,
  negdata,
  pival,
  alp1,
  beta1,
  alp2,
  beta2,
  maxite = 1000,
  eps = 1e-04
)
}
\arguments{
\item{posdata}{n by p feature matrix for the source (positive) data in PU data.}

\item{negdata}{m by p feature matrix for the target (unlabled) data in PU data.}

\item{pival}{initial value of the mixture proportion of positive data in the target dataset, is set for the EM algorithm.}

\item{alp1}{initial value of alpha_1 of DETM in EM algorithm}

\item{beta1}{A vector of initial value of beta_1 of DETM in EM algorithm}

\item{alp2}{initial value of alpha_2 of DETM in EM algorithm}

\item{beta2}{A vector of initial value of beta_2 of DETM in EM algorithm}

\item{maxite}{maximum iteration number of EM algorithm}

\item{eps}{tolerance value about stopping the algorithm when the increment in the log-EL after an iteration is no greater than, e.g. 1e-4.}
}
\value{
A list including mixture proportion estimation, parameter estimation, iterations, likelihood and penalized likelihood in DETM.
}
\description{
DETM fitting with full parameters
}
