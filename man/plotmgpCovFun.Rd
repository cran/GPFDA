% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mgp.functions.R
\name{plotmgpCovFun}
\alias{plotmgpCovFun}
\title{Plot auto- or cross-covariance function of a multivariate Gaussian process}
\usage{
plotmgpCovFun(
  type = "Cov",
  output,
  outputp,
  Data,
  hp,
  idx,
  ylim = NULL,
  xlim = NULL,
  mar = c(4.5, 5.1, 2.2, 0.8),
  oma = c(0, 0, 0, 0),
  cex.lab = 1.5,
  cex.axis = 1,
  cex.main = 1.5
)
}
\arguments{
\item{type}{Logical. It can be either 'Cov' (for covariance function) or
'Cor' (for corresponding correlation function).}

\item{output}{Integer identifying one element of the multivariate process.}

\item{outputp}{Integer identifying one element of the multivariate process.
If 'output' and 'outputp' are the same, the auto-covariance function will
be plotted. Otherwise, the cross-covariance function between 'output' and
'outputp' will be plotted.}

\item{Data}{List of two elements: 'input' and 'response'. The element 'input'
is a list of N vectors, where each vector represents the input covariate
values for a particular output. The element 'response' is the corresponding
list of N matrices (if there are multiple realisations) or vectors (for a
single realisation) representing the response variables.}

\item{hp}{Vector of hyperparameters}

\item{idx}{Index vector identifying to which output the elements of
concatenated vectors correspond to.}

\item{ylim}{Graphical parameter}

\item{xlim}{Graphical parameter}

\item{mar}{Graphical parameter passed to par().}

\item{oma}{Graphical parameter passed to par().}

\item{cex.lab}{Graphical parameter passed to par().}

\item{cex.axis}{Graphical parameter passed to par().}

\item{cex.main}{Graphical parameter passed to par().}
}
\value{
A plot
}
\description{
Plot auto- or cross-covariance function of a multivariate Gaussian process
}
\examples{
## See examples in vignette:
# vignette("mgpr", package = "GPFDA")
}
