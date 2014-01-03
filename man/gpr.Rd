\name{gpr}
\alias{gpr}

\title{
Gaussian Process regression for single curve
}
\description{
Gaussian Process regression for single curve with train data.
}
\usage{
gpr(Data, response, Cov, hyper=NULL, NewHyper=NULL, mean=0, gamma=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
The input data from train data.
}
  \item{hyper}{
The hyper parameters. Must be a list.
}
  \item{NewHyper}{
Vector of the names of the new hyper parameters from customized kernel function.
}
  \item{mean}{
Is the mean taken out when analysis? Default to be 0, which assumes the mean is zero. if assume mean is a constant, mean=1; if assume mean is a linear trend, mean='t'.
}
  \item{response}{
The response data from train data.
}
  \item{Cov}{
The kernel functions or covariance functions to use.
}
  \item{gamma}{
Power parameter used in powered exponential kernel function.
}
}
\value{
\item{CovFun}{Covariance function type}
\item{fitted}{Fitted value of training data}
\item{fitted.sd}{Standard deviation of the fitted value of training data}
\item{gamma}{Parameter used in powered exponential covariance function}
\item{hyper}{Hyper-parameter estimated from training data}
\item{I}{Variance of the estimated hyper-parameters}
\item{train.x}{Training covariates}
\item{train.y}{Training response}
\item{Q}{Covariance matrix}
\item{inv}{Inverse of the covariance matrix}
\item{mean}{The mean assumed in the analysis}
\item{lrm}{'lm' object if mean is a linear regression. NULL otherwise.}
}
\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}

\seealso{
  \code{\link{gpfr}}
}
\examples{
library(GPFDA)
library(MASS) ## used to generate data
hp <- list('pow.ex.w'=log(10),'linear.a'=log(10),'pow.ex.v'=log(5),'vv'=log(1))
c <- seq(0,1,len=40)
idx <- sort(sample(1:40,21))
X <- as.matrix(c[idx])
Y <- (mvrnorm(n=40,mu=c-c,Sigma=(cov.linear(hp,c)+cov.pow.ex(hp,c)))[1,])+sin(c*6)
Y <- as.matrix(Y[idx])
x <- as.matrix(seq(0,1,by=0.03))
a <- gpr(X,Y,c('linear','pow.ex'))
}
