\name{cov.pow.ex}
\alias{cov.pow.ex}

\title{
Covariance function. Powered exponential covariance function.
}
\description{
Stationary covariance function. 
}
\usage{
cov.pow.ex(hyper, Data, Data.new = NULL, gamma = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hyper}{
The hyper parameters. Must be a list with certain names.
}
  \item{Data}{
The input data. Must be a vector or a matrix.
}
  \item{Data.new}{
The data for prediction. Must be a vector or a matrix. Default to be NULL.
}
  \item{gamma}{
Power parameter that cannot be estimated by simple non-linear optimization.
}
}
\value{
\item{cov.pow.ex}{Covariance matrix}
}
\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}

\seealso{
  \code{\link{cov.linear}},
  \code{\link{cov.rat.qu}}
}

\examples{
library(GPFDA)
require(MASS)

set.seed(30)
hp <- list('pow.ex.w'=log(10),'linear.a'=log(10),'pow.ex.v'=log(5),'vv'=log(1))
c <- seq(0,1,len=40)
idx <- sort(sample(1:40,21))
X <- as.matrix(c[idx])
Y <- (mvrnorm(n=40,mu=c-c,Sigma=(cov.linear(hp,c)+cov.pow.ex(hp,c)))[1,])+sin(c*6)
Y <- as.matrix(Y[idx])
x <- as.matrix(seq(0,1,by=0.03))
a <- gpr(X,Y,c('pow.ex'),hp)
b <- gppredict(a,x)

upper=b$mu+1.96*b$sigma
lower=b$mu-1.96*b$sigma
plot(-100,-100,col=0,xlim=range(x[,1]),ylim=c(min(upper,lower,Y)-0.1*abs(min(upper,lower,Y)),
  max(upper,lower,Y)+0.1*abs(max(upper,lower,Y))),main="Prediction", 
  xlab="input ( x )",ylab="response")
polygon(c(x[,1], rev(x[,1])), c(upper, rev(lower)),col = "grey90", border = NA)
points(X[,1],Y,pch=4,col=2)
#points(x[,1],a$mu)
lines(X[,1],Y)
lines(x[,1],b$mu,col=3,lwd=2)
}

