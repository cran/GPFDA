\name{D2}
\alias{D2}

\title{
Second derivative of the likelihood
}
\description{
Computer the second derivative of the likelihood function with first and second derivative of the kernel function given.
}
\usage{
D2(d1, d2, inv.Q, Alpha.Q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d1}{
First derivative of the kernel function.
}
  \item{d2}{
Second derivative of the kernel function.
}
  \item{inv.Q}{
Inverse matrix of the covariance matrix
}
  \item{Alpha.Q}{
\deqn{Q^{-1}\%*\%Y\%*\%(Q^{-1}\%*\%Y)^T-Q^{-1}} where Q is the covariance matrix, Y is the response.
}
}
\value{
\item{out}{A number}
}
\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}