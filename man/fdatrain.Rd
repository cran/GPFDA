\name{fdatrain}
\alias{fdatrain}

\title{
Learn the training data using functional data analysis.
}
\description{
Use B-spline basis to fit the curves and use functional regression with functional response as the model to learn the data.
}
\usage{
fdatrain(data, nbasis = NULL, norder = 6, lambda1 = 1e-07, lambda2 = 1e-05, pen.order = 2)
}

\arguments{
  \item{data}{
Training data. Should be the data after wrap or with certain data.table format.
}
  \item{nbasis}{
Number of basis functions. Default is the integer of dimension of data over 5.
}
  \item{norder}{
Number of norder. Default to be 4.
}
  \item{lambda1}{
Smoothing parameter for curve fitting. Default to be 1e-07
}
  \item{lambda2}{
Smoothing paramter for functional regression. Default to be 1e-05
}
  \item{pen.order}{
Order of penalty term. Default to be 2.
}
}
\details{
Use penalized fitting and regression with large number of basis functions to avoid some computational costs. 
}
\value{
A list of
\item{betaestlist}{the list of coefficients basis}
\item{sigv2fd}{Functional object of the variance structure}
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}
\seealso{
create.bspline.basis, fdapred
}



