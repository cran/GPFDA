\name{gpfrpred}
\alias{gpfrpred}
\title{Prediction of the Gaussian Process using functional regression}
\description{
Predict the new points in Gaussian Process using the training results or manual input}

\usage{
gpfrpred(object,data=NULL,newtime=NULL,data.new = NULL,type=1,yregfd=NULL, 
  hyper.p=NULL, iuu=NULL, Cov=NULL,gamma=1)
}
\arguments{
  \item{object}{
The result from training. Default to be NULL, if not NULL, other arguments (except for Data.new) will replaced by NULL.}
  \item{data}{
The test data. Must have the format from wrap. For type 1 prediction. If `newtime' is missing, but type is 2, time in data will be used for type 2 prediction.
}
  \item{newtime}{
The new time point to predict. For type 2 prediction.}
  \item{data.new}{
The test data for predict. For type 2 prediction.}
  \item{type}{
Prediction types. See details below.
}
  \item{yregfd}{
Functional coefficients. 
}
  \item{hyper.p}{
Hyper-parameter estimated from training. Can use manual input. Default to be NULL.
}
  \item{iuu}{
Scaler variables for prediction. Always required for both type of predictions.
}
  \item{Cov}{
Names of covariance functions used. Default to be NULL. 
}
  \item{gamma}{
Parameter used in power exponential covariance function. Default to be NULL. 
}
}
\details{
Two types of prediction are supplied. Type one is the new batch has a few observations, type two is the new batch has no observations.}
\value{
\item{ypred}{matrix of predicted value with confidence interval}
\item{time}{new time points}
\item{y_gppred}{predicted value}
\item{object}{all items trained from gpfr if exists}
}
\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}

\seealso{
  \code{\link{gpr}}
}

