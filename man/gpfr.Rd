\name{gpfr}
\alias{gpfr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gaussian Process in functional data.
}
\description{
Use functional regression to be the mean structure and Gaussian Process to be the covariance structure.
}

\usage{
gpfr(TrainData,CovType = c("linear", "pow.ex"), hyper=NULL,NewHyper=NULL, gamma = 1,
  nbasis = NULL, norder = 6, lambda1 = 1e-07, lambda2 = 1e-05, pen.order = 2, 
  new.sample = NULL,accuracy=c('high','normal','low'),trace.iter=5,fitting=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TrainData}{
The train data. Must from "wrap" or with certain data.table format.
}
  \item{CovType}{
Kernel function or covariance function type(s).
}
  \item{hyper}{
Hyper parameter initial value. Default to be NULL.
}
  \item{NewHyper}{
Vector of the names of the new hyper parameters from customized kernel function.
}
  \item{gamma}{
Power parameter used in powered exponential.
}
  \item{nbasis}{
Number of basis to use in the functional regression. Default to be the integer of dimension of train data over 5.
}
  \item{norder}{
Order of basis function in functional regression. Default to be 4.
}
  \item{lambda1}{
Smoothing parameter for curve fitting. Default to be 1e-07
}
  \item{lambda2}{
Smoothing parameter for functional regression. Default to be 1e-05
}
  \item{pen.order}{
Order of penalty term. Default to be 2.
}
  \item{new.sample}{
The dimension of the new sampled matrix (if applicable). Default to be the maximum dimension in the data.
}
  \item{accuracy}{
Optimization accuracy. Default to be high. 
}
  \item{trace.iter}{
Print the processing of iterations of optimization. 
}
  \item{fitting}{
Is fitting required or not. Default to be F.
}
}

\value{
A list of
\item{hyper}{estimated hyper-parameters}
\item{I}{A vector of estimated standard deviation of hyper-parameters}
\item{betaestlist}{estimated functional coefficient}
\item{iuu}{scale part of the data}
\item{CovFun}{Covariance function}
}
\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}

\seealso{
  \code{\link{gpr}}
}

\examples{
library(GPFDA)

traindata=vector('list',20)
for(i in 1:20) traindata[[i]]=i
n=500
traindata=lapply(traindata,function(i) {
  x=seq(-3,3,len=n)
  y=sin(x^2)-x+0.2*rnorm(n,runif(1,-3,3),runif(1,0.5,3))
  x1=0.5*x^3+exp(x)+rnorm(n,runif(1,-3,3),runif(1,0.5,5))
  x2=cos(x^3)+0.2*rnorm(n,runif(1,-3,3),runif(1,0.5,5))
  mat=cbind(x,x1,x2,y)
  colnames(mat)=c('time','x1','x2','y')
  scale=t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
  i=list(mat,scale)
})

n=300 #test input
x=seq(-3,3,len=n)
y=sin(x^2)-x+0.2*rnorm(n,runif(1,-3,3),runif(1,0.5,3))
x1=0.5*x^3+exp(x)+rnorm(n,runif(1,-3,3),runif(1,0.5,5))
x2=cos(x^3)+0.2*rnorm(n,runif(1,-3,3),runif(1,0.5,5))
mat=cbind(x,x1,x2,y)
colnames(mat)=c('time','x1','x2','y')
scale=t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
# testdata[[1]]=vector('list',3)
n=100 # test new points
xt=seq(1,3,len=n)
yt=sin(xt^2)-xt+0.2*rnorm(n,runif(1,-3,3),runif(1,0.5,3))
xt1=0.5*xt^3+exp(xt)+rnorm(n,runif(1,-3,3),runif(1,0.5,5))
xt2=cos(xt^3)+0.2*rnorm(n,runif(1,-3,3),runif(1,0.5,5))
mat_t=cbind(xt,xt1,xt2,yt)
colnames(mat_t)=c('time','xt1','xt2','y')
td=list(mat,scale,mat_t)


trdata=wrap(functional=lapply(traindata,function(i)i[[1]]),
     do.call('rbind',lapply(traindata,function(i)i[[2]])),
     list='traindata',time='time',response='y')
tedata=wrap(functional=td[[1]],scale=td[[2]],testdata=td[[3]],
     list='testdata',time='time',response='y')

# system.time(result<-gpfr(trdata,tedata))
}
