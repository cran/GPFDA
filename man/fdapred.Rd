\name{fdapred}
\alias{fdapred}
\title{
Predict using result from learning the training data.
}
\description{
Use fuctional data analysis to do the prediction using result from learning training data.
}
\usage{
fdapred(yregfd, sigv2fd, data, TrainData)
}

\arguments{
  \item{yregfd}{
Object from "fdatrain". 
}
  \item{sigv2fd}{
Object from "fdatrain".
}
  \item{data}{
Test data. Must be a object from "wrap" or a data.table with certain format.
}
  \item{TrainData}{
Train data. Must be a object from "wrap" or a data.table with certain format.
}
}
\value{
A list of 
\item{ypred}{The matrix of predicted value and confidence interval}
\item{rcc}{Rooted mean squared error and correlation of predicted value and test data}
\item{time}{The time scale that the curves lie on}
\item{s2}{covariance structure}
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}
\examples{
library(GPFDA)

traindata=vector('list',50)
for(i in 1:50) traindata[[i]]=i
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


# fda_train=fdatrain(trdata)
# fda_result=fdapred(yregfd=fda_train$betaestlist,sigv2fd=fda_train$sigv2fd,tedata,trdata)
}
