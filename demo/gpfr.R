rm(list=ls())
library(GPFDA)

set.seed(100)
traindata <- vector('list',20)
for(i in 1:20) traindata[[i]]=i
n <- 500
traindata <- lapply(traindata,function(i) {
  x <- seq(-3,3,len=n)
  y <- sin(x^2)-x+0.2*rnorm(n,0,3)
  x1 <- 0.5*x^3+exp(x)+rnorm(n,0,3)
  x2 <- cos(x^3)+0.2*rnorm(n,0,3)
  mat <- cbind(x,x1,x2,y)
  colnames(mat) <- c('time','x1','x2','y')
  scale <- t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
  i <- list(mat,scale)
})

n <- 800 #test input
x <- seq(-3,3,len=n)
y <- sin(x^2)-x+0.2*rnorm(n,0,3)
x1 <- 0.5*x^3+exp(x)+rnorm(n,0,3)
x2 <- cos(x^3)+0.2*rnorm(n,0,3)
mat <- cbind(x,x1,x2,y)
colnames(mat) <- c('time','x1','x2','y')
scale <- t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
# testdata[[1]]=vector('list',3)
n <- 100 # test new points
xt <- seq(1,3,len=n)
yt <- sin(xt^2)-xt+0.2*rnorm(n,0,3)
xt1 <- 0.5*xt^3+exp(xt)+rnorm(n,0,3)
xt2 <- cos(xt^3)+0.2*rnorm(n,0,3)
mat_t <- cbind(xt,xt1,xt2)
colnames(mat_t) <- c('time','xt1','xt2')
td <- list(mat,scale,mat_t)


trdata <- wrap(functional=lapply(traindata,function(i)i[[1]]),
     do.call('rbind',lapply(traindata,function(i)i[[2]])),
     list='traindata',time='time',response='y')
tedata <- wrap(functional=td[[1]],scale=td[[2]],testdata=td[[3]],
     list='testdata',time='time',response='y')

library(ggplot2)
qplot(sample,input,data=trdata[trdata$type=='functional'&trdata$col=='response'& trdata$batch==1,],col=batch)


a<-gpfr(trdata)
result_type1<-gpfrpred(a,tedata)
result_type2<-gpfrpred(a,newtime=xt,data.new=mat_t[,-1],type=2)
# result<-gpfrpred(a,tedata,type=2)

plot(-1000,col=0,xlim=range(result_type1$time),ylim=range(result_type1$ypred),xlab='time',ylab='prediction',
     main='Prediction by GPFR: type I')
lines(result_type1$time,result_type1$ypred[,1])
lines(result_type1$time,result_type1$ypred[,2],lty=2,col=2)
lines(result_type1$time,result_type1$ypred[,3],lty=2,col=2)
points(result_type1$time,yt)

plot(-1000,col=0,xlim=range(result_type2$time),ylim=range(result_type2$ypred),xlab='time',ylab='prediction',
     main='Prediction by GPFR: type II')
lines(result_type2$time,result_type2$ypred[,1])
lines(result_type2$time,result_type2$ypred[,2],lty=2,col=2)
lines(result_type2$time,result_type2$ypred[,3],lty=2,col=2)
points(result_type2$time,yt)

fda_train <- fdatrain(trdata)
fda_result <- fdapred(yregfd=fda_train$betaestlist,sigv2fd=fda_train$sigv2fd,tedata,trdata)

plot(-1000,col=0,xlim=range(fda_result$time),ylim=range(fda_result$ypred),xlab='time',ylab='prediction',main='Prediction by FDR')
lines(fda_result$time,fda_result$ypred[,1])
lines(fda_result$time,fda_result$ypred[,2],lty=2)
lines(fda_result$time,fda_result$ypred[,3],lty=2)



