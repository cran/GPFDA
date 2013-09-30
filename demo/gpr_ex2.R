library(GPFDA)
require(MASS)

set.seed(60)
hp <- list('pow.ex.w'=rep(log(10),4),'linear.a'=rep(log(10),4),'pow.ex.v'=log(5),'vv'=log(1))
kernal <- c('linear','pow.ex')

nn=100; # training sample size
mm=1000; # testing sample size

p=dnorm((rnorm(800)))
idx=sort(sample(1:800,nn,prob=p/sum(p)))
X=matrix(0,ncol=4,nrow=800)
X[,1]=seq(-5,10,len=800)
X[,2]=seq(0,1,len=800)
X[,3]=seq(-15,-10,len=800)
X[,4]=seq(1,2,len=800)
# Y=mvrnorm(n=200,mu=as.matrix(X[,1]-X[,1]),Sigma=(cov.linear(hp,X)+cov.pow.ex(hp,X)))[1,]+0.2*sign(X[,1])*abs(X[,1])^(1/3)-4*sin(X[,2])
Y=(mvrnorm(n=800,mu=as.matrix(X[,1]-X[,1]),Sigma=(cov.linear(hp,X)+cov.pow.ex(hp,X)))[1,])+
  (0.2*sign(X[,1])*abs(X[,1])^(1/3)-4*sin(X[,2])+exp(X[,3])+log(X[,4]))*3
X=X[idx,];Y=as.matrix(Y[idx])
x=matrix(0,ncol=4,nrow=mm)
x[,1]=seq(-5,10,len=mm)
x[,2]=seq(0,1,len=mm)
x[,3]=seq(-15,-10,len=mm)
x[,4]=seq(1,2,len=mm)

a=gpr(X,Y,kernal,hp)
b=gppredict(a,x)

upper=b$mu+1.96*b$sigma
lower=b$mu-1.96*b$sigma
plot(-100,-100,col=0,xlim=range(x[,1]),ylim=c(min(upper,lower,Y)-0.1*abs(min(upper,lower,Y)),max(upper,lower,Y)+0.1*abs(max(upper,lower,Y))),main="Prediction", xlab="input ( x )",ylab="responce")
polygon(c(x[,1], rev(x[,1])), c(upper, rev(lower)),col = "grey90", border = NA)
points(X[,1],Y,pch=4,col=2)
#points(x[,1],a$mu)
lines(X[,1],Y)
lines(x[,1],b$mu,col=3,lwd=2)

