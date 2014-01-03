if(getRversion() >= "3.0")  utils::globalVariables(c("J"))
#### Main function ####
gpfr=function(TrainData,CovType=c('linear','pow.ex'),hyper=NULL,NewHyper=NULL,gamma=1,
              nbasis=NULL,norder=6,lambda1=1e-7,lambda2=1e-5,pen.order=2,
              new.sample=NULL,accuracy=c('high','normal','low'),trace.iter=5,
              fitting=F){
  col.no=length(unique(TrainData$col[TrainData$type=='functional']))-2
  if(is.null(hyper)){
    hyper=list()
    if(any(CovType=='linear'))
      hyper$linear.a=rnorm(col.no)
    if(any(CovType=='pow.ex')){
      hyper$pow.ex.v=runif(1,-1,1)
      hyper$pow.ex.w=(-abs(rnorm(col.no)))
    }
    if(any(CovType=='rat.qu')){
      hyper$rat.qu.w=rnorm(col.no)
      hyper$rat.qu.s=runif(1,0.01,1)
      hyper$rat.qu.a=runif(1,0.01,1)
    }
    hyper$vv=runif(1,-1,-0.01)
    hyper.nam=names(hyper)
    
    if(!is.null(NewHyper)){
      hyper.nam=c(hyper.nam,NewHyper)
      nh.length=length(NewHyper)
      for(i in 1:nh.length){
        hyper=c(hyper,runif(1,-1,1))
      }
      names(hyper)=hyper.nam
    }
  }
  
  
  
  TrainData=wrap2(TrainData,new.sample=new.sample)
  a1=gpfrtrain(TrainData,hyper,CovType,gamma,nbasis,norder,lambda1,lambda2,pen.order,accuracy,trace.iter,fitting)
  return(a1)
}


#### fda train ####
fdatrain=function(data,nbasis=NULL,norder=6,lambda1=1e-7,lambda2=1e-5,pen.order=2){
  if(is.null(nbasis)) nbasis=as.integer(length(data[J('functional',1,'time')]$input)/5)
  nRng=range(data[J('functional',1,'time')]$input)
  nbatch=max(data[J('functional')]$batch)
  trainBasis = create.bspline.basis(nRng, nbasis, norder)
  
  D2fdPar = fdPar(trainBasis, lambda=lambda1)
  t=replicate(nbatch,data[J('functional',1,'time')]$input)
  
  ys=reshape(data[J('functional',unique(data$batch),'response')],v.names='input',timevar='batch',
             idvar=c('sample'),direction='wide')
  trainfd = smooth.basis(as.matrix(data[J('functional',1,'time')]$input), 
                         as.matrix(ys[,!names(ys)%in%c('type','col','sample'),with=F]), D2fdPar)$fd 
  
  trainBetaBasis = create.bspline.basis(nRng, nbasis)
  
  trainBetaPar = fdPar(trainBetaBasis, pen.order, lambda2)
  
  p=unique(data[J('scale')]$col)
  np = length(p)
  xfdlist <- vector("list",np)
  for (j in 1:np) xfdlist[[j]] <- data[J('scale',unique(data$batch),p[j])]$input
  
  betalist=vector('list',np)
  for (j in 1:np) betalist[[j]] <- trainBetaPar
  
  fRegressList <- fRegress(trainfd, xfdlist, betalist)
  
  betaestlist <- fRegressList$betaestlist #retunr varlue 1
  yhatfdobj   <- fRegressList$yhatfdobj 
  
  yhatmat    <- predict(yhatfdobj, t[,1])
  ymat       <- eval.fd(t[,1], trainfd)
  resmat <- ymat - yhatmat
  SigmaE     <- var(t(resmat))
  
  sigv2 = apply((resmat)^2,1,sum)/(nbatch-1);
  sigv2fd = smooth.basis(t[,1],sigv2, trainBetaBasis); # return value 2
  
  return(list('betaestlist'=betaestlist,'sigv2fd'=sigv2fd))
}

#### fdapred ####
fdapred=function(yregfd, sigv2fd, data,TrainData){
  zmat=reshape(TrainData[J('scale')],v.names='input',timevar='col',idvar=c('sample','batch'),direction='wide')
  zmat=data.matrix(zmat[,!(names(zmat) %in% c('sample','batch','type')),with=F])
  
  iuu=mymatrix2(t(zmat)%*%zmat)$res
  ttest = data[J('testdata',1,'time')]$input
  yregtmp=unlist(lapply(yregfd,function(i) eval.fd(i$fd,ttest)))
  yregtmp=matrix(yregtmp,ncol=length(yregfd));#yyyy<<-yregtmp
  Ut=reshape(data[J('scale')],v.names='input',timevar='col',idvar='batch',direction='wide')
  Ut=t(data.matrix(Ut[,!names(Ut)%in%c('type','batch','sample'),with=F])); # uutt<<-Ut
  ypredfda = yregtmp%*%Ut
  
  
  sigv2 = eval.fd(sigv2fd$fd,ttest)
  s2 = sigv2%*%(1 + t(Ut)%*%iuu%*%Ut)
  ypredfdaup = ypredfda + 1.96*sqrt(s2);
  ypredfdalo = ypredfda - 1.96*sqrt(s2);
  CI=cbind(ypredfda,ypredfdaup,ypredfdalo)

  return(list(ypred=CI, time=ttest,s2=s2))
}

#### gpfrtrain #####
gpfrtrain=function(data,hyper=NULL,Cov,gamma=1,nbasis=NULL,norder=6,lambda1=1e-7,lambda2=1e-5,pen.order=2,accuracy=c('high','normal','low'),trace.iter=5,fitting=F){
  fRegList = fdatrain(data,nbasis=nbasis,norder=norder,lambda1=lambda1,lambda2=lambda2,pen.order=pen.order)
  
  stepsize = length(data[J('functional',1,'time')]$input)  #training data size
  tdsize = as.integer(stepsize/2)  #choose half data for training
  
  gptraindata = data[data$type=='functional'&data$sample%in%1:tdsize]
  bbat=unique(data[J('functional')]$batch)
  for(i in seq_along(bbat)){
    tdata = data[J('functional',bbat[i],'time')]$input
    yregtmp = do.call('cbind',lapply(fRegList$betaestlist,function(i) predict(i,tdata)));#yyyy<<-yregtmp
    U=reshape(data[J('scale')],v.names='input',timevar='col',idvar=c('sample','batch'),direction='wide',)
    U=U[,!names(U)%in%c('type','batch','sample'),with=F]
    meanvec = yregtmp%*%t(U[i,])
    
    stepsize = length(tdata);#sts<<-stepsize
    xindtrain=sort(sample(stepsize,tdsize))
    
    tempmat = data[J('functional',bbat[i])]
    
    tempmat[J('functional',bbat[i],'response')]$input= 
      tempmat[J('functional',bbat[i],'response')]$input - meanvec[,1]
    tempmat=tempmat[tempmat$sample%in%xindtrain]$input
    gptraindata[J('functional',bbat[i])]$input=tempmat
  }
  
  init0 = unlist(hyper)
  accuracy=accuracy[1]
  if(accuracy=='high') acc=1e-10
  if(accuracy=='normal') acc=1e-6
  if(accuracy=='low') acc=1e-2
  
  optm.idx=1:as.integer(max(gptraindata$sample)/4)*3
  init1=nlminb(init0,repgp.loglikelihood,repgp.Dloglikelihood,Data=gptraindata[gptraindata$sample%in%optm.idx],Cov=Cov,gamma=gamma,control=list(iter.max=5,rel.tol=1e-2))[[1]]
  
  pp=nlminb(init1,repgp.loglikelihood,repgp.Dloglikelihood,Data=gptraindata,Cov=Cov,gamma=gamma,control=list(trace=trace.iter,rel.tol=acc))#,Xprior=prior,Xprior2=NA)
  cat('optimization done','\n','\n')
  pp.cg=pp[[1]]
  names(pp.cg)=names(init0)
  pp.df=data.frame(pp.cg=pp.cg,pp.N=substr(names(init0),1,8))
  names(pp.df)=c('pp.cg','pp.N')
  pp.cg=split(pp.df$pp.cg,pp.df$pp.N)
  
  x.nam=unique(gptraindata[J('functional')]$col)[!unique(gptraindata[J('functional')]$col)%in%c('time','response')]
  da=gptraindata[gptraindata$col%in%x.nam & gptraindata$type=='functional']
  Data=reshape(da,v.names='input',timevar='batch',idvar=c('sample','col'),direction='wide',drop='type')
  Data$col=NULL;Data$sample=NULL
  Y=gptraindata[gptraindata$col=='response' & gptraindata$type=='functional']
  Y=reshape(Y,v.names='input',timevar='batch',idvar=c('sample'),direction='wide',drop='type')
  Y$col=NULL;Y$sample=NULL
  allbat=vector('list',length=length(bbat))
  for(i in 1:length(bbat)){
    X=matrix(data.matrix(Data)[,i],ncol=length(x.nam))
    allbat[[i]]=cbind(data.matrix(Y)[,i],X)
  }
  Qlist=lapply(allbat,function(l) fisherinfo(pp.cg=pp.cg,X=l[,-1],Y=l[,1],Cov=Cov,gamma=gamma))
  II=abs(-1/apply(do.call('cbind',Qlist),1,sum))
  cat("fisher's information done",'\n','\n')
  
  zmat=reshape(data[J('scale')],v.names='input',timevar='col',idvar=c('sample','batch'),direction='wide')
  zmat=data.matrix(zmat[,!(names(zmat) %in% c('sample','batch','type')),with=F])
  
  iuu=mymatrix2(t(zmat)%*%zmat)$res
  
  fitted=fitted.sd=NULL
  
  yregfd=fRegList[[1]]
  if(fitting==T){
    fitted=matrix(0,ncol=length(bbat),nrow=length(data[J('functional',bbat[1],'time')]$input))
    fitted.sd=fitted
    
    for(i in seq_along(bbat)){
      tdata = data[J('functional',bbat[i],'time')]$input
      yregtmp=do.call('cbind',lapply(yregfd,function(m) predict(object=m,newdata=tdata)))
      
      nam=unique(data[J('functional')]$col);x.nam=nam[which(!nam%in%c('response','time'))]
      xinput = reshape(data[J('functional',bbat[i],x.nam)],v.names='input',timevar='col',idvar=c('sample','batch'),direction='wide')
      setnames(xinput,names(xinput),c('type','batch','sample',x.nam))
      xinput = data.matrix(xinput[,names(xinput)%in%x.nam,with=F])
      
      yinput = data[J('functional',bbat[i],'response')]$input 
      U=as.matrix(data[J('scale',bbat[i])]$input)
      ygpobs = yinput-yregtmp%*%U
      y_gpr=gppredict(train=F,Data.new=tdata,Data=as.matrix(xinput),Y=as.matrix(ygpobs), hyper=pp.cg, Cov=Cov,gamma=gamma)
      
      ygppred = y_gpr$fitted
      s2 = y_gpr$fitted.sd^2%*%(1 + t(U)%*%iuu%*%U)
      ypred = yregtmp%*%U + ygppred ## fd rgression plus gp regression
      fitted[,i]=ypred
      fitted.sd[,i]=sqrt(s2)
      if(i%%5==0) cat('fitting',i,' th curve','\n')
    }
  }
  
  result=list('hyper'=pp.cg,'I'=II, 'betaestlist'=fRegList[[1]],'iuu'=iuu,'CovFun'=Cov,'gamma'=gamma,
              'fitted'=fitted,'fitted.sd'=fitted.sd,'bat'=bbat,'train'=data)
  class(result)='gpfda'
  return(result)
}

fisherinfo=function(pp.cg,X,Y,Cov,gamma){
  n=length(Cov)
  CovList=vector('list',n)
  for(i in 1:n) CovList[i]=list(paste0('cov.',Cov[i]))
  CovL=lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex')
      return(f(pp.cg,X,X,gamma=gamma))
    if(j!='cov.pow.ex')
      return(f(pp.cg,X,X))
  }  )
  if(length(CovL)==1)
    Q=CovL[[1]]
  if(length(CovL)>1)
    Q=Reduce('+',CovL)
  
  response=as.matrix(Y)
  X=as.matrix(X)
  Q=Q+diag(exp(pp.cg$vv),dim(Q)[1])
  QR=mymatrix2(Q,response)$res
  invQ=mymatrix2(Q)$res
  AlphaQ=QR%*%t(QR)-invQ
  
  D2=function(d1,d2,inv.Q,Alpha.Q){
    Aii=t(d1)%*%inv.Q%*%d1
    al=Alpha.Q+inv.Q
    return(0.5*(sum(Alpha.Q*(d2-Aii))-sum(al*Aii)))
  }
  
  D2fx=lapply(seq_along(pp.cg),function(i){
    Dp=pp.cg[i]
    name.Dp=names(Dp)
    f=get(paste0('D2',name.Dp))
    if(name.Dp%in%c('pow.ex.w','pow.ex.v') )
      D2para=f(pp.cg,X,gamma=gamma,inv.Q=invQ,Alpha.Q=AlphaQ)
    if(!name.Dp%in%c('pow.ex.w','pow.ex.v'))
      D2para=f(pp.cg,X,inv.Q=invQ,Alpha.Q=AlphaQ)
    return(D2para)
  })
  names(D2fx)=names(pp.cg)
  II=abs(-1/(unlist(D2fx)*dim(X)[1]))
  return(II)
}
#### gpfrpred ####
gpfrpred=function(object,data=NULL,newtime=NULL,data.new=NULL,type=1,yregfd=NULL, hyper.p=NULL, iuu=NULL, Cov=NULL,gamma=1){
  if(class(object)=='gpfda'){
    yregfd=object$betaestlist
    hyper.p=object$hyper
    iuu=object$iuu
    Cov=object$CovFun
    gamma=object$gamma
  }
  if(type==1){
    tinput = data[J('inputdata',1,'time')]$input
    nam=unique(data[J('inputdata')]$col);x.nam=nam[which(!nam%in%c('response','time'))]
    xinput = reshape(data[J('inputdata',1,x.nam)],v.names='input',timevar='col',idvar=c('sample','batch'),direction='wide')
    setnames(xinput,names(xinput),c('type','batch','sample',x.nam))
    xinput = data.matrix(xinput[,names(xinput)%in%x.nam,with=F])
    nam2=unique(data[J('testdata')]$col);x.nam2=nam2[which(!nam%in%c('response','time'))-1]
    ttest = data[J('testdata',1,'time')]$input
    xtest = reshape(data[J('testdata',1,x.nam2)],v.names='input',timevar='col',idvar=c('sample','batch'),direction='wide')
    setnames(xtest,names(xtest),c('type','batch','sample',x.nam2))
    xtest = data.matrix(xtest[,names(xtest)%in%x.nam2,with=F])
    
    yregtmp=do.call('cbind',lapply(yregfd,function(i) predict(object=i,newdata=tinput)))
    Ut=as.matrix(data[J('scale',1)]$input)
    
    yinput = data[J('inputdata',1,'response')]$input 
    ygpobs = yinput-yregtmp%*%Ut
    
    y_gppred=gppredict(train=F,hyper=hyper.p, Data=as.matrix(xinput), Y=as.matrix(ygpobs), Data.new=as.matrix(xtest), Cov=Cov,gamma=gamma)
    ygppred = y_gppred$mu
    s2 = y_gppred$sigma^2%*%(1 + t(Ut)%*%iuu%*%Ut)
    yregtmp=unlist(lapply(yregfd,function(i) predict(i,ttest)))
    yregtmp=matrix(yregtmp,ncol=length(yregfd))
    ypred = yregtmp%*%Ut + ygppred ## fd rgression plus gp regression
  }
  
  if(type==2){
    if(!is.null(newtime)){
      ttest=newtime
      xtest=data.new
    } 
    if(is.null(newtime) & !is.null(data)){
      ttest=data[J('testdata',1,'time')]$input
      nam2=unique(data[J('testdata')]$col);x.nam2=nam2[which(!nam%in%c('response','time'))-1]
      xtest = reshape(data[J('testdata',1,x.nam2)],v.names='input',timevar='col',idvar=c('sample','batch'),direction='wide')
      setnames(xtest,names(xtest),c('type','batch','sample',x.nam2))
      xtest = data.matrix(xtest[,names(xtest)%in%x.nam2,with=F])
    } 
    train=object$train;bbat=object$bat
    fitted=matrix(0,ncol=length(bbat),nrow=length(ttest))
    fitted.var=fitted
    tdata = train[J('functional',bbat[i],'time')]$input
    yregtmp=do.call('cbind',lapply(yregfd,function(i) predict(object=i,newdata=tdata)))
    yregpred=do.call('cbind',lapply(yregfd,function(i) predict(object=i,newdata=ttest)))
    for(i in seq_along(bbat)){
      nam=unique(train[J('functional')]$col);x.nam=nam[which(!nam%in%c('response','time'))]
      xinput = reshape(train[J('functional',bbat[i],x.nam)],v.names='input',timevar='col',idvar=c('sample','batch'),direction='wide')
      setnames(xinput,names(xinput),c('type','batch','sample',x.nam))
      xinput = data.matrix(xinput[,names(xinput)%in%x.nam,with=F])
      
      yinput = train[J('functional',bbat[i],'response')]$input 
      U=as.matrix(train[J('scale',bbat[i])]$input)
      ygpobs = yinput-yregtmp%*%U
      y_gpr=gppredict(train=F,Data.new=xtest,hyper=hyper.p,Data=as.matrix(xinput), Y=as.matrix(ygpobs),Cov=Cov,gamma=gamma)
      
      ygppred = y_gpr$mu
      s2 = y_gpr$sigma^2%*%(1 + t(U)%*%iuu%*%U)
      ypred = yregpred%*%U + ygppred ## fd rgression plus gp regression
      fitted[,i]=ypred
      fitted.var[,i]=s2
    }
    ypred=as.matrix(apply(fitted,1,mean))
    s2=as.matrix(apply(fitted.var,1,mean))+apply(fitted,2,function(j) mean((j-ypred)^2))
  }
  
  
  ypredup = ypred + 1.96*sqrt(s2)
  ypredlo = ypred - 1.96*sqrt(s2)
  
  CI=cbind(ypred, ypredup, ypredlo)
  
  result=c(list(ypred=CI, time=ttest),unclass(object))
  class(result)='gpfda'
  return(result)
}


#### repeated likelihood ####
repgp.loglikelihood=function(hyper.p,  Data,Cov,gamma=1,...){#,Xprior,Xprior2){
  r.col=unique(Data[J('functional')]$col)
  f_llh=function(bat){
    old_input_data=reshape(Data[J('functional',bat)],v.names='input',timevar='col',idvar='sample',direction='wide')
    setnames(old_input_data,names(old_input_data),c('type','batch','sample',unique(Data$col)))
    drop=c('type','batch','sample')
    old_input_data=old_input_data[,!(names(old_input_data) %in% drop),with=F]
    old_input=data.matrix(old_input_data[,r.col[!r.col%in%c('response','time')],with=F])
    response=as.matrix(Data[J('functional',bat,'response')]$input)
    gp.loglikelihood2(hyper.p,Data=old_input,response=response,Cov=Cov,gamma=gamma,...)
  }
  sum(apply(as.matrix(seq_along(unique(Data$batch))), 1, function(i) f_llh(i)))
}

repgp.Dloglikelihood=function(hyper.p, Data,Cov,gamma=1){#,Xprior=prior_D1_likelihood,Xprior2=prior_likelihood){
  r.col=unique(Data[J('functional')]$col)
  f_llh=function(bat){
    old_input_data=reshape(Data[J('functional',bat)],v.names='input',timevar='col',idvar='sample',direction='wide')
    setnames(old_input_data,names(old_input_data),c('type','batch','sample',unique(Data$col)))
    drop=c('sample','batch','type')
    old_input_data=old_input_data[,!(names(old_input_data) %in% drop),with=F]
    old_input=data.matrix(old_input_data[,r.col[!r.col%in%c('response','time')],with=F])
    response=as.matrix(Data[J('functional',bat,'response')]$input)
    gp.Dlikelihood2(hyper.p,Data=old_input,response=response,Cov=Cov,gamma=gamma)
  }
  out=apply(as.matrix(seq_along(unique(Data$batch))), 1, function(i) f_llh(i))
  apply(out,1,sum)
}

#### merge data to the proper form ####

wrap_train=function(functional,scale=NULL,time=NULL,response=NULL){
  if(is.list(functional) & !is.data.frame(functional)){
    nbatch=length(functional)
    ncol=dim(functional[[1]])[2]
    nsample=dim(functional[[1]])[1]
    input=unlist(functional)
    batch=rep(1:nbatch,each=nsample*ncol)
    sample=rep(1:nsample,nbatch*ncol)
    col=rep(rep(unique(colnames(functional[[1]])),each=nsample),nbatch)
    col[col==time]='time';col[col==response]='response'
    type=rep('functional',nbatch*nsample*ncol)
    data=data.table(type=type,batch=batch,col=col,sample=sample,input=input)
    setnames(data,names(data),c('type','batch','col','sample','input'))
  }
  
  if(!is.null(scale)){
    if(dim(scale)[1]!=nbatch) warning('number of batches in functional data does not match number of batches in scale data')
    n=dim(scale)[1]
    nscale=dim(scale)[2]
    input2=matrix(t(scale),ncol=1)[,1]
    batch2=rep(1:n,each=nscale)
    sample2=rep(1,n*nscale)
    col2=rep(names(data.table(scale)),n)
    type2=rep('scale',n*nscale)
    data2=data.table(type=type2,batch=batch2,col=col2,sample=sample2,input=input2)
    setnames(data2,names(data2),c('type','batch','col','sample','input'))
    data=rbind(data,data2)
  }
  setkey(data,type,batch,col,sample)
  return(data)
}

wrap_test=function(functional,scale=NULL,test,time=NULL,response=NULL){
  if(is.matrix(functional) | is.data.frame(functional)){
    nbatch=1
    ncol=dim(functional)[2]
    nsample=dim(functional)[1]
    input=matrix(functional,ncol=1)
    batch=rep(1:nbatch,each=nsample*ncol)
    sample=rep(1:nsample,nbatch*ncol)
    col=rep(rep(unique(colnames(functional)),each=nsample),nbatch)
    col[col==time]='time';col[col==response]='response'
    type=rep('inputdata',nbatch*nsample*ncol)
    data=data.table(type=type,batch=batch,col=col,sample=sample,input=input)
    setnames(data,names(data),c('type','batch','col','sample','input'))
    
    nbatch=1
    ncol=dim(test)[2]
    nsample=dim(test)[1]
    input=matrix(test,ncol=1)
    batch=rep(1:nbatch,each=nsample*ncol)
    sample=rep(1:nsample,nbatch*ncol)
    col=rep(rep(unique(colnames(test)),each=nsample),nbatch)
    col[col==time]='time';col[col==response]='response'
    type=rep('testdata',nbatch*nsample*ncol)
    data2=data.table(type=type,batch=batch,col=col,sample=sample,input=input)
    setnames(data2,names(data2),c('type','batch','col','sample','input'))
  }
  data=rbind(data,data2)
  if(!is.null(scale)){
    #     ds1<<-dim(scale);nb<<-nbatch
    if(dim(scale)[1]!=nbatch) warning('number of batches in functional data does not match number of batches in scale data')
    n=dim(scale)[1]
    nscale=dim(scale)[2]
    input2=matrix(scale,ncol=1)[,1]
    batch2=rep(1:n,each=nscale)
    sample2=rep(1,n*nscale)
    col2=rep(names(data.table(scale)),n)
    type2=rep('scale',n*nscale)
    data2=data.table(type=type2,batch=batch2,col=col2,sample=sample2,input=input2)
    data=rbind(data,data2)
    setnames(data,names(data),c('type','batch','col','sample','input'))
  }
  setkey(data,type,batch,col,sample)
  data=data[!is.na(data$col)]
  return(data)
}

wrap=function(functional,scale=NULL,testdata=NULL,list=c('traindata','testdata'),time=NULL,response=NULL){
  if(list=='traindata') data=wrap_train(functional=functional,scale=scale,time=time,response=response)
  if(list=='testdata') data=wrap_test(functional=functional,scale=scale,test=testdata,time=time,response=response)
  return(data)
}

mysmooth=function(data,time=NULL,nbasis=NULL,norder=4,NewTime=100,pen.o=2,lambda=.05,fd=F){
  if(is.vector(data)) data=as.matrix(data)
  if(is.data.frame(data)) data=data.matrix(data)
  if(is.null(nbasis)) nbasis=round(dim(data)[1]/3)
  if(is.null(time)){
    rng=c(0,1)
    time=seq(rng[1],rng[2],len=dim(data)[1])
    new.time=seq(rng[1],rng[2],len=NewTime)
  }
  if(!is.null(time)){
    rng=range(time)
    new.time=seq(rng[1],rng[2],len=NewTime)
  }
  basis=create.bspline.basis(rng,nbasis,norder)
  pen.fdPar=fdPar(basis, pen.o, lambda)
  xs=smooth.basis(time,as.matrix(data),pen.fdPar)
  xsfd=xs$fd
  xsy2c=xs$y2cMap
  c2rMap=eval.basis(time,basis)
  
  res=as.matrix(predict(xsfd,time))-as.matrix(data)
  #   a4<<-res
  res.var=apply(res^2,1,sum)/(dim(res)[2]-(dim(res)[2]>1))
  res=smooth.basis(time,res.var,pen.fdPar)
  resvec=predict(res$fd,time)
  Sigmayhat = c2rMap %*% xsy2c %*% diag(as.vector(resvec)) %*% t(xsy2c) %*% t(c2rMap)
  sqrt_sigma=sqrt(diag(Sigmayhat))
  pen.fdPar2=fdPar(create.bspline.basis(rng,length(sqrt_sigma)-1,norder=6), 2, lambda=1e-10)
  sd_fd=predict(smooth.basis(time,sqrt_sigma,pen.fdPar2)$fd,new.time)
  
  if(fd==F)  xsvec=predict(xsfd,new.time)
  if(fd==T)  xsvec=xsfd
  return(list(xsvec,sd_fd))
}

fda_smooth=function(dt,pre.rng,new.sam){
  dtdim=dim(dt)
  prebasis=create.bspline.basis(range(dt$input.time),nbasis=dtdim[1]-2,norder=4)
  D2fdPar = fdPar(prebasis, lambda=1e-10)
  prefd=smooth.basis(as.matrix(dt$input.time),as.matrix(dt[,!names(dt)%in%c('type','batch','sample','input.time'),with=F]),D2fdPar)$fd
  pre.mat=eval.fd(seq(0,1,len=new.sam), prefd)
  colnames(pre.mat)=names(dt[,!names(dt)%in%c('type','batch','sample','input.time'),with=F])
  pre.mat=data.table(type=dt$type[1],batch=dt$batch[1],sample=1:new.sam,pre.mat,input.time=seq(pre.rng[1],pre.rng[2],len=new.sam))
  return(pre.mat)
}

wrap2=function(data,new.sample=NULL){
  if(var(data[J('functional'),length(sample),by='batch']$V1)>0 | !is.null(new.sample)){
    fdata=reshape(data[J('functional')],v.names='input',timevar='col',idvar=c('sample','batch'),direction='wide')
    setkey(fdata,'type','batch','sample')
    if(is.null(new.sample)) new.sample=200
    rng.min=min(fdata$input.time)
    rng.max=max(fdata$input.time)
    fdata.new=NULL
    for(i in seq_along(unique(fdata$batch))){
      fdata.new=rbind(fdata.new,fda_smooth(fdata[J('functional',i)],c(rng.min,rng.max)))
    }
    fn.dim=dim(fdata.new)
    fdata.new=reshape(fdata.new,idvar=c('batch','sample'),varying=list(4:fn.dim[2]),v.names='input',direction='long')
    coll=fdata.new$time
    col=vector(length=dim(fdata.new)[1])
    col[coll==1]='response'
    col[coll==max(coll)]='time'
    col[coll>1&coll<max(coll)]=paste0('x',(coll[coll>1&coll<max(coll)]-1))
    fdata.new=data.table(fdata.new$type,fdata.new$batch,col,fdata.new$sample,fdata.new$input)
    setnames(fdata.new,names(fdata.new),c('type','batch','col','sample','input'))
    data=rbind(fdata.new,data[J('scale')])
  }
  else
    data=data
  return(data)
}
