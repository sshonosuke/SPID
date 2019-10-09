library(MASS)
library(MCMCpack)
library(SparseM)
library(statmod)



#######  PWD: Pair-wise Difference Prior  #######
# W is (non-scaled) adjacent matrix
PWD=function(Data,Z,W,mcmc=12000,burn=2000,family="LN",print=F){
  m=dim(Data)[1]
  N=length(Z)-1
  
  # Distributions
  if(family=="LN"){
    p=2
    Dist=function(x,U){ 
      mm=U[1]; ss=sqrt(exp(U[2]))
      plnorm(x,mm,ss)
    }
  }
  if(family=="SM"){
    p=3
    Dist=function(x,U){ 
      a=exp(U[1]); b=exp(U[2]); c=exp(U[3])
      y=(x/b)^a; 1-(1+y)^(-c)
    }
  }
  if(family=="DG"){
    p=3
    Dist=function(x,U){ 
      a=exp(U[1]); b=exp(U[2]); c=exp(U[3])
      y=(x/b)^a; (1+1/y)^(-c)
    }
  }
  
  # log-likelihood
  LogLike=function(data,U){
    val=0
    for(k in 1:N){
      dd=Dist(Z[k+1],U)-Dist(Z[k],U)
      dd[dd==0]=10^(-10)
      val=val+(data[k])*log(dd)
    }
    return(val)
  }
  
  # Gaissian approximation of likelihood
  Ut=matrix(NA,m,p); P=array(NA,c(m,p,p))
  for(i in 1:m){
    opt=function(u){ -LogLike(Data[i,],u) }
    ml=optim(par=rep(0,p),fn=opt,hessian=T)
    Ut[i,]=ml$par
    P[i,,]=ml$hessian
    if(min(eigen(P[i,,])$values)<0){ 
      P[i,,]=diag(diag(P[i,,])) 
      diag(P[i,,])[diag(P[i,,])<0]=0.01
    }
  }
  
  # prior
  b0=1; c0=1  # gamma priors for random effect precision and spatial precision
  a0=0.001   # precision parameter in priors for grand means
  
  # initial values
  U.pos=array(NA,c(mcmc,m,p))
  Mu.pos=matrix(NA,mcmc,p)
  Tau.pos=matrix(NA,mcmc,p)
  Lam.pos=matrix(NA,mcmc,p)
  U=Ut
  Mu=rep(0,p)
  Tau=rep(10,p)
  Lam=rep(10,p)
  
  # MCMC iterations
  W=as(W,"sparseMatrix")
  rW=apply(W,1,sum)
  WW=as(diag(rW)-W,"sparseMatrix")
  
  Q=function(tau,lam){ tau*diag(m)+lam*WW }
  tr=function(x){ sum(diag(x)) }
  
  for(r in 1:mcmc){
    # update mu
    mm=Tau*apply(U,2,sum)/(m*Tau+a0)
    ss=1/(m*Tau+a0)
    Mu=rnorm(p,mm,sqrt(ss))
    Mu.pos[r,]=Mu
    resid=t(t(U)-Mu)
    WD=c()
    for(l in 1:p){
      U.mat=matrix(rep(U[,l],m),m,m)
      Dif=(U.mat-t(U.mat))^2
      WD[l]=sum(W*Dif)/2
    }
    
    # update precision parameters (Langevin)
    for(l in 1:p){
      h=0.1
      cu=c(Tau[l],Lam[l])
      Q1=Q(cu[1],cu[2])
      invQ1=solve(Q1)
      aa=sum(resid[,l]^2); bb=WD[l]
      dU1=c()
      dU1[1]=-0.5*tr(invQ1)+0.5*aa
      dU1[2]=-0.5*tr(invQ1%*%WW)+0.5*bb
      prop=cu-h*dU1+sqrt(2*h)*rnorm(2)
      prop[prop<0.01]=0.01
      
      Q2=Q(prop[1],prop[2])
      invQ2=solve(Q2)
      val1=chol(as(Q1,"matrix.csr"))@log.det-0.5*cu[1]*aa-0.5*cu[2]*bb
      val2=chol(as(Q2,"matrix.csr"))@log.det-0.5*prop[1]*aa-0.5*prop[2]*bb
      dU2=c()
      dU2[1]=-0.5*tr(invQ2)+0.5*aa
      dU2[2]=-0.5*tr(invQ2%*%WW)+0.5*bb
      q2=-sum((prop-cu+h*dU1)^2)/(4*h)
      q1=-sum((cu-prop+h*dU2)^2)/(4*h)
      
      prob=min(1,exp(val2-q2-(val1-q1)))
      ch=rbinom(1,1,prob)
      Tau[l]=cu[1]+ch*(prop[1]-cu[1])
      Lam[l]=cu[2]+ch*(prop[2]-cu[2])
    }
    Tau.pos[r,]=Tau
    Lam.pos[r,]=Lam
    
    # update U (independent MH)
    for(i in 1:m){
      pp1=U[i,]
      ai=Lam*apply(W[i,]*U,2,sum)
      bi=P[i,,]%*%Ut[i,]
      ci=Tau*Mu
      mm=ai+bi+ci   # aproximated mean 
      vv=solve(P[i,,]+diag(Tau+Lam*rW[i])+diag(rep(0.0001,p)))   # approximated covariance
      pp2=mvrnorm(1,vv%*%mm,vv)   # proposal
      resid1=pp1-Ut[i,]
      L1=LogLike(Data[i,],pp1)+0.5*as.vector(t(resid1)%*%P[i,,]%*%resid1)
      resid2=pp2-Ut[i,]
      L2=LogLike(Data[i,],pp2)+0.5*as.vector(t(resid2)%*%P[i,,]%*%resid2)
      prob=min(1,exp(L2-L1))
      U[i,]=pp1+rbinom(1,1,prob)*(pp2-pp1)
    }
    U.pos[r,,]=U
    if(print==T & round(r/100)==r/100){ print(r) }
  }
  
  om=1:burn
  Res=list(U.pos[-om,,],Mu.pos[-om,],Tau.pos[-om,],Lam.pos[-om,],Ut,P)
  names(Res)=c("U","Mu","Tau","Lam","ML","Hessian")
  return(Res)
}









######  PWL: Pair-wise Difference Laplace Prior  ######
# W is (non-scaled) adjacent matrix
PWL=function(Data,Z,W,mcmc=12000,burn=2000,family="LN",Rn=10,print=F){
  m=dim(Data)[1]
  N=length(Z)-1
  if(family=="LN"){
    p=2
    Dist=function(x,U){ 
      mm=U[1]; ss=sqrt(exp(U[2]))
      plnorm(x,mm,ss)
    }
  }
  if(family=="SM"){
    p=3
    Dist=function(x,U){ 
      a=exp(U[1]); b=exp(U[2]); c=exp(U[3])
      y=(x/b)^a; 1-(1+y)^(-c)
    }
  }
  if(family=="DG"){
    p=3
    Dist=function(x,U){ 
      a=exp(U[1]); b=exp(U[2]); c=exp(U[3])
      y=(x/b)^a; (1+1/y)^(-c)
    }
  }
  
  # log-likelihood
  LogLike=function(data,U){
    val=0
    for(k in 1:N){
      dd=Dist(Z[k+1],U)-Dist(Z[k],U)
      dd[dd==0]=10^(-10)
      val=val+(data[k])*log(dd)
    }
    return(val)
  }
  
  # Gaissian approximation of likelihood
  Ut=matrix(NA,m,p); P=array(NA,c(m,p,p))
  for(i in 1:m){
    opt=function(u){ -LogLike(Data[i,],u) }
    ml=optim(par=rep(0,p),fn=opt,hessian=T)
    Ut[i,]=ml$par
    P[i,,]=ml$hessian
    if(min(eigen(P[i,,])$values)<0){ 
      P[i,,]=diag(diag(P[i,,])) 
      diag(P[i,,])[diag(P[i,,])<0]=0.01
    }
  }
  
  # prior
  b0=1; c0=1  # gamma prior for random effect precision and spatial precision
  a0=0.001   # precision parameter for grand mean
  
  # initial values
  U.pos=array(NA,c(mcmc,m,p))
  Mu.pos=matrix(NA,mcmc,p)
  Tau.pos=matrix(NA,mcmc,p)
  Lam.pos=c()
  U=Ut
  Mu=rep(0,p)
  Tau=rep(0.5*mean(P),p)
  Lam=0.5*mean(P)
  uL=2*mean(P)
  
  # MCMC iterations
  W=as(W,"sparseMatrix")
  S=as(matrix(1/mean(P),m,m)*W,"sparseMatrix")
  
  Q=function(tau,SS){ 
    sW=W/SS; sW[is.na(sW)]=0
    rW=apply(sW,1,sum)
    tau*diag(m)+as(diag(rW)-sW,"sparseMatrix")
  }
  
  tr=function(x){ sum(diag(x)) }
  delta=sum(W)/2
  
  # MCMC
  for(r in 1:mcmc){
    # update grand mean
    mm=Tau*apply(U,2,sum)/(m*Tau+a0)
    ss=1/(m*Tau+a0)
    Mu=rnorm(p,mm,sqrt(ss))
    Mu.pos[r,]=Mu
    resid=t(t(U)-Mu)
    aa=apply(resid^2,2,sum)
    
    Dif=0
    for(l in 1:p){
      U.mat=matrix(rep(U[,l],m),m,m)
      Dif=Dif+W*(U.mat-t(U.mat))^2
    }
    bb=sum(sqrt(Dif))/2
    
    # update precision parameters
    band=0.5
    prop.lam=Lam+band*rnorm(1)
    prop.tau=Tau+band*rnorm(p)
    prop.lam[prop.lam<0.01]=0.01
    prop.tau[prop.tau<0.01]=0.01
    prop.lam[prop.lam>uL]=uL
    prop.tau[prop.tau>uL]=uL
    
    LD1=c(); LD2=c()
    for(k in 1:Rn){
      rS1=matrix(0,m,m); rS2=matrix(0,m,m)
      for(i in 1:m){ 
        rS1[i,1:i]=rgamma(i,1,0.5*prop.lam^2) 
        rS2[i,1:i]=rgamma(i,1,0.5*Lam^2) 
      }
      rS1=rS1*W; rS2=rS2*W
      SS1=rS1+t(rS1); SS2=rS2+t(rS2)
      val1=c(); val2=c()
      for(l in 1:p){
        Q1=Q(prop.tau[l],SS1)
        val1[l]=-chol(as(Q1,"matrix.csr"))@log.det
        Q2=Q(Tau[l],SS2)
        val2[l]=-chol(as(Q2,"matrix.csr"))@log.det
      }
      LD1[k]=sum(val1)-sum(log(attr(SS1,"x")))/4-delta*log(prop.lam)
      LD2[k]=sum(val2)-sum(log(attr(SS2,"x")))/4-delta*log(Lam)
    }
    
    mv1=max(LD1); mv2=max(LD2)
    dd1=mv1+log(sum(exp((LD1-mv1))))   
    dd2=mv2+log(sum(exp((LD2-mv2))))   
    
    val1=-dd1-0.5*sum(aa*prop.tau)-prop.lam*bb
    val2=-dd2-0.5*sum(aa*Tau)-Lam*bb
    prob=min(1,exp(val1-val2))
    ch=rbinom(1,1,prob)
    Tau=Tau+ch*(prop.tau-Tau)
    Lam=Lam+ch*(prop.lam-Lam)
    
    Tau.pos[r,]=Tau
    Lam.pos[r]=Lam
    
    # update S
    newS=matrix(0,m,m)
    for(i in 1:m){
      newS[i,(1:i)]=1/rinvgauss(i,sqrt(Lam^2/Dif[i,1:i]),Lam^2)
    }
    newS=newS*W
    S=(newS+t(newS))
    
    # update U
    for(i in 1:m){
      pp1=U[i,]
      Si=diag(Tau+sum(1/S[i,W[i,]>0]))
      Ai=U/S[i,]; Ai[abs(Ai)==Inf]=0
      ai=apply(W[i,]*Ai,2,sum)
      bi=P[i,,]%*%Ut[i,]
      ci=Tau*Mu
      mm=ai+bi+ci   # aproximated mean 
      vv=solve(P[i,,]+Si+diag(rep(0.0001,p)))   # approximated covariance
      pp2=mvrnorm(1,vv%*%mm,vv)   # proposal
      
      resid1=pp1-Ut[i,]
      L1=LogLike(Data[i,],pp1)+0.5*as.vector(t(resid1)%*%P[i,,]%*%resid1)
      resid2=pp2-Ut[i,]
      L2=LogLike(Data[i,],pp2)+0.5*as.vector(t(resid2)%*%P[i,,]%*%resid2)
      prob=min(1,exp(L2-L1))
      U[i,]=pp1+rbinom(1,1,prob)*(pp2-pp1)
    }
    U.pos[r,,]=U
    if(print & round(r/100)==r/100){ print(r) }
  }
  om=1:burn
  
  Res=list(U.pos[-om,,],Mu.pos[-om,],Tau.pos[-om,],Lam.pos[-om],Ut,P)
  names(Res)=c("U","Mu","Tau","Lam","ML","Hessian")
  return(Res)
}





#######  Independent Random Effect Model  #######
IRE=function(Data,Z,mcmc=12000,burn=2000,family="LN",print=F){
  m=dim(Data)[1]
  N=length(Z)-1
  
  if(family=="LN"){
    p=2
    Dist=function(x,U){ 
      mm=U[1]; ss=sqrt(exp(U[2]))
      plnorm(x,mm,ss)
    }
  }
  if(family=="SM"){
    p=3
    Dist=function(x,U){ 
      a=exp(U[1]); b=exp(U[2]); c=exp(U[3])
      y=(x/b)^a; 1-(1+y)^(-c)
    }
  }
  if(family=="DG"){
    p=3
    Dist=function(x,U){ 
      a=exp(U[1]); b=exp(U[2]); c=exp(U[3])
      y=(x/b)^a; (1+1/y)^(-c)
    }
  }
  
  # log-likelihood
  LogLike=function(data,U){
    val=0
    for(k in 1:N){
      dd=Dist(Z[k+1],U)-Dist(Z[k],U)
      dd[dd==0]=10^(-10)
      val=val+(data[k])*log(dd)
    }
    return(val)
  }
  
  # Gaissian approximation of likelihood
  Ut=matrix(NA,m,p); P=array(NA,c(m,p,p))
  for(i in 1:m){
    opt=function(u){ -LogLike(Data[i,],u) }
    ml=optim(par=rep(0,p),fn=opt,hessian=T)
    Ut[i,]=ml$par
    P[i,,]=ml$hessian
    if(min(eigen(P[i,,])$values)<0){ 
      P[i,,]=diag(diag(P[i,,])) 
      diag(P[i,,])[diag(P[i,,])<0]=0.01
    }
  }
  
  # prior
  b0=1; c0=1  # gamma prior for random effect precision 
  a0=0.001   # precision parameter for grand mean
  
  # initial values
  U.pos=array(NA,c(mcmc,m,p))
  Mu.pos=matrix(NA,mcmc,p)
  Tau.pos=matrix(NA,mcmc,p)
  U=Ut
  M=rep(0,p)
  Tau=rep(10,p)
  
  # MCMC iterations
  for(r in 1:mcmc){
    # update mu
    mm=Tau*apply(U,2,sum)/(m*Tau+a0)
    ss=1/(m*Tau+a0)
    Mu=rnorm(p,mm,sqrt(ss))
    Mu.pos[r,]=Mu
    
    # update tau 
    resid=t(t(U)-Mu)
    sq=apply(resid^2,2,sum)
    Tau=rgamma(p,m/2+b0,sq/2+c0)
    Tau.pos[r,]=Tau
    
    # update U
    for(i in 1:m){
      pp1=U[i,]
      bi=as.vector(P[i,,]%*%Ut[i,])
      ci=Tau*Mu
      mm=bi+ci   # aproximated mean 
      vv=solve(P[i,,]+diag(Tau+rep(0.0001,p)))   # approximated covariance
      pp2=mvrnorm(1,vv%*%mm,vv)   # proposal
      resid1=pp1-Ut[i,]
      L1=LogLike(Data[i,],pp1)+0.5*as.vector(t(resid1)%*%P[i,,]%*%resid1)
      resid2=pp2-Ut[i,]
      L2=LogLike(Data[i,],pp2)+0.5*as.vector(t(resid2)%*%P[i,,]%*%resid2)
      prob=min(1,exp(L2-L1))
      U[i,]=pp1+rbinom(1,1,prob)*(pp2-pp1)
    }
    U.pos[r,,]=U
    if(print==T & round(r/100)==r/100){ print(r) }
  }
  
  om=1:burn
  Res=list(U.pos[-om,,],Mu.pos[-om,],Tau.pos[-om,],Ut,P)
  names(Res)=c("U","Mu","Tau","ML","Hessian")
  return(Res)
}






