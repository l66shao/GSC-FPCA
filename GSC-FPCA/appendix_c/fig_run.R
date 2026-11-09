#!/usr/bin/env Rscript
# Curves and pre-exposure residuals for Figures 1-3: scenario B, N_tr=40, N_co=100 (cell 5).
# Same seeds and random-number use as run_simulation.R, so the ATT columns must match the main results exactly.
# Usage: Rscript fig_run.R <worker 1..K> <K> <first iter> <last iter>
a=commandArgs(TRUE); w=as.integer(a[1]); K=as.integer(a[2]); i0=as.integer(a[3]); i1=as.integer(a[4])
library(Rcpp);library(MASS);library(splines);pdf(NULL)
Sys.setenv(COMEGA="1")
Rcpp::sourceCpp("prior_mean_zero/sampler.cpp"); source("prior_mean_zero/mfpca.R"); source("prior_mean_zero/sim_functions.R")
params=make_params()
tg=seq(0,1,.01)

run_rep_fig=function(Ntr,Nco,lambda,n_min,params,nloop,burnin,thin=1,bins=seq(.1,1,.1),effect=5){
  N=Ntr+Nco; d=simulate_cell(N,lambda,n_min,params)
  b_control=c(rep(TRUE,Nco),rep(FALSE,Ntr))
  Y=T_=X=Xc=vector("list",N); post=list(); pre=list()
  for(i in 1:N){
    n=length(d[[i]]$Y); k=n
    if(!b_control[i]){
      k=max(floor(2*n/3),1)
      s=(k+1):n; T0=d[[i]]$T[k]; dl=d[[i]]$T[s]-T0
      post[[length(post)+1]]=list(i=i,t=d[[i]]$T[s],delta=dl,Y0=d[[i]]$Y[s],e=d[[i]]$e[s],Yobs=d[[i]]$Y[s]+effect*dl,X=d[[i]]$X[s,,drop=FALSE],Xc=d[[i]]$Xc[s,,drop=FALSE])
      s=1:k
      pre[[length(pre)+1]]=list(i=i,t=d[[i]]$T[s],delta=d[[i]]$T[s]-T0,Yobs=d[[i]]$Y[s],X=d[[i]]$X[s,,drop=FALSE],Xc=d[[i]]$Xc[s,,drop=FALSE])
    }
    s=1:k
    Y[[i]]=list(matrix(d[[i]]$Y[s])); T_[[i]]=list(matrix(d[[i]]$T[s]))
    X[[i]]=list(d[[i]]$X[s,,drop=FALSE]); Xc[[i]]=list(d[[i]]$Xc[s,,drop=FALSE])
  }
  fit=mfpca_1d(Y_sparse=Y,time_sparse=T_,X_sparse=X,X_c_sparse=Xc,group=rep(1,N),Q=2,nknots=3,
               nloop=nloop,burnin=burnin,thin=thin,orth=FALSE,sim=FALSE,PSLOO=0,b_control=as.numeric(b_control))
  A=fit[[1]]; Mu=fit[[2]]; Dd=fit[[3]]; Sg=unlist(fit[[5]]); Th=fit[[6]]; Be=fit[[7]]; Bc=fit[[8]]; pc=fit[[9]][[5]][[1]]
  S=length(A)
  draws=function(L){
    po=do.call(rbind,lapply(L,function(p) data.frame(i=p$i,t=p$t,delta=p$delta,Yobs=p$Yobs)))
    Xm=do.call(rbind,lapply(L,function(p) p$X)); Xcm=do.call(rbind,lapply(L,function(p) p$Xc))
    Pm=pc[,round(po$t,3)*1000+1,drop=FALSE]
    Fm=matrix(NA,nrow(po),S)
    for(s in 1:S){
      f=crossprod(Pm,Mu[[s]])+rowSums(t(Pm)*t(Th[[s]]%*%A[[s]][,po$i,drop=FALSE]))+Xcm%*%Bc[[s]]
      for(j in 1:2) f=f+Xm[,j]*crossprod(Pm,Be[[s]][,j])
      Fm[,s]=f
    }
    list(po=po,Fm=Fm)
  }
  dp=draws(post); po=dp$po; Fm=dp$Fm
  po$Y0=unlist(lapply(post,function(p) p$Y0)); po$e=unlist(lapply(post,function(p) p$e))
  Nz=matrix(rnorm(length(Fm)),nrow(Fm),S)*rep(sqrt(Sg),each=nrow(Fm))
  po$bin=round(po$delta*10)/10
  one=function(w,lab){
    e=po$Yobs[w]-Fm[w,,drop=FALSE]; a=colMeans(e); ap=a-colMeans(Nz[w,,drop=FALSE])
    data.frame(bin=lab,n=sum(w),truth=mean(effect*po$delta[w]),est=mean(a),pi_lo=quantile(ap,.025,names=FALSE),pi_hi=quantile(ap,.975,names=FALSE))
  }
  r=list()
  for(b in bins){w=abs(po$bin-b)<1e-8; if(any(w)) r[[length(r)+1]]=one(w,sprintf("%.1f",b))}
  r[[length(r)+1]]=one(rep(TRUE,nrow(po)),"all")
  att=do.call(rbind,r)
  # in-sample pre-exposure residuals of treated units, by nearest 0.1 before the last pre-exposure time
  dq=draws(pre); q=dq$po; q$res=q$Yobs-rowMeans(dq$Fm); q$bin=round(q$delta*10)/10
  res=aggregate(res~bin,q,mean); res$n=aggregate(res~bin,q,length)$res
  # posterior means of the identified functionals on a common grid
  idx=round(tg*1000)+1; Bg=t(pc[,idx,drop=FALSE])
  mu=as.vector(Bg%*%Reduce(`+`,Mu)/S)
  M=Reduce(`+`,lapply(1:S,function(s) Th[[s]]%*%Dd[[s]]%*%t(Th[[s]])))/S
  list(att=att,res=res,mu=mu,cov_coef=M,Bg=Bg,sigma2=mean(Sg))
}

dir.create("fig",showWarnings=FALSE)
cl=5; s=c(40,100); p=c(4,3)
for(iter in i0:i1){ if(iter%%K!=w%%K) next
  f=sprintf("fig/rep_%04d.rds",iter); if(file.exists(f)) next
  set.seed(iter*100+cl)
  r=run_rep_fig(s[1],s[2],p[1],p[2],params,1000,500)
  saveRDS(r,f,compress=TRUE)
}
cat("fig worker",w,"finished\n")
