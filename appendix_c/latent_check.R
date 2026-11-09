#!/usr/bin/env Rscript
# Corrected sampler under latent selection: treated scores ~ N(nu, D). Compares empirical ATT error with B_N,
# for controls-only and pooled fitting of mu, Gamma, beta. Usage: Rscript latent_check.R <controls|pooled>
variant=commandArgs(TRUE)[1]
library(Rcpp);library(MASS);library(splines);pdf(NULL); Sys.setenv(COMEGA="1")
Rcpp::sourceCpp(if(variant=="controls") "prior_mean_zero/sampler.cpp" else "prior_mean_zero/sampler_pooled.cpp")
source("prior_mean_zero/mfpca.R"); source("prior_mean_zero/sim_functions.R")
params=make_params(); grid=seq(0,1,1/39); nu=c(1.5,1.5); z=qnorm(.975)
one=function(Ntr,Nco,m){ N=Ntr+Nco
  tm=lapply(1:N,function(i){ if(i<=Nco) list(grid[sort(sample(40)[1:min(max(rpois(1,4),3),10)])])
    else list(c(sort(sample(grid[grid<=.6])[1:m]),sort(sample(grid[grid>.6])[1:2])))})
  phi=basis_setup_sparse(tm,3,delta=1/1000)[[3]][[1]]; Th=params[[6]]; Y=T_=X=Xc=vector("list",N); po=list(); Bsum=0; vc0=0; nO=0
  for(i in 1:N){ P=t(phi[[i]]); n=nrow(P); A=rnorm(2)+(i>Nco)*nu; e=rnorm(n)
    Xci=cbind(rnorm(n),rbinom(n,1,.5)); Xi=cbind(rep(rnorm(1),n),rep(rbinom(1,1,.5),n)); colnames(Xci)=colnames(Xi)=c("X1","X2")
    y=as.vector(P%*%params[[7]]+P%*%Th%*%A+e+Xci%*%params[[10]]); for(j in 1:2) y=y+Xi[,j]*as.vector(P%*%params[[9]][,j])
    k=if(i<=Nco) n else m; s=1:k
    if(i>Nco){ F_=P[s,]%*%Th; Mi=solve(crossprod(F_)+diag(2)); fp=P[-s,,drop=FALSE]%*%Th; ci=colSums(fp)
      Bsum=Bsum+sum(ci*(Mi%*%nu)); vc0=vc0+as.numeric(t(ci)%*%Mi%*%ci); nO=nO+nrow(fp)
      po[[length(po)+1]]=list(i=i,t=tm[[i]][[1]][-s],Y0=y[-s],e=e[-s],X=Xi[-s,,drop=FALSE],Xc=Xci[-s,,drop=FALSE])}
    Y[[i]]=list(matrix(y[s])); T_[[i]]=list(matrix(tm[[i]][[1]][s])); X[[i]]=list(Xi[s,,drop=FALSE]); Xc[[i]]=list(Xci[s,,drop=FALSE])}
  fit=mfpca_1d(Y_sparse=Y,time_sparse=T_,X_sparse=X,X_c_sparse=Xc,group=rep(1,N),Q=2,nknots=3,nloop=600,burnin=200,thin=1,orth=FALSE,sim=FALSE,PSLOO=0,b_control=as.numeric((1:N)<=Nco))
  AA=fit[[1]];Mu=fit[[2]];Sg=unlist(fit[[5]]);TT=fit[[6]];Be=fit[[7]];Bc=fit[[8]];pc=fit[[9]][[5]][[1]]
  ids=unlist(lapply(po,function(p) rep(p$i,length(p$t)))); ts=unlist(lapply(po,function(p) p$t)); Y0=unlist(lapply(po,function(p) p$Y0)); eb=mean(unlist(lapply(po,function(p) p$e)))
  Xm=do.call(rbind,lapply(po,function(p) p$X)); Xcm=do.call(rbind,lapply(po,function(p) p$Xc)); Pm=pc[,round(ts,3)*1000+1,drop=FALSE]; S=length(AA); dr=numeric(S)
  for(s in 1:S){ fv=crossprod(Pm,Mu[[s]])+rowSums(t(Pm)*t(TT[[s]]%*%AA[[s]][,ids,drop=FALSE]))+Xcm%*%Bc[[s]]; for(j in 1:2) fv=fv+Xm[,j]*crossprod(Pm,Be[[s]][,j]); dr[s]=mean(Y0-fv)}
  dp=dr-rnorm(S,0,sqrt(Sg/nO))
  data.frame(est=mean(dr),ebar=eb,B=Bsum/nO,v0=vc0/nO^2+1/nO,sigma2=mean(Sg),post_sd=sd(dr),pred_sd=sd(dp),
    flo=quantile(dr,.025,names=FALSE),fhi=quantile(dr,.975,names=FALSE),plo=quantile(dp,.025,names=FALSE),phi=quantile(dp,.975,names=FALSE))
}
f=sprintf("latent_%s.csv",variant)
cells=rbind(data.frame(m=c(2,4,8,16),Ntr=40,Nco=100,reps=100),data.frame(m=c(2,8,20),Ntr=160,Nco=400,reps=50),data.frame(m=8,Ntr=160,Nco=1600,reps=25))
done=if(file.exists(f)) read.csv(f) else NULL
for(r in 1:100) for(cl in 1:nrow(cells)){ if(r>cells$reps[cl]) next
  if(!is.null(done)&&any(done$r==r&done$cl==cl)) next
  set.seed(r*100+cl); o=one(cells$Ntr[cl],cells$Nco[cl],cells$m[cl])
  write.table(cbind(r=r,cl=cl,cells[cl,1:3],o),f,sep=",",append=file.exists(f),col.names=!file.exists(f),row.names=FALSE)}
cat("done",variant,"\n")
