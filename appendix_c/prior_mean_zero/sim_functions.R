# GSC-FPCA simulation rerun: helpers. mfpca_1d is Lucy's apart from c_omega being set from COMEGA.
# The C++ changes are listed in run_simulation.R.
library(splines); library(MASS)

# Same calculations as basis_setup_sparse (basis_setup_sparse_Lucy.R) for the elements the
# simulation uses, with the time lookup vectorized. It also accepts subjects with a single
# pre-exposure observation, which the original's Phi_D_T loop cannot (2:1 indexing).
basis_setup_sparse=function(time_sparse,nknots=3,orth=FALSE,delta=1/1000,plotit=FALSE){
  N=length(time_sparse)
  for(i in 1:N) time_sparse[[i]][[1]]=round(time_sparse[[i]][[1]]/delta)*delta
  tu=round(sort(unique(unlist(lapply(time_sparse,function(x) x[[1]]))))/delta)*delta
  time_cont=round(seq(0,1/delta)*delta/delta)*delta
  qs=1/(nknots+1); kn=quantile(tu,qs)
  if(nknots>1) for(q in 2:nknots) kn=c(kn,q*quantile(tu,qs))
  kn=as.vector(kn)
  temp=bs(time_cont,knots=kn,degree=3,intercept=TRUE)
  for(k in 1:(nknots+4)) temp[,k]=temp[,k]/sqrt(sum(temp[,k]*temp[,k]))
  phi_cont=t(sqrt(1/delta)*temp)
  phi_t=lapply(1:N,function(i){
    idx=vapply(as.vector(time_sparse[[i]][[1]]),function(t) which.min(abs(time_cont-t)),1L)
    phi_cont[,idx,drop=FALSE]})
  list(list(kn),time_sparse,list(phi_t),time_cont,list(phi_cont))
}

make_params=function(){
  set.seed(123)                       # same truth as the original script, in every cell
  nknots=3; Q=2
  D=diag(rep(1,Q))
  MU=mvrnorm(1,mu=rep(0,nknots+4),Sigma=diag(nknots+4))
  beta=matrix(rep(1,2*(nknots+4)),ncol=2); beta_c=c(1,1)
  x=seq(0,2*pi,length.out=7); s=sin(x); co=cos(x)
  THETA=cbind(s/norm(s,"2"),co/norm(co,"2"))
  list(nknots,Q,D,D,1,THETA,MU,0,beta,beta_c)
}

# n_i = max(Poisson(lambda), n_min), capped at n_cap, times drawn without replacement from a
# grid of 20 points on [0,1]; truth generated on a separate basis as in the original.
simulate_cell=function(N,lambda,n_min,params,n_cap=10,N_T_max=20,design="paper"){
  grid=seq(0,1,1/(N_T_max-1))
  times=lapply(1:N,function(i){
    n=if(design=="paper") min(max(rpois(1,lambda),n_min),n_cap)
      else max(rpois(1,max(0,lambda-3))+3,n_min)
    list(grid[sort(sample(N_T_max)[1:n])])})
  bs_=basis_setup_sparse(times,params[[1]],delta=1/600)
  phi=bs_[[3]][[1]]
  A=mvrnorm(N,c(0,0),params[[3]]%*%params[[4]]%*%params[[3]])
  out=vector("list",N)
  for(i in 1:N){
    n=ncol(phi[[i]]); P=t(phi[[i]])
    Xc=cbind(rnorm(n),rbinom(n,1,.5)); colnames(Xc)=c("X1","X2")
    X=cbind(rep(rnorm(1),n),rep(rbinom(1,1,.5),n)); colnames(X)=c("X1","X2")
    e=rnorm(n,0,sqrt(params[[5]]))
    Y=P%*%params[[7]]+P%*%params[[6]]%*%A[i,]+e+Xc%*%params[[10]]
    for(j in 1:2) Y=Y+X[,j]*(P%*%params[[9]][,j])
    out[[i]]=list(T=bs_[[2]][[i]][[1]],Y=as.vector(Y),e=e,X=X,Xc=Xc)
  }
  out
}

run_rep=function(Ntr,Nco,lambda,n_min,params,nloop,burnin,thin=1,bins=seq(.1,1,.1),design="paper",pre_rule="two_thirds",effect=5){
  N=Ntr+Nco; d=simulate_cell(N,lambda,n_min,params,design=design)
  b_control=c(rep(TRUE,Nco),rep(FALSE,Ntr))
  Y=T_=X=Xc=vector("list",N); post=list()
  for(i in 1:N){
    n=length(d[[i]]$Y); k=n
    if(!b_control[i]){
      k=if(pre_rule=="two_thirds") max(floor(2*n/3),1) else ceiling(n/2)
      s=(k+1):n; T0=d[[i]]$T[k]; dl=d[[i]]$T[s]-T0
      post[[length(post)+1]]=list(i=i,t=d[[i]]$T[s],delta=dl,Y0=d[[i]]$Y[s],e=d[[i]]$e[s],Yobs=d[[i]]$Y[s]+effect*dl,X=d[[i]]$X[s,,drop=FALSE],Xc=d[[i]]$Xc[s,,drop=FALSE])
    }
    s=1:k
    Y[[i]]=list(matrix(d[[i]]$Y[s])); T_[[i]]=list(matrix(d[[i]]$T[s]))
    X[[i]]=list(d[[i]]$X[s,,drop=FALSE]); Xc[[i]]=list(d[[i]]$Xc[s,,drop=FALSE])
  }
  fit=mfpca_1d(Y_sparse=Y,time_sparse=T_,X_sparse=X,X_c_sparse=Xc,group=rep(1,N),Q=2,nknots=3,
               nloop=nloop,burnin=burnin,thin=thin,orth=FALSE,sim=FALSE,PSLOO=0,b_control=as.numeric(b_control))
  A=fit[[1]]; Mu=fit[[2]]; Sg=unlist(fit[[5]]); Th=fit[[6]]; Be=fit[[7]]; Bc=fit[[8]]; pc=fit[[9]][[5]][[1]]
  S=length(A)
  # counterfactual draws for each treated post-exposure observation: obs x draws
  po=do.call(rbind,lapply(post,function(p) data.frame(i=p$i,t=p$t,delta=p$delta,Y0=p$Y0,e=p$e,Yobs=p$Yobs)))
  Xm=do.call(rbind,lapply(post,function(p) p$X)); Xcm=do.call(rbind,lapply(post,function(p) p$Xc))
  Pm=pc[,round(po$t,3)*1000+1,drop=FALSE]
  Fm=matrix(NA,nrow(po),S)
  for(s in 1:S){
    f=crossprod(Pm,Mu[[s]])+rowSums(t(Pm)*t(Th[[s]]%*%A[[s]][,po$i,drop=FALSE]))+Xcm%*%Bc[[s]]
    for(j in 1:2) f=f+Xm[,j]*crossprod(Pm,Be[[s]][,j])
    Fm[,s]=f
  }
  Nz=matrix(rnorm(length(Fm)),nrow(Fm),S)*rep(sqrt(Sg),each=nrow(Fm))
  po$bin=round(po$delta*10)/10
  one=function(w,lab){
    e=po$Yobs[w]-Fm[w,,drop=FALSE]; a=colMeans(e); ap=a-colMeans(Nz[w,,drop=FALSE])
    data.frame(bin=lab,n=sum(w),n_units=length(unique(po$i[w])),truth=mean(effect*po$delta[w]),
      noise_mean=mean(po$e[w]),
      est=mean(a),est_median=median(a),post_sd=sd(a),cr_lo=quantile(a,.025,names=FALSE),cr_hi=quantile(a,.975,names=FALSE),
      pi_lo=quantile(ap,.025,names=FALSE),pi_hi=quantile(ap,.975,names=FALSE),
      pred_sd=sd(ap),cf_mse=mean((po$Y0[w]-po$e[w]-rowMeans(Fm[w,,drop=FALSE]))^2))
  }
  r=list()
  for(b in bins){w=abs(po$bin-b)<1e-8; if(any(w)) r[[length(r)+1]]=one(w,sprintf("%.1f",b))}
  r[[length(r)+1]]=one(rep(TRUE,nrow(po)),"all")
  out=do.call(rbind,r); out$sigma2=mean(Sg); out$n_obs_total=sum(sapply(d,function(x) length(x$Y)))
  attr(out,"draws")=NULL; out
}
