#!/usr/bin/env Rscript
# Figures 1 and 2 of the simulation section from fig/rep_*.rds (scenario B, N_tr=40, N_co=100).
# Figure 3 of the paper is made by simulation/make_figure3.R.
library(splines); library(MASS)
source("prior_mean_zero/sim_functions.R")
params=make_params()
fs=sort(Sys.glob("fig/rep_*.rds")); R=length(fs); cat("replicates:",R,"\n")
tg=seq(0,1,.01)
# truth on the same grid; the 20 grid times fix the knots
gr=seq(0,1,1/19); bs_=basis_setup_sparse(list(list(gr)),3,delta=1/1000)
Bt=t(bs_[[5]][[1]][,round(tg*1000)+1,drop=FALSE])
mu0=as.vector(Bt%*%params[[7]]); F0=Bt%*%params[[6]]
cat("int f1^2, f2^2, f1f2:",round(c(mean(F0[,1]^2),mean(F0[,2]^2),mean(F0[,1]*F0[,2])),3),"\n")
MU=matrix(NA,R,101); F1=F2=matrix(NA,R,101); sig=numeric(R); att=res=vector("list",R)
for(r in 1:R){ x=readRDS(fs[r]); MU[r,]=x$mu; sig[r]=x$sigma2
  C=x$Bg%*%x$cov_coef%*%t(x$Bg); e=eigen((C+t(C))/2,symmetric=TRUE)
  Fh=e$vectors[,1:2]%*%diag(sqrt(pmax(e$values[1:2],0)))
  sv=svd(crossprod(Fh,F0)); Fa=Fh%*%(sv$u%*%t(sv$v))
  F1[r,]=Fa[,1]; F2[r,]=Fa[,2]
  att[[r]]=cbind(iter=r,x$att[x$att$bin!="all",c("bin","n","truth","est")]); res[[r]]=cbind(iter=r,x$res)
}
cat("rms error mean curve:",round(sqrt(mean((colMeans(MU)-mu0)^2)),3)," f1:",round(sqrt(mean((colMeans(F1)-F0[,1])^2)),3)," f2:",round(sqrt(mean((colMeans(F2)-F0[,2])^2)),3),"\n")
saveRDS(list(tg=tg,mu0=mu0,F0=F0,MU=MU,F1=F1,F2=F2),"fig_curves.rds")
col_run=rgb(.45,.45,.45,.14); col_true="#D62728"; col_mean="#000000"
ink="#333333"
png("Figure1.png",width=2000,height=1250,res=300,type="cairo")
par(mar=c(4,4,2.2,1),mgp=c(2.4,.7,0),col.axis=ink,col.lab=ink,col.main=ink,fg=ink,cex.axis=.9,cex.main=1,font.main=1,las=1)
plot(NA,xlim=c(0,1),ylim=range(MU,mu0),xlab="Time",ylab=expression(mu(t)),main="Mean function",bty="l")
matlines(tg,t(MU),col=col_run,lty=1,lwd=.6)
lines(tg,colMeans(MU),col=col_mean,lwd=2); lines(tg,mu0,col=col_true,lwd=2,lty=2)
legend("topright",c("True value","Average of fits","Individual fits"),col=c(col_true,col_mean,"grey55"),lwd=c(2,2,1),lty=c(2,1,1),bty="n",cex=.85,seg.len=2.4)
dev.off()
png("Figure2.png",width=2600,height=1250,res=300,type="cairo")
par(mfrow=c(1,2),mar=c(4,4,2.2,1),mgp=c(2.4,.7,0),col.axis=ink,col.lab=ink,col.main=ink,fg=ink,cex.axis=.9,cex.main=1,font.main=1,las=1)
yl=range(F1,F2,F0); yl=c(max(yl[1],-3.5),min(yl[2],3.5))
for(j in 1:2){ Fm=if(j==1) F1 else F2
  plot(NA,xlim=c(0,1),ylim=yl,xlab="Time",ylab=bquote(f[.(j)](t)),main=paste0("FPC function ",j),bty="l")
  matlines(tg,t(Fm),col=col_run,lty=1,lwd=.6)
  lines(tg,colMeans(Fm),col=col_mean,lwd=2); lines(tg,F0[,j],col=col_true,lwd=2,lty=2)
  if(j==1) legend("top",c("True value","Average of fits","Individual fits"),col=c(col_true,col_mean,"grey55"),lwd=c(2,2,1),lty=c(2,1,1),bty="n",cex=.72,seg.len=2.4)
}
dev.off()
cat("mean sigma2",round(mean(sig),3),"\n")
