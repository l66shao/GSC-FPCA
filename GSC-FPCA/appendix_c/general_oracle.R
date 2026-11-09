# Known components, treated scores N(nu_i, Sigma_i) with Sigma_i != D: check variance and coverage formulas for all three intervals.
set.seed(11); f=function(t) cbind(sqrt(2)*sin(2*pi*t),sqrt(2)*cos(2*pi*t)); D=diag(c(2,1)); s2=1; z=qnorm(.975)
pi_=function(x,s) pnorm(s*z-x)-pnorm(-s*z-x)
for(cfg in list(list(N=60,ms=c(2,3,6),nu=c(0,0),Sig=diag(c(.5,.25))),list(N=60,ms=c(1,2,4),nu=c(1,-.5),Sig=diag(c(4,3))),list(N=60,ms=c(12,20),nu=c(1,1),Sig=diag(c(.5,3))))){
 N=cfg$N; m=sample(cfg$ms,N,TRUE); pre=lapply(m,function(k) sort(runif(k,0,.6))); post=lapply(1:N,function(i) runif(1+rbinom(1,1,.5),.6,1)); nO=sum(lengths(post))
 W=R=Om=cs=list(); B=vc0=vcS=0; Li=solve(D)
 for(i in 1:N){F_=f(pre[[i]]); Mi=solve(crossprod(F_)+s2*Li); W[[i]]=Mi%*%t(F_); R[[i]]=Mi%*%(s2*Li); Om[[i]]=s2*Mi; cs[[i]]=colSums(f(post[[i]]))
  B=B+sum(cs[[i]]*(R[[i]]%*%cfg$nu)); vc0=vc0+as.numeric(t(cs[[i]])%*%Om[[i]]%*%cs[[i]]); vcS=vcS+as.numeric(t(cs[[i]])%*%(s2*W[[i]]%*%t(W[[i]])+R[[i]]%*%cfg$Sig%*%t(R[[i]]))%*%cs[[i]])}
 B=B/nO; vc0=vc0/nO^2; vcS=vcS/nO^2; v0=vc0+s2/nO; vS=vcS+s2/nO; reps=20000; a=b=c_=e=numeric(reps); L=chol(cfg$Sig)
 for(r in 1:reps){tot=0;ebar=0; for(i in 1:N){al=cfg$nu+as.vector(t(L)%*%rnorm(2)); ep=rnorm(m[i]); ah=W[[i]]%*%(f(pre[[i]])%*%al+ep); ej=rnorm(length(post[[i]])); tot=tot+sum(f(post[[i]])%*%(al-ah))+sum(ej); ebar=ebar+sum(ej)}
  err=tot/nO; eb=ebar/nO; e[r]=err; a[r]=abs(err)<=z*sqrt(v0); b[r]=abs(err-eb)<=z*sqrt(vc0); c_[r]=abs(err)<=z*sqrt(vc0)}
 cat(sprintf("B %.3f emp %.3f | var %.4f emp %.4f | pred/theta %.3f th %.3f | fit/realized %.3f th %.3f | fit/theta %.3f th %.3f\n",B,mean(e),vS,var(e),mean(a),pi_(B/sqrt(vS),sqrt(v0/vS)),mean(b),pi_(B/sqrt(vcS),sqrt(vc0/vcS)),mean(c_),pi_(B/sqrt(vS),sqrt(vc0/vS))))}
