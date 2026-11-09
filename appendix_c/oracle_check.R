# Known components: check bias formula B_N, variance identity, and exact coverage claims.
set.seed(20261104)
f=function(t) cbind(sqrt(2)*sin(2*pi*t),sqrt(2)*cos(2*pi*t))
D=diag(c(2,1)); s2=1; z=qnorm(.975)
run=function(Ntr,m,nu,reps=4000){
  out=matrix(NA,reps,7)
  # design fixed across replicates, so results are conditional on it as in the theorem
  pre=lapply(1:Ntr,function(i) sort(runif(m,0,.6)))
  post=lapply(1:Ntr,function(i) sort(runif(1+rbinom(1,1,.5),.6,1)))
  nO=sum(lengths(post)); Binv=numeric(Ntr); vc=0; B=0; Ws=list(); Rs=list(); cs=list()
  for(i in 1:Ntr){F_=f(pre[[i]]); M=crossprod(F_)+s2*solve(D); Mi=solve(M)
    Ws[[i]]=Mi%*%t(F_); Rs[[i]]=Mi%*%(s2*solve(D)); cs[[i]]=colSums(f(post[[i]]))
    vc=vc+as.numeric(t(cs[[i]])%*%(s2*Mi)%*%cs[[i]]); B=B+sum(cs[[i]]*(Rs[[i]]%*%nu))}
  B=B/nO; vc=vc/nO^2; v=vc+s2/nO
  # identity check: s2 W W' + R D R' = s2 M^{-1}
  i=1;Mi=Rs[[i]]%*%D/s2; id=max(abs(s2*Ws[[i]]%*%t(Ws[[i]])+Rs[[i]]%*%D%*%t(Rs[[i]])-s2*Mi))
  for(r in 1:reps){ tot=0; epsbar=0
    for(i in 1:Ntr){a=nu+sqrt(diag(D))*rnorm(2); ep=rnorm(m,0,sqrt(s2)); Fp=f(pre[[i]])
      ahat=Ws[[i]]%*%(Fp%*%a+ep); ej=rnorm(length(post[[i]]),0,sqrt(s2))
      tot=tot+sum(f(post[[i]])%*%(a-ahat))+sum(ej); epsbar=epsbar+sum(ej)}
    err=tot/nO; epsbar=epsbar/nO
    out[r,]=c(err,err-epsbar,abs(err)<=z*sqrt(v),abs(err)<=z*sqrt(vc),abs(err-epsbar)<=z*sqrt(vc),0,0)}
  cat(sprintf("Ntr=%3d m=%2d nu=%.1f |O|=%3d N*-ish | B formula %.4f emp %.4f (se %.4f) | var formula %.5f emp %.5f | cov pred/theta %.3f [theory %.3f] | cov cred/theta %.3f [theory %.3f] | cov cred/realized %.3f [theory %.3f] | id %.1e\n",
   Ntr,m,nu[1],nO,B,mean(out[,1]),sd(out[,1])/sqrt(reps),v,var(out[,1]),mean(out[,3]),pnorm(z-B/sqrt(v))-pnorm(-z-B/sqrt(v)),
   mean(out[,4]),pnorm((z*sqrt(vc)-B)/sqrt(v))-pnorm((-z*sqrt(vc)-B)/sqrt(v)),mean(out[,5]),pnorm(z-B/sqrt(vc))-pnorm(-z-B/sqrt(vc)),id))
}
for(nu in list(c(0,0),c(1,1))) for(N in c(20,80,320)) for(m in c(1,2,4,8,16,32)) run(N,m,nu,reps=if(N==320) 1500 else 3000)
