mfpca_1d=function(Y_sparse,time_sparse,X_sparse,X_c_sparse,group,Q,nknots,nloop,burnin,thin,orth=FALSE,PSLOO=1,b_control,sim=TRUE,path="sim_data_rho_75.Rdata"){


  N=length(Y_sparse)
  P=1
  
  basis_stuff=basis_setup_sparse(time_sparse,nknots,delta=1/1000,plotit=FALSE,orth=FALSE)
  knots=basis_stuff[[1]]#number of knots 
  phi_t=basis_stuff[[3]]# basis at each observed time point
  time_cont=basis_stuff[[4]]#consecutive times
  phi_t_cont=basis_stuff[[5]]# basis at consecutive times
  K=nknots+4#number of spline basis for each dimension of Y, so dimension 1*P
  
  P=length(Y_sparse[[1]])#dimension of Y, always 1 here
  N=length(Y_sparse)#number of subjects
  id=1:N
  G=length(unique(group))#not using group here
  N_g=sum(group==unique(group)[1])
  if(G>1){
    for(g in 2:G){
      N_g=c(N_g,sum(group==unique(group)[g]))
    }
  }
  
  ## initialize parameters   
  
  Alpha=array(0,dim=c(sum(Q),N)) # FPCA scores for all N subjects
  Mu=cbind(array(0,dim=c(sum(K),1))) # Basis coefficients for mean curves
  Theta=array(0,dim=c(sum(K),sum(Q))) # FPCA loadings
  
  
  #building a block of FPCA loadings with orthonormal columns via a manual Gram–Schmidt 
  #then placing it into a big theta matrix
  for(p in 1){
    Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),sum(Q[min(1,p-1):(p-1)])+1]=mvrnorm(n=1,rep(0,K[p]),diag(rep(1,K[p])))
    Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),sum(Q[min(1,p-1):(p-1)])+1]=
      Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),sum(Q[min(1,p-1):(p-1)])+1]/
      sqrt(sum(Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),sum(Q[min(1,p-1):(p-1)])+1]^2))
    if(Q[p]>=2) for(q in (sum(Q[min(1,p-1):(p-1)])+2):sum(Q[1:p])){  # guard added so that k = 1 runs
      Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]=mvrnorm(n=1,rep(0,K[p]),diag(rep(1,K[p])))
      for(r in 1:(q-sum(Q[min(1,p-1):(p-1)])-1)){
        Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]=Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]-
          sum(Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]*Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q-r])/
          sum(Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q-r]^2)*Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q-r]
        Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]=Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]/
          sqrt(sum(Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]^2))
      }        
    }
  }	
  
  Beta=matrix(0,ncol=ncol(X_sparse[[1]][[1]]),nrow =( nknots+4))
  
  D=diag(rep(1,sum(Q))) #Variance of FPCA scores
  Sigma_alpha=D #Variance-Covariance of FPCA scores

  Sigma_omega=0.1*rep(1,P) # Variances of errors
  #eta_eps=1
  #S_eps=.1*diag(rep(1,P))
  #MCMC parameters?
  c_omega=as.numeric(Sys.getenv("COMEGA","1"))
  psi=0
  nu=1000^.5
  delta=1000^.5
  eta1_eps=.01
  eta2_eps=.01
  array_ind=0
  #The following are initiaizations
  Alpha_mean=Alpha
  Mu_mean=Mu
  Theta_mean=Theta
  D_mean=D
  Sigma_omega_mean=Sigma_omega
  Beta_mean=Beta
  
  
  Y_sparse2 <- lapply(Y_sparse, function(x) x[[1]])# Extract the matrix from each nested list
  time_sparse2=lapply(time_sparse, function(x) x[[1]])  # Extract the matrix from each nested list
  X_sparse2=lapply(X_sparse, function(x) x[[1]])
  X_c_sparse2=lapply(X_c_sparse, function(x) x[[1]])
  phi_t2=phi_t[[1]]  # Extract the matrix from each nested list
  P_X_c=ncol(X_c_sparse[[1]][[1]])
  Beta_c=array(0,P_X_c)
  
  # set up arrays
  ALPHA_array=NULL
  MU_array=NULL
  D_array=NULL
  R_array=NULL
  THETA_array=NULL
  SIGMA_OMEGA_array=NULL
  BETA_array=NULL
  BETA_c_array=NULL
  
  est_mu=1
  est_alpha=1
  est_d=1
  est_r=1
   est_sigma_omega=1
   est_theta=1
   est_beta=1
 
   est_beta_c=1
  
  #Rcpp::sourceCpp("FPCA_with_both_cov_with_diag.cpp")
  mcmc_results=run_mcmc(nloop, burnin,thin,Y_sparse2, time_sparse2, X_sparse2,X_c_sparse2,
                        phi_t2, Mu, Beta,Beta_c, Alpha, Theta, 
                        D, Sigma_omega, Sigma_alpha, 
                        c_omega,est_mu=1,est_alpha, est_d, 
                        est_sigma_omega,est_theta,est_beta=1,est_beta_c,PSLOO,b_control)
  #mcmc_results[[7]]
 mcmc_results[[12]]=mcmc_results[[9]]  # keep the pointwise log likelihood for LOOIC
 mcmc_results[[9]]=basis_stuff
 if(sim==TRUE){
   mcmc_results[[10]]=params
 }	
  return (mcmc_results)
}
