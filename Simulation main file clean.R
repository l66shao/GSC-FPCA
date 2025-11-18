rm(list=ls())
#list.files()

##############################
##  SIMULATION PARAMETERS   ##
##############################
#library(Rcpp)
#library(RcppArmadillo)

library(Rcpp)
library(RcppArmadillo)
library(splines)
library(MASS)
library(Matrix)
#library(MCMCpack)
setwd("/Users/lucyshao/Downloads/folder for Wes FPCA-SC updated/")
source("mfpca_with_both_cov.R")
source("simulate_data_sparse_both_cov.R")
source("basis_setup_sparse_Lucy.R")
Rcpp::sourceCpp("FPCA_with_both_cov_with_diag.cpp")


Ns=c(200,300,150)
N_T_means=c(2,4,7)
N=100;N_T_mean=7;
for(N in Ns){
  for(N_T_mean in N_T_means){
          
          set.seed(123)
          T_min=0;T_max=1;T_range=c(T_min,T_max)
          N_T_max=20
          P=1#y dimension, always 1 herre
          nknots=c(3)
          params=list()
          orth=FALSE
          
          P_true=2#number of time varying covariates

          Q_true=c(2) # number of modes of variation(FPC)
          D_true=diag(rep(1,sum(Q_true)));D_true[1,1]=1; #correlations of FPCA scores
          SIGMA_ALPHA_true=D_true # Variance-covariance matrix of FPCA scores
          SIGMA_OMEGA_true=1*rep(1,P) # Error variances
          #MU_true=rep(0,sum(nknots+4)) # Coefficients of basis functions for mean curves across subjects
        
          Nt=100 #number of treated
          
          MU_true=mvrnorm(1,mu=rep(0,nknots+4),Sigma =diag(nknots+4))
          beta=matrix(rep(1,P_true*(nknots+4)),ncol = P_true)#this is fixed
          beta_c=c(1,1)#beta_c is non-time varying coefficients for X_c
          
   #generate cos sin based FPC function       
          
          x <- seq(0, 2*pi, length.out = 7)
          
          # Sin and cosine vectors
          sin_vec <- sin(x)
          cos_vec <- cos(x)
          x=x/(2*pi)
          # Print
          sin_vec/5;cos_vec/5
          norm(matrix(cos_vec))
          THETA_true=cbind(sin_vec/norm(sin_vec,"2"),cos_vec/(norm(cos_vec,"2")))
          (t(THETA_true[,1:2])%*%THETA_true[,1:2])
          
          time=seq(T_range[1],T_range[2],(T_range[2]-T_range[1])/(N_T_max-1))
          time_sparse=list()
          for(i in 1:N){
            time_sparse[[i]]=list()
            for(p in 1:P){
              time_sparse[[i]][[p]]=time
            }
          }
          basis_stuff=basis_setup_sparse(time_sparse,nknots,delta=1/1000,orth=TRUE)
          knots=basis_stuff[[1]]
          time_sparse=basis_stuff[[2]]
          time_sparse_combined=basis_stuff[[3]]
          phi_t=basis_stuff[[4]]
          time_cont=basis_stuff[[5]]
          phi_t_cont=	basis_stuff[[6]]
          Phi_D=basis_stuff[[7]]
          Phi_D_T=basis_stuff[[8]]
          
          FPC_true=array(0,dim=c(N_T_max,sum(Q_true)))
          for(t in 1:N_T_max){
            ind=0
            for(q in 1:length(Q_true)){
              FPC_true[t,(ind+1):(ind+Q_true[q])]=(Phi_D[[1]][[t]]%*%THETA_true)[q,(ind+1):(ind+Q_true[q])]
              ind=ind+Q_true[q]
            }	
          }
          par(mfrow=c(1,1))
          
          params[[1]]=nknots
          params[[2]]=Q_true
          params[[3]]=D_true
          params[[4]]=D_true
          params[[5]]=SIGMA_OMEGA_true
          params[[6]]=THETA_true
          params[[7]]=MU_true
          params[[8]]=FPC_true
          params[[9]]=beta
          params[[10]]=beta_c
          
          
          save(file=paste("true_params.Rdata"),params)
          
          #######################
          #######################
          ##   GENERATE DATA   ##
          #######################
          #######################
          
          
          results=list()
          
          for( iters in 1:1000){
                set.seed(iters)
                
                TIME_SPARSE=list()
                Y_SPARSE=list()
                MU_SPARSE=list()
                F_SPARSE=list()
                OMEGA_SPARSE=list()
                ALPHA=list()
                X_SPARSE=list()
                X_c_SPARSE=list()
                id_SPARSE=list()
                nsims=1
                for(ii in 1:nsims){
                    simdata=simulate_data_sparse(T_range=T_range,N_T_max=N_T_max,N_T_mean=N_T_mean,P,N,params,orth=orth)
                    TIME_SPARSE[[ii]]=simdata[[1]]
                    Y_SPARSE[[ii]]=simdata[[2]]
                    MU_SPARSE[[ii]]=simdata[[3]]
                    F_SPARSE[[ii]]=simdata[[4]]
                    OMEGA_SPARSE[[ii]]=simdata[[5]]
                    ALPHA[[ii]]=simdata[[6]]
                    X_SPARSE[[ii]]=simdata[[9]]
                    X_c_SPARSE[[ii]]=simdata[[10]]
                    id_SPARSE[[ii]]=simdata[[11]]
                }
                
                #####################################################
                
                ii=sample(N,1)
                ss=sample(nsims,1)
                
                ## Plot FPC_trues and individual mean Bold responses
                time_sparse=TIME_SPARSE[[ss]]
                Mu_sparse=MU_SPARSE[[ss]]
                F_sparse=F_SPARSE[[ss]]
                Y_sparse=Y_SPARSE[[ss]]
                
                time_sparse_i=time_sparse[[ii]]
                Mu_sparse_i=Mu_sparse[[ii]]
                F_sparse_i=F_sparse[[ii]]
                Y_i=Y_sparse[[ii]]
                
                Y1=as.vector(Y_sparse[[1]][[1]])
                id1=rep(1,length(Y1))
                times1=time_sparse[[1]][[1]]

                for(i in 2:N){
                    Y1=c(Y1,as.vector(Y_sparse[[i]][[1]]))
                    id1=c(id1,rep(i,length(as.vector(Y_sparse[[i]][[1]]))))
                    times1=c(times1,time_sparse[[i]][[1]])

                }
                data1=cbind.data.frame(id=id1,time=times1,Y=Y1)
            
                
                par(mfrow=c(2,2))
            
                #############################
                #############################
                ##  RUN SFPCA MODEL   ##
                #############################
                #############################

                
              
           
                
                #----train test split--------
              
    
                
                #split into train and test set for the treated units
                train_Y=Y_SPARSE;train_T=TIME_SPARSE;train_X=X_SPARSE;train_X_c=X_c_SPARSE
                train_id=id_SPARSE
                
                test_Y=Y_SPARSE;test_T=TIME_SPARSE;test_X=X_SPARSE;test_X_c=X_c_SPARSE;
                test_id=id_SPARSE
                
                T_start=ceiling(N_T_mean/2)# treatment starting time
            
                ids=1:N
                b_control=!(ids%in%((N-Nt+1):N))
                #b_control=rep(1,N)
                for(i in (N-Nt+1):N){
                      T_start=length(TIME_SPARSE[[1]][[i]][[1]])-2
                      tmp=Y_SPARSE[[1]][[i]][[1]]
                      train_Y[[1]][[i]][[1]]=matrix(tmp[1:T_start])
                      test_Y[[1]][[i]][[1]]=(tmp[(T_start+1):length(tmp)])
                      
                      tmp2=TIME_SPARSE[[1]][[i]][[1]]
                      train_T[[1]][[i]][[1]]=matrix(tmp2[1:(T_start)])
                      test_T[[1]][[i]][[1]]=(tmp2[(T_start+1):length(tmp2)])
                      
                      tmp5=id_SPARSE[[1]][[i]][[1]]
                      train_id[[1]][[i]][[1]]=matrix(tmp5[1:(T_start)])
                      test_id[[1]][[i]][[1]]=(tmp5[(T_start+1):length(tmp2)])
                      
                      tmp3=X_SPARSE[[1]][[i]][[1]]
                      train_X[[1]][[i]][[1]]=(tmp3[1:(T_start),, drop = FALSE])
                      test_X[[1]][[i]][[1]]=(tmp3[(T_start+1):length(tmp2),, drop = FALSE])
                      
                      tmp4=X_c_SPARSE[[1]][[i]][[1]]
                      train_X_c[[1]][[i]][[1]]=(tmp4[1:(T_start),, drop = FALSE])
                      test_X_c[[1]][[i]][[1]]=(tmp4[(T_start+1):length(tmp2),, drop = FALSE])
                      
                }
                
                source("mfpca_with_both_cov.R")
                
                
                #-------starting simulation----------
                nsims=length(train_Y)
                sim_results=list()
                ii=1
                for(ii in 1:nsims){
                      print(paste("iter:",iters))
                      Y_sparse=train_Y[[ii]]
                      time_sparse=train_T[[ii]]
                      X_sparse=train_X[[ii]]
                      X_c_sparse=train_X_c[[ii]]
                      N=length(Y_sparse)
                      group=rep(1,N)
                      Q=c(2)
                      nknots=c(3)
                      nloop=100
                      burnin=10
                      thin=1
                      orth=FALSE
                      sim=TRUE
                      
                      sim_results[[ii]]=mfpca_1d(Y_sparse=Y_sparse,time_sparse=time_sparse,X_sparse=X_sparse,X_c_sparse=X_c_sparse,
                                                 group=group,Q,nknots=nknots,nloop=nloop,burnin=burnin,thin=thin,orth=FALSE,sim=TRUE,PSLOO=1,b_control)
                }
                #save.image(file="result_train_time_cov.Rdata")
                
                
                ######################
                ######################
                ##  CHECK RESULTS   ##
                ######################
                ######################
                
            
                
                ALPHA_array=sim_results[[1]][[1]]
                MU_array=sim_results[[1]][[2]]
                D_array=sim_results[[1]][[3]]
                #R_array=sim_results[[1]][[4]]
                SIGMA_OMEGA_array=sim_results[[1]][[5]]
                THETA_array=sim_results[[1]][[6]]
                basis_stuff=sim_results[[1]][[9]]
                BETA_array=sim_results[[1]][[7]]
                BETA_c_array=sim_results[[1]][[8]]
                knots=basis_stuff[[1]]
                time_sparse=basis_stuff[[2]]
                phi_t=basis_stuff[[3]]
                time_cont=basis_stuff[[4]]
                phi_t_cont=basis_stuff[[5]]
                #params=sim_results[[1]][[8]]
                nknots=params[[1]]
                Q_true=params[[2]]
                #D_true=params[[3]]
                D_true=params[[4]]
                SIGMA_OMEGA_true=params[[5]]
                THETA_true=params[[6]]
                MU_true=params[[7]]
                FPC_true=params[[8]]
                BETA_true=params[[9]]
                BETA_c_true=params[[10]]
                nloop=length(ALPHA_array)
                first=1
                last=nloop
                
                #plot mean functions
                
             
                MU_mean=MU_array[[1]]
                Beta_mean=BETA_array[[1]]
                for(iter in 2:length(MU_array)){
                MU_mean=MU_mean+MU_array[[iter]]
                Beta_mean=Beta_mean+BETA_array[[iter]]
                }
                MU_mean=cbind(MU_mean/(last-first+1))
                Beta_mean=cbind(Beta_mean/(last-first+1))
                Mu_true_functions=t(cbind(phi_t_cont[[1]]))%*%MU_true
                Mu_functions=t(cbind(phi_t_cont[[1]]))%*%MU_mean
                Mu1_true=Mu_true_functions[1:length(time_cont)]
                
                Mu1=Mu_functions[1:length(time_cont)]
                
                par(mfrow=c(1,1))

                
                #----alpha mean----
                #plot mean functions
                ALPHA_mean=ALPHA_array[[1]]
                for(iter in 2:length(MU_array)){
                ALPHA_mean=ALPHA_mean+ALPHA_array[[iter]]
                }
                ALPHA_mean=cbind(ALPHA_mean/(last-first+1))
                
                #plot PC functions
                
                THETA_mean=THETA_array[[first]]
                for(iter in 2:nloop){
                THETA_mean=THETA_mean+THETA_array[[iter]]
                }
                THETA_mean=cbind(THETA_mean/(last-first+1))
                
                FPC1_mean=t(phi_t_cont[[1]])%*%THETA_mean
                
     
                #-------get parameter matrix---------
                
                df_tmp <- do.call(rbind, lapply(MU_array, t))
                df_mu=df_tmp
                df_MU_summary <- rbind(
                  mean=apply(df_tmp, 2, mean),        # Row 1: mean
                  median=apply(df_tmp, 2, median),      # Row 2: median
                  Q0.025=apply(df_tmp, 2, quantile, probs = 0.025),  # Row 3: 2.5% quantile
                  Q0.975=apply(df_tmp, 2, quantile, probs = 0.975) ,  # Row 4: 97.5% quantile
                  sd=apply(df_tmp, 2, sd)
                )
                
                THETA_array2=lapply(THETA_array, as.vector)
                df_tmp <- do.call(rbind, lapply(THETA_array2, t))
                df_theta=df_tmp
                df_theta_summary <- rbind(
                  mean=apply(df_tmp, 2, mean),        # Row 1: mean
                  median=apply(df_tmp, 2, median),      # Row 2: median
                  quantile0.025=apply(df_tmp, 2, quantile, probs = 0.025),  # Row 3: 2.5% quantile
                  quantile0.975=apply(df_tmp, 2, quantile, probs = 0.975)  ,  # Row 4: 97.5% quantile
                  sd=apply(df_tmp, 2, sd)
                )
                
                ALPHA_array2=lapply(ALPHA_array, t)
                ALPHA_array2=lapply(ALPHA_array2, as.vector)
                df_tmp <- do.call(rbind, lapply(ALPHA_array2, t))
                df_alpha=df_tmp
                df_alpha_summary <- rbind(
                  mean=apply(df_tmp, 2, mean),        # Row 1: mean
                  median=apply(df_tmp, 2, median),      # Row 2: median
                  quantile0.025=apply(df_tmp, 2, quantile, probs = 0.025),  # Row 3: 2.5% quantile
                  quantile0.975=apply(df_tmp, 2, quantile, probs = 0.975)   ,  # Row 4: 97.5% quantile
                  sd=apply(df_tmp, 2, sd)
                )
                
                
                D_array2=lapply(D_array, as.vector)
                df_tmp <- do.call(rbind, lapply(D_array2, t))
                df_D=df_tmp
                df_D_summary <- rbind(
                  mean=apply(df_tmp, 2, mean),        # Row 1: mean
                  median=apply(df_tmp, 2, median),      # Row 2: median
                  quantile0.025=apply(df_tmp, 2, quantile, probs = 0.025),  # Row 3: 2.5% quantile
                  quantile0.975=apply(df_tmp, 2, quantile, probs = 0.975)   ,  # Row 4: 97.5% quantile
                  sd=apply(df_tmp, 2, sd)
                )
                
                
                Beta_array2=lapply(BETA_array, as.vector)
                
                df_tmp <- do.call(rbind, lapply(Beta_array2, t))
                df_beta=df_tmp
                df_beta_summary <- rbind(
                  mean=apply(df_tmp, 2, mean),        # Row 1: mean
                  median=apply(df_tmp, 2, median),      # Row 2: median
                  quantile0.025=apply(df_tmp, 2, quantile, probs = 0.025),  # Row 3: 2.5% quantile
                  quantile0.975=apply(df_tmp, 2, quantile, probs = 0.975)   ,  # Row 4: 97.5% quantile
                  sd=apply(df_tmp, 2, sd)
                )
                
                
                Beta_c_array2=lapply(BETA_c_array, as.vector)
                
                df_tmp <- do.call(rbind, lapply(Beta_c_array2, t))
                df_beta_c=df_tmp
                df_beta_c_summary <- rbind(
                  mean=apply(df_tmp, 2, mean),        # Row 1: mean
                  median=apply(df_tmp, 2, median),      # Row 2: median
                  quantile0.025=apply(df_tmp, 2, quantile, probs = 0.025),  # Row 3: 2.5% quantile
                  quantile0.975=apply(df_tmp, 2, quantile, probs = 0.975)   ,  # Row 4: 97.5% quantile
                  sd=apply(df_tmp, 2, sd)
                )
                
                
                
                df_tmp=data.frame(unlist(SIGMA_OMEGA_array))
                df_sigma=df_tmp
                df_sigma_summary <- rbind(
                  mean=apply(df_tmp, 2, mean),        # Row 1: mean
                  median=apply(df_tmp, 2, median),      # Row 2: median
                  quantile0.025=apply(df_tmp, 2, quantile, probs = 0.025),  # Row 3: 2.5% quantile
                  quantile0.975=apply(df_tmp, 2, quantile, probs = 0.975)   ,  # Row 4: 97.5% quantile
                  sd=apply(df_tmp, 2, sd)
                )
                
                
                param_est_list=list(mu=df_MU_summary,beta=df_beta_summary,beta_c=df_beta_c_summary,theta=df_theta_summary,ALPHA=df_alpha_summary,D=df_D_summary,sigma=df_sigma_summary,phi_cont=phi_t_cont)
                
                #-----predict att-----------
                #time_cont
                
                
             
                
                #--------unlist all data to one df--------
              
                
                #train first
                df_train=data.frame(unlist(train_Y),unlist(train_T),unlist(train_id))
                df_train$treat=0
                df_train$X1_c=unlist(lapply(train_X_c[[1]], function(x) x[[1]][,1]))
                df_train$X2_c=unlist(lapply(train_X_c[[1]], function(x) x[[1]][,2]))
                df_train$X1=unlist(lapply(train_X[[1]], function(x) x[[1]][,1]))
                df_train$X2=unlist(lapply(train_X[[1]], function(x) x[[1]][,2]))
                colnames(df_train)[1:4]=c("Y","T","id","treat")
                df_test=data.frame(unlist(test_Y[[1]][(N-Nt+1):N]),unlist(test_T[[1]][(N-Nt+1):N]),unlist(test_id[[1]][(N-Nt+1):N]))
                df_test$treat=1
                df_test$X1_c=unlist(lapply(test_X_c[[1]][(N-Nt+1):N], function(x) x[[1]][,1]))
                df_test$X2_c=unlist(lapply(test_X_c[[1]][(N-Nt+1):N], function(x) x[[1]][,2]))
                df_test$X1=unlist(lapply(test_X[[1]][(N-Nt+1):N], function(x) (x[[1]])[,1]))
                df_test$X2=unlist(lapply(test_X[[1]][(N-Nt+1):N], function(x) x[[1]][,2]))
                
                colnames(df_test)[1:4]=c("Y","T","id","treat")
                
                df=rbind(df_train,df_test)
                
                cov_t_df=df[,c("X1","X2")]
                cov_c_df=df[,c("X1_c","X2_c")]
                
                
                df$fitted_y=NA
                df$fitted_y_upper=NA
                df$fitted_y_lower=NA
                df$fitted_y_median=NA
                fit_summary=matrix(NA,nrow=nrow(df),ncol=length(ALPHA_array))
                phi_cont=phi_t_cont
                
                for(i in 1:nrow(df)){
                  id_i=as.numeric(df$id[i])
                  alpha_i=cbind(df_alpha[, id_i],df_alpha[, id_i+N])
                  #alpha_i=cbind(df_alpha[,id_i])
                  fitted1=t(phi_cont[[1]][,round(df$T[i],3)*1000+1])%*%t(df_theta[,1:(nknots+4)])*alpha_i[,1]
                  fitted1=fitted1+t(phi_cont[[1]][,round(df$T[i],3)*1000+1])%*%t(df_theta[,8:14])*alpha_i[,2]
                  
                  fitted1=fitted1+t(phi_cont[[1]][,round(df$T[i],3)*1000+1])%*%t(df_mu[,1:(nknots+4)])
                  
                  for(j in 1:ncol(cov_t_df)){
                    fitted1=fitted1+as.numeric(cov_t_df[i,j])*t(phi_cont[[1]][,round(df$T[i],3)*1000+1])%*%t(df_beta[,(1+(j-1)*(nknots+4)):(j*(nknots+4))])#1*500
                    
                  } 
                  
                  for(k in 1:ncol(cov_c_df)){
                    fitted1=fitted1+t(as.matrix(as.numeric(cov_c_df[i,k])*(df_beta_c[,k])))
                  }
                  

                  fit_summary[i,]=fitted1
                  
                }
                
                
             
                df$fitted_y=rowMeans(fit_summary)
                df$fitted_y_lower=apply(fit_summary, 1, quantile, probs = 0.025)  # Row 3: 2.5% quantile
                df$fitted_y_upper=apply(fit_summary, 1, quantile, probs = 0.975) 
                df$residual=rowMeans(fit_summary-df$Y)
                df$residual_lower=apply(fit_summary-df$Y, 1, quantile, probs = 0.025)  # Row 3: 2.5% quantile
                df$residual_upper=apply(fit_summary-df$Y
                                        , 1, quantile, probs = 0.975) 
                df$fitted_y_median=apply(fit_summary, 1, quantile, probs = 0.5) 
                df$fitted_y_sd=apply(fit_summary, 1, sd) 
                df$iter=iters
                
                if(iters==1){
                  result_mat=as.matrix(df)
                }else{
                  result_mat=rbind(result_mat,as.matrix(df))
                }
                
                
                results_current=list()
                results_current=list(parameters=param_est_list)
                results[[iters]]= results_current
          
          }
          
          #results[[2]]
          save(file=paste("N_",N,"T_mean_",N_T_mean,"100 test switched.Rdata"),results)
          write.table(result_mat, file = paste("N_",N,"T_mean_",N_T_mean,"100 test switched.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
          

  }
}
