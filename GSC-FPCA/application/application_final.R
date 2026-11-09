#Binge drinking application clean

rm(list=ls())
list.files()

library(Matrix)
library(MASS)
library(RcppArmadillo)
library(Rcpp)
library(splines)
library(lme4)
#library(lmerexposed) 
library(ggplot2)
library(nlme)
pdf(NULL)
APP_Q=as.numeric(Sys.getenv("APP_Q","1")); APP_NLOOP=as.numeric(Sys.getenv("APP_NLOOP","1000")); APP_BURN=as.numeric(Sys.getenv("APP_BURN","500")); APP_THIN=as.numeric(Sys.getenv("APP_THIN","1")); APP_SEED=as.numeric(Sys.getenv("APP_SEED","123"))
APP_S=floor((APP_NLOOP-APP_BURN-1)/APP_THIN)+1
APP_ROOT=Sys.getenv("APP_ROOT","."); APP_OUT=Sys.getenv("APP_OUT","out"); setwd(APP_ROOT)
source("basis_setup_sparse_Lucy.R")



library(dplyr)
set.seed(APP_SEED)
#----------read data------------

dat <- read.csv(Sys.getenv("NCANDA_DAT", "data/ncanda_dat.csv"))
# Outcome: true left/right average of superior frontal gray volume, taken from FullData_Y9.csv via build_ext.R. In the
# version of ncanda_dat.csv used for the first posting, the superiorfrontal_grayvol column averaged the left volume with
# the left frontal pole vertex count, because the extraction script paired columns by position.
ext <- readRDS("ext_all.rds")
new_y <- ext$superiorfrontal_grayvol.new[match(paste(dat$id,dat$visit),paste(ext$id,ext$visit))]
stopifnot(all(is.na(new_y)==is.na(dat$superiorfrontal_grayvol)))
dat$superiorfrontal_grayvol <- new_y
dim(dat)
#dat=dat[dat$sexM==1,]
length(unique(dat$id))

## keep only rows with complete cases
dim(dat)
dat <- dat[complete.cases(dat),]
dim(dat)
n <- length(unique(dat$id)); n

## get rid of abnormal cases
dim(dat)
dat <- dat[!dat$abnormal == 1,]
dim(dat)
n <- length(unique(dat$id)); n


## keep only youngest age cohort(muted)
#dat=dat[dat$base_age_cat15_17 == 1,]
#dim(dat)
length(unique(dat$id))

## keep only observations <= 21 years of age
#quantile(dat$visit_age,0.95)
dat=dat[dat$visit_age >=12,]
dat=dat[dat$visit_age <=21,]
dim(dat)
length(unique(dat$id))

## keep only the exposed that were initially exposed from ages 16 and older
#dat$keep <- 1
#for(i in 1:length(unique(dat$id))){
#  dat_i = dat[dat$id == unique(dat$id)[i],]
#  if(max(dat_i$binge12)>0){
#    if(min(dat_i$visit_age[dat_i$binge12 > 0]) < 16) dat_i$keep = 0
#  }
#  dat[dat$id == unique(dat$id)[i],] = dat_i
#}
#dat <- dat[dat$keep == 1, ]
#dim(dat)
#length(unique(dat$id))

## keep only those exposed continuously after exposure start
dat$keep <- 1
for(i in 1:length(unique(dat$id))){
  dat_i = dat[dat$id == unique(dat$id)[i],]
  if(max(dat_i$binge12)>0){
    if(min(dat_i$binge[dat_i$binge12 > 0])<12) dat_i$keep = 0
  }
  dat[dat$id == unique(dat$id)[i],] = dat_i
}
dat <- dat[dat$keep == 1, ]
dat$keep <- NULL
dim(dat)
n = length(unique(dat$id)); show(n)

## compute age at exposure start

dat$age_start <- NA
for(i in 1:n){
  dat_i = dat[dat$id == unique(dat$id)[i],]
  if(max(dat_i$binge12)>0){
    dat_i$age_start = min(dat_i$visit_age[dat_i$binge12>0])
  }
  dat[dat$id == unique(dat$id)[i],] = dat_i
}
summary(dat$age_start)
hist(dat$age_start)

n <- length(unique(dat$id)); n

# get baseline age, baseline icv volume and baseline visit
dat <- dat %>%
  arrange(id, visit) %>%  # or use visit number if available
  group_by(id) %>%
  mutate(icv_baseline = first(sri24_icv_brain_volume),base_age= first(visit_age),base_visit=first(visit)) %>%
  ungroup()

dat <- dat %>%group_by(id) %>%mutate(visit_count = n()) %>%ungroup()
dim(dat)
names(dat)

dat_all <- dat

dat_all$group <- "control"
dat_all$group[dat_all$binge12==1] <- "exposed"
if(Sys.getenv("APP_PLACEBO")=="1"){
  for(i in unique(dat_all$id[dat_all$binge12==1])){
    w=which(dat_all$id==i&dat_all$binge12==0)
    if(length(w)>=2) dat_all$group[w[which.max(dat_all$visit_age[w])]]="exposed"
  }
}

# get baseline age, baseline icv volume and baseline visit
dat_all <- dat_all %>%
  arrange(id, visit) %>%  # or use visit number if available
  group_by(id) %>%
  mutate(icv_mean = mean(sri24_icv_brain_volume,na.rm=TRUE),icv_baseline = first(sri24_icv_brain_volume),base_age= first(visit_age),base_visit=first(visit)) %>%
  ungroup()

#dat_all$binge[dat_all$group=="control"]=0 # we are controlling the counterfactual to be 0 vs. 12+

dat_all$race_labelWhite <- as.numeric(dat_all$race_label == "Caucasian/White")
dat_all$race_labelAsian <- as.numeric(dat_all$race_label == "Asian")
dat_all$race_labelBlack <- as.numeric(dat_all$race_label == "African-American/Black")
dat_all$race_labelOther <- as.numeric(dat_all$race_label == "Other")

dat_all$siteB <- as.numeric(dat_all$site == "B")
dat_all$siteC <- as.numeric(dat_all$site == "C")
dat_all$siteD <- as.numeric(dat_all$site == "D")
dat_all$siteE <- as.numeric(dat_all$site == "E")

dat_all$icv_baseline <- scale(dat_all$icv_baseline)
dat_all$superiorfrontal_grayvol_original=dat_all$superiorfrontal_grayvol 
dat_all$superiorfrontal_grayvol <- scale(dat_all$superiorfrontal_grayvol)

z <- dat_all$superiorfrontal_grayvol  # currently scaled

center <- attr(z, "scaled:center")
sc     <- attr(z, "scaled:scale")



dat_all$icv_baseline_F=dat_all$icv_baseline*dat_all$sexF
dat_all$icv_baseline_M=dat_all$icv_baseline*dat_all$sexM
dat_all$zeros=0
dat_all$ones=0
control=dat_all[dat_all$group=="control",]
exposed=dat_all[dat_all$group=="exposed",]

dat_all$demo_group="control"
dat_all$demo_group[dat_all$id%in%unique(exposed$id)]="exposed"

dat_all$race_label <- relevel(factor(dat_all$race_label), ref = "Other")

#-----------lme here------------------
#check if race confounds: 
table(dat_all$race_label, dat_all$binge12)
prop.table(table(dat_all$race_label, dat_all$binge12), 1)
chisq.test(table(dat_all$race_label, dat_all$binge12))
fit=lme(superiorfrontal_grayvol~binge12+icv_baseline+sex+visit+base_age+race_label,random=~visit|id,data=dat_all)

sm=summary(fit)
sm$tTable[,c(1,2,5)]
library(dplyr)

#--------get some population feature here------------
## 1. Number of visits per id (using distinct visit values)
visits_per_id <- dat_all %>%
  group_by(id) %>%
  summarise(n_visit = n_distinct(visit), .groups = "drop")

## Overall mean and sd of visits per id
visits_per_id %>%
  summarise(
    median=median(n_visit),
    min=min(n_visit),
    max=max(n_visit),
    mean_visits = mean(n_visit),
    sd_visits   = sd(n_visit)
  )


visits_per_id_exposed <- dat_all[dat_all$id%in%exposed$id,] %>%
  group_by(id) %>%
  summarise(n_visit = n_distinct(visit), .groups = "drop")

visits_per_id_exposed %>%
  summarise(
    median=median(n_visit),
    min=min(n_visit),
    max=max(n_visit),
    mean_visits_exposed = mean(n_visit),
    sd_visits_exposed   = sd(n_visit)
  )

#-------Demographics------------

dat_all_first <- dat_all %>%
  group_by(id) %>%          # one group per subject
  slice_head(n = 1) %>%     # keep the first row in each id
  ungroup()
table_demo=dat_all_first %>%
  group_by(demo_group) %>%
  summarise(
    icv_baseline_mean = mean(icv_baseline, na.rm = TRUE),
    icv_baseline_sd   = sd(icv_baseline,   na.rm = TRUE),
    base_age_mean     = mean(base_age,     na.rm = TRUE),
    base_age_sd       = sd(base_age,       na.rm = TRUE),
    
    female_percent    = mean(sexF, na.rm = TRUE),
    female_count      = sum(sexF,  na.rm = TRUE),
    
    asian_percent     = mean(race_label == "Asian",           na.rm = TRUE),
    asian_count       = sum(race_label == "Asian",            na.rm = TRUE),
    
    black_percent     = mean(race_label == "African-American/Black",           na.rm = TRUE),
    black_count       = sum(race_label == "African-American/Black",            na.rm = TRUE),
    
    other_percent     = mean(race_label == "Other",           na.rm = TRUE),
    other_count       = sum(race_label == "Other",            na.rm = TRUE),
    
    white_percent     = mean(race_label == "Caucasian/White", na.rm = TRUE),
    white_count       = sum(race_label == "Caucasian/White",  na.rm = TRUE),
    
    .groups = "drop"
  )
t(table_demo)

#-------Select Covariates----------
#covs_t=c("icv_baseline","base_age","sexF")
#covs_c=c("race_labelAsian","race_labelBlack","race_labelOther","siteB","siteC","siteD","siteE")

covs_t=c("icv_baseline","sexF","race_labelBlack","race_labelAsian","race_labelWhite")#this is now non-time varying covariate
covs_c=c("zeros")#this is now the time varying covariate, if you don't want any, just include the variable zeros in the data
est_beta=1 #,"siteB"    ,               "siteC"    ,               "siteD"    ,             "siteE" 
est_beta_c=0 #if you want any time varying,set this to 1 ,"race_labelAsian","race_labelBlack","race_labelOther"
#
est_mu=1

control <- control %>% arrange(id, visit)
control$visit_age2=(control$visit_age-min(dat_all$visit_age))/(max(dat_all$visit_age)-min(dat_all$visit_age))
min(control$visit_age)
max(control$visit_age)

exposed <- exposed %>% arrange(id, visit)
exposed$visit_age2=(exposed$visit_age-min(dat_all$visit_age))/(max(dat_all$visit_age)-min(dat_all$visit_age))
min(exposed$visit_age)
max(exposed$visit_age)

id_names <- control %>% group_by(id) %>% summarise(n = n(), .groups = "drop") %>% pull(id)

id_list <- control %>%
  group_by(id) %>%
  summarise(id_vec = list(rep(first(id), n())), .groups = "drop") %>%
  pull(id_vec)

control_ids=sapply(id_list, `[`, 1)
grayvol_list <- control %>%
  group_by(id) %>%
  summarise(grayvol = list(matrix(superiorfrontal_grayvol,ncol=1)), .groups = "drop") %>%
  pull(grayvol)

exposed_id=unique(exposed$id)

b_control=as.numeric(!id_names%in%exposed_id) #whether a control:yes 1
Y_sparse2_b=grayvol_list

Y_sparse_b = control %>%
  group_by(id) %>%
  summarise(grayvol = list(list(matrix(superiorfrontal_grayvol,ncol=1))), .groups = "drop") %>%
  pull(grayvol)
names(Y_sparse2_b)=id_names
Y_sparse_b_control=list()


#------basis set up for age-----------------

visit_list_control <- control %>%
  group_by(id) %>%
  summarise(time = list(list(visit_age2)), .groups = "drop") %>%
  pull(time)
time_sparse <- visit_list_control
visit_list_exposed <- exposed %>%
  group_by(id) %>%
  summarise(time = list(list(visit_age2)), .groups = "drop") %>%
  pull(time)

#visit_list=append(visit_list_control,visit_list_exposed)#1:785 control, 786:1031 exposed
visit_list=visit_list_control
N=length(Y_sparse2_b)
nknots=3

basis_stuff=basis_setup_sparse(visit_list,nknots,delta=1/1000,plotit=FALSE,orth=FALSE)

knots=basis_stuff[[1]]#number of knots is 1 for Y1 for now
phi_t=basis_stuff[[3]]# basis at each observed time point?
time_cont=basis_stuff[[4]]#consecutive times
phi_t_cont=basis_stuff[[5]]# basis at consecutive times
K=nknots+4#number of spline basis for each dimension of Y, so dimension 1*P

visit_list2 <- control %>%
  group_by(id) %>%
  summarise(time = list((visit_age2)), .groups = "drop") %>%
  pull(time)
time_sparse2_b=visit_list2

#sri24_icv_brain_volume,race_labelAsian,race_labelCaucasian.White
cov_list <- control %>%
  group_by(id) %>%
  summarise(time = list(as.matrix(across((covs_t)))), .groups = "drop") %>%
  pull(time)


control$intercept_x=1

#sri24_icv_brain_volume
cov_list_c <- control %>%
  group_by(id) %>%
  summarise(time = list(as.matrix(across((covs_c)))), .groups = "drop") %>%
  pull(time)

X_sparse_b=cov_list 
X_sparse_b_c=cov_list_c# non-time varying coefficient
P_X=ncol(X_sparse_b[[1]])
P_X_c=ncol(X_sparse_b_c[[1]])
Q=APP_Q#number of FPC




#--------- initialize parameters  ------------------

P=1#dimension of Y is 1
N=length(Y_sparse2_b)#number of subjects

#initialize parameters

K=nknots+4#number of spline basis for each dimension of Y, so dimension 1*P

P=length(Y_sparse_b[[1]])#dimension of Y
N=length(Y_sparse_b)#number of subjects
id=1:N

Alpha=array(0,dim=c(sum(Q),N)) # FPCA scores for all N subjects
Mu=cbind(array(0,dim=c(sum(K),1))) # Basis coefficients for mean curves
Theta=array(0,dim=c(sum(K),sum(Q))) # FPCA loadings
#this is trying to do the theta constraint, however, not zero instead replace by really small number, may be can make very small number zero

Beta_c=array(0,P_X_c)
Beta=matrix(0,nrow = (nknots+4),ncol=P_X)

D=diag(rep(1,sum(Q))) #Variance of FPCA scores
Sigma_alpha=D #Variance-Covariance of FPCA scores

Sigma_omega=0.1*rep(1,P) # Variances of errors

c_omega=1 
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
#R_mean=R
D_mean=D
Sigma_omega_mean=Sigma_omega
Beta_mean=Beta


Y_sparse2_b <- lapply(Y_sparse_b, function(x) x[[1]])# Extract the matrix from each nested list
time_sparse2=lapply(time_sparse, function(x) x[[1]])  # Extract the matrix from each nested list


#X_sparse2=lapply(X_sparse_b, function(x) x[[1]])
phi_t2=phi_t[[1]]  # Extract the matrix from each nested list

# set up arrays
ALPHA_array=NULL
MU_array=NULL
D_array=NULL
R_array=NULL
THETA_array=NULL
SIGMA_OMEGA_array=NULL
BETA_array=NULL


est_alpha=1
est_d=1
est_r=1
est_sigma_omega=1
est_theta=1

PSLOO=1


nloop=APP_NLOOP; burnin=APP_BURN;thin=APP_THIN

#-------control model----------------

Rcpp::sourceCpp("../method/sampler_tm.cpp",cacheDir="cpp_cache")
mcmc_results=run_mcmc(nloop, burnin,thin,Y_sparse2_b, time_sparse2_b, X_sparse_b,X_sparse_b_c,
                      phi_t2, Mu, Beta,Beta_c, Alpha, Theta, 
                      D, Sigma_omega, Sigma_alpha, 
                      c_omega,est_mu,est_alpha, est_d, 
                      est_sigma_omega,est_theta,est_beta,est_beta_c,PSLOO,b_control)

sim_results=list(mcmc_results)
library(loo)
loglike_mat=mcmc_results[[9]]
source("../method/loo_cut.R")
loo_result=elpd_cut(loglike_mat,b_control==1)
print(loo_result$psis)
saveRDS(loo_result,file.path(APP_OUT,"loo.rds"))
if(Sys.getenv("APP_STOP")=="loo") quit(save="no")




sim_results=list(mcmc_results)
ALPHA_array=sim_results[[1]][[10]]  # treated scores under the separate prior mean; controls unchanged
MU_array=sim_results[[1]][[2]]
D_array=sim_results[[1]][[3]]
SIGMA_OMEGA_array=sim_results[[1]][[5]]
THETA_array=sim_results[[1]][[6]]
BETA_array=sim_results[[1]][[7]]
BETA_c_array=sim_results[[1]][[8]]
basis_stuff=basis_setup_sparse(visit_list,nknots,delta=1/1000,plotit=FALSE,orth=FALSE)
knots=basis_stuff[[1]]
time_sparse=basis_stuff[[2]]
phi_t=basis_stuff[[3]]

time_cont=basis_stuff[[4]]
phi_t_cont=basis_stuff[[5]]
#params=sim_results[[1]][[8]]


nloop=length(ALPHA_array)
first=1
last=nloop


#-----plot the overall prediction--------


# Convert list into a matrix of samples: each row is a sample, each column is a component of mu
mu_matrix <- do.call(cbind, lapply(MU_array, as.vector))
mu_matrix <- t(mu_matrix)  # shape: iterations × mu dimension (3 × 7)


# Compute 95% credible intervals (2.5% and 97.5% quantiles)
CI_mu <- apply(mu_matrix, 2, quantile, probs = c(0.025, 0.975))
CI_mu

# Convert list into a matrix of samples: each row is a sample, each column is a component of mu
beta_matrix <- do.call(cbind, lapply(BETA_array, as.vector))
beta_matrix <- t(beta_matrix)  # shape: iterations × mu dimension (3 × 7)


# Convert list into a matrix of samples: each row is a sample, each column is a component of mu
beta_c_matrix <- do.call(cbind, lapply(BETA_c_array, as.vector))
beta_c_matrix <- t(beta_c_matrix)  # shape: iterations × mu dimension (3 × 7)
Beta_c_mean=colMeans(beta_c_matrix)
# Compute 95% credible intervals (2.5% and 97.5% quantiles)
CI_beta <- apply(beta_matrix, 2, quantile, probs = c(0.025, 0.975))
CI_beta





# Convert list into a matrix of samples: each row is a sample, each column is a component of mu
theta_matrix <- do.call(cbind, lapply(THETA_array, as.vector))
theta_matrix <- t(theta_matrix)  # shape: iterations × mu dimension (3 × 7)

# Compute 95% credible intervals (2.5% and 97.5% quantiles)
CI_theta <- apply(theta_matrix, 2, quantile, probs = c(0.025, 0.975))
CI_theta


# Convert list into a matrix of samples: each row is a sample, each column is a component of mu
alpha_matrix <- do.call(cbind, lapply(ALPHA_array, as.vector))
alpha_matrix <- t(alpha_matrix)  # shape: iterations × mu dimension (3 × 7)

# Compute 95% credible intervals (2.5% and 97.5% quantiles)
CI_alpha <- apply(alpha_matrix, 2, quantile, probs = c(0.025, 0.975))
CI_alpha

apply(alpha_matrix, 2, sd)#the first pc have very large SEs

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




param_est_list=list(mu=df_MU_summary,beta=df_beta_summary,theta=df_theta_summary,ALPHA=df_alpha_summary,D=df_D_summary,sigma=df_sigma_summary,phi_cont=phi_t_cont)

#--------get fitted/predicted value for all data-----------

phi_cont=phi_t_cont
theta_summary=df_theta_summary

dat_all$binge=0

#get initial treated time
dat_all$T_scaled=NA #get scaled time
dat_all$visit_scaled=NA


train_id=unique(dat_all$id[!(dat_all$id%in%exposed$id)])
for( i in unique(train_id)){
  dat_all$T_scaled[dat_all$id==i]=dat_all$visit_age[dat_all$id==i]-max(dat_all$visit_age[dat_all$id==i])
  dat_all$visit_scaled[dat_all$id==i]=dat_all$visit[dat_all$id==i]-max(dat_all$visit[dat_all$id==i])
}
table(dat_all$visit_scaled)
counter=0


#align by previous year where the exposure actually starts

table(dat_all$T_scaled)
counter=0
exposed_id=unique(dat_all$id[(dat_all$id%in%exposed$id)])

for( i in unique(exposed_id)){
  
  if(sum(dat_all$id==i&dat_all$group=="exposed")>0){
    dat_all$T_scaled[dat_all$id==i]=dat_all$visit_age[dat_all$id==i]-max(dat_all$visit_age[dat_all$id==i&dat_all$group!="exposed"])
    dat_all$visit_scaled[dat_all$id==i]=dat_all$visit[dat_all$id==i]-max(dat_all$visit[dat_all$id==i&dat_all$group!="exposed"])
  }else{
    dat_all$T_scaled[dat_all$id==i]=dat_all$visit_age[dat_all$id==i]-max(dat_all$visit_age[dat_all$id==i])
    dat_all$visit_scaled[dat_all$id==i]=dat_all$visit[dat_all$id==i]-max(dat_all$visit[dat_all$id==i])
    
  }
  
  counter=max(counter,sum(dat_all$id==i&dat_all$group=="exposed"))
}

table(dat_all$visit_scaled)
dat_all$T_scaled_group=floor((dat_all$T_scaled)/0.5)
table(dat_all$T_scaled_group)
# Replace `visit_date` with your actual visit variable if different
#dat_all$LifeBinge=0
df=dat_all

df=df[df$id%in%exposed$id,]

#----get fitted Y first?---
alpha=param_est_list$ALPHA
mu_summary=param_est_list$mu
ncol(alpha)
N=length(unique(control$id))  
df$visit_age2=(df$visit_age-min(control$visit_age))/(max(control$visit_age)-min(control$visit_age))
df$T=(df$visit_age2)
df$Binge=0
cov_t_df=df[,covs_t]
cov_c_df=df[,covs_c]


df$fitted_y=NA
df$fitted_y_upper=NA
df$fitted_y_lower=NA
df$fitted_y_median=NA
fit_summary=matrix(NA,nrow=nrow(df),ncol=APP_S)
id_map=unique(control$id)

for(i in 1:nrow(df)){
  id_i=as.numeric(df$id[i])
  id_pos=(1:length(id_map))[id_i==id_map]
  #alpha_i=cbind(df_alpha[,id_pos],df_alpha[,id_pos+ncol(alpha)/2])
  alpha_i=cbind(df_alpha[,id_pos])
  fitted1=t(phi_cont[[1]][,round(df$T[i],3)*1000+1])%*%t(df_theta[,1:(nknots+4)])*alpha_i[,1]
  #fitted1=fitted1+t(phi_cont[[1]][,round(df$T[i],3)*1000+1])%*%t(df_theta[,8:14])*alpha_i[,2]
  if(Q>=2) for(qq in 2:Q){
    fitted1=fitted1+t(phi_cont[[1]][,round(df$T[i],3)*1000+1])%*%t(df_theta[,(1+(qq-1)*(nknots+4)):(qq*(nknots+4))])*df_alpha[,id_pos+(qq-1)*ncol(alpha)/Q]
  }
  
  fitted1=fitted1+t(phi_cont[[1]][,round(df$T[i],3)*1000+1])%*%t(df_mu[,1:(nknots+4)])
if(length(covs_t)>0)  {
  for(j in 1:length(covs_t)){
    fitted1=fitted1+as.numeric(cov_t_df[i,j])*t(phi_cont[[1]][,round(df$T[i],3)*1000+1])%*%t(df_beta[,(1+(j-1)*(nknots+4)):(j*(nknots+4))])#1*500
    
  }
}
 
  
  for(k in 1:length(covs_c)){
    fitted1=fitted1+t(as.matrix(as.numeric(cov_c_df[i,k])*(df_beta_c[,k])))
  }
  
  #df$fitted_y[i]=mean(fitted1)
  #df$fitted_y_lower[i]=quantile(fitted1,0.025)
  #df$fitted_y_upper[i]=quantile(fitted1,0.975)
  #df$fitted_y_median[i]=median(fitted1)
  fit_summary[i,]=fitted1

}


nrows=length(unique(exposed_id))*dim(phi_t_cont[[1]])[2]
fit_cts=matrix(NA,nrow=nrows,ncol=5)
colnames(fit_cts)=c("id",'est',"lower","upper","T")
id_map=unique(control$id)

for(i in 1:length(unique(exposed_id))){
  id_i=unique(exposed_id)[i]
  id_pos=(1:length(id_map))[id_i==id_map]
  #alpha_i=cbind(df_alpha[,id_pos],df_alpha[,id_pos+ncol(alpha)/2])
  alpha_i=cbind(df_alpha[,id_pos])
  
  fitted1=sweep(t(phi_cont[[1]])%*%t(df_theta[,1:(nknots+4)]),2,alpha_i[,1],"*")
  if(Q>=2) for(qq in 2:Q){
    fitted1=fitted1+sweep(t(phi_cont[[1]])%*%t(df_theta[,(1+(qq-1)*(nknots+4)):(qq*(nknots+4))]),2,df_alpha[,id_pos+(qq-1)*ncol(alpha)/Q],"*")
  }
  
  fitted1=fitted1+t(phi_cont[[1]])%*%t(df_mu[,1:(nknots+4)])
  if(length(covs_t)>0)  {
    for(j in 1:length(covs_t)){
      fitted1=fitted1+as.numeric(cov_t_df[df$id==id_i,j][1,])*t(phi_cont[[1]])%*%t(df_beta[,(1+(j-1)*(nknots+4)):(j*(nknots+4))])#1*500
      
    }
  }
  
 #no non-time var for now 
#  for(k in 1:length(covs_c)){
#    fitted1=fitted1+t(as.matrix(as.numeric(cov_c_df[df$id==id_i,k])*(df_beta_c[,k])))
#  }
  
  #df$fitted_y[i]=mean(fitted1)
  #df$fitted_y_lower[i]=quantile(fitted1,0.025)
  #df$fitted_y_upper[i]=quantile(fitted1,0.975)
  #df$fitted_y_median[i]=median(fitted1)
  
  fit_cts[((i-1)*1001+1):(i*1001),1]=id_i
  fit_cts[((i-1)*1001+1):(i*1001),2]=rowMeans(fitted1)
  fit_cts[((i-1)*1001+1):(i*1001),3:4]=t(apply(fitted1, 1, quantile, probs = c(0.025,  0.975)))
  fit_cts[((i-1)*1001+1):(i*1001),5]=(0:1000)/1000
}

att_array=c()
att_upper=c()
att_lower=c()
#-4:2
df$T_rounded=floor(df$T_scaled)

table(df$visit_scaled)

#  floor(df$T_scaled)
t=-3:3
for(i in t){
  b=df$T_rounded==i
  #ATT_median= apply(fit_summary[b,]-as.data.frame(replicate(APP_S, df$superiorfrontal_grayvol[b])), 1, median)
  
  ATT_ests= -colMeans(fit_summary[b,]-as.data.frame(replicate(APP_S, df$superiorfrontal_grayvol[b])),na.rm=TRUE)
  att_array=c(att_array,mean(ATT_ests,na.rm=TRUE))
  att_lower=c(att_lower,quantile(ATT_ests,0.025,na.rm=TRUE))
  att_upper=c(att_upper,quantile(ATT_ests,0.975,na.rm=TRUE))
}





t;att_array;att_upper;att_lower
plot(t,att_array,type="l")
abline(h=0,col="red")
abline(v=0,col="red")
lines(t,att_lower,col="grey")
lines(t,att_upper,col="grey")


# Create a data frame for ggplot
plot_data_att <- data.frame(
  Time = t,
  ATT = att_array,
  Upper = att_upper,
  Lower = att_lower,
  n=table(df$T_rounded[df$T_rounded%in%t])
)

# Plot with confidence bands
p1=ggplot(plot_data_att, aes(x = Time, y = ATT)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Average Treatment Effect on Treated (ATT)",
    x = "Years before and after Binge Drinking",
    y = "ATT Estimate of Superior Frontal Gray Volume(standardized)",
    caption = "Dashed red lines indicate t=0 and ATT=0"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = t) 
#ylim(-150, 110)  # Adjust y-axis limits if needed





#--------ATT rerun: post-exposure observations only, with noise draws----------
sigma_draws=unlist(SIGMA_OMEGA_array)
stopifnot(length(sigma_draws)==ncol(fit_summary))
set.seed(APP_SEED+1)
att_row=function(b,label,xpos,kind){
  Fm=fit_summary[b,,drop=FALSE]; y=as.numeric(df$superiorfrontal_grayvol[b]); nO=length(y)
  dr=colMeans(y-Fm); dp=dr-rnorm(length(dr),0,sqrt(sigma_draws/nO))
  data.frame(label=label,x=xpos,kind=kind,n=nO,units=length(unique(df$id[b])),est=mean(dr),median=median(dr),
    fit_lo=quantile(dr,.025,names=FALSE),fit_hi=quantile(dr,.975,names=FALSE),
    pred_lo=quantile(dp,.025,names=FALSE),pred_hi=quantile(dp,.975,names=FALSE),post_sd=sd(dr),pred_sd=sd(dp),post_prob_neg=mean(dp<0))
}
placebo=Sys.getenv("APP_PLACEBO")=="1"
if(placebo){
  bp=df$group=="exposed"&df$binge12==0
  att_table=att_row(bp,"Placebo: held-out last pre-exposure visit",NA,"placebo")
}else{
  ispost=df$binge12==1
  rows=list()
  yr=floor(df$T_scaled+0.5)  # nearest year since the last pre-exposure visit
  for(k in -3:0){b=!ispost&yr==k; rows[[length(rows)+1]]=att_row(b,sprintf("year %d",k),k,"in-sample")}
  for(k in 1:5){b=ispost&yr==k; rows[[length(rows)+1]]=att_row(b,sprintf("year %d",k),k,if(k<=3) "post" else "post-sparse")}
  for(k in 0:3){b=ispost&df$T_scaled>=k&df$T_scaled<k+1; rows[[length(rows)+1]]=att_row(b,sprintf("floor [%d,%d)",k,k+1),NA,"post-floor")}
  b=ispost&df$T_scaled<2; rows[[length(rows)+1]]=att_row(b,"floor [0,2)",NA,"post-floor")
  rows[[length(rows)+1]]=att_row(ispost,"all post-exposure",NA,"post-all")
  att_table=do.call(rbind,rows)
}
att_table$sc=sc; att_table$sigma2=mean(sigma_draws)
write.csv(att_table,file.path(APP_OUT,"att_table.csv"),row.names=FALSE)
print(att_table)
if(FALSE){
  # Shrinkage bias B_N from the proof, with the treated mean score estimated by GLS at each draw
  S=length(MU_array); ids_t=unique(df$id); bins=list("year 1"=1,"year 2"=2,"year 3"=3)
  yr=floor(df$T_scaled+0.5); ispost=df$binge12==1
  rowid=match(ids_t,id_map)
  nu_s=d_s=numeric(S); Bmat=matrix(NA,S,4,dimnames=list(NULL,c("year 1","year 2","year 3","all"))); 
  pc=phi_t_cont[[1]]
  for(s in 1:S){
    mu_=MU_array[[s]]; be=BETA_array[[s]]; th=as.vector(THETA_array[[s]]); d_=D_array[[s]][1,1]; s2=SIGMA_OMEGA_array[[s]][1]
    num=den=0; Rv=numeric(length(ids_t))
    for(u in seq_along(ids_t)){
      pos=rowid[u]; P=t(phi_t2[[pos]]); y=as.vector(Y_sparse2_b[[pos]]); X=X_sparse_b[[pos]]
      g=as.vector(P%*%mu_); for(j in 1:ncol(X)) g=g+X[,j]*as.vector(P%*%be[,j])
      F_=as.vector(P%*%th); FF=sum(F_^2); R=(s2/d_)/(FF+s2/d_); Rv[u]=R
      num=num+R*sum(F_*(y-g))/s2; den=den+FF*R/s2
    }
    nu=num/den; nu_s[s]=nu; d_s[s]=d_
    fp=as.vector(th%*%pc[,round(df$T,3)*1000+1,drop=FALSE])
    Ri=Rv[match(df$id,ids_t)]
    for(nm in colnames(Bmat)){ b=if(nm=="all") ispost else ispost&yr==as.numeric(sub("year ","",nm)); Bmat[s,nm]=mean(fp[b]*Ri[b])*nu }
  }
  sb=data.frame(bin=colnames(Bmat),B_mean=colMeans(Bmat),B_lo=apply(Bmat,2,quantile,.025),B_hi=apply(Bmat,2,quantile,.975))
  sb$nu_over_sd_mean=mean(nu_s/sqrt(d_s)); sb$nu_over_sd_lo=quantile(nu_s/sqrt(d_s),.025); sb$nu_over_sd_hi=quantile(nu_s/sqrt(d_s),.975)
  sb$mean_R=mean(Rv); sb$fsign=sign(mean(fp))
  adj=list(); for(nm in colnames(Bmat)){ b=if(nm=="all") ispost else ispost&yr==as.numeric(sub("year ","",nm)); dr=colMeans(as.numeric(df$superiorfrontal_grayvol[b])-fit_summary[b,,drop=FALSE])-Bmat[,nm]; dp=dr-rnorm(S,0,sqrt(sigma_draws/sum(b))); adj[[nm]]=c(adj_est=mean(dr),adj_lo=quantile(dp,.025,names=FALSE),adj_hi=quantile(dp,.975,names=FALSE)) }
  sb=cbind(sb,do.call(rbind,adj)); write.csv(sb,file.path(APP_OUT,"shrinkage.csv"),row.names=FALSE); print(sb)
}
if(!placebo){
pd=att_table[att_table$kind%in%c("in-sample","post"),]; pd$lo=ifelse(pd$kind=="post",pd$pred_lo,pd$fit_lo); pd$hi=ifelse(pd$kind=="post",pd$pred_hi,pd$fit_hi)
p1=ggplot(pd,aes(x=x,y=est,color=kind,fill=kind))+
  geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.2,color=NA)+geom_line(linewidth=1)+geom_point(size=2)+
  geom_hline(yintercept=0,color="red",linetype="dashed")+geom_vline(xintercept=0,color="red",linetype="dashed")+
  scale_color_manual(values=c("in-sample"="grey40","post"="blue"),labels=c("in-sample"="Pre-exposure residual (in sample)","post"="ATT (post-exposure)"))+
  scale_fill_manual(values=c("in-sample"="grey40","post"="blue"),labels=c("in-sample"="Pre-exposure residual (in sample)","post"="ATT (post-exposure)"))+
  scale_x_continuous(breaks=-3:3)+
  labs(title="Average treatment effect on the treated (ATT)",x="Years since last pre-exposure visit (nearest year)",y="Superior frontal gray volume (standardized)",color=NULL,fill=NULL)+
  theme_minimal()+theme(plot.title=element_text(hjust=.5,face="bold"),legend.position="top",panel.grid.major=element_line(color="grey90"),panel.grid.minor=element_blank())
}else{ saveRDS(att_table,file.path(APP_OUT,"placebo.rds")); quit(save="no") }

#----------fit plot--------------



att_array=c()
att_upper=c()
att_lower=c()
y_actual=c()
y_fitted=c()
y_median=c()
y_lower=c()
y_upper=c()
t=c()
#-4:2
for(i in 12:21){
  b=round(df$visit_age)==i &df$binge12==0
  
  #ATT_median= apply(fit_summary[b,]-as.data.frame(replicate(APP_S, df$superiorfrontal_grayvol[b])), 1, median)
  if(sum(b)>1){
    ATT_ests= colMeans(-fit_summary[b,]+as.data.frame(replicate(APP_S, df$superiorfrontal_grayvol[b])),na.rm=TRUE)
    att_array=c(att_array,mean(ATT_ests,na.rm=TRUE))
    att_lower=c(att_lower,quantile(ATT_ests,0.025,na.rm=TRUE))
    att_upper=c(att_upper,quantile(ATT_ests,0.975,na.rm=TRUE))
    y_actual=c(y_actual,mean(df$superiorfrontal_grayvol[b]))
    y_fitted=c(y_fitted,mean(colMeans(fit_summary[b,],na.rm=TRUE)))
    y_median=c(y_median,median(colMeans(fit_summary[b,],na.rm=TRUE),na.rm=TRUE))
    y_lower=c(y_lower,quantile(colMeans(fit_summary[b,],na.rm=TRUE),0.025,na.rm=TRUE))
    y_upper=c(y_upper,quantile(colMeans(fit_summary[b,],na.rm=TRUE),0.975,na.rm=TRUE))
    t=c(t,i)
  }
  
}


#t=-7:0
t;att_array;att_upper;att_lower
plot(t,att_array,ylim = c(-150,50),type="l")
abline(h=0,col="red")
abline(v=0,col="red")
lines(t,att_lower,col="grey")
lines(t,att_upper,col="grey")


# Create a data frame for ggplot
plot_data <- data.frame(
  Time = t,
  ATT = att_array,
  Upper = att_upper,
  Lower = att_lower,
  y_actual,
  y_fitted,
  y_upper,
  y_lower
  
)



p2=ggplot(plot_data, aes(x = Time)) +
  # Confidence interval band
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), fill = "red", alpha = 0.2) +
  
  # Observed line (blue, solid)
  geom_line(aes(y = y_actual, color = "observed", linetype = "observed"), linewidth = 1) +
  
  # Fitted line (red, dotted)
  geom_line(aes(y = y_fitted, color = "fitted", linetype = "fitted"), linewidth = 1) +
  
  # Manual scale mappings
  scale_color_manual(values = c("observed" = "blue", "fitted" = "red")) +
  scale_linetype_manual(values = c("observed" = "solid", "fitted" = "dashed")) +
  
  # Labels and theme
  labs(
    title = "Pretreatment Fit for Treated group by age",
    x = "Age at visit",
    y = "Superior frontal gray volume (standardized)",
    color = "Line Type",
    linetype = "Line Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = t)

df$fitted_y=rowMeans(fit_summary)
df$fitted_y_upper=apply(as.matrix(fit_summary), 1, quantile, probs = 0.975, na.rm = TRUE)
df$fitted_y_lower=apply(as.matrix(fit_summary), 1, quantile, probs = 0.275, na.rm = TRUE)
plot(df$visit_age,df$superiorfrontal_grayvol)
points(df$visit_age,df$fitted_y,col="red")

plot(df$visit_age,df$superiorfrontal_grayvol-df$fitted_y)

library(gridExtra)
 # side by side


png(file.path(APP_OUT,"ATT and pretreat.png"), width = 1400, height = 600, res = 150)
grid.arrange(p1, p2, ncol = 2) 
dev.off()



png(file.path(APP_OUT,"diagnostic plot.png"), width = 1000, height = 600, res = 150)
plot(loo_result$psis)
dev.off()

#------------------------FPC and mean function visualization----------
BANDS=list()
coef_band=function(draws,name){ fd=t(phi_t_cont[[1]])%*%t(draws); out=cbind(mean=rowMeans(fd),median=apply(fd,1,median),lower=apply(fd,1,quantile,.025),upper=apply(fd,1,quantile,.975),sd=apply(fd,1,sd))
  a=time_cont*(max(control$visit_age)-min(control$visit_age))+min(control$visit_age); BANDS[[name]]<<-data.frame(name=name,share_excl0=mean(out[,"lower"]>0|out[,"upper"]<0),share_pos=mean(out[,"lower"]>0),share_neg=mean(out[,"upper"]<0),mean_min=min(out[,"mean"]),mean_max=max(out[,"mean"]),mean12=out[which.min(abs(a-12)),"mean"],mean16=out[which.min(abs(a-16)),"mean"],mean20=out[which.min(abs(a-20)),"mean"]); out }

#------a plot visualizing FPC for each iteration---------
#with CI
phi_cont=phi_t_cont
#FPC1
theta_summary=df_theta_summary

theta_mat <- theta_summary[, 1:(nknots + 4)]
ref=df_theta[1,1:(nknots+4)]
th_n=t(apply(df_theta[,1:(nknots+4),drop=FALSE],1,function(v){v=v/sqrt(sum(v^2)); if(sum(v*ref)<0) -v else v}))
fd=t(phi_cont[[1]])%*%t(th_n)
if(mean(fd)<0) fd=-fd
FPC1=cbind(mean=rowMeans(fd),median=apply(fd,1,median),lower=apply(fd,1,quantile,.025),upper=apply(fd,1,quantile,.975),sd=apply(fd,1,sd))

#FPC2_t=t(FPC2)
colnames(FPC1)=c("mean","median","lower","upper","sd")
# Plot
age_cont=time_cont*(max(control$visit_age)-min(control$visit_age))+min(control$visit_age)
p_fpc=ggplot(as.data.frame(FPC1), aes(x = age_cont)) +
  geom_line(aes(y = mean, color = "Mean"), linewidth = 1.5) +
  geom_line(aes(y = median, color = "Median"), linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  labs(
    title = "Functional Principal Component",
    x = "Age",
    y = "Normalized first FPC",
    color = "Legend"
  ) +
  scale_color_manual(
    values = c("Mean" = "darkblue", "Median" = "red"),
    breaks = c("Mean", "Median")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) 


MU=coef_band(df_mu[,1:(nknots+4),drop=FALSE],'mu')

#FPC2_t=t(FPC2)
colnames(MU)=c("mean","median","lower","upper","sd")
# Plot
age_cont=time_cont*(max(control$visit_age)-min(control$visit_age))+min(control$visit_age)
p_mu=ggplot(as.data.frame(MU), aes(x = age_cont)) +
  geom_line(aes(y = mean, color = "Mean"), linewidth = 1.5) +
  geom_line(aes(y = median, color = "Median"), linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  labs(
    title = "Mean Function (without Covariates)",
    x = "Age",
    y = "Superior frontal gray volume (standardized)",
    color = "Legend"
  ) +
  scale_color_manual(
    values = c("Mean" = "darkblue", "Median" = "red"),
    breaks = c("Mean", "Median")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) 



if(Q>=2){
  # Several FPCs: the sampler leaves Theta unconstrained, so report eigenfunctions of the covariance function f(s)'Df(t).
  # Each draw's eigenfunctions are rotated toward the eigenfunctions of the posterior mean covariance function.
  gi=seq(1,1001,by=5); P_g=phi_cont[[1]][,gi,drop=FALSE]; dx=5/1000; ag=age_cont[gi]; S_=length(THETA_array)
  Cbar=matrix(0,length(gi),length(gi))
  for(s in 1:S_){ F_=t(P_g)%*%THETA_array[[s]]; Cbar=Cbar+F_%*%D_array[[s]]%*%t(F_) }
  Cbar=Cbar/S_
  e0=eigen(Cbar,symmetric=TRUE); v0=e0$vectors[,1:Q,drop=FALSE]/sqrt(dx)
  if(mean(v0[,1])<0) v0[,1]=-v0[,1]
  for(j in 2:Q) if(cor(v0[,j],ag)<0) v0[,j]=-v0[,j]
  share=e0$values[1:Q]/sum(e0$values)
  arr=array(NA,c(length(gi),Q,S_)); lam_s=matrix(NA,S_,Q)
  for(s in 1:S_){
    F_=t(P_g)%*%THETA_array[[s]]; q_=qr(F_); Q_=qr.Q(q_); R_=qr.R(q_)
    e_=eigen(R_%*%D_array[[s]]%*%t(R_),symmetric=TRUE); v_=Q_%*%e_$vectors/sqrt(dx)
    sv=svd(t(v_)%*%v0*dx); arr[,,s]=v_%*%sv$u%*%t(sv$v); lam_s[s,]=e_$values*dx
  }
  bnd=function(m) data.frame(age=ag,mean=rowMeans(m),median=apply(m,1,median),lower=apply(m,1,quantile,.025),upper=apply(m,1,quantile,.975))
  BB=lapply(1:Q,function(j) bnd(arr[,j,]))
  write.csv(do.call(rbind,lapply(1:Q,function(j) cbind(fpc=j,BB[[j]]))),file.path(APP_OUT,"fpc.csv"),row.names=FALSE)
  write.csv(data.frame(fpc=1:Q,share=share,lambda=e0$values[1:Q]*dx,lambda_mean=colMeans(lam_s),lambda_lo=apply(lam_s,2,quantile,.025),lambda_hi=apply(lam_s,2,quantile,.975)),file.path(APP_OUT,"fpc_eigen.csv"),row.names=FALSE)
  pf=function(B,ttl,yl) ggplot(B,aes(x=age))+
    geom_line(aes(y=mean,color="Mean"),linewidth=1.5)+geom_line(aes(y=median,color="Median"),linewidth=1,linetype="dashed")+
    geom_ribbon(aes(ymin=lower,ymax=upper),fill="blue",alpha=.2)+
    labs(title=ttl,x="Age",y=yl,color="Legend")+
    scale_color_manual(values=c("Mean"="darkblue","Median"="red"),breaks=c("Mean","Median"))+
    theme_minimal()+theme(plot.title=element_text(hjust=.5,face="bold"),legend.position="top",panel.grid.major=element_line(color="grey90"),panel.grid.minor=element_blank())
  ords=c("First","Second","Third","Fourth","Fifth")
  pl=lapply(1:Q,function(j) pf(BB[[j]],sprintf("%s FPC (%.1f%% of variance)",ords[j],100*share[j]),paste(ords[j],"FPC")))
  png(file.path(APP_OUT,"FPC and mean.png"), width = 1500, height = 1200, res = 150)
  do.call(grid.arrange,c(list(p_mu),pl,list(ncol=2)))
  dev.off()
}else{
png(file.path(APP_OUT,"FPC and mean.png"), width = 1200, height = 600, res = 150)
grid.arrange( p_mu,p_fpc, ncol = 2)
dev.off()
}


png(file.path(APP_OUT,"ATT.png"), width = 1000, height = 800, res = 150)
p1
dev.off()

#------a plot visualizing sex and icv baseline---------
#with CI
phi_cont=phi_t_cont
#FPC1
beta=coef_band(df_beta[,1:7,drop=FALSE],'icv')


colnames(beta)=c("mean","median","lower","upper","sd")
# Plot
age_cont=time_cont*(max(control$visit_age)-min(control$visit_age))+min(control$visit_age)
p_icv=ggplot(as.data.frame(beta), aes(x = age_cont)) +
  geom_line(aes(y = mean, color = "Mean"), linewidth = 1.5) +
  geom_line(aes(y = median, color = "Median"), linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  labs(
    title = "Baseline ICV",
    x = "Age",
    y = "Superior frontal gray volume (standardized)",
    color = "Legend"
  ) +
  scale_color_manual(
    values = c("Mean" = "darkblue", "Median" = "red"),
    breaks = c("Mean", "Median")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) 

#FPC1
beta=coef_band(df_beta[,8:14,drop=FALSE],'sexF')

#FPC2_t=t(FPC2)
colnames(beta)=c("mean","median","lower","upper","sd")
# Plot
age_cont=time_cont*(max(control$visit_age)-min(control$visit_age))+min(control$visit_age)
p_sexF=ggplot(as.data.frame(beta), aes(x = age_cont)) +
  geom_line(aes(y = mean, color = "Mean"), linewidth = 1.5) +
  geom_line(aes(y = median, color = "Median"), linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  labs(
    title = "Sex",
    x = "Age",
    y = "Superior frontal gray volume (standardized)",
    color = "Legend"
  ) +
  scale_color_manual(
    values = c("Mean" = "darkblue", "Median" = "red"),
    breaks = c("Mean", "Median")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) 


#png("baseline icv and age.png", width = 1200, height = 600, res = 150)
#grid.arrange( p_icv,p_age, ncol = 2) 
#dev.off()

#------a plot visualizing baseline sex and race---------
#with CI
phi_cont=phi_t_cont
#FPC1
beta=coef_band(df_beta[,15:21,drop=FALSE],'black')


colnames(beta)=c("mean","median","lower","upper","sd")
# Plot
age_cont=time_cont*(max(control$visit_age)-min(control$visit_age))+min(control$visit_age)
p_race_black=ggplot(as.data.frame(beta), aes(x = age_cont)) +
  geom_line(aes(y = mean, color = "Mean"), linewidth = 1.5) +
  geom_line(aes(y = median, color = "Median"), linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  labs(
    title = "Race Label Black ",
    x = "Age",
    y = "Superior frontal gray volume (standardized)",
    color = "Legend"
  ) +
  scale_color_manual(
    values = c("Mean" = "darkblue", "Median" = "red"),
    breaks = c("Mean", "Median")
  ) + scale_x_continuous(limits = c(min(age_cont, na.rm = TRUE), 20.5))+
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) 


phi_cont=phi_t_cont
#FPC1
beta=coef_band(df_beta[,22:28,drop=FALSE],'asian')


colnames(beta)=c("mean","median","lower","upper","sd")
# Plot
age_cont=time_cont*(max(control$visit_age)-min(control$visit_age))+min(control$visit_age)
p_race_asian=ggplot(as.data.frame(beta), aes(x = age_cont)) +
  geom_line(aes(y = mean, color = "Mean"), linewidth = 1.5) +
  geom_line(aes(y = median, color = "Median"), linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  labs(
    title = "Race Label Asian ",
    x = "Age",
    y = "Superior frontal gray volume (standardized)",
    color = "Legend"
  ) + scale_x_continuous(limits = c(min(age_cont, na.rm = TRUE), 20.5))+
  scale_color_manual(
    values = c("Mean" = "darkblue", "Median" = "red"),
    breaks = c("Mean", "Median")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) 
phi_cont=phi_t_cont
#FPC1
beta=coef_band(df_beta[,29:35,drop=FALSE],'white')


colnames(beta)=c("mean","median","lower","upper","sd")
# Plot
age_cont=time_cont*(max(control$visit_age)-min(control$visit_age))+min(control$visit_age)
p_race_white=ggplot(as.data.frame(beta), aes(x = age_cont)) +
  geom_line(aes(y = mean, color = "Mean"), linewidth = 1.5) +
  geom_line(aes(y = median, color = "Median"), linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  labs(
    title = "Race Label White",
    x = "Age",
    y = "Superior frontal gray volume (standardized)",
    color = "Legend"
  ) + scale_x_continuous(limits = c(min(age_cont, na.rm = TRUE), 20.5))+
  scale_color_manual(
    values = c("Mean" = "darkblue", "Median" = "red"),
    breaks = c("Mean", "Median")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) 


png(file.path(APP_OUT,"covariate plots.png"), width = 1000, height = 600, res = 150)
grid.arrange( p_icv,p_sexF, ncol = 2) 
dev.off()


png(file.path(APP_OUT,"race plot.png"), width = 1500, height = 600, res = 150)
grid.arrange(p_race_black,p_race_asian,p_race_white, ncol = 3) 
dev.off()


#------------individual fit-------------------------
df$T=df$T*(max(dat_all$visit_age)-min(dat_all$visit_age))+min(dat_all$visit_age)
fit_cts[,"T"]=fit_cts[,"T"]*(max(dat_all$visit_age)-min(dat_all$visit_age))+min(dat_all$visit_age)
df<- df %>%
  group_by(id) %>%
  mutate(
    first_treat_time = min(T[group == "exposed"], na.rm = TRUE)
#    T_scaled = T - first_treat_time
  ) %>%
  ungroup()


ids <- unique(exposed$id)
df$y=df$superiorfrontal_grayvol


set.seed(APP_SEED+2); pick <- sample(ids,16 )
# named vector of means for those ids

mean_by_id <- setNames(
  sapply(pick, function(k) mean(df$y[df$id == k], na.rm = TRUE)),
  pick
)
# ascending sort (low → high)
mean_by_id_sorted <- sort(mean_by_id)

ids_sorted_asc <- names(mean_by_id_sorted)

# ascending sort (low → high)
mean_by_id_sorted <- sort(mean_by_id)
ids_sorted_asc <- names(mean_by_id_sorted)

subd <- subset(df, id %in% ids_sorted_asc)
subd$id <- factor(subd$id, levels = ids_sorted_asc)
df_cts= data.frame(fit_cts[fit_cts[,1] %in% ids_sorted_asc,])
df_cts$fitted=1
subd$fitted=0
#colnames(df_cts)=c("id","y","fitted_y_lower","fitted_y_upper","T","fitted")
#rbind(df_cts,subd[,c("id","y","fitted_y_lower","fitted_y_upper","T","fitted")])

library(ggplot2)
library(dplyr)

## One row per id for the treatment time (in case it's repeated)
vline_df <- subd %>%
  distinct(id, first_treat_time)

p_ind=ggplot() +
  # detection-limit band: lower/upper vs T (per id)
  geom_ribbon(
    data = df_cts,
    aes(x = T, ymin = lower, ymax = upper, group = id),
    alpha = 0.2
  ) +
  # fitted curve est vs T (per id)
  geom_line(
    data = df_cts,
    aes(x = T, y = est, group = id)
  ) +
  # observed points y vs T (per id)
  geom_point(
    data = subd,
    aes(x = T, y = y)
  ) +
  # vertical red line at first_treat_time (per id)
  geom_vline(
    data = vline_df,
    aes(xintercept = first_treat_time),
    colour = "red"
  ) +
  # 3x3 grid, one panel per id
  facet_wrap(~ id, nrow = 4, ncol = 4) +
  labs(
    x = "Age (years)",
    y = "Superior frontal gray volume (standardized)",
    title = "Fitted counterfactuals (curves) and observed values (dots) by participant"
  ) +
  theme_bw() + theme(strip.background = element_blank(), strip.text = element_blank())


png(file.path(APP_OUT,"individual fit.png"), width = 3000, height = 3000, res = 300)
p_ind 
dev.off()


first_binge <- dat_all %>%
  group_by(id) %>%
  summarise(
    first_binge_age =
      if (any(binge12 == 1, na.rm = TRUE)) {
        min(visit_age[binge12 == 1], na.rm = TRUE)
      } else {
        NA_real_
      },
    .groups = "drop"
  )


png(
  filename = file.path(APP_OUT,"age_distribution_heavy_binge_drinking.png"),
  width = 6, height = 4, units = "in",
  res = 300
)

ggplot(first_binge, aes(x = first_binge_age)) +
  geom_histogram(binwidth = 0.5, na.rm = TRUE) +
  labs(x = "Age at heavy binge onset", y = "Number of subjects")
dev.off()

sm$tTable[,c(1,2,5)]
plot_data_att
plot_data_att[,2:4]= plot_data_att[,2:4]* sc 
plot_data_att
write.csv(plot_data_att,file.path(APP_OUT,"att_original_bins.csv"))


print(loo_result$psis)

saveRDS(list(lme=sm$tTable,sc=sc,center=center,sigma2=mean(unlist(SIGMA_OMEGA_array)),visits=list(all=visits_per_id,exposed=visits_per_id_exposed),demo=table_demo),file.path(APP_OUT,"summary.rds"))
write.csv(do.call(rbind,BANDS),file.path(APP_OUT,"coef_bands.csv"),row.names=FALSE)
write.csv(as.data.frame(t(table_demo)),file.path(APP_OUT,"demo.csv"))
write.csv(rbind(cbind(group="all",as.data.frame(visits_per_id %>% summarise(median=median(n_visit),min=min(n_visit),max=max(n_visit),mean=mean(n_visit)))),cbind(group="exposed",as.data.frame(visits_per_id_exposed %>% summarise(median=median(n_visit),min=min(n_visit),max=max(n_visit),mean=mean(n_visit)))),cbind(group="control",as.data.frame(visits_per_id[!visits_per_id$id%in%visits_per_id_exposed$id,] %>% summarise(median=median(n_visit),min=min(n_visit),max=max(n_visit),mean=mean(n_visit))))),file.path(APP_OUT,"visits.csv"),row.names=FALSE)
write.csv(as.data.frame(sm$tTable),file.path(APP_OUT,"lme.csv"))
write.csv(data.frame(sc=sc,center=center,sigma2=mean(unlist(SIGMA_OMEGA_array)),d=mean(unlist(D_array)),theta_norm=mean(sapply(THETA_array,function(x) sqrt(sum(x^2))))),file.path(APP_OUT,"scale.csv"),row.names=FALSE)
cat("done\n")
