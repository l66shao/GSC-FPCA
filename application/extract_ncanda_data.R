############################################################################################## 
############################################################################################## 
############################################################################################## 
##                   Extract NCANDA data for GSC-FPCA analyses                      	        ##
##              				  Wesley Thompson                      ##
##                                October 20, 2025    						                ##
############################################################################################## 
##############################################################################################
############################################################################################## 

rm(list=ls())

######################### Packages #####################################

library(readr)

######################### Functions #####################################

spaghetti <- function(y,times,ids,ylab="y",xlab="times", main=NA,
                      xlim=c(min(times,na.rm=T),max(times,na.rm=T)),
                      ylim=c(min(y,na.rm=T),max(y,na.rm=T)),
                      cols=rep(alpha("black", 0.7),length(unique(ids))),lwd=1){
  plot(times,y,type="n",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab, main=main)
  for(i in 1:length(unique(ids))){
    if(length(y[ids==unique(ids)[i]])==1) lines(times[ids==unique(ids)[i]],
                                                y[ids==unique(ids)[i]],type="p",col=cols[i],pch=15,cex=.25)
    if(length(y[ids==unique(ids)[i]])>1) lines(times[ids==unique(ids)[i]],
                                               y[ids==unique(ids)[i]],type="l",col=cols[i],lwd=lwd)
  }
}


LifeUseRunningTotal <- function(data){	
	data %>% group_by(id) %>% arrange(id, visit) %>% 
		mutate(
			LifeDrnk = coalesce(youthreport1_cddr19[visit == 1],0) + cumsum(coalesce(youthreport1_cddr20,0)),
			LifeMJ = coalesce(youthreport1_cddr55a[visit == 1],0) + cumsum(coalesce(youthreport1_cddr55b,0)),
			LifeBinge = cumsum(coalesce(youthreport1_cddr30,0)),
			LifeTob = coalesce(youthreport1_cddr4a[visit == 1],0) + cumsum(coalesce(youthreport1_cddr4b,0))
		)
}

######################### Dataset Assembly #####################################

dat <- read.csv(Sys.getenv("NCANDA_FULLDATA", "data/FullData_Y9.csv"))

## exposure variable-binge drinking
dat$binge_YN = dat$youthreport1_cddr27
dat$binge = dat$youthreport1_cddr30
dat$binge[dat$binge_YN==0 & !is.na(dat$binge_YN)] = 0
table(dat$binge,useNA="always")

length(unique(dat$subject))

## delete any participant with exceeds = "Y"
dat = dat[dat$exceeds_bl_drinking == "N",]
length(unique(dat$subject))

## Fill in 0 for any subject with dat$binge=NA at baseline
dat$binge[is.na(dat$binge & dat$visit == "baseline")] = 0

## delete any row with dat$binge=NA
dat = dat[!is.na(dat$binge),]
length(unique(dat$subject))


dat$site = as.factor(dat$site)

## recode visit
visit = rep(0, length(dat$visit))
visit[dat$visit == "baseline"] = 1
visit[dat$visit == "followup_1y"] = 2
visit[dat$visit == "followup_2y"] = 3
visit[dat$visit == "followup_3y"] = 4
visit[dat$visit == "followup_4y"] = 5
visit[dat$visit == "followup_5y"] = 6
visit[dat$visit == "followup_6y"] = 7
visit[dat$visit == "followup_7y"] = 8
visit[dat$visit == "followup_8y"] = 9
visit[dat$visit == "followup_9y"] = 10
dat$visit <- visit; rm(visit)
table(dat$visit)

## recode subject ID
id <- rep(NA, length(dat$subject))
N <- length(unique(dat$subject))
for(i in 1:N){
	id[dat$subject == unique(dat$subject)[i]]	<- i
}
dat$id <- id; 
  
## recode visit age
dat$base_age <- NA
for(i in 1:N){
	dat$base_age[dat$id == i] <- dat$visit_age[dat$id == i & dat$visit == 1]	
}
dat$ageD <- dat$visit_age - dat$base_age

## clump small race categories
dat$race_label[which(dat$race_label!="Caucasian/White" & dat$race_label!="Asian" 
                         & dat$race_label!="African-American/Black" 
                         & !is.na(dat$race_label))] = "Other"
dat$race_label = factor(dat$race_label)

dat$sex = factor(dat$sex)
dat$sexM = as.numeric(dat$sex == "M")
dat$sexF = as.numeric(dat$sex == "F")

## add scanner
dat$scanner = "siemans"
dat$scanner[dat$site=="A"] = "siemans"
dat$scanner[dat$site=="B"] = "ge"
dat$scanner[dat$site=="C"] = "ge"
dat$scanner[dat$site=="D"] = "siemans"
dat$scanner[dat$site=="E"] = "ge"
table(dat$scanner)
  
## response variables 
# gray matter volumes - averaging left (.x) and right (.y) volumes, paired by name. The version of this
# script used for the first posting paired each left volume with the column 34 positions over, which for
# superior frontal was the left frontal pole vertex count. The analyses recompute the outcome from
# FullData_Y9.csv in build_ext.R, so they do not depend on this column.
for(lx in grep("_grayvol\\.x$", colnames(dat), value = TRUE)){
  temp_name = sub("\\.x$", "", lx)
  dat[,temp_name] = (dat[,lx]+dat[,paste0(temp_name, ".y")])/2
}
colnames(dat)
out_gm_idx = c(which(colnames(dat) == "sri24_icv_brain_gm"),
               which( colnames(dat) == "bankssts_grayvol"):which( colnames(dat) == "insula_grayvol"))
out_gm = colnames(dat)[out_gm_idx]
  
# white matter volumes
dat$cerebellum_white_matter_volume_mm3 = NA
dat$cerebellum_white_matter_volume_mm3 = 0.5*(dat$left_cerebellum_white_matter_volume_mm3 +
                                                    dat$right_cerebellum_white_matter_volume_mm3)
out_wm_idx = c(which(names(dat)=="sri24_icv_brain_wm"), #total WM
               which(names(dat)=="sri24_parc116_pons_wm"),  #pons
               which(names(dat)=="sri24_parc116_corpus_callosum_wm"), #corpus callous
               which(names(dat)=="cc_posterior_volume_mm3"), #corpus callous segment
               which(names(dat)=="cc_mid_posterior_volume_mm3"), 
               which(names(dat)=="cc_central_volume_mm3"), 
               which(names(dat)=="cc_mid_anterior_volume_mm3"), 
               which(names(dat)=="cc_anterior_volume_mm3"), 
               which(names(dat)=="measure_corticalwhitemattervol"), # total cortical wm
               which(names(dat)=="cerebellum_white_matter_volume_mm3") #cerebellum
)  
out_wm = colnames(dat)[out_wm_idx]
  
## final processing of data
dim(dat)
length(unique(dat$id))
id = unique(dat$id); n = length(id)

dat$maxvisit = 0
dat$age_m = NA
dat$binge_ever = NA
dat$binge_to_date = NA
for(i in 1:n){
  	dat$maxvisit[dat$id == id[i]] = length(dat$maxvisit[dat$id == id[i]])
  	dat$age_m[dat$id == id[i]] = mean(dat$visit_age[dat$id == id[i]])
  	if(!length(dat$binge_YN[dat$id == id[i]]) == sum(is.na(dat$binge_YN[dat$id == id[i]]))){
		dat$binge_ever[dat$id == id[i]] = max(dat$binge_YN[dat$id == id[i]], na.rm=TRUE)
    		dat$binge_to_date[dat$id == id[i]] = 0
    		for(j in 1:length(dat$binge_to_date[dat$id == id[i]])){
    			if(!sum(is.na(dat$binge[dat$id == id[i]][1:j]))==j){
    				dat$binge_to_date[dat$id == id[i]][j] = 
    					max(dat$binge[dat$id == id[i]][1:j]>0, na.rm=TRUE)
    			}			
    		}
  	}	
}	

dat$binge12 = as.numeric(dat$binge>=12)
for(i in 1:n){
	dat_i = dat[dat$id == id[i],]
	n_i = nrow(dat_i)
	if(n_i > 1){
		for(j in 2:n_i){
			dat_i$binge12[j] = max(dat_i$binge12[1:j], na.rm=TRUE)
		}
	}
	dat[dat$id == id[i],] = dat_i	
}

dat$first_binge_visit = 11
for(i in id){
	dat_i = dat[dat$id == id[i],]
	if(max(dat_i$binge, na.rm=TRUE)>0) dat_i$first_binge_visit = min(dat_i$visit[dat_i$binge>0], na.rm=TRUE)
	dat[dat$id == id[i],] = dat_i
}	

dat$base_age_cat12_14 = as.numeric(dat$base_age < 15)
dat$base_age_cat15_17 = as.numeric(dat$base_age >= 15 & dat$base_age < 18)
dat$base_age_cat18_21 = as.numeric(dat$base_age >= 18)

# brain abnormal 
abnormal <- read.csv(Sys.getenv("NCANDA_ABNORMAL", "data/brain-abnormalities.csv"))
dat$abnormal <- 0
for(i in unique(dat$subject)){
	if(i %in% abnormal$mri_xnat_sid) dat$abnormal[dat$subject == i] = 1
}


nms = c("id","visit","superiorfrontal_grayvol","sri24_icv_brain_volume","LifeBinge",
	"base_age","visit_age","sex","sexF","sexM","race_label","site","base_age_cat12_14","base_age_cat15_17",
	"base_age_cat18_21","binge","binge12","first_binge_visit","abnormal")
ncanda_dat = dat[,names(dat) %in% nms]

### write out to csv
write_csv(ncanda_dat, Sys.getenv("NCANDA_DAT", "data/ncanda_dat.csv"))
