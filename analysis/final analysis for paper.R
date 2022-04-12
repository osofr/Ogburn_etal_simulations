
library(tmlenet)
library(igraph)
library(boot)

############################################################################################### 
#
# PREPROCESS THE DATA using original FHS data available through dbGap
# ( data can be requested here: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000007.v32.p13 )
# The data creation is based on work by Youjin Lee for the paper Lee & Ogburn (2020; JASA);
# her original code is available at https://github.com/youjin1207
#
############################################################################################### 

# READ IN NETWORK DATA, CREATE THE NETWORK MATRIX THAT TMLENET NEEDS. 
# remove the first colummn, which is "shareid"

	network_c1 = read.table("phs000153.v9.pht000836.v8.p8.c1.vr_sntwk_2008_m_0641s.DS-SN-IRB-MDS.txt", sep = "\t", header = TRUE)
	network_c2 = read.table("phs000153_c2.txt", sep = "\t", header = TRUE)
	network = rbind(network_c1, network_c2)
	head(network)

	wave1_network = network[(network$SPELLBEGIN <= 12*4) & (network$SPELLEND > 0*12),  ]
	#wave2_network = network[(network$SPELLBEGIN <= 12*12) & (network$SPELLEND > 8*12),  ]

# *** Fix this network creation to 
	shareid = unique(wave1_network$shareid)
	n = length(shareid)
	network = matrix(0, nrow=n,ncol=89)
	network[,1]<-unique(wave1_network$shareid)
	for (i in c(1:n)){
		friends<-unique(wave1_network[which(wave1_network$shareid==shareid[i]),]$sharealterid)
		network[i,2:(1+length(friends))]<-friends
		remove(friends)
	}
	colnames(network)<-c("shareid",rep("friends",88))


# CREATE THE DATASET NEEDED FOR TMLE
# each subject has covariates (age1, sex, education, smoke), exposure (bmi1), outcome (bmi2)

	#BMI for waves 1 and 2
	wave1_c1 = read.table("phs000007.v29.pht000030.v7.p10.c1.ex1_1s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
	wave2_c1 = read.table("phs000007.v29.pht000031.v7.p10.c1.ex1_2s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
	wave1_c2 = read.table("phs000007.v29.pht000030.v7.p10.c2.ex1_1s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
	wave2_c2 = read.table("phs000007.v29.pht000031.v7.p10.c2.ex1_2s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE) 

	BMI1_c1 = (wave1_c1$A50)*703 / (wave1_c1$A51)^2
	BMI2_c1 = (wave2_c1$B13) / (wave2_c1$B14 / 100)^2
	BMI1_c2 = (wave1_c2$A50)*703 / (wave1_c2$A51)^2
	BMI2_c2 = (wave2_c2$B13) / (wave2_c2$B14 / 100)^2

	BMI_wave1 =  as.data.frame(cbind(c(wave1_c1$shareid,wave1_c2$shareid), c(BMI1_c1,BMI1_c2)))
	BMI_wave2 =  as.data.frame(cbind(c(wave2_c1$shareid,wave2_c2$shareid), c(BMI2_c1,BMI2_c2)))
	colnames(BMI_wave1) = c("shareid", "BMI 1")
	colnames(BMI_wave2) = c("shareid", "BMI 2")

	## year of education
	edu_c1 = cbind(wave2_c1$shareid, wave2_c1$B43)
	edu_c2 = cbind(wave2_c2$shareid, wave2_c2$B43)
	edu.dat = as.data.frame(rbind(edu_c1, edu_c2))
	colnames(edu.dat) = c("shareid", "eduyear")
	

	## age and sex
	age.whole_c1 = read.table("phs000007.v29.pht003099.v4.p10.c1.vr_dates_2014_a_0912s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
	age.off_c1 = age.whole_c1[age.whole_c1$idtype == 1,] ## meaning "Offspring cohort"
	agesex_c1 = as.data.frame(cbind(age.off_c1$shareid, age.off_c1$sex, age.off_c1$age1))
	colnames(agesex_c1) = c("shareid", "sex", "age1")

	age.whole_c2 = read.table("phs000007.v29.pht003099.v4.p10.c2.vr_dates_2014_a_0912s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
	age.off_c2 = age.whole_c2[age.whole_c2$idtype == 1,]
	agesex_c2 = as.data.frame(cbind(age.off_c2$shareid, age.off_c2$sex, age.off_c2$age1))
	colnames(agesex_c2) = c("shareid", "sex", "age1")

	agesex.dat = rbind(agesex_c1, agesex_c2)


	# SUBSET FOR SUBJECTS WHO ARE IN BOTH WAVES 1 AND 2     
        
	shareid_network<-shareid
	shareid_wave1<-BMI_wave1$shareid
	shareid_wave2<-BMI_wave2$shareid
	final_shareid<-intersect(shareid_wave2,shareid_network)
	final_shareid<-final_shareid[order(final_shareid)]

	agesex.dat = subset(agesex.dat, shareid%in%final_shareid)
	edu.dat = subset(edu.dat, shareid%in%final_shareid)
	BMIwave1 = subset(BMI_wave1, shareid%in%final_shareid)
	BMIwave2 = subset(BMI_wave2, shareid%in%final_shareid)

	agesex.dat = agesex.dat[order(agesex.dat$shareid),]
	edu.dat = edu.dat[order(edu.dat$shareid),]
	BMIwave1 = BMIwave1[order(BMIwave1$shareid),]
	BMIwave2 = BMIwave2[order(BMIwave2$shareid),]
	
	#need outcome to be normalized to be between 0 and 1
	tmledata<-as.data.frame(cbind(final_shareid, agesex.dat$sex-1, agesex.dat$age1, 
							(agesex.dat$age1-min(agesex.dat$age1))/max(agesex.dat$age1)	,
							edu.dat$eduyear, BMIwave1$"BMI 1",
							(BMIwave1$"BMI 1"-min(BMIwave1$"BMI 1",na.rm=T))/max(BMIwave1$"BMI 1",na.rm=T), 
							BMIwave2$"BMI 2",
							(BMIwave2$"BMI 2"-min(BMIwave2$"BMI 2",na.rm=T))/max(BMIwave2$"BMI 2",na.rm=T)))
	colnames(tmledata) = c("shareid","sex","age","agenorm","edu","bmi1","bmi1norm","bmi2","bmi2norm")


# *** Fix network creation to solve the problem Oleg found, that friends may not be in the data 
# *** subset everything by final_shareid
	n = length(final_shareid)
	network = matrix(0, nrow=n,ncol=89)
	network[,1]<-final_shareid
	for (i in c(1:n)){
		friends<-unique(wave1_network[which(wave1_network$shareid==final_shareid[i]),]$sharealterid)
		friends_update<-friends[which(friends%in%final_shareid)]
		if(length(friends_update)>0){
		network[i,2:(1+length(friends_update))]<-friends_update}
		remove(friends)
		remove(friends_update)
	}
	colnames(network)<-c("ID",rep("friends",88))
	
write.table(network,"updatednetwork.txt")	

# *** fix ends here


#subset the network data for the subjects in the tmledata
finalnetwork = subset(network, network[,1]%in%final_shareid)
finalnetwork = finalnetwork[order(finalnetwork[,1]),]
finalnetwork[finalnetwork == 0] <- NA

#reducedfinalnetwork=(finalnetwork[,2:11])
reducedfinalnetwork=(finalnetwork[,2:4])

# save data for Oleg
colnames(finalnetwork)<-c("ID",rep("friends",88))
colnames(tmledata) = c("ID","sex","age","agenorm","edu","bmi1","bmi1norm","bmi2","bmi2norm")
write.table(finalnetwork, "finalnetwork.txt")
save(tmledata,file="tmledata.Rda")


############################################################################################### 
#
# ANALYZE THE DATA
#
############################################################################################### 


load("tmledata.Rda")
finalnetwork<-read.table("updatednetwork.txt")
reducedfinalnetwork=as.matrix((finalnetwork[,2:11]))

tmledata$obese1<-as.numeric(tmledata$bmi1>30)
tmledata$obese2<-as.numeric(tmledata$bmi2>30)

## OS: converting network matrix into string of IDs, sep by " ", removing NAs
friends_str = apply(reducedfinalnetwork, 1, function(x) paste(x[!is.na(x)], collapse = " "))
head(friends_str)
class(friends_str)
tmledata2 = cbind(tmledata, friends_str = friends_str, stringsAsFactors=FALSE)
head(tmledata2)
sapply(tmledata2, class)



############################################################################################### 
#
# first let's look at whether an increase in friends' bmi causes an increase in ego's bmi
#
############################################################################################### 

#options(tmlenet.verbose=TRUE)
sW=def_sW(sex="sex", age="agenorm", edu="edu", meanBMI1=sum(bmi1norm[[1:Kmax]])/nF)
sA<-def_sA(bmi1="bmi1norm", meanBMI1=sum(bmi1norm[[1:Kmax]])/nF, meanBMI2 = sum(bmi2norm[[1:Kmax]])/nF, sex="sex",age="agenorm", edu="edu", replaceNAw0 = TRUE) 
res <- eval.summaries(sW = sW, sA = sA,  Kmax = 10, data = tmledata2,
                      IDnode = "ID", NETIDnode = "friends_str", verbose = F)


tmle_output <- tmlenet(data = tmledata2, 
		sW=sW,
		sA=sA,
		Ynode = "bmi2norm",
		Kmax=10 ,
		# NETIDmat=reducedfinalnetwork,
		# IDnode="shareid",
		IDnode = "ID", NETIDnode = "friends_str",
		intervene1.sA = def_new_sA(meanBMI2 = meanBMI2 + .04), #adding half a std dev to the mean
		Qform = "bmi2norm ~ bmi1 + meanBMI1 + meanBMI2 +age+sex+edu", #this is essentially CF's outcome model
		hform.g0 = "edu ~ age", #innocuous h model to mimic CF's outcome regression only analysis
		optPars=list(n_MCsims=50,bootstrap.var = TRUE, n.bootstrap = 100),
		#optPars = list(bootstrap.var = bootstrap.var, n.bootstrap = n.bootstrap, boot.nodes = boot.nodes, boot.form = boot.form)
		#verbose = TRUE
		)


# to back transform from normalized BMI to BMI, multiply by max(bmi2)=52.65356 and add min(bmi2)=15.41289
# mean(bmi2) = 25.51257
# stdev(bmi2) = 4.418731

## point estimates		
tmle_output$EY_gstar1$estimates*52.65356+15.41289
## bootstrap CIs
tmle_output$EY_gstar1$boot.CIs*52.65356+15.41289
##IC.CIs
tmle_output$EY_gstar1$IC.CIs*52.65356+15.41289
##iid.CIs
tmle_output$EY_gstar1$iid.CIs*52.65356+15.41289
remove(tmle_output)


######################################################################################################### 
#
# now let's look at whether an increase in number of obese friends causes an increase in P(ego is obese)
#
#########################################################################################################

#options(tmlenet.verbose=TRUE)
sW=def_sW(sex="sex", age="agenorm", edu="edu",obese1="obese1", meanobese1=sum(obese1[[1:Kmax]])/nF)
sA<-def_sA(obese1="obese1", meanobese1=sum(obese1[[1:Kmax]])/nF, meanobese2 = sum(obese2[[1:Kmax]])/nF, sex="sex",age="agenorm", edu="edu", replaceNAw0 = TRUE) 
res <- eval.summaries(sW = sW, sA = sA,  Kmax = 10, data = tmledata2,
                      IDnode = "ID", NETIDnode = "friends_str", verbose = F)


tmle_output <- tmlenet(data = tmledata2, 
		sW=sW,
		sA=sA,
		Ynode = "obese2",
		Kmax=10 ,
		# NETIDmat=reducedfinalnetwork,
		# IDnode="shareid",
		IDnode = "ID", NETIDnode = "friends_str",
		#intervene1.sA = def_new_sA(meanobese2 = meanobese2 + .1), #adding half a std dev to the mean
		intervene1.sA = def_new_sA(meanobese2 = (sum(obese2[[1:Kmax]])+1)/nF), #adding one obese friend
		Qform = "obese2 ~ obese1 + meanobese1 +age+sex+edu", #this is essentially CF's outcome model
		hform.g0 = "edu~age", #innocuous h model to mimic CF's outcome regression only analysis
		optPars=list(n_MCsims=50,bootstrap.var = TRUE, n.bootstrap = 100),
		#optPars = list(bootstrap.var = bootstrap.var, n.bootstrap = n.bootstrap, boot.nodes = boot.nodes, boot.form = boot.form)
		#verbose = TRUE
		)
		

## point estimates		
tmle_output$EY_gstar1$estimates
## bootstrap CIs
tmle_output$EY_gstar1$boot.CIs
##IC.CIs
tmle_output$EY_gstar1$IC.CIs
##iid.CIs
tmle_output$EY_gstar1$iid.CIs
remove(tmle_output)


######################################################################################################### 
#
# in response to reviewers, we will run a version of the CF analysis that more closely tracks our own analysis
#    - use the same data: wave 2 outcome data with wave 1 for control
#    - use the pairwise CF logistic regression to predict E[Y(1)]
#    - bootstrap the confidence interval without accounting for dependence
#########################################################################################################

###############################################
# create pairwise data needed for CF analysis #
# (with multiple observations per ego)        #
###############################################

## c1 data
wave1_c1 = read.table("phs000007.v29.pht000030.v7.p10.c1.ex1_1s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
wave2_c1 = read.table("phs000007.v29.pht000031.v7.p10.c1.ex1_2s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
wave3_c1 = read.table("phs000007.v29.pht000032.v6.p10.c1.ex1_3s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
wave4_c1 = read.table("phs000007.v29.pht000033.v8.p10.c1.ex1_4s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
wave5_c1 = read.table("phs000007.v29.pht000034.v7.p10.c1.ex1_5s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
wave6_c1 = read.table("phs000007.v29.pht000035.v8.p10.c1.ex1_6s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
wave7_c1 = read.table("phs000007.v29.pht000036.v8.p10.c1.ex1_7s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)


BMI1_c1 = (wave1_c1$A50)*703 / (wave1_c1$A51)^2
BMI2_c1 = (wave2_c1$B13) / (wave2_c1$B14 / 100)^2
BMI3_c1 = (wave3_c1$C416)*703 / (wave3_c1$C417)^2
BMI4_c1 = (wave4_c1$D401)*703 / (wave4_c1$D402)^2
BMI5_c1 = (wave5_c1$E024)*703 / (wave5_c1$E025)^2
BMI6_c1 = (wave6_c1$F007)*703 / (wave6_c1$F008)^2
BMI7_c1 = (wave7_c1$G440)*703 / (wave7_c1$G441)^2

BMI_c1 = rbind( cbind(wave1_c1$shareid, BMI1_c1, rep(1, nrow(wave1_c1))),
                cbind(wave2_c1$shareid, BMI2_c1, rep(2, nrow(wave2_c1))),
                cbind(wave3_c1$shareid, BMI3_c1, rep(3, nrow(wave3_c1))),
                cbind(wave4_c1$shareid, BMI4_c1, rep(4, nrow(wave4_c1))),
                cbind(wave5_c1$shareid, BMI5_c1, rep(5, nrow(wave5_c1))),
                cbind(wave6_c1$shareid, BMI6_c1, rep(6, nrow(wave6_c1))),
                cbind(wave7_c1$shareid, BMI7_c1, rep(7, nrow(wave7_c1))))
obesity.ind = ifelse(BMI_c1[,2] > 30, 1, 0)                
obesity_c1 = cbind(BMI_c1, obesity.ind) 
colnames(obesity_c1) = c("shareid", "BMI", "wave", "obesity")
obesity_c1 = as.data.frame(obesity_c1)
## c2 data
wave1_c2 = read.table("phs000007.v29.pht000030.v7.p10.c2.ex1_1s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
wave2_c2 = read.table("phs000007.v29.pht000031.v7.p10.c2.ex1_2s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
wave3_c2 = read.table("phs000007.v29.pht000032.v6.p10.c2.ex1_3s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
wave4_c2 = read.table("phs000007.v29.pht000033.v8.p10.c2.ex1_4s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
wave5_c2 = read.table("phs000007.v29.pht000034.v7.p10.c2.ex1_5s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
wave6_c2 = read.table("phs000007.v29.pht000035.v8.p10.c2.ex1_6s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
wave7_c2 = read.table("phs000007.v29.pht000036.v8.p10.c2.ex1_7s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)

BMI1_c2 = (wave1_c2$A50)*703 / (wave1_c2$A51)^2
BMI2_c2 = (wave2_c2$B13) / (wave2_c2$B14 / 100)^2
BMI3_c2 = (wave3_c2$C416)*703 / (wave3_c2$C417)^2
BMI4_c2 = (wave4_c2$D401)*703 / (wave4_c2$D402)^2
BMI5_c2 = (wave5_c2$E024)*703 / (wave5_c2$E025)^2
BMI6_c2 = (wave6_c2$F007)*703 / (wave6_c2$F008)^2
BMI7_c2 = (wave7_c2$G440)*703 / (wave7_c2$G441)^2

BMI_c2 = rbind( cbind(wave1_c2$shareid, BMI1_c2, rep(1, nrow(wave1_c2))),
                cbind(wave2_c2$shareid, BMI2_c2, rep(2, nrow(wave2_c2))),
                cbind(wave3_c2$shareid, BMI3_c2, rep(3, nrow(wave3_c2))),
                cbind(wave4_c2$shareid, BMI4_c2, rep(4, nrow(wave4_c2))),
                cbind(wave5_c2$shareid, BMI5_c2, rep(5, nrow(wave5_c2))),
                cbind(wave6_c2$shareid, BMI6_c2, rep(6, nrow(wave6_c2))),
                cbind(wave7_c2$shareid, BMI7_c2, rep(7, nrow(wave7_c2))))
obesity.ind = ifelse(BMI_c2[,2] > 30, 1, 0)                
obesity_c2 = cbind(BMI_c2, obesity.ind) 
colnames(obesity_c2) = c("shareid", "BMI", "wave", "obesity")
obesity_c2 = as.data.frame(obesity_c2)


obesity.dat = rbind(obesity_c1, obesity_c2)
############# EDUCATION, AGE, SEX ###############
## year of education
edu_c1 = cbind(wave2_c1$shareid, wave2_c1$B43)
edu_c2 = cbind(wave2_c2$shareid, wave2_c2$B43)
edu.dat = rbind(edu_c1, edu_c2)
colnames(edu.dat) = c("shareid", "eduyear"); edu.dat = as.data.frame(edu.dat)
edu = rep(NA, nrow(tie.data))
for(i in 1:nrow(tie.data)){
  ind = which(edu.dat$shareid %in% tie.data$shareid[i])
  if(length(ind)!=0){
    edu[i] = edu.dat$eduyear[ind]
  }
}
## age and sex
#age.whole_c1 = read.table("data/age/phs000007.v29.pht003099.v4.p10.c1.vr_dates_2014_a_0912s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
age.whole_c1 = read.table("phs000007.v29.pht003099.v4.p10.c1.vr_dates_2014_a_0912s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
age.off_c1 = age.whole_c1[age.whole_c1$idtype == 1,]
agesex_c1 = as.data.frame(cbind(age.off_c1$shareid, age.off_c1$sex, age.off_c1$age1,
                                age.off_c1$age2, age.off_c1$age3, age.off_c1$age4, 
                                age.off_c1$age5, age.off_c1$age6, age.off_c1$age7))
colnames(agesex_c1) = c("shareid", "sex", "age1", "age2", "age3", "age4", "age5", "age6", "age7")


## age and sex
age.whole_c2 = read.table("phs000007.v29.pht003099.v4.p10.c2.vr_dates_2014_a_0912s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
age.off_c2 = age.whole_c2[age.whole_c2$idtype == 1,]
agesex_c2 = as.data.frame(cbind(age.off_c2$shareid, age.off_c2$sex, age.off_c2$age1,
                                age.off_c2$age2, age.off_c2$age3, age.off_c2$age4, 
                                age.off_c2$age5, age.off_c2$age6, age.off_c2$age7))
colnames(agesex_c2) = c("shareid", "sex", "age1", "age2", "age3", "age4", "age5", "age6", "age7")

## match education,age and sex
## age data at exam 1 are available for everyone so if ages are missing at other exams,
## we 'extrapolate' ages based on what we have at exam 1 for each subjects. 
agesex.dat = rbind(agesex_c1, agesex_c2)

for(i in 1:nrow(agesex.dat)){
  if(is.na(agesex.dat$age2[i])){
    agesex.dat$age2[i] = agesex.dat$age1[i] + 8
  }
  if(is.na(agesex.dat$age3[i])){
    agesex.dat$age3[i] = agesex.dat$age2[i] + 4
  }
  if(is.na(agesex.dat$age4[i])){
    agesex.dat$age4[i] = agesex.dat$age3[i] + 4
  }
  if(is.na(agesex.dat$age5[i])){
    agesex.dat$age5[i] = agesex.dat$age4[i] + 4
  }
  if(is.na(agesex.dat$age6[i])){
    agesex.dat$age6[i] = agesex.dat$age5[i] + 4
  }
  if(is.na(agesex.dat$age7[i])){
    agesex.dat$age7[i] = agesex.dat$age6[i] + 3
  }
}

age.current = rep(NA, nrow(tie.data)); sex = rep(NA, nrow(tie.data));
for(i in 1:nrow(tie.data)){
  ind = which(agesex.dat$shareid %in% tie.data$shareid[i])
  if(length(ind)!=0){
    age.current[i] = agesex.dat[ind,(tie.data$wave[i]+2)]
    sex[i] = agesex.dat$sex[ind]
  }
}

## make Y_{i,t}, Y_{j,t}, Y_{i,t-1}, and Y_{j, t-1} at each wave t.
ego.obesity.current = rep(NA, nrow(tie.data)); ego.obesity.previous = rep(NA, nrow(tie.data));
alter.obesity.previous = rep(NA, nrow(tie.data)); alter.obesity.current=rep(NA,nrow(tie.data))

for(i in 1:nrow(tie.data)){
  
  ego.current = which(obesity.dat$shareid %in% tie.data$shareid[i] & obesity.dat$wave == tie.data$wave[i])
  if(length(ego.current) == 1) ego.obesity.current[i] = obesity.dat$obesity[ego.current]
  
  ego.previous = which(obesity.dat$shareid %in% tie.data$shareid[i] & obesity.dat$wave == (tie.data$wave[i]-1))
  if(length(ego.previous) == 1) ego.obesity.previous[i] = obesity.dat$obesity[ego.previous]
  
  alter.current = which(obesity.dat$shareid %in% tie.data$sharealterid[i] & obesity.dat$wave == tie.data$wave[i])
  if(length(alter.current) == 1) alter.obesity.current[i] = obesity.dat$obesity[alter.current]
  
  alter.previous = which(obesity.dat$shareid %in% tie.data$sharealterid[i] & obesity.dat$wave == (tie.data$wave[i]-1))
  if(length(alter.previous) == 1) alter.obesity.previous[i] = obesity.dat$obesity[alter.previous]
    
}


## merge every data for each pair of siblings.
all.data = cbind(tie.data, sex, age.current,edu,
                 ego.obesity.current, ego.obesity.previous,
                 alter.obesity.previous)
                 
wave2.data = all.data[all.data$wave == 2,]

wave2.data<-read.csv("wave2.data1.csv")
wave2.data.na<-wave2.data1[complete.cases(wave2.data),]

###############################################
# function to obtain E[Y(1)] from the data  ###
###############################################

predictY1 <- function(data, indices) {
  d <- data[indices,] 
  fit <- glm(ego.obesity.current ~ age.current + sex + edu + ego.obesity.previous + alter.obesity.previous + alter.obesity.current, 
			data=d, na.action=na.exclude,family="binomial")
  newdata.txt<-data.frame(age.current=d$age.current, sex=d$sex, edu=d$edu, 		
					ego.obesity.previous=d$ego.obesity.previous,  
					alter.obesity.previous=d$alter.obesity.previous, alter.obesity.current=rep(1,nrow(d)))				
  
  predict1<-mean(predict(fit, newdata.txt, type="response"), na.rm=T)
  return(predict1)
}


#################
# bootstrap   ###
#################

set.seed(1)
resultsY1 <- boot(data=wave2.data.na, statistic=predictY1, R=500)
resultsY1$t0
print(c(sort(resultsY1$t)[13], sort(resultsY1$t)[487]))

