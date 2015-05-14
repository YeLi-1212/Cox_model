coreData <- read.csv('CoreTable_training.csv',header = T)

library(rms)
library(survival)

#store the survival time and status
origin_survival_time <- coreData$LKADT_P
origin_death <- coreData$DEATH
#transfrom the death factor into binary vector
origin_death <- ifelse(origin_death == 'YES',1,0)

#ALP ~ CREACLCA
commonfeatures <- coreData[30:54]
commonfeaturesnames <- colnames(coreData)[30:54]

mulfeature <- matrix(nrow = 1600)
for(i in 1:length(commonfeatures)){
  mulfeature <- cbind(mulfeature,as.numeric(as.matrix(commonfeatures[i],ncol = 1,nrow = 1600)))  
}

mulfeature <- mulfeature[,-1]

colnames(mulfeature,do.NULL=F)
colnames(mulfeature) <- commonfeaturesnames


#delete those variables have less than 1000 records
index <- c()

for(i in 1:ncol(mulfeature)){
  if(sum(!is.na(mulfeature[,i])) < 1000)
    index <- c(index,i)
}

mulfeature <- mulfeature[,-index]


#add age variable

age <- as.numeric(as.matrix(coreData[14],ncol = 1,nrow = 1600))

#first delete age_group variable
#add age_group variable
#age_group <- ifelse(coreData[15] == '18-64',1,ifelse(coreData[15] == '65-74',2,3))

#add bmi variable
bmi <- as.numeric(as.matrix(coreData[17],ncol = 1,nrow = 1600))

#add weight
weight <- as.numeric(as.matrix(coreData[20],ncol = 1,nrow = 1600))

mulfeature <- cbind(mulfeature,age,bmi,weight)

#process those binary features
binary_variable <- coreData[55:length(coreData)]

index <- c()
for(i in 1:length(binary_variable)){
  if(all(is.na(binary_variable[i])))
    index <- c(index,i)
}

binary_variable <- binary_variable[-index]

for(i in 1:length(binary_variable))
  binary_variable[i] <- ifelse(binary_variable[i] == 'Y',1,ifelse(binary_variable[i] == 'YES',1,0))

#all the variables are stores in the maritx called mulfeature
mulfeature <- as.matrix(cbind(mulfeature,binary_variable))

#analysis this matrix

#1
#first select all those rows with all of the variables accessible
#998 record


#since i delete ARTTHROM variable i have to delete this one from 
#the selection as well

#colnames(first_selection) == 'ARTTHROM'
mulfeature <- mulfeature[,-57]

index <- c()
for(i in 1:1600){
  if(any(is.na(mulfeature[i,]))) index <- c(index,i)
}

first_selection <- data.frame(mulfeature[-index,])



response <- Surv(origin_survival_time[-index],origin_death[-index],type='right')

#store survival_time and status to use for timeROC functions
survival_time <- origin_survival_time[-index]
death <- origin_death[-index]

#implement 10-fold check upon the first selection rows

folds <- split(sample(nrow(first_selection),size = nrow(first_selection)),1:10)

for(i in 1:10)
{
  print(i)
  
  index_chose <- folds[[i]]
  training_set <- first_selection[-index_chose,]
  test_set <- first_selection[index_chose,]
  
  
  #try to solve 
  #Error in Design(eval.parent(m)) : 
  #dataset ddist not found for options(datadist=)
  dd <- datadist(training_set)
  options(datadist = 'dd')
  
  cox_model <- cph(response[-index_chose] ~ ALP+ALT+AST+CA+CREAT+HB+NEU+PLT+PSA+TBILI+WBC+
                     NA.+MG+PHOS+ALB+TPRO+GLU+age+
                     #AGEGRP2
                     bmi+weight+
                     NON_TARGET+TARGET+BONE+RECTAL+LYMPH_NODES+KIDNEYS+
                     LUNGS+LIVER+PLEURA+OTHER+PROSTATE+ADRENAL+BLADDER+
                     PERITONEUM+COLON+SOFT_TISSUE+ABDOMINAL+ORCHIDECTOMY+
                     PROSTATECTOMY+TURP+LYMPHADENECTOMY+SPINAL_CORD_SURGERY+
                     BILATERAL_ORCHIDECTOMY+PRIOR_RADIOTHERAPY+ANALGESICS+
                     ANTI_ANDROGENS+GLUCOCORTICOID+GONADOTROPIN+BISPHOSPHONATE+
                     CORTICOSTEROID+IMIDAZOLE+ACE_INHIBITORS+BETA_BLOCKING+
                     HMG_COA_REDUCT+ESTROGENS+ANTI_ESTROGENS+
                     #ARTTHROM+
                     CEREBACC+
                     CHF+DVT+DIAB+GASTREFL+GIBLEED+MI+PUD+PULMEMB+PATHFRAC+
                     SPINCOMP+COPD+MHBLOOD+MHCARD+MHCONGEN+MHEAR+MHENDO+MHEYE+
                     MHGASTRO+MHGEN+MHHEPATO+MHIMMUNE+MHINFECT+MHINJURY+MHINVEST
                   +MHMETAB+MHMUSCLE+MHNEOPLA+MHNERV+MHPSYCH+MHRENAL+MHRESP+
                     MHSKIN+MHSOCIAL+MHSURG+MHVASC
                   
                   ,data = training_set,x = T ,y = T, surv = T)
  #PAN PAN PAN
  #note 8th iteration the ARTTHROM variables causes a failure
  #since the contingency table of this guy is this 
  #ARTTHROM 0   1
  #death0 555   0
       #1 442   1
  #i decided to ture off this guy
  
  #get survival ratio at specify time points
  
  g <- Survival(cox_model)
  time_points <- c(360,540,720)
  pr <- g(times=time_points)
  
  #there are two different ways to calculate the HR for test set
  #1 using predict function and use exp
  x_beta <- predict(cox_model,test_set)
  
  #2 using coefficitents and variables and exp
  coef <- cox_model$coef
  x_beta_not_centered <- as.matrix(test_set) %*% coef
  
  
  #evalute the three AUC
  #360
  #print('using x_beta')
  print('using x_beta')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=exp(x_beta)*pr[1],cause=1,times=360,iid=T)
  print(risk)
  print('using x_beta_not_centered')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=(exp(x_beta_not_centered)*pr[1]),cause=1,times=360,iid=T)
  print(risk)
  #540
  print('using x_beta')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=exp(x_beta)*pr[2],cause=1,times=540,iid=T)
  print(risk)
  print('using x_beta_not_centered')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=(exp(x_beta_not_centered)*pr[2]),cause=1,times=540,iid=T)
  print(risk)
  #720
  print('using x_beta')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=exp(x_beta)*pr[3],cause=1,times=720,iid=T)
  print(risk)
  print('using x_beta_not_centered')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=(exp(x_beta_not_centered)*pr[3]),cause=1,times=720,iid=T)
  print(risk)
  
}


#2
#then delete the last 6 variables which cause 433 missing rows
#NA. GLU ALB PHOS TPRO  MG
#    uni-HR  mul-HR
#NA.  0.952  0.984
#GLU     1   0.960
#ALB  0.95   0.977
#PHOS 0.756  0.829
#TPRO 0.993  1.01
#MG   0.959  0.910

delete_index <- c()
delete_variables <- c('NA.','GLU','ALB','PHOS','TPRO','MG')

for(i in 1:length(delete_variables)){
  delete_index <- c(delete_index,which(colnames(mulfeature) == delete_variables[i]))
}
second_selection <- mulfeature[,-delete_index]

#since following ERROR occurs:
#Warning message:
#  In fitter(X, Y, strata = Strata, offset = offset, weights = weights,  :
#             Loglik converged before variable  35 ; beta may be infinite. 
# i try to delete to the 35th variable to find out what will happen
#"LYMPHADENECTOMY"
#                          coef exp(coef) se(coef)      z    p
#second_selection[[35]] -0.0967     0.908    0.126 -0.768 0.44
#it seems this variable has little contribution to the risk score

second_selection <- second_selection[,-35]

index <- c()

for(i in 1:1600){
  if(any(is.na(second_selection[i,])))
     index <- c(index,i)
}

second_selection <- data.frame(second_selection[-index,])

#recreate response 

response <- Surv(origin_survival_time[-index],origin_death[-index],type='right')
#save survival time and status for timeROC functions

survival_time <- origin_survival_time[-index]
death <- origin_death[-index]

#implement 10-fold check upon the first selection rows

folds <- split(sample(nrow(second_selection),size = nrow(second_selection)),1:10)

for(i in 1:10)
{
  print(i)
  
  index_chose <- folds[[i]]
  training_set <- second_selection[-index_chose,]
  test_set <- second_selection[index_chose,]
  
  
  #try to solve 
  #Error in Design(eval.parent(m)) : 
  #dataset ddist not found for options(datadist=)
  dd <- datadist(training_set)
  options(datadist = 'dd')
  
  cox_model <- cph(response[-index_chose] ~ ALP+ALT+AST+CA+CREAT+HB+NEU+PLT+PSA+TBILI+WBC+
                       age
                     #+AGEGRP2
                     +bmi+weight+
                       NON_TARGET+TARGET+BONE+RECTAL+LYMPH_NODES+KIDNEYS+
                       LUNGS+LIVER+PLEURA+OTHER+PROSTATE+ADRENAL+BLADDER+
                       PERITONEUM+COLON+SOFT_TISSUE+ABDOMINAL+ORCHIDECTOMY+
                       PROSTATECTOMY+TURP
                       #+LYMPHADENECTOMY
                       +SPINAL_CORD_SURGERY+
                       BILATERAL_ORCHIDECTOMY+PRIOR_RADIOTHERAPY+ANALGESICS+
                       ANTI_ANDROGENS+GLUCOCORTICOID+GONADOTROPIN+BISPHOSPHONATE+
                       CORTICOSTEROID+IMIDAZOLE+ACE_INHIBITORS+BETA_BLOCKING+
                       HMG_COA_REDUCT+ESTROGENS+ANTI_ESTROGENS+
                       #ARTTHROM+
                       CEREBACC+
                       CHF+DVT+DIAB+GASTREFL+GIBLEED+MI+PUD+PULMEMB+PATHFRAC+
                       SPINCOMP+COPD+MHBLOOD+MHCARD+MHCONGEN+MHEAR+MHENDO+MHEYE+
                       MHGASTRO+MHGEN+MHHEPATO+MHIMMUNE+MHINFECT+MHINJURY+MHINVEST
                     +MHMETAB+MHMUSCLE+MHNEOPLA+MHNERV+MHPSYCH+MHRENAL+MHRESP+
                       MHSKIN+MHSOCIAL+MHSURG+MHVASC
                     
                     ,data = training_set,x = T ,y = T, surv = T)
  
  
  
  #get survival ratio at specify time points
  
  g <- Survival(cox_model)
  time_points <- c(360,540,720)
  pr <- g(times=time_points)
  
  #there are two different ways to calculate the HR for test set
  #1 using predict function and use exp
  x_beta <- predict(cox_model,test_set)
  
  #2 using coefficitents and variables and exp
  coef <- cox_model$coef
  x_beta_not_centered <- as.matrix(test_set) %*% coef
  
  
  #evalute the three AUC
  #360
  #print('using x_beta')
  print('using x_beta')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=exp(x_beta)*pr[1],cause=1,times=360,iid=T)
  print(risk)
  print('using x_beta_not_centered')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=(exp(x_beta_not_centered)*pr[1]),cause=1,times=360,iid=T)
  print(risk)
  #540
  print('using x_beta')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=exp(x_beta)*pr[2],cause=1,times=540,iid=T)
  print(risk)
  print('using x_beta_not_centered')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=(exp(x_beta_not_centered)*pr[2]),cause=1,times=540,iid=T)
  print(risk)
  #720
  print('using x_beta')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=exp(x_beta)*pr[3],cause=1,times=720,iid=T)
  print(risk)
  print('using x_beta_not_centered')
  risk <- timeROC(T=survival_time[index_chose],delta=death[index_chose],marker=(exp(x_beta_not_centered)*pr[3]),cause=1,times=720,iid=T)
  print(risk)
  
}












