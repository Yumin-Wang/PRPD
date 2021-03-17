###########################################################################################################################

# Simulation code for Partially Randomized Preference Design with Propensity Score Stratification

###########################################################################################################################

###########################################################################################################################
#Function: search

#Author: Yumin Wang, Fan Li

#Creation Date: July 26, 2020 (version 3.6.1)

#Purpose: This function is used to search alpha0 in the logistic model (true propensity score model)
#         which corresponds to a marginal choice arm and random arm allocation ratio: theta.
#         The function takes in alpha0 and returns marginal theta. The user is expected to vary 
#         their input alpha0 to try to get their desired marginal theta (output). The value of alpha0 
#         is used in the following function: gendata.
#
#         Further explanation: The true propensity score in the simulation population is assumed to follow a 
#         logistic model:e(X)=(1+exp[-(alpha0+alpha1*X1+alpha2*X2)])^(-1), where X1~N(0,1), X2~Bern(0.5). 
#         To control for the marginal choice arm and random arm alllocation ratio theta and find the 
#         corresponding alpha0, first setting (alpha1,alpha2)=(0.35,-0.3). Then try to vary input 
#         alpha0 to get one's desired marginal theta (output) for the simulation population.

# Required Parameters: 
#         a0: alpha0 in the logistic model

#Output:  The function takes in alpha0 and returns marginal theta for the simulation population.

#Example: The following gives the marginal theta=0.5 which corresponds to the input alpha0=0.15
#         search(0.15)
#         We have found the corresponding alpha0 when varying theta from 0.1-0.9 with a step of 0.05
#         alpha0=(-2.1,-1.63,-1.28,-0.98,-0.73,-0.49,-0.26,-0.06,0.15,0.36,0.57,0.79,1.03,1.29,1.58,1.94,2.4)
#         which corresponds to marginal theta=(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9)

###########################################################################################################################
search=function(a0){
  n=50000                    # set number of people to be generated to 50000
  Z1<-rnorm(n,mean=0,sd=1)   # generate X1
  Z2<-rbinom(n,1,0.5)        # generate X2
  W <- cbind(1, Z1, Z2)      # Design matrix 
  alpha1<-0.35               # set alpha1
  alpha2<- -0.3              # set alpha2
  a <- c(a0, alpha1, alpha2) # set vector (alpha0,alpha1,alpha2)
  mean(plogis(c(W%*%a)))     # calculate marginal theta
}

###########################################################################################################################
#Function: gendata

#Author: Yumin Wang, Fan Li

#Creation Date: September 25, 2020 (version 3.6.1)

#Purpose: This function generates one simulation dataset and one vector which verifies the true treatment, selection ,preferece 
#         effects of the dataset. 
#       The dataset has the following variables:
#       (1)Preference: Binary variable; Individuals choice and random arm assigned. 1 denotes choice arm, 0 denotes random arm
#       (2)Var1: A continuous variable generated from N(0,1)
#       (3)Var2: A binary variable generated from Bern(0.5)
#       (4)PS: True propensity score for each individual, 0<PS<1
#       (5)Arm: Factor version for (1), Random Arm, Choice Arm
#       (6)PS.hat: Estimated propensity score that used both X1 and X2 in estimating propensity score
#       (7)PS.hat.X1: Estimated propensity score that used only X1 in estimating propensity score
#       (8)PS.hat.X2: Estimated propensity score that used only X2 in estimating propensity score
#       (9)Trt: Binary variable; Individuals treatment assigned. 1 denotes Treatment A, 0 denotes Treatment B
#       (10)Treatment: Factor version for Trt: Treatment A, Treatment B
#       (11)pretruth: Binary variable. Individual's true preference for treatment A or treatment B. 1 denotes preferring Treatment A, 0 denotes preferring Treatment B
#       (12)Four counterfactuals: mu11, mu12, mu21, mu22 (Y11, Y12, Y21, Y22). The first subscript letter means the treatment patient received, 
#           the second subscript letter means the treatment patient preferred
#       (13)Y: Individual's actual response
#       (14)str.hat.X1: Stratum number the individual belongs to when using only X1 in estimating propensity score
#       (15)str.hat: Stratum number the individual belongs to when using both X1 and X2 in estimating propensity score
#       (16)str: Stratum number the individual belongs to when using true propensity score stratification
#       (17)str.hat.X2: Stratum number the individual belongs to when using only X2 in estimating propensity score
#       (18)count: Auxilliary varible used to count individuals, count=1.
#
#       The vector verifies the true treatment, selection ,preferece effects of the dataset: truth=(dtrt,dsel,dpre)

# Required Parameters: 
#          phi: This is the preference rate for treatment A, we assume it is constant in the entire study population. Range from 0 to 1.
#        theta: Marginal proportion of individual having a preference (marginal choice arm and random arm allocation ratio). Range from 0.1 to 0.9, step=0.05.
#       strata: Number of propensity score strata for PRPD. Range starting from 1.
#            n: Number of patients in the generated dataset
#          trt: True treatment effect of the simulation population
#          sel: True selection effect of the simulation population
#          pre: True preference effect of the simulation population
#        model: The model used to generate counterfacutals. Homo: Homogeneous effects model. Hetero: Heterogeneous effects model.

#Output:   
#       A list which contains one dataset and one vector.

#Example: The following code gives the dataset and the vector when phi=0.3,theta=0.3,strata=5,n=1000,trearment effect=0.2,selection effect=0.3,
#         preference effect=0.4 under homogeneous effects model.
#         dataset<-gendata(phi=0.3,theta=0.3,strata=5,n=1000,trt=0.2,sel=0.3,pre=0.4,model="Homo")
#         dataset$table1, dataset$truth

#############################################################################################################################################
gendata<-function(phi,theta,strata,n,trt,sel,pre,model){
  Z1<-rnorm(n,mean=0,sd=1)  # generate var1 (X1)
  p_Z2=0.5 #set the probability of binominal distribution
  Z2<-rbinom(n,1,p_Z2)  # generate var2 (X2)
  
  #Find the corresponding alpha0 in the logistic model (propensity score model) based on theta
  ALP<-c(-2.1,-1.63,-1.28,-0.98,-0.73,-0.49,-0.26,-0.06,0.15,0.36,0.57,0.79,1.03,1.29,1.58,1.94,2.4)
  a0=ALP[ceiling(theta/0.05-1)]
  #User could delete the above two lines and use the line below if they want to explore theta other than theta=0.1-0.9 (step=0.05)
  #a0=?
  #
  
  alpha1<-0.35 #set alpha1
  alpha2<- -0.3 # set alpha2
  W <- cbind(1, Z1, Z2) # Design matrix 
  a <- c(a0, alpha1, alpha2)
  e<-plogis(c(W%*%a)) # generate true PS
  #generate preference indicator from ps
  z <- rbinom(n, 1, e)
  
  #generate table1
  table1<-data.frame(Preference=z,Var1=Z1,Var2=Z2,PS=e)
  #define another column Arm which is the factor version of column Preference
  table1$Arm<-factor(table1$Preference,levels = c(0,1),labels = c("Random Arm","Choice Arm"))
  
  ##################################################################
  #calculate quantiles based on true PS and number of strata needed
  quan<-quantile(table1$PS,probs=seq(0,1,by=1/strata)) 
  quan[1]=0
  quan[strata+1]=1
  
  #estimate PS using both X1 and X2, the correct specification of the propensity score model
  emodel=glm(Preference~Var1+Var2,family = binomial,data=table1)
  #define another column PS.hat which is estimated propensity score which used both X1 and X2 in estimating propensity score
  table1$PS.hat=emodel$fitted.values
  #calculate quantiles based on estimated PS (both X1 and X2) and number of strata needed
  quan.hat<-quantile(table1$PS.hat,probs=seq(0,1,by=1/strata))
  quan.hat[1]=0
  quan.hat[strata+1]=1
  
  #estimate PS only X1, the incorrect specification of the propensity score model
  emodel_X1=glm(Preference~Var1,family = binomial,data=table1)
  #define another column PS.hat.X1 which is estimated propensity score which used only X1 in estimating propensity score
  table1$PS.hat.X1=emodel_X1$fitted.values
  #calculate quantiles based on estimated PS (only X1) and number of strata needed
  quan.hat.X1<-quantile(table1$PS.hat.X1,probs=seq(0,1,by=1/strata))
  quan.hat.X1[1]=0
  quan.hat.X1[strata+1]=1
  
  #estimate PS only X2, the incorrect specification of the propensity score model
  emodel_X2=glm(Preference~Var2,family = binomial,data=table1)
  #define another column PS.hat.X2 which is estimated propensity score which used only X2 in estimating propensity score
  table1$PS.hat.X2=emodel_X2$fitted.values
  ####################################################################
  
  table1$Trt<-9 #initiate column Trt
  # for people having a preference(choice arm),their Trt is generated from Bern(phi).
  table1$Trt[table1$Preference==1]<-rbinom(sum(table1$Preference==1),1,phi)
  # for people in the random arm, flip a coin (50%)
  table1$Trt[table1$Preference==0]<-rbinom(sum(table1$Preference==0),1,1/2)
  # create factor version for Trt
  table1$Treatment<-factor(table1$Trt,levels=c(0,1),labels = c("TreatmentB","TreatmentA")) 
  
  # create their true preference for A or B
  # for people having a preference(choice arm),their treatment assignment and their true preference are the same
  table1$pretruth[table1$Preference==1]<-table1$Trt[table1$Preference==1]
  # for people not having a preference(random arm),their true preference is generated from Bern(phi).
  table1$pretruth[table1$Preference==0]<-rbinom(sum(table1$Preference==0),1,phi)
  
  #generate counterfactuals using homogeneous effects model if user specify it
  if(model=="Homo"){
  #calculate parameters in homogeneous effects model based on true effects
  tao1<-trt/2
  tao2<--trt/2
  sigma1<-(1-phi)*sel
  sigma2<--phi*sel
  pi11<-(1-phi)*pre
  pi12<--phi*pre
  pi21<--(1-phi)*pre
  pi22<-phi*pre
  #check the truth (if these values are the same as input, then the above equation solver is correct)
  dtrt=tao1-tao2
  dsel=sigma1-sigma2
  dpre=(pi11+pi22-pi12-pi21)/2
  
  # generate four counterfactuals Y11,Y12,Y21,Y22
  table1$mu11<-tao1+sigma1+pi11+0.1*table1$Var1-0.15*table1$Var2+rnorm(n,mean=0,sd=1)
  table1$mu12<-tao1+sigma2+pi12+0.1*table1$Var1-0.15*table1$Var2+rnorm(n,mean=0,sd=1)
  table1$mu21<-tao2+sigma1+pi21+0.1*table1$Var1-0.15*table1$Var2+rnorm(n,mean=0,sd=1)
  table1$mu22<-tao2+sigma2+pi22+0.1*table1$Var1-0.15*table1$Var2+rnorm(n,mean=0,sd=1)
  }
  
  #generate counterfactuals using heterogeneous effects model if user specify it
  if(model=="Hetero"){
    #specify gamma110,gamma111,gamma112,gamma120,gamma121,gamma122,gamma210,gamma211,gamma212,gamma220,gamma221,gamma222
    gamma111<-0.8; gamma112<-0.2
    gamma121<-0.4; gamma122<-0.4
    gamma211<-0.7; gamma212<-0.6
    gamma221<-0.2; gamma222<-0.8
    
    # solve the other parameters
    gamma220<-0.1
    gamma210<- sel-pre+(gamma220+gamma222*p_Z2)-gamma212*p_Z2
    mu21=gamma210+gamma212/2
    mu22=gamma220+gamma222/2
    
    gamma110<-trt+phi*mu21+(1-phi)*mu22+(1-phi)*(sel+pre)-gamma112*p_Z2
    mu11=gamma110+gamma112/2
    
    gamma120<-mu11-(sel+pre)-gamma122*p_Z2
    mu12=gamma120+gamma122/2
    
    # check the truth (if these values are the same as input, then the above equation solver is correct)
    dtrt=phi*(mu11-mu21)+(1-phi)*(mu12-mu22)
    dsel=(mu11+mu21-mu12-mu22)/2
    dpre=(mu11-mu21-mu12+mu22)/2   
    # generate four counterfactuals Y11,Y12,Y21,Y22
    table1$mu11<-gamma110+gamma111*table1$Var1+gamma112*table1$Var2+rnorm(n,mean=0,sd=1)
    table1$mu12<-gamma120+gamma121*table1$Var1+gamma122*table1$Var2+rnorm(n,mean=0,sd=1)
    table1$mu21<-gamma210+gamma211*table1$Var1+gamma212*table1$Var2+rnorm(n,mean=0,sd=1)
    table1$mu22<-gamma220+gamma221*table1$Var1+gamma222*table1$Var2+rnorm(n,mean=0,sd=1)
  }
  
  #initiate true response Y
  table1$Y<-9
  #Get the true response Y based on patient's actual treatment assigned and preferred
  table1$Y[table1$Preference==1&table1$Trt==1]<-table1$mu11[table1$Preference==1&table1$Trt==1]
  table1$Y[table1$Preference==1&table1$Trt==0]<-table1$mu22[table1$Preference==1&table1$Trt==0]
  table1$Y[table1$Preference==0&table1$Trt==1&table1$pretruth==1]<-
    table1$mu11[table1$Preference==0&table1$Trt==1&table1$pretruth==1]
  table1$Y[table1$Preference==0&table1$Trt==1&table1$pretruth==0]<-
    table1$mu12[table1$Preference==0&table1$Trt==1&table1$pretruth==0]
  table1$Y[table1$Preference==0&table1$Trt==0&table1$pretruth==1]<-
    table1$mu21[table1$Preference==0&table1$Trt==0&table1$pretruth==1]
  table1$Y[table1$Preference==0&table1$Trt==0&table1$pretruth==0]<-
    table1$mu22[table1$Preference==0&table1$Trt==0&table1$pretruth==0]
  
  # label each patient's stratum number based on true PS, estimated PS (both X1 and X2), estimated PS (only X1), estimated PS (only X2) stratification
  table1$str=table1$str.hat=table1$str.hat.X1=99
  table1$str.hat.X2=strata
  for(s in 1:strata){
    table1$str[(table1$PS<=quan[s+1]) & (table1$PS>quan[s])]=s
    table1$str.hat[(table1$PS.hat<=quan.hat[s+1]) & (table1$PS.hat>quan.hat[s])]=s
    table1$str.hat.X1[(table1$PS.hat.X1<=quan.hat.X1[s+1]) & (table1$PS.hat.X1>quan.hat.X1[s])]=s
  }
  table1$count<-1
  table1<-table1[order(table1$PS.hat.X2),]
  width<-ceiling(n/strata)
  for(s in 1:(strata-1)){
    table1$str.hat.X2[(width*(s-1)+1):(width*s)]=s
  }
  
  
  # return output: the dataset and the vector of true effects
  return(list(table1=table1,truth=c(dtrt,dsel,dpre)))
}

#################################################################################################################################################
#Function: genresult

#Author: Yumin Wang, Fan Li

#Creation Date: October 10, 2020 (version 3.6.1)

#Purpose: This function returns the simulation results for empirical estimated effects, empirical variance, 
#         empirical power (or empirical Type I error) and the corrsponding true values for treatment, selection, preference effects under
#         different parameter settings (phi, theta, strata, n, trt, sel, pre), stratification strategies (true PS, estimated PS stratification) 
#         and data generating methods (Homogeneous or Heterogeneous effects model).


# Required Parameters: 
#          phi: This is the preference rate for treatment A, we assume it is constant in the entire study population. Range from 0 to 1.
#        theta: Marginal proportion of individual having a preference (marginal choice arm and random arm allocation ratio). Range from 0.1 to 0.9, step=0.05.
#       strata: Number of propensity score strata for PRPD. Range starting from 1.
#            n: Number of patients in the generated dataset
#          trt: True treatment effect of the simulation population
#          sel: True selection effect of the simulation population
#          pre: True preference effect of the simulation population
#        model: The model used to generate counterfacutals. Homo: Homogeneous effects model. Hetero: Heterogeneous effects model.     
#         seed: Specify to reproduce results


#Optional Parameters: 
#        alpha: Nominal type I error rate, default=0.05
#         nism: Number of simulation to be conducted, default=1000
     
#Output:  The output first displays four lines of text which includes: 
#         (1) The true treatment, selection, preference effects of the simulation population
#         (2) Empirical estimated treatment effect without stratification
#         (3) Empirical estimated selection effect without stratification
#         (4) Empirical estimated preference effect without stratification
#         Then the output displays a summary table, some basic information of the summary table:
#         Rows (1)-(5) are for treatment effects, rows (6)-(10) are for selection effects, and rows (11)-(15) are for preference effects
#         Column (1) is for effect estimates, Column (2) is for variance, Column (3) is for power (or Type I error rate)

#Other Functions: gendata()

#Example:
#        The following gives the summary table when phi=0.3,theta=0.3,strata=5,n=1000,trearment effect=0.4,selection effect=0.4,
#        preference effect=0.4,alpha=0.05,nsim=1000,seed=98431 under homogeneous effects model.
#        genresults(phi=0.3,theta=0.3,strata=5,n=1000,trt=0.4,sel=0.4,pre=0.4,model="Homo",alpha=0.05,nsim=1000,seed=98431)

#################################################################################################################################################
genresult=function(phi,theta,strata,n,trt,sel,pre,model,alpha=0.05,nsim=1000,seed){
  #set seed to reproduce results
  set.seed(seed)
  #compute predicted variance, power for treatment, selection, preference effects
  #first generate a dataset with 100,000 individuals in order to get true values for stratum-specfic theta, stratum-specific sigma
  #stratum-specific phi and stratum-specific selection, preference effects
  restrue=gendata(phi=phi,theta=theta,strata=strata,n=100000,trt=trt,sel=sel,pre=pre,model=model)
  tabletrue=restrue$table1
  #extract stratum-specific theta
  Theta=tapply(tabletrue$Preference,tabletrue$str,mean)
  #extract stratum-specific sigma
  SigmaR=tapply(tabletrue$Y[tabletrue$Preference==0],tabletrue$str[tabletrue$Preference==0],var)
  SigmaCA=tapply(tabletrue$Y[tabletrue$Preference==1&tabletrue$Trt==1],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==1],var)
  SigmaCB=tapply(tabletrue$Y[tabletrue$Preference==1&tabletrue$Trt==0],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==0],var)
  #extract stratum-specific phi
  Phi<-tapply(tabletrue$Trt[tabletrue$Preference==1],tabletrue$str[tabletrue$Preference==1],mean)
  #extract stratum-specific selection,preference effects
  Z1<-tapply(tabletrue$Y[tabletrue$Preference==1&tabletrue$Trt==1],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==1],sum)-tapply(tabletrue$count[tabletrue$Preference==1&tabletrue$Trt==1],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==1],sum)*tapply(tabletrue$Y[tabletrue$Preference==0&tabletrue$Trt==1],tabletrue$str[tabletrue$Preference==0&tabletrue$Trt==1],mean)
  Z2<-tapply(tabletrue$Y[tabletrue$Preference==1&tabletrue$Trt==0],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==0],sum)-tapply(tabletrue$count[tabletrue$Preference==1&tabletrue$Trt==0],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==0],sum)*tapply(tabletrue$Y[tabletrue$Preference==0&tabletrue$Trt==0],tabletrue$str[tabletrue$Preference==0&tabletrue$Trt==0],mean)
  M<-tapply(tabletrue$count[tabletrue$Preference==1],tabletrue$str[tabletrue$Preference==1],sum)
  DeltaNu<-(Z1-Z2)/(2*Phi*(1-Phi)*M)
  DeltaPi<-(Z1+Z2)/(2*Phi*(1-Phi)*M)
  #calculate predicted variances for treatment, selection, preference effects based on formula
  pred.var.trt=4/(n*strata)*sum(SigmaR/(1-Theta))
  pred.var.sel=1/(4*n*strata)*sum((Phi*SigmaCA+(1-Phi)*SigmaCB+Phi*(1-Phi)*((2*Phi-1)*DeltaNu+DeltaPi)^2+2*SigmaR*(Phi^2+(1-Phi)^2)*(Theta)/(1-Theta))/(Theta*Phi^2*(1-Phi)^2))
  pred.var.pre=1/(4*n*strata)*sum((Phi*SigmaCA+(1-Phi)*SigmaCB+Phi*(1-Phi)*((2*Phi-1)*DeltaPi+DeltaNu)^2+2*SigmaR*(Phi^2+(1-Phi)^2)*(Theta)/(1-Theta))/(Theta*Phi^2*(1-Phi)^2))
  #calculate predicted power for treatment, selection, preference effects based on formula
  pred.power.trt=pnorm(abs(trt)/sqrt(pred.var.trt)-qnorm(1-alpha/2))
  pred.power.sel=pnorm(abs(sel)/sqrt(pred.var.sel)-qnorm(1-alpha/2))
  pred.power.pre=pnorm(abs(pre)/sqrt(pred.var.pre)-qnorm(1-alpha/2))
  
  #initiate vectors to store values
  #store estimated treatment, selection, preference effects without subclassification
  EST_0_trt=NULL 
  EST_0_sel=NULL 
  EST_0_pre=NULL
  # store overall point estimates for treatment, selection, preference effects with different stratification strategies
  EST_tPS_trt=EST_ePS_trt=EST_ePS_X1_trt=EST_ePS_X2_trt=NULL   
  EST_tPS_sel=EST_ePS_sel=EST_ePS_X1_sel=EST_ePS_X2_sel=NULL   
  EST_tPS_pre=EST_ePS_pre=EST_ePS_X1_pre=EST_ePS_X2_pre=NULL 
  # store p-values for treatment, selection, preference effects with different stratification strategies
  PVAL_tPS_trt=PVAL_ePS_trt=PVAL_ePS_X1_trt=PVAL_ePS_X2_trt=NULL 
  PVAL_tPS_sel=PVAL_ePS_sel=PVAL_ePS_X1_sel=PVAL_ePS_X2_sel=NULL 
  PVAL_tPS_pre=PVAL_ePS_pre=PVAL_ePS_X1_pre=PVAL_ePS_X2_pre=NULL 
  
  # carry out nsim simulations
  for(i in 1:nsim){
    #generate one table
    table1=gendata(phi=phi,theta=theta,strata=strata,n=n,trt=trt,sel=sel,pre=pre,model=model)$table1
    
    # stratum-specific effect estimators and variances for treatment effect with different stratification strategies
    est.stratum.trt=est.stratum.hat.trt=est.stratum.hat.X1.trt=est.stratum.hat.X2.trt=rep(NA,strata)
    var.stratum.trt=var.stratum.hat.trt=var.stratum.hat.X1.trt=var.stratum.hat.X2.trt=rep(NA,strata)
    # stratum-specific effect estimators and variances for selection effect with different stratification strategies
    est.stratum.sel=est.stratum.hat.sel=est.stratum.hat.X1.sel=est.stratum.hat.X2.sel=rep(NA,strata)
    var.stratum.sel=var.stratum.hat.sel=var.stratum.hat.X1.sel=var.stratum.hat.X2.sel=rep(NA,strata)
    # stratum-specific effect estimators and variances for preference effect with different stratification strategies
    est.stratum.pre=est.stratum.hat.pre=est.stratum.hat.X1.pre=est.stratum.hat.X2.pre=rep(NA,strata)
    var.stratum.pre=var.stratum.hat.pre=var.stratum.hat.X1.pre=var.stratum.hat.X2.pre=rep(NA,strata)
    #consider strata one by one
    for(s in 1:strata){
      #subset for random arm in stratum s with different stratification strategies
      temp.r=subset(table1,str==s & Arm=="Random Arm",drop=T)
      temp.hat.r=subset(table1,str.hat==s & Arm=="Random Arm",drop=T)
      temp.hat.X1.r=subset(table1,str.hat.X1==s & Arm=="Random Arm",drop=T)
      temp.hat.X2.r=subset(table1,str.hat.X2==s & Arm=="Random Arm",drop=T)
      
      #subset for choice arm in stratum s with different stratification strategies
      temp.c=subset(table1,str==s & Arm=="Choice Arm",drop=T)
      temp.hat.c=subset(table1,str.hat==s & Arm=="Choice Arm",drop=T)
      temp.hat.X1.c=subset(table1,str.hat.X1==s & Arm=="Choice Arm",drop=T)
      temp.hat.X2.c=subset(table1,str.hat.X2==s & Arm=="Choice Arm",drop=T)
      
      #calculate phi in stratum s with different stratification strategies
      temp.phi<-mean(temp.c$Trt)
      temp.hat.phi<-mean(temp.hat.c$Trt)
      temp.hat.X1.phi<-mean(temp.hat.X1.c$Trt)
      temp.hat.X2.phi<-mean(temp.hat.X2.c$Trt)
      
      #calculate d1,d2 in stratum s with different stratification strategies
      d1<-mean(temp.c$Y[temp.c$Trt==1])-mean(temp.r$Y[temp.r$Trt==1])
      d2<-mean(temp.c$Y[temp.c$Trt==0])-mean(temp.r$Y[temp.r$Trt==0])
      
      d1.hat<-mean(temp.hat.c$Y[temp.hat.c$Trt==1])-mean(temp.hat.r$Y[temp.hat.r$Trt==1])
      d2.hat<-mean(temp.hat.c$Y[temp.hat.c$Trt==0])-mean(temp.hat.r$Y[temp.hat.r$Trt==0])
      
      d1.hat.X1<-mean(temp.hat.X1.c$Y[temp.hat.X1.c$Trt==1])-mean(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==1])
      d2.hat.X1<-mean(temp.hat.X1.c$Y[temp.hat.X1.c$Trt==0])-mean(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==0])
      
      d1.hat.X2<-mean(temp.hat.X2.c$Y[temp.hat.X2.c$Trt==1])-mean(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==1])
      d2.hat.X2<-mean(temp.hat.X2.c$Y[temp.hat.X2.c$Trt==0])-mean(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==0])
      
      #calculate number of people in treatment A or B in choice or random arm in stratum s with different stratification strategies
      m1l=sum(temp.c$count[temp.c$Trt==1])
      m2l=sum(temp.c$count[temp.c$Trt==0])
      n1l=sum(temp.r$count[temp.r$Trt==1])
      n2l=sum(temp.r$count[temp.r$Trt==0])
      ml=m1l+m2l
      
      m1l.hat=sum(temp.hat.c$count[temp.hat.c$Trt==1])
      m2l.hat=sum(temp.hat.c$count[temp.hat.c$Trt==0])
      n1l.hat=sum(temp.hat.r$count[temp.hat.r$Trt==1])
      n2l.hat=sum(temp.hat.r$count[temp.hat.r$Trt==0])
      ml.hat=m1l.hat+m2l.hat
      
      m1l.hat.X1=sum(temp.hat.X1.c$count[temp.hat.X1.c$Trt==1])
      m2l.hat.X1=sum(temp.hat.X1.c$count[temp.hat.X1.c$Trt==0])
      n1l.hat.X1=sum(temp.hat.X1.r$count[temp.hat.X1.r$Trt==1])
      n2l.hat.X1=sum(temp.hat.X1.r$count[temp.hat.X1.r$Trt==0])
      ml.hat.X1=m1l.hat.X1+m2l.hat.X1
      
      m1l.hat.X2=sum(temp.hat.X2.c$count[temp.hat.X2.c$Trt==1])
      m2l.hat.X2=sum(temp.hat.X2.c$count[temp.hat.X2.c$Trt==0])
      n1l.hat.X2=sum(temp.hat.X2.r$count[temp.hat.X2.r$Trt==1])
      n2l.hat.X2=sum(temp.hat.X2.r$count[temp.hat.X2.r$Trt==0])
      ml.hat.X2=m1l.hat.X2+m2l.hat.X2
      
      
      #calculate Z1,Z2, var(Z1+Z2), var(Z1-Z2) in stratum s with different stratification strategies
      temp.Z1=sum(temp.c$Y[temp.c$Trt==1])-m1l*mean(temp.r$Y[temp.r$Trt==1])
      temp.Z2=sum(temp.c$Y[temp.c$Trt==0])-m2l*mean(temp.r$Y[temp.r$Trt==0])
      temp.var.Z1=m1l*var(temp.c$Y[temp.c$Trt==1])+(1+(ml-1)*m1l/ml)*m1l*var(temp.r$Y[temp.r$Trt==1])/n1l+(m1l*m2l)*(mean(temp.c$Y[temp.c$Trt==1])-mean(temp.r$Y[temp.r$Trt==1]))^2/ml
      temp.var.Z2=m2l*var(temp.c$Y[temp.c$Trt==0])+(1+(ml-1)*m2l/ml)*m2l*var(temp.r$Y[temp.r$Trt==0])/n2l+(m1l*m2l)*(mean(temp.c$Y[temp.c$Trt==0])-mean(temp.r$Y[temp.r$Trt==0]))^2/ml
      temp.covZ1Z2=-(m1l*m2l)*(mean(temp.c$Y[temp.c$Trt==1])-mean(temp.r$Y[temp.r$Trt==1]))*(mean(temp.c$Y[temp.c$Trt==0])-mean(temp.r$Y[temp.r$Trt==0]))/ml
      temp.var.Z1plusZ2=temp.var.Z1+temp.var.Z2+2*temp.covZ1Z2
      temp.var.Z1minusZ2=temp.var.Z1+temp.var.Z2-2*temp.covZ1Z2
      
      temp.hat.Z1=sum(temp.hat.c$Y[temp.hat.c$Trt==1])-m1l.hat*mean(temp.hat.r$Y[temp.hat.r$Trt==1])
      temp.hat.Z2=sum(temp.hat.c$Y[temp.hat.c$Trt==0])-m2l.hat*mean(temp.hat.r$Y[temp.hat.r$Trt==0])
      temp.hat.var.Z1=m1l.hat*var(temp.hat.c$Y[temp.hat.c$Trt==1])+(1+(ml.hat-1)*m1l.hat/ml.hat)*m1l.hat*var(temp.hat.r$Y[temp.hat.r$Trt==1])/n1l.hat+(m1l.hat*m2l.hat)*(mean(temp.hat.c$Y[temp.hat.c$Trt==1])-mean(temp.hat.r$Y[temp.hat.r$Trt==1]))^2/ml.hat
      temp.hat.var.Z2=m2l.hat*var(temp.hat.c$Y[temp.hat.c$Trt==0])+(1+(ml.hat-1)*m2l.hat/ml.hat)*m2l.hat*var(temp.hat.r$Y[temp.hat.r$Trt==0])/n2l.hat+(m1l.hat*m2l.hat)*(mean(temp.hat.c$Y[temp.hat.c$Trt==0])-mean(temp.hat.r$Y[temp.hat.r$Trt==0]))^2/ml.hat
      temp.hat.covZ1Z2=-(m1l.hat*m2l.hat)*(mean(temp.hat.c$Y[temp.hat.c$Trt==1])-mean(temp.hat.r$Y[temp.hat.r$Trt==1]))*(mean(temp.hat.c$Y[temp.hat.c$Trt==0])-mean(temp.hat.r$Y[temp.hat.r$Trt==0]))/ml.hat
      temp.hat.var.Z1plusZ2=temp.hat.var.Z1+temp.hat.var.Z2+2*temp.hat.covZ1Z2
      temp.hat.var.Z1minusZ2=temp.hat.var.Z1+temp.hat.var.Z2-2*temp.hat.covZ1Z2
      
      temp.hat.X1.Z1=sum(temp.hat.X1.c$Y[temp.hat.X1.c$Trt==1])-m1l.hat.X1*mean(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==1])
      temp.hat.X1.Z2=sum(temp.hat.X1.c$Y[temp.hat.X1.c$Trt==0])-m2l.hat.X1*mean(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==0])
      temp.hat.X1.var.Z1=m1l.hat.X1*var(temp.hat.X1.c$Y[temp.hat.X1.c$Trt==1])+(1+(ml.hat.X1-1)*m1l.hat.X1/ml.hat.X1)*m1l.hat.X1*var(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==1])/n1l.hat.X1+(m1l.hat.X1*m2l.hat.X1)*(mean(temp.hat.X1.c$Y[temp.hat.X1.c$Trt==1])-mean(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==1]))^2/ml.hat.X1
      temp.hat.X1.var.Z2=m2l.hat.X1*var(temp.hat.X1.c$Y[temp.hat.X1.c$Trt==0])+(1+(ml.hat.X1-1)*m2l.hat.X1/ml.hat.X1)*m2l.hat.X1*var(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==0])/n2l.hat.X1+(m1l.hat.X1*m2l.hat.X1)*(mean(temp.hat.X1.c$Y[temp.hat.X1.c$Trt==0])-mean(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==0]))^2/ml.hat.X1
      temp.hat.X1.covZ1Z2=-(m1l.hat.X1*m2l.hat.X1)*(mean(temp.hat.X1.c$Y[temp.hat.X1.c$Trt==1])-mean(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==1]))*(mean(temp.hat.X1.c$Y[temp.hat.X1.c$Trt==0])-mean(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==0]))/ml.hat.X1
      temp.hat.X1.var.Z1plusZ2=temp.hat.X1.var.Z1+temp.hat.X1.var.Z2+2*temp.hat.X1.covZ1Z2
      temp.hat.X1.var.Z1minusZ2=temp.hat.X1.var.Z1+temp.hat.X1.var.Z2-2*temp.hat.X1.covZ1Z2
      
      temp.hat.X2.Z1=sum(temp.hat.X2.c$Y[temp.hat.X2.c$Trt==1])-m1l.hat.X2*mean(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==1])
      temp.hat.X2.Z2=sum(temp.hat.X2.c$Y[temp.hat.X2.c$Trt==0])-m2l.hat.X2*mean(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==0])
      temp.hat.X2.var.Z1=m1l.hat.X2*var(temp.hat.X2.c$Y[temp.hat.X2.c$Trt==1])+(1+(ml.hat.X2-1)*m1l.hat.X2/ml.hat.X2)*m1l.hat.X2*var(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==1])/n1l.hat.X2+(m1l.hat.X2*m2l.hat.X2)*(mean(temp.hat.X2.c$Y[temp.hat.X2.c$Trt==1])-mean(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==1]))^2/ml.hat.X2
      temp.hat.X2.var.Z2=m2l.hat.X2*var(temp.hat.X2.c$Y[temp.hat.X2.c$Trt==0])+(1+(ml.hat.X2-1)*m2l.hat.X2/ml.hat.X2)*m2l.hat.X2*var(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==0])/n2l.hat.X2+(m1l.hat.X2*m2l.hat.X2)*(mean(temp.hat.X2.c$Y[temp.hat.X2.c$Trt==0])-mean(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==0]))^2/ml.hat.X2
      temp.hat.X2.covZ1Z2=-(m1l.hat.X2*m2l.hat.X2)*(mean(temp.hat.X2.c$Y[temp.hat.X2.c$Trt==1])-mean(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==1]))*(mean(temp.hat.X2.c$Y[temp.hat.X2.c$Trt==0])-mean(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==0]))/ml.hat.X2
      temp.hat.X2.var.Z1plusZ2=temp.hat.X2.var.Z1+temp.hat.X2.var.Z2+2*temp.hat.X2.covZ1Z2
      temp.hat.X2.var.Z1minusZ2=temp.hat.X2.var.Z1+temp.hat.X2.var.Z2-2*temp.hat.X2.covZ1Z2
      
      # stratum-specific effect estimators and variances for treatment effect with true PS stratification
      est.stratum.trt[s]<-mean(temp.r$Y[temp.r$Trt==1])-mean(temp.r$Y[temp.r$Trt==0])
      var.stratum.trt[s]<-var(temp.r$Y)*(1/sum(temp.r$Trt==1)+1/sum(temp.r$Trt==0))
      # stratum-specific effect estimators and variances for treatment effect with estimated PS stratification (both X1 and X2)
      est.stratum.hat.trt[s]<-mean(temp.hat.r$Y[temp.hat.r$Trt==1])-mean(temp.hat.r$Y[temp.hat.r$Trt==0])
      var.stratum.hat.trt[s]<-var(temp.hat.r$Y)*(1/sum(temp.hat.r$Trt==1)+1/sum(temp.hat.r$Trt==0))
      
      # stratum-specific effect estimators and variances for treatment effect with estimated PS stratification (only X1)
      est.stratum.hat.X1.trt[s]<-mean(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==1])-mean(temp.hat.X1.r$Y[temp.hat.X1.r$Trt==0])
      var.stratum.hat.X1.trt[s]<-var(temp.hat.X1.r$Y)*(1/sum(temp.hat.X1.r$Trt==1)+1/sum(temp.hat.X1.r$Trt==0))
      
      # stratum-specific effect estimators and variances for treatment effect with estimated PS stratification (only X2)
      est.stratum.hat.X2.trt[s]<-mean(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==1])-mean(temp.hat.X2.r$Y[temp.hat.X2.r$Trt==0])
      var.stratum.hat.X2.trt[s]<-var(temp.hat.X2.r$Y)*(1/sum(temp.hat.X2.r$Trt==1)+1/sum(temp.hat.X2.r$Trt==0))
      
      # stratum-specific effect estimators and variances for selection effect with true PS stratification
      est.stratum.sel[s]<-(temp.Z1-temp.Z2)/(2*temp.phi*(1-temp.phi)*ml)
      var.stratum.sel[s]<-(temp.var.Z1minusZ2+2*(m1l*d1-m2l*d2)*(2*temp.phi-1)*(d1+d2)+(2*temp.phi-1)^2*(m1l*d1-m2l*d2)^2/(ml*temp.phi*(1-temp.phi)))/(4*temp.phi^2*(1-temp.phi)^2*ml^2)
      # stratum-specific effect estimators and variances for selection effect with estimated PS stratification (both X1 and X2)
      est.stratum.hat.sel[s]<-(temp.hat.Z1-temp.hat.Z2)/(2*temp.hat.phi*(1-temp.hat.phi)*ml.hat)
      var.stratum.hat.sel[s]<-(temp.hat.var.Z1minusZ2+2*(m1l.hat*d1.hat-m2l.hat*d2.hat)*(2*temp.hat.phi-1)*(d1.hat+d2.hat)+(2*temp.hat.phi-1)^2*(m1l.hat*d1.hat-m2l.hat*d2.hat)^2/(ml.hat*temp.hat.phi*(1-temp.hat.phi)))/(4*temp.hat.phi^2*(1-temp.hat.phi)^2*ml.hat^2)
      
      # stratum-specific effect estimators and variances for selection effect with estimated PS stratification (only X1)
      est.stratum.hat.X1.sel[s]<-(temp.hat.X1.Z1-temp.hat.X1.Z2)/(2*temp.hat.X1.phi*(1-temp.hat.X1.phi)*ml.hat.X1)
      var.stratum.hat.X1.sel[s]<-(temp.hat.X1.var.Z1minusZ2+2*(m1l.hat.X1*d1.hat.X1-m2l.hat.X1*d2.hat.X1)*(2*temp.hat.X1.phi-1)*(d1.hat.X1+d2.hat.X1)+(2*temp.hat.X1.phi-1)^2*(m1l.hat.X1*d1.hat.X1-m2l.hat.X1*d2.hat.X1)^2/(ml.hat.X1*temp.hat.X1.phi*(1-temp.hat.X1.phi)))/(4*temp.hat.X1.phi^2*(1-temp.hat.X1.phi)^2*ml.hat.X1^2)
      
      # stratum-specific effect estimators and variances for selection effect with estimated PS stratification (only X2)
      est.stratum.hat.X2.sel[s]<-(temp.hat.X2.Z1-temp.hat.X2.Z2)/(2*temp.hat.X2.phi*(1-temp.hat.X2.phi)*ml.hat.X2)
      var.stratum.hat.X2.sel[s]<-(temp.hat.X2.var.Z1minusZ2+2*(m1l.hat.X2*d1.hat.X2-m2l.hat.X2*d2.hat.X2)*(2*temp.hat.X2.phi-1)*(d1.hat.X2+d2.hat.X2)+(2*temp.hat.X2.phi-1)^2*(m1l.hat.X2*d1.hat.X2-m2l.hat.X2*d2.hat.X2)^2/(ml.hat.X2*temp.hat.X2.phi*(1-temp.hat.X2.phi)))/(4*temp.hat.X2.phi^2*(1-temp.hat.X2.phi)^2*ml.hat.X2^2)
      
      
      # stratum-specific effect estimators and variances for preference effect with true PS stratification
      est.stratum.pre[s]<-(temp.Z1+temp.Z2)/(2*temp.phi*(1-temp.phi)*ml)
      var.stratum.pre[s]<-(temp.var.Z1plusZ2+2*(m1l*d1+m2l*d2)*(2*temp.phi-1)*(d1-d2)+(2*temp.phi-1)^2*(m1l*d1+m2l*d2)^2/(ml*temp.phi*(1-temp.phi)))/(4*temp.phi^2*(1-temp.phi)^2*ml^2)
      # stratum-specific effect estimators and variances for preference effect with estimated PS stratification (both X1 and X2)
      est.stratum.hat.pre[s]<-(temp.hat.Z1+temp.hat.Z2)/(2*temp.hat.phi*(1-temp.hat.phi)*ml.hat)
      var.stratum.hat.pre[s]<-(temp.hat.var.Z1plusZ2+2*(m1l.hat*d1.hat+m2l.hat*d2.hat)*(2*temp.hat.phi-1)*(d1.hat-d2.hat)+(2*temp.hat.phi-1)^2*(m1l.hat*d1.hat+m2l.hat*d2.hat)^2/(ml.hat*temp.hat.phi*(1-temp.hat.phi)))/(4*temp.hat.phi^2*(1-temp.hat.phi)^2*ml.hat^2)
      
      # stratum-specific effect estimators and variances for preference effect with estimated PS stratification (only X1)
      est.stratum.hat.X1.pre[s]<-(temp.hat.X1.Z1+temp.hat.X1.Z2)/(2*temp.hat.X1.phi*(1-temp.hat.X1.phi)*ml.hat.X1)
      var.stratum.hat.X1.pre[s]<-(temp.hat.X1.var.Z1plusZ2+2*(m1l.hat.X1*d1.hat.X1+m2l.hat.X1*d2.hat.X1)*(2*temp.hat.X1.phi-1)*(d1.hat.X1-d2.hat.X1)+(2*temp.hat.X1.phi-1)^2*(m1l.hat.X1*d1.hat.X1+m2l.hat.X1*d2.hat.X1)^2/(ml.hat.X1*temp.hat.X1.phi*(1-temp.hat.X1.phi)))/(4*temp.hat.X1.phi^2*(1-temp.hat.X1.phi)^2*ml.hat.X1^2)
      
      # stratum-specific effect estimators and variances for preference effect with estimated PS stratification (only X2)
      est.stratum.hat.X2.pre[s]<-(temp.hat.X2.Z1+temp.hat.X2.Z2)/(2*temp.hat.X2.phi*(1-temp.hat.X2.phi)*ml.hat.X2)
      var.stratum.hat.X2.pre[s]<-(temp.hat.X2.var.Z1plusZ2+2*(m1l.hat.X2*d1.hat.X2+m2l.hat.X2*d2.hat.X2)*(2*temp.hat.X2.phi-1)*(d1.hat.X2-d2.hat.X2)+(2*temp.hat.X2.phi-1)^2*(m1l.hat.X2*d1.hat.X2+m2l.hat.X2*d2.hat.X2)^2/(ml.hat.X2*temp.hat.X2.phi*(1-temp.hat.X2.phi)))/(4*temp.hat.X2.phi^2*(1-temp.hat.X2.phi)^2*ml.hat.X2^2)
      
    }
    #store empirical estimated treatment,selection, preference effects without stratification
    EST_0_trt=c(EST_0_trt, mean(table1$Y[(table1$Arm=="Random Arm") & (table1$Trt==1)])-mean(table1$Y[(table1$Arm=="Random Arm") & (table1$Trt==0)]))
    PHI<-mean(table1$Trt[table1$Arm=="Choice Arm"])
    M1L<-sum(table1$count[(table1$Arm=="Choice Arm")&(table1$Trt==1)])
    M2L<-sum(table1$count[(table1$Arm=="Choice Arm")&(table1$Trt==0)])
    ML<-M1L+M2L
    EST_0_sel=c(EST_0_sel,(sum(table1$Y[(table1$Arm=="Choice Arm")&(table1$Trt==1)])-M1L*mean(table1$Y[(table1$Arm=="Random Arm")&(table1$Trt==1)])-(sum(table1$Y[(table1$Arm=="Choice Arm")&(table1$Trt==0)])-M2L*mean(table1$Y[(table1$Arm=="Random Arm")&(table1$Trt==0)])))/(2*(PHI)*(1-PHI)*ML))
    EST_0_pre=c(EST_0_pre,(sum(table1$Y[(table1$Arm=="Choice Arm")&(table1$Trt==1)])-M1L*mean(table1$Y[(table1$Arm=="Random Arm")&(table1$Trt==1)])+(sum(table1$Y[(table1$Arm=="Choice Arm")&(table1$Trt==0)])-M2L*mean(table1$Y[(table1$Arm=="Random Arm")&(table1$Trt==0)])))/(2*(PHI)*(1-PHI)*ML))
    
    #store overall treatment effect by averaging stratum-specific one with different stratification strategies
    EST_tPS_trt=c(EST_tPS_trt, mean(est.stratum.trt))
    EST_ePS_trt=c(EST_ePS_trt, mean(est.stratum.hat.trt))
    EST_ePS_X1_trt=c(EST_ePS_X1_trt, mean(est.stratum.hat.X1.trt))
    EST_ePS_X2_trt=c(EST_ePS_X2_trt, mean(est.stratum.hat.X2.trt))
    #store p values for treatment effect with different stratification strategies
    PVAL_tPS_trt=c(PVAL_tPS_trt, 2*(1-pnorm(abs(mean(est.stratum.trt))/sqrt(sum(var.stratum.trt)/strata^2))))
    PVAL_ePS_trt=c(PVAL_ePS_trt, 2*(1-pnorm(abs(mean(est.stratum.hat.trt))/sqrt(sum(var.stratum.hat.trt)/strata^2))))
    PVAL_ePS_X1_trt=c(PVAL_ePS_X1_trt, 2*(1-pnorm(abs(mean(est.stratum.hat.X1.trt))/sqrt(sum(var.stratum.hat.X1.trt)/strata^2))))
    PVAL_ePS_X2_trt=c(PVAL_ePS_X2_trt, 2*(1-pnorm(abs(mean(est.stratum.hat.X2.trt))/sqrt(sum(var.stratum.hat.X2.trt)/strata^2))))
    
    #store overall selection effect by averaging stratum-specific one with different stratification strategies
    EST_tPS_sel=c(EST_tPS_sel, mean(est.stratum.sel))
    EST_ePS_sel=c(EST_ePS_sel, mean(est.stratum.hat.sel))
    EST_ePS_X1_sel=c(EST_ePS_X1_sel, mean(est.stratum.hat.X1.sel))
    EST_ePS_X2_sel=c(EST_ePS_X2_sel, mean(est.stratum.hat.X2.sel))
    #store p values for selection effect with different stratification strategies
    PVAL_tPS_sel=c(PVAL_tPS_sel, 2*(1-pnorm(abs(mean(est.stratum.sel))/sqrt(sum(var.stratum.sel)/strata^2))))
    PVAL_ePS_sel=c(PVAL_ePS_sel, 2*(1-pnorm(abs(mean(est.stratum.hat.sel))/sqrt(sum(var.stratum.hat.sel)/strata^2))))
    PVAL_ePS_X1_sel=c(PVAL_ePS_X1_sel, 2*(1-pnorm(abs(mean(est.stratum.hat.X1.sel))/sqrt(sum(var.stratum.hat.X1.sel)/strata^2))))
    PVAL_ePS_X2_sel=c(PVAL_ePS_X2_sel, 2*(1-pnorm(abs(mean(est.stratum.hat.X2.sel))/sqrt(sum(var.stratum.hat.X2.sel)/strata^2))))
    
    #store overall preference effect by averaging stratum-specific one with different stratification strategies
    EST_tPS_pre=c(EST_tPS_pre, mean(est.stratum.pre))
    EST_ePS_pre=c(EST_ePS_pre, mean(est.stratum.hat.pre))
    EST_ePS_X1_pre=c(EST_ePS_X1_pre, mean(est.stratum.hat.X1.pre))
    EST_ePS_X2_pre=c(EST_ePS_X2_pre, mean(est.stratum.hat.X2.pre))
    #store p values for preference effect with different stratification strategies
    PVAL_tPS_pre=c(PVAL_tPS_pre, 2*(1-pnorm(abs(mean(est.stratum.pre))/sqrt(sum(var.stratum.pre)/strata^2))))
    PVAL_ePS_pre=c(PVAL_ePS_pre, 2*(1-pnorm(abs(mean(est.stratum.hat.pre))/sqrt(sum(var.stratum.hat.pre)/strata^2))))
    PVAL_ePS_X1_pre=c(PVAL_ePS_X1_pre, 2*(1-pnorm(abs(mean(est.stratum.hat.X1.pre))/sqrt(sum(var.stratum.hat.X1.pre)/strata^2))))
    PVAL_ePS_X2_pre=c(PVAL_ePS_X2_pre, 2*(1-pnorm(abs(mean(est.stratum.hat.X2.pre))/sqrt(sum(var.stratum.hat.X2.pre)/strata^2))))
    
    # loop control
    if(i%%2000==0) print(i)
  }
  # construct summary table
  # the true effect, true variance and true power for treatment, selection, preference effects
  gold.trt=c(trt,pred.var.trt,pred.power.trt)
  gold.sel=c(sel,pred.var.sel,pred.power.sel)
  gold.pre=c(pre,pred.var.pre,pred.power.pre)
  
  # set infinite values to NA for calculation purposes
  EST_tPS_sel[is.infinite(EST_tPS_sel)]<-NA
  EST_ePS_sel[is.infinite(EST_ePS_sel)]<-NA
  EST_ePS_X1_sel[is.infinite(EST_ePS_X1_sel)]<-NA
  EST_ePS_X2_sel[is.infinite(EST_ePS_X2_sel)]<-NA
  
  EST_tPS_pre[is.infinite(EST_tPS_pre)]<-NA
  EST_ePS_pre[is.infinite(EST_ePS_pre)]<-NA
  EST_ePS_X1_pre[is.infinite(EST_ePS_X1_pre)]<-NA
  EST_ePS_X2_pre[is.infinite(EST_ePS_X2_pre)]<-NA
  
  
  EST_tPS_trt[is.infinite(EST_tPS_trt)]<-NA
  EST_ePS_trt[is.infinite(EST_ePS_trt)]<-NA
  EST_ePS_X1_trt[is.infinite(EST_ePS_X1_trt)]<-NA
  EST_ePS_X2_trt[is.infinite(EST_ePS_X2_trt)]<-NA
  
  #results for treatment effect with different stratification strategies (empirial estimated effects, empirical variance, empirical power (alpha))
  res_tPS_trt=c(mean(EST_tPS_trt,na.rm=TRUE),var(EST_tPS_trt,na.rm=TRUE),mean(PVAL_tPS_trt<alpha,na.rm=TRUE))
  res_ePS_trt=c(mean(EST_ePS_trt,na.rm=TRUE),var(EST_ePS_trt,na.rm=TRUE),mean(PVAL_ePS_trt<alpha,na.rm=TRUE))
  res_ePS_X1_trt=c(mean(EST_ePS_X1_trt,na.rm=TRUE),var(EST_ePS_X1_trt,na.rm=TRUE),mean(PVAL_ePS_X1_trt<alpha,na.rm=TRUE))
  res_ePS_X2_trt=c(mean(EST_ePS_X2_trt,na.rm=TRUE),var(EST_ePS_X2_trt,na.rm=TRUE),mean(PVAL_ePS_X2_trt<alpha,na.rm=TRUE))
  
  #results for selection effect with different stratification strategies (empirial estimated effects, empirical variance, empirical power (alpha))
  res_tPS_sel=c(mean(EST_tPS_sel,na.rm=TRUE),var(EST_tPS_sel,na.rm=TRUE),mean(PVAL_tPS_sel<alpha,na.rm=TRUE))
  res_ePS_sel=c(mean(EST_ePS_sel,na.rm=TRUE),var(EST_ePS_sel,na.rm=TRUE),mean(PVAL_ePS_sel<alpha,na.rm=TRUE))
  res_ePS_X1_sel=c(mean(EST_ePS_X1_sel,na.rm=TRUE),var(EST_ePS_X1_sel,na.rm=TRUE),mean(PVAL_ePS_X1_sel<alpha,na.rm=TRUE))
  res_ePS_X2_sel=c(mean(EST_ePS_X2_sel,na.rm=TRUE),var(EST_ePS_X2_sel,na.rm=TRUE),mean(PVAL_ePS_X2_sel<alpha,na.rm=TRUE))
  
  #results for preference effect with different stratification strategies (empirial estimated effects, empirical variance, empirical power (alpha))
  res_tPS_pre=c(mean(EST_tPS_pre,na.rm=TRUE),var(EST_tPS_pre,na.rm=TRUE),mean(PVAL_tPS_pre<alpha,na.rm=TRUE))
  res_ePS_pre=c(mean(EST_ePS_pre,na.rm=TRUE),var(EST_ePS_pre,na.rm=TRUE),mean(PVAL_ePS_pre<alpha,na.rm=TRUE))
  res_ePS_X1_pre=c(mean(EST_ePS_X1_pre,na.rm=TRUE),var(EST_ePS_X1_pre,na.rm=TRUE),mean(PVAL_ePS_X1_pre<alpha,na.rm=TRUE))
  res_ePS_X2_pre=c(mean(EST_ePS_X2_pre,na.rm=TRUE),var(EST_ePS_X2_pre,na.rm=TRUE),mean(PVAL_ePS_X2_pre<alpha,na.rm=TRUE))
  
  #summary table
  sumtable=rbind(gold.trt,res_tPS_trt,res_ePS_trt,res_ePS_X1_trt,res_ePS_X2_trt,gold.sel,res_tPS_sel,res_ePS_sel,res_ePS_X1_sel,res_ePS_X2_sel,gold.pre,res_tPS_pre,res_ePS_pre,res_ePS_X1_pre,res_ePS_X2_pre)
  rownames(sumtable)=c("predicted","truePS","estimatedPS","estimatedPS-onlyX1","estimatedPS-onlyX2","predicted","truePS","estimatedPS","estimatedPS-onlyX1","estimatedPS-onlyX2","predicted","truePS","estimatedPS","estimatedPS-onlyX1","estimatedPS-onlyX2")
  colnames(sumtable)=c("dtau/dnu/dpi","var","power")
  
  # messages
  cat("true trt, sel, and pre effect verified as ", c(restrue$truth), "\n")
  cat("dtau with one strata has mean ", mean(EST_0_trt), "\n")
  cat("dnu with one strata has mean ", mean(EST_0_sel), "\n")
  cat("dpi with one strata has mean ", mean(EST_0_pre), "\n")
  return(sumtable)
}
###############################################################################################################################################
#Function: pred.power

#Author: Yumin Wang

#Creation Date: October 10, 2020 (version 3.6.1)

#Purpose: This function returns the simulated predicted power for treatment, selection, preference effects under different parameter settings
#         (phi, theta, strata, n, trt, sel, pre) and different data generating methods (homogeneous or heterogeneous effects model). 
#         This function could also be used to find n (number of people) needed to make power 80% or 90% under specific parameter settings.
#         One could use binary search to identify this n. For example, if target power is 0.8 and n=1000 makes power smaller than 0.8 and 
#         n=2000 makes power bigger than 0.8, then try (1000+2000)/2=1500 to see if power is bigger or smaller than 0.8. If it's bigger, then
#         try (1500+1000)/2=1250. If it's smaller, then try (1500+2000)/2=1750. Repeat the process until one is comfortable with the closeness 
#         to 0.8.

# Required Parameters: 
#          phi: This is the preference rate for treatment A, we assume it is constant in the entire study population. Range from 0 to 1.
#        theta: Marginal proportion of individual having a preference (marginal choice arm and random arm allocation ratio). Range from 0.1 to 0.9, step=0.05.
#       strata: Number of propensity score strata for PRPD. Range starting from 1.
#            n: Number of patients in the generated dataset
#          trt: True treatment effect of the simulation population
#          sel: True selection effect of the simulation population
#          pre: True preference effect of the simulation population
#        model: The model used to generate counterfacutals. Homo: Homogeneous effects model. Hetero: Heterogeneous effects model.     
#         seed: Specify to reproduce results


#Optional Parameters: 
#        alpha: Nominal type I error rate, default=0.05

#Other Functions: gendata()

#Output:  
#         The output displays a summary table, some basic information of the summary table:
#         Column (1) is power for treatment effect, Column (2) is power for selection effect, Column (3) is power for preference effect

#Example:
#        The following gives the summary table when phi=0.3,theta=0.3,strata=5,n=1000,trearment effect=0.4,selection effect=0.4,
#        preference effect=0.4,alpha=0.05,seed=98431 under homogeneous effects model.
#        pred.power(phi=0.3,theta=0.3,strata=5,n=1000,trt=0.4,sel=0.4,pre=0.4,model="Homo",alpha=0.05,seed=98431)

################################################################################################################################################
pred.power<-function(phi,theta,strata,n,trt,sel,pre,model,alpha=0.05,seed){
  #set seed to reproduce results
  set.seed(seed)
  #compute predicted variance, power for treatment, selection, preference effects
  #first generate a dataset with 100,000 individuals in order to get true values for stratum-specfic theta, stratum-specific sigma
  #stratum-specific phi and stratum-specific selection, preference effects
  restrue=gendata(phi=phi,theta=theta,strata=strata,n=100000,trt=trt,sel=sel,pre=pre,model)
  tabletrue=restrue$table1
  #extract stratum-specific theta
  Theta=tapply(tabletrue$Preference,tabletrue$str,mean)
  #extract stratum-specific sigma
  SigmaR=tapply(tabletrue$Y[tabletrue$Preference==0],tabletrue$str[tabletrue$Preference==0],var)
  SigmaCA=tapply(tabletrue$Y[tabletrue$Preference==1&tabletrue$Trt==1],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==1],var)
  SigmaCB=tapply(tabletrue$Y[tabletrue$Preference==1&tabletrue$Trt==0],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==0],var)
  #extract stratum-specific phi
  Phi<-tapply(tabletrue$Trt[tabletrue$Preference==1],tabletrue$str[tabletrue$Preference==1],mean)
  #extract stratum-specific selection, preference effects
  Z1<-tapply(tabletrue$Y[tabletrue$Preference==1&tabletrue$Trt==1],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==1],sum)-tapply(tabletrue$count[tabletrue$Preference==1&tabletrue$Trt==1],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==1],sum)*tapply(tabletrue$Y[tabletrue$Preference==0&tabletrue$Trt==1],tabletrue$str[tabletrue$Preference==0&tabletrue$Trt==1],mean)
  Z2<-tapply(tabletrue$Y[tabletrue$Preference==1&tabletrue$Trt==0],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==0],sum)-tapply(tabletrue$count[tabletrue$Preference==1&tabletrue$Trt==0],tabletrue$str[tabletrue$Preference==1&tabletrue$Trt==0],sum)*tapply(tabletrue$Y[tabletrue$Preference==0&tabletrue$Trt==0],tabletrue$str[tabletrue$Preference==0&tabletrue$Trt==0],mean)
  M<-tapply(tabletrue$count[tabletrue$Preference==1],tabletrue$str[tabletrue$Preference==1],sum)
  DeltaNu<-(Z1-Z2)/(2*Phi*(1-Phi)*M)
  DeltaPi<-(Z1+Z2)/(2*Phi*(1-Phi)*M)
  #calculate predicted variances for treatment, selection, preference effects based on formula
  pred.var.trt=4/(n*strata)*sum(SigmaR/(1-Theta))
  pred.var.sel=1/(4*n*strata)*sum((Phi*SigmaCA+(1-Phi)*SigmaCB+Phi*(1-Phi)*((2*Phi-1)*DeltaNu+DeltaPi)^2+2*SigmaR*(Phi^2+(1-Phi)^2)*(Theta)/(1-Theta))/(Theta*Phi^2*(1-Phi)^2))
  pred.var.pre=1/(4*n*strata)*sum((Phi*SigmaCA+(1-Phi)*SigmaCB+Phi*(1-Phi)*((2*Phi-1)*DeltaPi+DeltaNu)^2+2*SigmaR*(Phi^2+(1-Phi)^2)*(Theta)/(1-Theta))/(Theta*Phi^2*(1-Phi)^2))
  #calculate predicted power for treatment, selection, preference effects based on formula
  pred.power.trt=pnorm(abs(trt)/sqrt(pred.var.trt)-qnorm(1-alpha/2))
  pred.power.sel=pnorm(abs(sel)/sqrt(pred.var.sel)-qnorm(1-alpha/2))
  pred.power.pre=pnorm(abs(pre)/sqrt(pred.var.pre)-qnorm(1-alpha/2))
  #construct power table
  powertable<-rbind(c(pred.power.trt,pred.power.sel,pred.power.pre))
  colnames(powertable)=c("Treatment effect","Selection effect","Preference effect")
  rownames(powertable)=c("Power")
  #return the table
  return(powertable)
}
##################################################################################################################################################
#Example:
genresult(phi=0.3,theta=0.3,strata=5,n=1000,trt=0.4,sel=0.4,pre=0.4,model="Homo",nsim=1000,alpha=0.05,seed=98431)
# true trt, sel, and pre effect verified as  0.4 0.4 0.4
# dtau with one strata has mean  0.3994181
# dnu with one strata has mean  0.3577221
# dpi with one strata has mean  0.5030281
# dtau/dnu/dpi         var     power
# predicted             0.4000000 0.006406651 0.9988070
# truePS                0.4002429 0.006129680 0.9980000
# estimatedPS           0.4005905 0.006115718 0.9990000
# estimatedPS-onlyX1    0.4004650 0.006123528 0.9990000
# estimatedPS-onlyX2    0.3996770 0.006125746 1.0000000
# predicted             0.4000000 0.030953113 0.6230888
# truePS                0.3986247 0.034505889 0.5860000
# estimatedPS           0.3984036 0.034786156 0.5760000
# estimatedPS-onlyX1    0.3885694 0.033729578 0.5600000
# estimatedPS-onlyX2    0.3656822 0.033565855 0.5140000
# predicted             0.4000000 0.030955270 0.6230588
# truePS                0.4089253 0.033636118 0.6040000
# estimatedPS           0.4080864 0.034015821 0.6070000
# estimatedPS-onlyX1    0.4330290 0.033114933 0.6500000
# estimatedPS-onlyX2    0.4868206 0.032738743 0.7610000

genresult(phi=0.7,theta=0.5,strata=8,n=1000,trt=0.4,sel=0.4,pre=0.4,model="Hetero",nsim=1000,alpha=0.05,seed=12345)
# true trt, sel, and pre effect verified as  0.4 0.4 0.4 
# dtau with one strata has mean  0.3623873 
# dnu with one strata has mean  0.7477148 
# dpi with one strata has mean  0.8529129 
# dtau/dnu/dpi        var     power
# predicted             0.4000000 0.01175959 0.9580646
# truePS                0.3979720 0.01102713 0.9580000
# estimatedPS           0.3984485 0.01097185 0.9540000
# estimatedPS-onlyX1    0.3846382 0.00987605 0.9580000
# estimatedPS-onlyX2    0.3758930 0.01112791 0.9330000
# predicted             0.4000000 0.03372243 0.5863829
# truePS                0.4418738 0.04092961 0.6060000
# estimatedPS           0.4422029 0.04102358 0.6060000
# estimatedPS-onlyX1    0.4507744 0.03654262 0.6610000
# estimatedPS-onlyX2    0.7729633 0.05041585 0.9450000
# predicted             0.4000000 0.03391573 0.5839597
# truePS                0.4519342 0.04015967 0.6220000
# estimatedPS           0.4538862 0.03634995 0.6280000
# estimatedPS-onlyX1    0.3858847 0.03656867 0.5440000
# estimatedPS-onlyX2    0.9512783 0.05002090 0.9950000

