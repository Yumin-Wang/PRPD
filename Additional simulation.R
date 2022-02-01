
#Degree of positivity near violation

#load packages
library(ggplot2)
#change the parameter gamma in the PS model to control degree of positivity violation. gamma can be 1,4,6,8
gamma=8
###########################################################################################################################
#Function: gendata

#Author: Yumin Wang

#Creation Date: NOV 25, 2021 (version 3.6.1)

#Purpose: This function generates one simulation dataset and one vector which verifies the true treatment, selection ,preferece 
#         effects of the dataset and one plot for PS overlap between choice arm and random arm varying degree of positivity near violation under homogeneous effect model.
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
#       The plot verifies the overlap for PS scores between choice arm and random arm.

# Required Parameters: 
#          phi: This is the preference rate for treatment A, we assume it is constant in the entire study population. Range from 0 to 1.
#        theta: Marginal proportion of individual having a preference (marginal choice arm and random arm allocation ratio). Could only be 0.5.
#       strata: Number of propensity score strata for PRPD. Range starting from 1.
#            n: Number of patients in the generated dataset
#          trt: True treatment effect of the simulation population
#          sel: True selection effect of the simulation population
#          pre: True preference effect of the simulation population
#       

#Output:   
#       A list which contains one dataset and one vector and one plot.

#Example: The following code gives the dataset and the vector and the plot when phi=0.3,theta=0.5,strata=5,n=1000,trearment effect=0.2,selection effect=0.3,
#         preference effect=0.4 under homogeneous effects model.
#         dataset<-gendata(phi=0.3,theta=0.5,strata=5,n=1000,trt=0.2,sel=0.3,pre=0.4)
#         dataset$table1, dataset$truth, dataset$plot

#############################################################################################################################################
gendata<-function(phi,theta,strata,n,trt,sel,pre){
  Z1<-rnorm(n,mean=0,sd=1)  # generate var1 (X1)
  p_Z2=0.5 #set the probability of binominal distribution
  Z2<-rbinom(n,1,p_Z2)  # generate var2 (X2)
  
  #Find the corresponding alpha0 in the logistic model (propensity score model) based on gamma when theta=0.5
  if(gamma==1){a0=0.15}
  if(gamma==4){a0=0.61}
  if(gamma==6){a0=0.90}
  if(gamma==8){a0=1.22}
  #
  
  alpha1<-0.35*gamma #set alpha1
  alpha2<- -0.3*gamma # set alpha2
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
  
  #generate counterfactuals using homogeneous effects model 
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
  
  #generate PS plot between choice arm and random arm.
  plot<-ggplot(table1, aes(x=PS, color=Arm, fill=Arm))+
    geom_density(alpha=0.6)+scale_color_grey()+scale_fill_grey() +theme_classic()+
    labs(x="Propensity Score", y = "Density")+theme_classic()+xlim(0,1)+theme(text = element_text(size = 20))
  
  
  # return output: the dataset and the vector of true effects and the PS overlap plot
  return(list(table1=table1,truth=c(dtrt,dsel,dpre),plot=plot))
}

#################################################################################################################################################
#Function: genresult

#Author: Yumin Wang

#Creation Date: NOV 25, 2021 (version 3.6.1)

#Purpose: This function returns the simulation results for empirical estimated effects, empirical variance, 
#         empirical power (or empirical Type I error) and the corrsponding true values for treatment, selection, preference effects under
#         different parameter settings (phi, theta=0.5, strata, n, trt, sel, pre), stratification strategies (true PS, estimated PS stratification) 
#         under homogeneous effects model


# Required Parameters: 
#          phi: This is the preference rate for treatment A, we assume it is constant in the entire study population. Range from 0 to 1.
#        theta: Marginal proportion of individual having a preference (marginal choice arm and random arm allocation ratio). Could only be 0.5
#       strata: Number of propensity score strata for PRPD. Range starting from 1.
#            n: Number of patients in the generated dataset
#          trt: True treatment effect of the simulation population
#          sel: True selection effect of the simulation population
#          pre: True preference effect of the simulation population
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
#        The following gives the summary table when phi=0.3,theta=0.5,strata=5,n=1000,trearment effect=0.4,selection effect=0.4,
#        preference effect=0.4,alpha=0.05,nsim=1000,seed=98431 under homogeneous effects model.
#        genresults(phi=0.3,theta=0.5,strata=5,n=1000,trt=0.4,sel=0.4,pre=0.4,alpha=0.05,nsim=1000,seed=98431)

#################################################################################################################################################




genresult=function(phi,theta,strata,n,trt,sel,pre,alpha=0.05,nsim=1000,seed){
  #set seed to reproduce results
  set.seed(seed)
  #compute predicted variance, power for treatment, selection, preference effects
  #first generate a dataset with 100,000 individuals in order to get true values for stratum-specfic theta, stratum-specific sigma
  #stratum-specific phi and stratum-specific selection, preference effects
  restrue=gendata(phi=phi,theta=theta,strata=strata,n=100000,trt=trt,sel=sel,pre=pre)
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
    table1=gendata(phi=phi,theta=theta,strata=strata,n=n,trt=trt,sel=sel,pre=pre)$table1
    
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

#Creation Date: NOV 25, 2021 (version 3.6.1)

#Purpose: This function returns the simulated predicted power for treatment, selection, preference effects under different parameter settings
#         (phi, theta=0.5, strata, n, trt, sel, pre) homogeneous effect model 
#         This function could also be used to find n (number of people) needed to make power 80% or 90% under specific parameter settings.
#         One could use binary search to identify this n. For example, if target power is 0.8 and n=1000 makes power smaller than 0.8 and 
#         n=2000 makes power bigger than 0.8, then try (1000+2000)/2=1500 to see if power is bigger or smaller than 0.8. If it's bigger, then
#         try (1500+1000)/2=1250. If it's smaller, then try (1500+2000)/2=1750. Repeat the process until one is comfortable with the closeness 
#         to 0.8.

# Required Parameters: 
#          phi: This is the preference rate for treatment A, we assume it is constant in the entire study population. Range from 0 to 1.
#        theta: Marginal proportion of individual having a preference (marginal choice arm and random arm allocation ratio). Could only be 0.5
#       strata: Number of propensity score strata for PRPD. Range starting from 1.
#            n: Number of patients in the generated dataset
#          trt: True treatment effect of the simulation population
#          sel: True selection effect of the simulation population
#          pre: True preference effect of the simulation population
#         seed: Specify to reproduce results


#Optional Parameters: 
#        alpha: Nominal type I error rate, default=0.05

#Other Functions: gendata()

#Output:  
#         The output displays a summary table, some basic information of the summary table:
#         Column (1) is power for treatment effect, Column (2) is power for selection effect, Column (3) is power for preference effect

#Example:
#        The following gives the summary table when phi=0.3,theta=0.5,strata=5,n=1000,trearment effect=0.4,selection effect=0.4,
#        preference effect=0.4,alpha=0.05,seed=98431 under homogeneous effects model.
#        pred.power(phi=0.3,theta=0.5,strata=5,n=1000,trt=0.4,sel=0.4,pre=0.4,alpha=0.05,seed=98431)

################################################################################################################################################
pred.power<-function(phi,theta,strata,n,trt,sel,pre,alpha=0.05,seed){
  #set seed to reproduce results
  set.seed(seed)
  #compute predicted variance, power for treatment, selection, preference effects
  #first generate a dataset with 100,000 individuals in order to get true values for stratum-specfic theta, stratum-specific sigma
  #stratum-specific phi and stratum-specific selection, preference effects
  restrue=gendata(phi=phi,theta=theta,strata=strata,n=100000,trt=trt,sel=sel,pre=pre)
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

#gradual impact confounders
###########################################################################################################################
#Function: search

#Author: Yumin Wang

#Creation Date: NOV 26, 2021 (version 3.6.1)

#Purpose: This function is used to search alpha0 in the logistic model (true propensity score model)
#         which corresponds to a marginal choice arm and random arm allocation ratio: theta.
#         The function takes in alpha0 and returns marginal theta. The user is expected to vary 
#         their input alpha0 to try to get their desired marginal theta (output). The value of alpha0 
#         is used in the following function: gendata.
#
#         Further explanation: The true propensity score in the simulation population is assumed to follow a 
#         logistic model:e(X)=(1+exp[-(alpha0+alpha1*X1+alpha2*X2+alpha3*X3+alpha4*X4+alpha5*X5)])^(-1), where X1,X3,X5~N(0,1), X2,X4~Bern(0.5). 
#         To control for the marginal choice arm and random arm alllocation ratio theta and find the 
#         corresponding alpha0, first setting (alpha1,alpha2,alpha3,alpha4,alpha5)=(0.35,-0.3,0.25,-0.4,0.15). Then try to vary input 
#         alpha0 to get one's desired marginal theta (output) for the simulation population.

# Required Parameters: 
#         a0: alpha0 in the logistic model

#Output:  The function takes in alpha0 and returns marginal theta for the simulation population.

#Example: The following gives the marginal theta=0.5 which corresponds to the input alpha0=0.35
#         search(0.35)
#         We have found the corresponding alpha0 when varying theta from 0.1-0.9 with a step of 0.05
#alpha0=(-1.95,-1.47,-1.1,-0.81,-0.54,-0.31,-0.08,0.14,0.35,0.56,0.78,1.01,1.25,1.52,1.81,2.18,2.66)
#which corresponds to marginal theta=(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9)

###########################################################################################################################
search=function(a0){
  n=50000                    # set number of people to be generated to 50000
  Z1<-rnorm(n,mean=0,sd=1)   # generate X1
  Z2<-rbinom(n,1,0.5)        # generate X2
  Z3<-rnorm(n,mean=0,sd=1)   # generate X3
  Z4<-rbinom(n,1,0.5)        # generate X4
  Z5<-rnorm(n,mean=0,sd=1)   # generate X5
  W <- cbind(1, Z1, Z2, Z3, Z4, Z5)      # Design matrix 
  alpha1<-0.35               # set alpha1
  alpha2<- -0.3              # set alpha2
  alpha3<- 0.25              # set alpha3
  alpha4<- -0.4              # set alpha4
  alpha5<- 0.15              # set alpha5
  a <- c(a0, alpha1, alpha2, alpha3, alpha4, alpha5) # set vector (alpha0,alpha1,alpha2,alpha3,alpha4,alpha5)
  mean(plogis(c(W%*%a)))     # calculate marginal theta
}

#gendata follows similar description above
gendata<-function(phi,theta,strata,n,trt,sel,pre){
  Z1<-rnorm(n,mean=0,sd=1)   # generate X1
  Z2<-rbinom(n,1,0.5)        # generate X2
  Z3<-rnorm(n,mean=0,sd=1)   # generate X3
  Z4<-rbinom(n,1,0.5)        # generate X4
  Z5<-rnorm(n,mean=0,sd=1)   # generate X5
  #set alpha0 based on theta
  ALP<-c(-1.95,-1.47,-1.1,-0.81,-0.54,-0.31,-0.08,0.14,0.35,0.56,0.78,1.01,1.25,1.52,1.81,2.18,2.66)
  a0=ALP[ceiling(theta/0.05-1)]
  #set alpha1-5
  alpha1<-0.35               
  alpha2<- -0.3              
  alpha3<- 0.25
  alpha4<- -0.4
  alpha5<- 0.15
  W <- cbind(1, Z1, Z2, Z3, Z4, Z5)      # Design matrix 
  a <- c(a0, alpha1, alpha2, alpha3, alpha4, alpha5) # set vector (alpha0,alpha1,alpha2,alpha3,alpha4,alpha5)
  e<-plogis(c(W%*%a)) # generate true PS
  z <- rbinom(n, 1, e)
  
  #generate table1
  table1<-data.frame(Preference=z,Var1=Z1,Var2=Z2,Var3=Z3,Var4=Z4,Var5=Z5,PS=e)
  #define another column Arm which is the factor version of column Preference
  table1$Arm<-factor(table1$Preference,levels = c(0,1),labels = c("Random Arm","Choice Arm"))
  
  
  #calculate quantiles based on true PS and number of strata needed
  quan<-quantile(table1$PS,probs=seq(0,1,by=1/strata)) 
  quan[1]=0
  quan[strata+1]=1
  
  #estimate PS using both X1-X5, the correct specification of the propensity score model
  emodel.X1_X5=glm(Preference~Var1+Var2+Var3+Var4+Var5,family = binomial,data=table1)
  #define another column PS.X1_X5 which is estimated propensity score which used X1-X5 in estimating propensity score
  table1$PS.X1_X5=emodel.X1_X5$fitted.values
  quan.X1_X5<-quantile(table1$PS.X1_X5,probs=seq(0,1,by=1/strata))
  quan.X1_X5[1]=0
  quan.X1_X5[strata+1]=1
  
  #estimate PS using X1-X4, the incorrect specification of the propensity score model
  emodel.X1_X4=glm(Preference~Var1+Var2+Var3+Var4,family = binomial,data=table1)
  #define another column PS.X1_X4 which is estimated propensity score which X1-X4 in estimating propensity score
  table1$PS.X1_X4=emodel.X1_X4$fitted.values
  quan.X1_X4<-quantile(table1$PS.X1_X4,probs=seq(0,1,by=1/strata))
  quan.X1_X4[1]=0
  quan.X1_X4[strata+1]=1
  
  #estimate PS using X1-X3, the incorrect specification of the propensity score model
  emodel.X1_X3=glm(Preference~Var1+Var2+Var3,family = binomial,data=table1)
  #define another column PS.X1_X3 which is estimated propensity score which X1-X3 in estimating propensity score
  table1$PS.X1_X3=emodel.X1_X3$fitted.values
  quan.X1_X3<-quantile(table1$PS.X1_X3,probs=seq(0,1,by=1/strata))
  quan.X1_X3[1]=0
  quan.X1_X3[strata+1]=1
  
  #estimate PS using X1-X2, the incorrect specification of the propensity score model
  emodel.X1_X2=glm(Preference~Var1+Var2,family = binomial,data=table1)
  #define another column PS.X1_X2 which is estimated propensity score which used X1-X2 in estimating propensity score
  table1$PS.X1_X2=emodel.X1_X2$fitted.values
  quan.X1_X2<-quantile(table1$PS.X1_X2,probs=seq(0,1,by=1/strata))
  quan.X1_X2[1]=0
  quan.X1_X2[strata+1]=1
  
  #estimate PS using X1, the incorrect specification of the propensity score model
  emodel.X1=glm(Preference~Var1,family = binomial,data=table1)
  #define another column PS.X1 which is estimated propensity score which used only X1 in estimating propensity score
  table1$PS.X1=emodel.X1$fitted.values
  quan.X1<-quantile(table1$PS.X1,probs=seq(0,1,by=1/strata))
  quan.X1[1]=0
  quan.X1[strata+1]=1
  
  #estimate PS using nothing, the incorrect specification of the propensity score model
  emodel.none=glm(Preference~1,family = binomial,data=table1)
  #define another column PS.none which is estimated propensity score which used nothing in estimating propensity score
  table1$PS.none=emodel.none$fitted.values
  
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
  table1$mu11<-tao1+sigma1+pi11+0.1*table1$Var1-0.15*table1$Var2+0.1*table1$Var3-0.15*table1$Var4+0.15*table1$Var5+rnorm(n,mean=0,sd=1)
  table1$mu12<-tao1+sigma2+pi12+0.1*table1$Var1-0.15*table1$Var2+0.1*table1$Var3-0.15*table1$Var4+0.15*table1$Var5+rnorm(n,mean=0,sd=1)
  table1$mu21<-tao2+sigma1+pi21+0.1*table1$Var1-0.15*table1$Var2+0.1*table1$Var3-0.15*table1$Var4+0.15*table1$Var5+rnorm(n,mean=0,sd=1)
  table1$mu22<-tao2+sigma2+pi22+0.1*table1$Var1-0.15*table1$Var2+0.1*table1$Var3-0.15*table1$Var4+0.15*table1$Var5+rnorm(n,mean=0,sd=1)
  
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
  
  # label each patient's stratum number 
  table1$str=table1$str.X1_X5=table1$str.X1_X4=table1$str.X1_X3=table1$str.X1_X2=table1$str.X1=99
  table1$str.none=strata
  
  for(s in 1:strata){
    table1$str.X1_X5[(table1$PS.X1_X5<=quan.X1_X5[s+1]) & (table1$PS.X1_X5>quan.X1_X5[s])]=s
    table1$str.X1_X4[(table1$PS.X1_X4<=quan.X1_X4[s+1]) & (table1$PS.X1_X4>quan.X1_X4[s])]=s
    table1$str.X1_X3[(table1$PS.X1_X3<=quan.X1_X3[s+1]) & (table1$PS.X1_X3>quan.X1_X3[s])]=s
    table1$str.X1_X2[(table1$PS.X1_X2<=quan.X1_X2[s+1]) & (table1$PS.X1_X2>quan.X1_X2[s])]=s
    table1$str.X1[(table1$PS.X1<=quan.X1[s+1]) & (table1$PS.X1>quan.X1[s])]=s
    table1$str[(table1$PS<=quan[s+1]) & (table1$PS>quan[s])]=s
    
  }
  table1$count<-1
  width<-ceiling(n/strata)
  for(s in 1:(strata-1)){
    table1$str.none[(width*(s-1)+1):(width*s)]=s
  }

  
  return(list(table1=table1,truth=c(dtrt,dsel,dpre)))
}

#genresult follows similar description above, the difference is that we now have 5 confounders.
genresult=function(phi,theta,strata,n,trt,sel,pre,alpha=0.05,nsim=1000,seed){
  #set seed to reproduce results
  set.seed(seed)
  restrue=gendata(phi=phi,theta=theta,strata=strata,n=100000,trt=trt,sel=sel,pre=pre)
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
  EST_tPS_trt=EST_ePS_X1_X5_trt=EST_ePS_X1_X4_trt=EST_ePS_X1_X3_trt=EST_ePS_X1_X2_trt=EST_ePS_X1_trt=EST_ePS_none_trt=NULL   
  EST_tPS_sel=EST_ePS_X1_X5_sel=EST_ePS_X1_X4_sel=EST_ePS_X1_X3_sel=EST_ePS_X1_X2_sel=EST_ePS_X1_sel=EST_ePS_none_sel=NULL   
  EST_tPS_pre=EST_ePS_X1_X5_pre=EST_ePS_X1_X4_pre=EST_ePS_X1_X3_pre=EST_ePS_X1_X2_pre=EST_ePS_X1_pre=EST_ePS_none_pre=NULL 
  # store p-values for treatment, selection, preference effects with different stratification strategies
  PVAL_tPS_trt=PVAL_ePS_X1_X5_trt=PVAL_ePS_X1_X4_trt=PVAL_ePS_X1_X3_trt=PVAL_ePS_X1_X2_trt=PVAL_ePS_X1_trt=PVAL_ePS_none_trt=NULL 
  PVAL_tPS_sel=PVAL_ePS_X1_X5_sel=PVAL_ePS_X1_X4_sel=PVAL_ePS_X1_X3_sel=PVAL_ePS_X1_X2_sel=PVAL_ePS_X1_sel=PVAL_ePS_none_sel=NULL 
  PVAL_tPS_pre=PVAL_ePS_X1_X5_pre=PVAL_ePS_X1_X4_pre=PVAL_ePS_X1_X3_pre=PVAL_ePS_X1_X2_pre=PVAL_ePS_X1_pre=PVAL_ePS_none_pre=NULL 
  
  # carry out nsim simulations
  for(i in 1:nsim){
    #generate one table
    table1=gendata(phi=phi,theta=theta,strata=strata,n=n,trt=trt,sel=sel,pre=pre)$table1
    
    # stratum-specific effect estimators and variances for treatment effect with different stratification strategies
    est.stratum.trt=est.stratum.X1_X5.trt=est.stratum.X1_X4.trt=est.stratum.X1_X3.trt=est.stratum.X1_X2.trt=est.stratum.X1.trt=est.stratum.none.trt=rep(NA,strata)
    var.stratum.trt=var.stratum.X1_X5.trt=var.stratum.X1_X4.trt=var.stratum.X1_X3.trt=var.stratum.X1_X2.trt=var.stratum.X1.trt=var.stratum.none.trt=rep(NA,strata)
    # stratum-specific effect estimators and variances for selection effect with different stratification strategies
    est.stratum.sel=est.stratum.X1_X5.sel=est.stratum.X1_X4.sel=est.stratum.X1_X3.sel=est.stratum.X1_X2.sel=est.stratum.X1.sel=est.stratum.none.sel=rep(NA,strata)
    var.stratum.sel=var.stratum.X1_X5.sel=var.stratum.X1_X4.sel=var.stratum.X1_X3.sel=var.stratum.X1_X2.sel=var.stratum.X1.sel=var.stratum.none.sel=rep(NA,strata)
    # stratum-specific effect estimators and variances for preference effect with different stratification strategies
    est.stratum.pre=est.stratum.X1_X5.pre=est.stratum.X1_X4.pre=est.stratum.X1_X3.pre=est.stratum.X1_X2.pre=est.stratum.X1.pre=est.stratum.none.pre=rep(NA,strata)
    var.stratum.pre=var.stratum.X1_X5.pre=var.stratum.X1_X4.pre=var.stratum.X1_X3.pre=var.stratum.X1_X2.pre=var.stratum.X1.pre=var.stratum.none.pre=rep(NA,strata)
    
    #consider strata one by one
    for(s in 1:strata){
      #subset for random arm in stratum s with different stratification strategies
      temp.r=subset(table1,str==s & Arm=="Random Arm",drop=T)
      temp.X1_X5.r=subset(table1,str.X1_X5==s & Arm=="Random Arm",drop=T)
      temp.X1_X4.r=subset(table1,str.X1_X4==s & Arm=="Random Arm",drop=T)
      temp.X1_X3.r=subset(table1,str.X1_X3==s & Arm=="Random Arm",drop=T)
      temp.X1_X2.r=subset(table1,str.X1_X2==s & Arm=="Random Arm",drop=T)
      temp.X1.r=subset(table1,str.X1==s & Arm=="Random Arm",drop=T)
      temp.none.r=subset(table1,str.none==s & Arm=="Random Arm",drop=T)
      
      #subset for choice arm in stratum s with different stratification strategies
      temp.c=subset(table1,str==s & Arm=="Choice Arm",drop=T)
      temp.X1_X5.c=subset(table1,str.X1_X5==s & Arm=="Choice Arm",drop=T)
      temp.X1_X4.c=subset(table1,str.X1_X4==s & Arm=="Choice Arm",drop=T)
      temp.X1_X3.c=subset(table1,str.X1_X3==s & Arm=="Choice Arm",drop=T)
      temp.X1_X2.c=subset(table1,str.X1_X2==s & Arm=="Choice Arm",drop=T)
      temp.X1.c=subset(table1,str.X1==s & Arm=="Choice Arm",drop=T)
      temp.none.c=subset(table1,str.none==s & Arm=="Choice Arm",drop=T)
      
      #calculate phi in stratum s with different stratification strategies
      temp.phi<-mean(temp.c$Trt)
      temp.X1_X5.phi<-mean(temp.X1_X5.c$Trt)
      temp.X1_X4.phi<-mean(temp.X1_X4.c$Trt)
      temp.X1_X3.phi<-mean(temp.X1_X3.c$Trt)
      temp.X1_X2.phi<-mean(temp.X1_X2.c$Trt)
      temp.X1.phi<-mean(temp.X1.c$Trt)
      temp.none.phi<-mean(temp.none.c$Trt)
      
      #calculate d1,d2 in stratum s with different stratification strategies
      d1<-mean(temp.c$Y[temp.c$Trt==1])-mean(temp.r$Y[temp.r$Trt==1])
      d2<-mean(temp.c$Y[temp.c$Trt==0])-mean(temp.r$Y[temp.r$Trt==0])
      
      d1.X1_X5<-mean(temp.X1_X5.c$Y[temp.X1_X5.c$Trt==1])-mean(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==1])
      d2.X1_X5<-mean(temp.X1_X5.c$Y[temp.X1_X5.c$Trt==0])-mean(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==0])
      
      d1.X1_X4<-mean(temp.X1_X4.c$Y[temp.X1_X4.c$Trt==1])-mean(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==1])
      d2.X1_X4<-mean(temp.X1_X4.c$Y[temp.X1_X4.c$Trt==0])-mean(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==0])
      
      d1.X1_X3<-mean(temp.X1_X3.c$Y[temp.X1_X3.c$Trt==1])-mean(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==1])
      d2.X1_X3<-mean(temp.X1_X3.c$Y[temp.X1_X3.c$Trt==0])-mean(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==0])
      
      d1.X1_X2<-mean(temp.X1_X2.c$Y[temp.X1_X2.c$Trt==1])-mean(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==1])
      d2.X1_X2<-mean(temp.X1_X2.c$Y[temp.X1_X2.c$Trt==0])-mean(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==0])
      
      d1.X1<-mean(temp.X1.c$Y[temp.X1.c$Trt==1])-mean(temp.X1.r$Y[temp.X1.r$Trt==1])
      d2.X1<-mean(temp.X1.c$Y[temp.X1.c$Trt==0])-mean(temp.X1.r$Y[temp.X1.r$Trt==0])
      
      d1.none<-mean(temp.none.c$Y[temp.none.c$Trt==1])-mean(temp.none.r$Y[temp.none.r$Trt==1])
      d2.none<-mean(temp.none.c$Y[temp.none.c$Trt==0])-mean(temp.none.r$Y[temp.none.r$Trt==0])
      
      #calculate number of people in treatment A or B in choice or random arm in stratum s with different stratification strategies
      m1l=sum(temp.c$count[temp.c$Trt==1])
      m2l=sum(temp.c$count[temp.c$Trt==0])
      n1l=sum(temp.r$count[temp.r$Trt==1])
      n2l=sum(temp.r$count[temp.r$Trt==0])
      ml=m1l+m2l
      
      m1l.X1_X5=sum(temp.X1_X5.c$count[temp.X1_X5.c$Trt==1])
      m2l.X1_X5=sum(temp.X1_X5.c$count[temp.X1_X5.c$Trt==0])
      n1l.X1_X5=sum(temp.X1_X5.r$count[temp.X1_X5.r$Trt==1])
      n2l.X1_X5=sum(temp.X1_X5.r$count[temp.X1_X5.r$Trt==0])
      ml.X1_X5=m1l.X1_X5+m2l.X1_X5
      
      m1l.X1_X4=sum(temp.X1_X4.c$count[temp.X1_X4.c$Trt==1])
      m2l.X1_X4=sum(temp.X1_X4.c$count[temp.X1_X4.c$Trt==0])
      n1l.X1_X4=sum(temp.X1_X4.r$count[temp.X1_X4.r$Trt==1])
      n2l.X1_X4=sum(temp.X1_X4.r$count[temp.X1_X4.r$Trt==0])
      ml.X1_X4=m1l.X1_X4+m2l.X1_X4
      
      m1l.X1_X3=sum(temp.X1_X3.c$count[temp.X1_X3.c$Trt==1])
      m2l.X1_X3=sum(temp.X1_X3.c$count[temp.X1_X3.c$Trt==0])
      n1l.X1_X3=sum(temp.X1_X3.r$count[temp.X1_X3.r$Trt==1])
      n2l.X1_X3=sum(temp.X1_X3.r$count[temp.X1_X3.r$Trt==0])
      ml.X1_X3=m1l.X1_X3+m2l.X1_X3
      
      m1l.X1_X2=sum(temp.X1_X2.c$count[temp.X1_X2.c$Trt==1])
      m2l.X1_X2=sum(temp.X1_X2.c$count[temp.X1_X2.c$Trt==0])
      n1l.X1_X2=sum(temp.X1_X2.r$count[temp.X1_X2.r$Trt==1])
      n2l.X1_X2=sum(temp.X1_X2.r$count[temp.X1_X2.r$Trt==0])
      ml.X1_X2=m1l.X1_X2+m2l.X1_X2
      
      m1l.X1=sum(temp.X1.c$count[temp.X1.c$Trt==1])
      m2l.X1=sum(temp.X1.c$count[temp.X1.c$Trt==0])
      n1l.X1=sum(temp.X1.r$count[temp.X1.r$Trt==1])
      n2l.X1=sum(temp.X1.r$count[temp.X1.r$Trt==0])
      ml.X1=m1l.X1+m2l.X1
      
      m1l.none=sum(temp.none.c$count[temp.none.c$Trt==1])
      m2l.none=sum(temp.none.c$count[temp.none.c$Trt==0])
      n1l.none=sum(temp.none.r$count[temp.none.r$Trt==1])
      n2l.none=sum(temp.none.r$count[temp.none.r$Trt==0])
      ml.none=m1l.none+m2l.none
      
      #calculate Z1,Z2, var(Z1+Z2), var(Z1-Z2) in stratum s with different stratification strategies
      temp.Z1=sum(temp.c$Y[temp.c$Trt==1])-m1l*mean(temp.r$Y[temp.r$Trt==1])
      temp.Z2=sum(temp.c$Y[temp.c$Trt==0])-m2l*mean(temp.r$Y[temp.r$Trt==0])
      temp.var.Z1=m1l*var(temp.c$Y[temp.c$Trt==1])+(1+(ml-1)*m1l/ml)*m1l*var(temp.r$Y[temp.r$Trt==1])/n1l+(m1l*m2l)*(mean(temp.c$Y[temp.c$Trt==1])-mean(temp.r$Y[temp.r$Trt==1]))^2/ml
      temp.var.Z2=m2l*var(temp.c$Y[temp.c$Trt==0])+(1+(ml-1)*m2l/ml)*m2l*var(temp.r$Y[temp.r$Trt==0])/n2l+(m1l*m2l)*(mean(temp.c$Y[temp.c$Trt==0])-mean(temp.r$Y[temp.r$Trt==0]))^2/ml
      temp.covZ1Z2=-(m1l*m2l)*(mean(temp.c$Y[temp.c$Trt==1])-mean(temp.r$Y[temp.r$Trt==1]))*(mean(temp.c$Y[temp.c$Trt==0])-mean(temp.r$Y[temp.r$Trt==0]))/ml
      temp.var.Z1plusZ2=temp.var.Z1+temp.var.Z2+2*temp.covZ1Z2
      temp.var.Z1minusZ2=temp.var.Z1+temp.var.Z2-2*temp.covZ1Z2
      
      temp.X1_X5.Z1=sum(temp.X1_X5.c$Y[temp.X1_X5.c$Trt==1])-m1l.X1_X5*mean(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==1])
      temp.X1_X5.Z2=sum(temp.X1_X5.c$Y[temp.X1_X5.c$Trt==0])-m2l.X1_X5*mean(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==0])
      temp.X1_X5.var.Z1=m1l.X1_X5*var(temp.X1_X5.c$Y[temp.X1_X5.c$Trt==1])+(1+(ml.X1_X5-1)*m1l.X1_X5/ml.X1_X5)*m1l.X1_X5*var(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==1])/n1l.X1_X5+(m1l.X1_X5*m2l.X1_X5)*(mean(temp.X1_X5.c$Y[temp.X1_X5.c$Trt==1])-mean(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==1]))^2/ml.X1_X5
      temp.X1_X5.var.Z2=m2l.X1_X5*var(temp.X1_X5.c$Y[temp.X1_X5.c$Trt==0])+(1+(ml.X1_X5-1)*m2l.X1_X5/ml.X1_X5)*m2l.X1_X5*var(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==0])/n2l.X1_X5+(m1l.X1_X5*m2l.X1_X5)*(mean(temp.X1_X5.c$Y[temp.X1_X5.c$Trt==0])-mean(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==0]))^2/ml.X1_X5
      temp.X1_X5.covZ1Z2=-(m1l.X1_X5*m2l.X1_X5)*(mean(temp.X1_X5.c$Y[temp.X1_X5.c$Trt==1])-mean(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==1]))*(mean(temp.X1_X5.c$Y[temp.X1_X5.c$Trt==0])-mean(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==0]))/ml.X1_X5
      temp.X1_X5.var.Z1plusZ2=temp.X1_X5.var.Z1+temp.X1_X5.var.Z2+2*temp.X1_X5.covZ1Z2
      temp.X1_X5.var.Z1minusZ2=temp.X1_X5.var.Z1+temp.X1_X5.var.Z2-2*temp.X1_X5.covZ1Z2
      
      temp.X1_X4.Z1=sum(temp.X1_X4.c$Y[temp.X1_X4.c$Trt==1])-m1l.X1_X4*mean(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==1])
      temp.X1_X4.Z2=sum(temp.X1_X4.c$Y[temp.X1_X4.c$Trt==0])-m2l.X1_X4*mean(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==0])
      temp.X1_X4.var.Z1=m1l.X1_X4*var(temp.X1_X4.c$Y[temp.X1_X4.c$Trt==1])+(1+(ml.X1_X4-1)*m1l.X1_X4/ml.X1_X4)*m1l.X1_X4*var(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==1])/n1l.X1_X4+(m1l.X1_X4*m2l.X1_X4)*(mean(temp.X1_X4.c$Y[temp.X1_X4.c$Trt==1])-mean(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==1]))^2/ml.X1_X4
      temp.X1_X4.var.Z2=m2l.X1_X4*var(temp.X1_X4.c$Y[temp.X1_X4.c$Trt==0])+(1+(ml.X1_X4-1)*m2l.X1_X4/ml.X1_X4)*m2l.X1_X4*var(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==0])/n2l.X1_X4+(m1l.X1_X4*m2l.X1_X4)*(mean(temp.X1_X4.c$Y[temp.X1_X4.c$Trt==0])-mean(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==0]))^2/ml.X1_X4
      temp.X1_X4.covZ1Z2=-(m1l.X1_X4*m2l.X1_X4)*(mean(temp.X1_X4.c$Y[temp.X1_X4.c$Trt==1])-mean(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==1]))*(mean(temp.X1_X4.c$Y[temp.X1_X4.c$Trt==0])-mean(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==0]))/ml.X1_X4
      temp.X1_X4.var.Z1plusZ2=temp.X1_X4.var.Z1+temp.X1_X4.var.Z2+2*temp.X1_X4.covZ1Z2
      temp.X1_X4.var.Z1minusZ2=temp.X1_X4.var.Z1+temp.X1_X4.var.Z2-2*temp.X1_X4.covZ1Z2
      
      temp.X1_X3.Z1=sum(temp.X1_X3.c$Y[temp.X1_X3.c$Trt==1])-m1l.X1_X3*mean(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==1])
      temp.X1_X3.Z2=sum(temp.X1_X3.c$Y[temp.X1_X3.c$Trt==0])-m2l.X1_X3*mean(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==0])
      temp.X1_X3.var.Z1=m1l.X1_X3*var(temp.X1_X3.c$Y[temp.X1_X3.c$Trt==1])+(1+(ml.X1_X3-1)*m1l.X1_X3/ml.X1_X3)*m1l.X1_X3*var(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==1])/n1l.X1_X3+(m1l.X1_X3*m2l.X1_X3)*(mean(temp.X1_X3.c$Y[temp.X1_X3.c$Trt==1])-mean(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==1]))^2/ml.X1_X3
      temp.X1_X3.var.Z2=m2l.X1_X3*var(temp.X1_X3.c$Y[temp.X1_X3.c$Trt==0])+(1+(ml.X1_X3-1)*m2l.X1_X3/ml.X1_X3)*m2l.X1_X3*var(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==0])/n2l.X1_X3+(m1l.X1_X3*m2l.X1_X3)*(mean(temp.X1_X3.c$Y[temp.X1_X3.c$Trt==0])-mean(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==0]))^2/ml.X1_X3
      temp.X1_X3.covZ1Z2=-(m1l.X1_X3*m2l.X1_X3)*(mean(temp.X1_X3.c$Y[temp.X1_X3.c$Trt==1])-mean(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==1]))*(mean(temp.X1_X3.c$Y[temp.X1_X3.c$Trt==0])-mean(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==0]))/ml.X1_X3
      temp.X1_X3.var.Z1plusZ2=temp.X1_X3.var.Z1+temp.X1_X3.var.Z2+2*temp.X1_X3.covZ1Z2
      temp.X1_X3.var.Z1minusZ2=temp.X1_X3.var.Z1+temp.X1_X3.var.Z2-2*temp.X1_X3.covZ1Z2
      
      temp.X1_X2.Z1=sum(temp.X1_X2.c$Y[temp.X1_X2.c$Trt==1])-m1l.X1_X2*mean(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==1])
      temp.X1_X2.Z2=sum(temp.X1_X2.c$Y[temp.X1_X2.c$Trt==0])-m2l.X1_X2*mean(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==0])
      temp.X1_X2.var.Z1=m1l.X1_X2*var(temp.X1_X2.c$Y[temp.X1_X2.c$Trt==1])+(1+(ml.X1_X2-1)*m1l.X1_X2/ml.X1_X2)*m1l.X1_X2*var(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==1])/n1l.X1_X2+(m1l.X1_X2*m2l.X1_X2)*(mean(temp.X1_X2.c$Y[temp.X1_X2.c$Trt==1])-mean(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==1]))^2/ml.X1_X2
      temp.X1_X2.var.Z2=m2l.X1_X2*var(temp.X1_X2.c$Y[temp.X1_X2.c$Trt==0])+(1+(ml.X1_X2-1)*m2l.X1_X2/ml.X1_X2)*m2l.X1_X2*var(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==0])/n2l.X1_X2+(m1l.X1_X2*m2l.X1_X2)*(mean(temp.X1_X2.c$Y[temp.X1_X2.c$Trt==0])-mean(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==0]))^2/ml.X1_X2
      temp.X1_X2.covZ1Z2=-(m1l.X1_X2*m2l.X1_X2)*(mean(temp.X1_X2.c$Y[temp.X1_X2.c$Trt==1])-mean(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==1]))*(mean(temp.X1_X2.c$Y[temp.X1_X2.c$Trt==0])-mean(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==0]))/ml.X1_X2
      temp.X1_X2.var.Z1plusZ2=temp.X1_X2.var.Z1+temp.X1_X2.var.Z2+2*temp.X1_X2.covZ1Z2
      temp.X1_X2.var.Z1minusZ2=temp.X1_X2.var.Z1+temp.X1_X2.var.Z2-2*temp.X1_X2.covZ1Z2
      
      temp.X1.Z1=sum(temp.X1.c$Y[temp.X1.c$Trt==1])-m1l.X1*mean(temp.X1.r$Y[temp.X1.r$Trt==1])
      temp.X1.Z2=sum(temp.X1.c$Y[temp.X1.c$Trt==0])-m2l.X1*mean(temp.X1.r$Y[temp.X1.r$Trt==0])
      temp.X1.var.Z1=m1l.X1*var(temp.X1.c$Y[temp.X1.c$Trt==1])+(1+(ml.X1-1)*m1l.X1/ml.X1)*m1l.X1*var(temp.X1.r$Y[temp.X1.r$Trt==1])/n1l.X1+(m1l.X1*m2l.X1)*(mean(temp.X1.c$Y[temp.X1.c$Trt==1])-mean(temp.X1.r$Y[temp.X1.r$Trt==1]))^2/ml.X1
      temp.X1.var.Z2=m2l.X1*var(temp.X1.c$Y[temp.X1.c$Trt==0])+(1+(ml.X1-1)*m2l.X1/ml.X1)*m2l.X1*var(temp.X1.r$Y[temp.X1.r$Trt==0])/n2l.X1+(m1l.X1*m2l.X1)*(mean(temp.X1.c$Y[temp.X1.c$Trt==0])-mean(temp.X1.r$Y[temp.X1.r$Trt==0]))^2/ml.X1
      temp.X1.covZ1Z2=-(m1l.X1*m2l.X1)*(mean(temp.X1.c$Y[temp.X1.c$Trt==1])-mean(temp.X1.r$Y[temp.X1.r$Trt==1]))*(mean(temp.X1.c$Y[temp.X1.c$Trt==0])-mean(temp.X1.r$Y[temp.X1.r$Trt==0]))/ml.X1
      temp.X1.var.Z1plusZ2=temp.X1.var.Z1+temp.X1.var.Z2+2*temp.X1.covZ1Z2
      temp.X1.var.Z1minusZ2=temp.X1.var.Z1+temp.X1.var.Z2-2*temp.X1.covZ1Z2
      
      temp.none.Z1=sum(temp.none.c$Y[temp.none.c$Trt==1])-m1l.none*mean(temp.none.r$Y[temp.none.r$Trt==1])
      temp.none.Z2=sum(temp.none.c$Y[temp.none.c$Trt==0])-m2l.none*mean(temp.none.r$Y[temp.none.r$Trt==0])
      temp.none.var.Z1=m1l.none*var(temp.none.c$Y[temp.none.c$Trt==1])+(1+(ml.none-1)*m1l.none/ml.none)*m1l.none*var(temp.none.r$Y[temp.none.r$Trt==1])/n1l.none+(m1l.none*m2l.none)*(mean(temp.none.c$Y[temp.none.c$Trt==1])-mean(temp.none.r$Y[temp.none.r$Trt==1]))^2/ml.none
      temp.none.var.Z2=m2l.none*var(temp.none.c$Y[temp.none.c$Trt==0])+(1+(ml.none-1)*m2l.none/ml.none)*m2l.none*var(temp.none.r$Y[temp.none.r$Trt==0])/n2l.none+(m1l.none*m2l.none)*(mean(temp.none.c$Y[temp.none.c$Trt==0])-mean(temp.none.r$Y[temp.none.r$Trt==0]))^2/ml.none
      temp.none.covZ1Z2=-(m1l.none*m2l.none)*(mean(temp.none.c$Y[temp.none.c$Trt==1])-mean(temp.none.r$Y[temp.none.r$Trt==1]))*(mean(temp.none.c$Y[temp.none.c$Trt==0])-mean(temp.none.r$Y[temp.none.r$Trt==0]))/ml.none
      temp.none.var.Z1plusZ2=temp.none.var.Z1+temp.none.var.Z2+2*temp.none.covZ1Z2
      temp.none.var.Z1minusZ2=temp.none.var.Z1+temp.none.var.Z2-2*temp.none.covZ1Z2
      
      # stratum-specific effect estimators and variances for treatment effect with true PS stratification
      est.stratum.trt[s]<-mean(temp.r$Y[temp.r$Trt==1])-mean(temp.r$Y[temp.r$Trt==0])
      var.stratum.trt[s]<-var(temp.r$Y)*(1/sum(temp.r$Trt==1)+1/sum(temp.r$Trt==0))
      # stratum-specific effect estimators and variances for selection effect with true PS stratification
      est.stratum.sel[s]<-(temp.Z1-temp.Z2)/(2*temp.phi*(1-temp.phi)*ml)
      var.stratum.sel[s]<-(temp.var.Z1minusZ2+2*(m1l*d1-m2l*d2)*(2*temp.phi-1)*(d1+d2)+(2*temp.phi-1)^2*(m1l*d1-m2l*d2)^2/(ml*temp.phi*(1-temp.phi)))/(4*temp.phi^2*(1-temp.phi)^2*ml^2)
      # stratum-specific effect estimators and variances for preference effect with true PS stratification
      est.stratum.pre[s]<-(temp.Z1+temp.Z2)/(2*temp.phi*(1-temp.phi)*ml)
      var.stratum.pre[s]<-(temp.var.Z1plusZ2+2*(m1l*d1+m2l*d2)*(2*temp.phi-1)*(d1-d2)+(2*temp.phi-1)^2*(m1l*d1+m2l*d2)^2/(ml*temp.phi*(1-temp.phi)))/(4*temp.phi^2*(1-temp.phi)^2*ml^2)
      
      
      
      est.stratum.X1_X5.trt[s]<-mean(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==1])-mean(temp.X1_X5.r$Y[temp.X1_X5.r$Trt==0])
      var.stratum.X1_X5.trt[s]<-var(temp.X1_X5.r$Y)*(1/sum(temp.X1_X5.r$Trt==1)+1/sum(temp.X1_X5.r$Trt==0))
      est.stratum.X1_X5.sel[s]<-(temp.X1_X5.Z1-temp.X1_X5.Z2)/(2*temp.X1_X5.phi*(1-temp.X1_X5.phi)*ml.X1_X5)
      var.stratum.X1_X5.sel[s]<-(temp.X1_X5.var.Z1minusZ2+2*(m1l.X1_X5*d1.X1_X5-m2l.X1_X5*d2.X1_X5)*(2*temp.X1_X5.phi-1)*(d1.X1_X5+d2.X1_X5)+(2*temp.X1_X5.phi-1)^2*(m1l.X1_X5*d1.X1_X5-m2l.X1_X5*d2.X1_X5)^2/(ml.X1_X5*temp.X1_X5.phi*(1-temp.X1_X5.phi)))/(4*temp.X1_X5.phi^2*(1-temp.X1_X5.phi)^2*ml.X1_X5^2)
      est.stratum.X1_X5.pre[s]<-(temp.X1_X5.Z1+temp.X1_X5.Z2)/(2*temp.X1_X5.phi*(1-temp.X1_X5.phi)*ml.X1_X5)
      var.stratum.X1_X5.pre[s]<-(temp.X1_X5.var.Z1plusZ2+2*(m1l.X1_X5*d1.X1_X5+m2l.X1_X5*d2.X1_X5)*(2*temp.X1_X5.phi-1)*(d1.X1_X5-d2.X1_X5)+(2*temp.X1_X5.phi-1)^2*(m1l.X1_X5*d1.X1_X5+m2l.X1_X5*d2.X1_X5)^2/(ml.X1_X5*temp.X1_X5.phi*(1-temp.X1_X5.phi)))/(4*temp.X1_X5.phi^2*(1-temp.X1_X5.phi)^2*ml.X1_X5^2)
      
      est.stratum.X1_X4.trt[s]<-mean(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==1])-mean(temp.X1_X4.r$Y[temp.X1_X4.r$Trt==0])
      var.stratum.X1_X4.trt[s]<-var(temp.X1_X4.r$Y)*(1/sum(temp.X1_X4.r$Trt==1)+1/sum(temp.X1_X4.r$Trt==0))
      est.stratum.X1_X4.sel[s]<-(temp.X1_X4.Z1-temp.X1_X4.Z2)/(2*temp.X1_X4.phi*(1-temp.X1_X4.phi)*ml.X1_X4)
      var.stratum.X1_X4.sel[s]<-(temp.X1_X4.var.Z1minusZ2+2*(m1l.X1_X4*d1.X1_X4-m2l.X1_X4*d2.X1_X4)*(2*temp.X1_X4.phi-1)*(d1.X1_X4+d2.X1_X4)+(2*temp.X1_X4.phi-1)^2*(m1l.X1_X4*d1.X1_X4-m2l.X1_X4*d2.X1_X4)^2/(ml.X1_X4*temp.X1_X4.phi*(1-temp.X1_X4.phi)))/(4*temp.X1_X4.phi^2*(1-temp.X1_X4.phi)^2*ml.X1_X4^2)
      est.stratum.X1_X4.pre[s]<-(temp.X1_X4.Z1+temp.X1_X4.Z2)/(2*temp.X1_X4.phi*(1-temp.X1_X4.phi)*ml.X1_X4)
      var.stratum.X1_X4.pre[s]<-(temp.X1_X4.var.Z1plusZ2+2*(m1l.X1_X4*d1.X1_X4+m2l.X1_X4*d2.X1_X4)*(2*temp.X1_X4.phi-1)*(d1.X1_X4-d2.X1_X4)+(2*temp.X1_X4.phi-1)^2*(m1l.X1_X4*d1.X1_X4+m2l.X1_X4*d2.X1_X4)^2/(ml.X1_X4*temp.X1_X4.phi*(1-temp.X1_X4.phi)))/(4*temp.X1_X4.phi^2*(1-temp.X1_X4.phi)^2*ml.X1_X4^2)
      
      est.stratum.X1_X3.trt[s]<-mean(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==1])-mean(temp.X1_X3.r$Y[temp.X1_X3.r$Trt==0])
      var.stratum.X1_X3.trt[s]<-var(temp.X1_X3.r$Y)*(1/sum(temp.X1_X3.r$Trt==1)+1/sum(temp.X1_X3.r$Trt==0))
      est.stratum.X1_X3.sel[s]<-(temp.X1_X3.Z1-temp.X1_X3.Z2)/(2*temp.X1_X3.phi*(1-temp.X1_X3.phi)*ml.X1_X3)
      var.stratum.X1_X3.sel[s]<-(temp.X1_X3.var.Z1minusZ2+2*(m1l.X1_X3*d1.X1_X3-m2l.X1_X3*d2.X1_X3)*(2*temp.X1_X3.phi-1)*(d1.X1_X3+d2.X1_X3)+(2*temp.X1_X3.phi-1)^2*(m1l.X1_X3*d1.X1_X3-m2l.X1_X3*d2.X1_X3)^2/(ml.X1_X3*temp.X1_X3.phi*(1-temp.X1_X3.phi)))/(4*temp.X1_X3.phi^2*(1-temp.X1_X3.phi)^2*ml.X1_X3^2)
      est.stratum.X1_X3.pre[s]<-(temp.X1_X3.Z1+temp.X1_X3.Z2)/(2*temp.X1_X3.phi*(1-temp.X1_X3.phi)*ml.X1_X3)
      var.stratum.X1_X3.pre[s]<-(temp.X1_X3.var.Z1plusZ2+2*(m1l.X1_X3*d1.X1_X3+m2l.X1_X3*d2.X1_X3)*(2*temp.X1_X3.phi-1)*(d1.X1_X3-d2.X1_X3)+(2*temp.X1_X3.phi-1)^2*(m1l.X1_X3*d1.X1_X3+m2l.X1_X3*d2.X1_X3)^2/(ml.X1_X3*temp.X1_X3.phi*(1-temp.X1_X3.phi)))/(4*temp.X1_X3.phi^2*(1-temp.X1_X3.phi)^2*ml.X1_X3^2)
      
      est.stratum.X1_X2.trt[s]<-mean(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==1])-mean(temp.X1_X2.r$Y[temp.X1_X2.r$Trt==0])
      var.stratum.X1_X2.trt[s]<-var(temp.X1_X2.r$Y)*(1/sum(temp.X1_X2.r$Trt==1)+1/sum(temp.X1_X2.r$Trt==0))
      est.stratum.X1_X2.sel[s]<-(temp.X1_X2.Z1-temp.X1_X2.Z2)/(2*temp.X1_X2.phi*(1-temp.X1_X2.phi)*ml.X1_X2)
      var.stratum.X1_X2.sel[s]<-(temp.X1_X2.var.Z1minusZ2+2*(m1l.X1_X2*d1.X1_X2-m2l.X1_X2*d2.X1_X2)*(2*temp.X1_X2.phi-1)*(d1.X1_X2+d2.X1_X2)+(2*temp.X1_X2.phi-1)^2*(m1l.X1_X2*d1.X1_X2-m2l.X1_X2*d2.X1_X2)^2/(ml.X1_X2*temp.X1_X2.phi*(1-temp.X1_X2.phi)))/(4*temp.X1_X2.phi^2*(1-temp.X1_X2.phi)^2*ml.X1_X2^2)
      est.stratum.X1_X2.pre[s]<-(temp.X1_X2.Z1+temp.X1_X2.Z2)/(2*temp.X1_X2.phi*(1-temp.X1_X2.phi)*ml.X1_X2)
      var.stratum.X1_X2.pre[s]<-(temp.X1_X2.var.Z1plusZ2+2*(m1l.X1_X2*d1.X1_X2+m2l.X1_X2*d2.X1_X2)*(2*temp.X1_X2.phi-1)*(d1.X1_X2-d2.X1_X2)+(2*temp.X1_X2.phi-1)^2*(m1l.X1_X2*d1.X1_X2+m2l.X1_X2*d2.X1_X2)^2/(ml.X1_X2*temp.X1_X2.phi*(1-temp.X1_X2.phi)))/(4*temp.X1_X2.phi^2*(1-temp.X1_X2.phi)^2*ml.X1_X2^2)
      
      est.stratum.X1.trt[s]<-mean(temp.X1.r$Y[temp.X1.r$Trt==1])-mean(temp.X1.r$Y[temp.X1.r$Trt==0])
      var.stratum.X1.trt[s]<-var(temp.X1.r$Y)*(1/sum(temp.X1.r$Trt==1)+1/sum(temp.X1.r$Trt==0))
      est.stratum.X1.sel[s]<-(temp.X1.Z1-temp.X1.Z2)/(2*temp.X1.phi*(1-temp.X1.phi)*ml.X1)
      var.stratum.X1.sel[s]<-(temp.X1.var.Z1minusZ2+2*(m1l.X1*d1.X1-m2l.X1*d2.X1)*(2*temp.X1.phi-1)*(d1.X1+d2.X1)+(2*temp.X1.phi-1)^2*(m1l.X1*d1.X1-m2l.X1*d2.X1)^2/(ml.X1*temp.X1.phi*(1-temp.X1.phi)))/(4*temp.X1.phi^2*(1-temp.X1.phi)^2*ml.X1^2)
      est.stratum.X1.pre[s]<-(temp.X1.Z1+temp.X1.Z2)/(2*temp.X1.phi*(1-temp.X1.phi)*ml.X1)
      var.stratum.X1.pre[s]<-(temp.X1.var.Z1plusZ2+2*(m1l.X1*d1.X1+m2l.X1*d2.X1)*(2*temp.X1.phi-1)*(d1.X1-d2.X1)+(2*temp.X1.phi-1)^2*(m1l.X1*d1.X1+m2l.X1*d2.X1)^2/(ml.X1*temp.X1.phi*(1-temp.X1.phi)))/(4*temp.X1.phi^2*(1-temp.X1.phi)^2*ml.X1^2)
      
      est.stratum.none.trt[s]<-mean(temp.none.r$Y[temp.none.r$Trt==1])-mean(temp.none.r$Y[temp.none.r$Trt==0])
      var.stratum.none.trt[s]<-var(temp.none.r$Y)*(1/sum(temp.none.r$Trt==1)+1/sum(temp.none.r$Trt==0))
      est.stratum.none.sel[s]<-(temp.none.Z1-temp.none.Z2)/(2*temp.none.phi*(1-temp.none.phi)*ml.none)
      var.stratum.none.sel[s]<-(temp.none.var.Z1minusZ2+2*(m1l.none*d1.none-m2l.none*d2.none)*(2*temp.none.phi-1)*(d1.none+d2.none)+(2*temp.none.phi-1)^2*(m1l.none*d1.none-m2l.none*d2.none)^2/(ml.none*temp.none.phi*(1-temp.none.phi)))/(4*temp.none.phi^2*(1-temp.none.phi)^2*ml.none^2)
      est.stratum.none.pre[s]<-(temp.none.Z1+temp.none.Z2)/(2*temp.none.phi*(1-temp.none.phi)*ml.none)
      var.stratum.none.pre[s]<-(temp.none.var.Z1plusZ2+2*(m1l.none*d1.none+m2l.none*d2.none)*(2*temp.none.phi-1)*(d1.none-d2.none)+(2*temp.none.phi-1)^2*(m1l.none*d1.none+m2l.none*d2.none)^2/(ml.none*temp.none.phi*(1-temp.none.phi)))/(4*temp.none.phi^2*(1-temp.none.phi)^2*ml.none^2)
      
      
    }
    
    #store empirical estimated treatment,selection, preference effects without stratification
    EST_0_trt=c(EST_0_trt, mean(table1$Y[(table1$Arm=="Random Arm") & (table1$Trt==1)])-mean(table1$Y[(table1$Arm=="Random Arm") & (table1$Trt==0)]))
    PHI<-mean(table1$Trt[table1$Arm=="Choice Arm"])
    M1L<-sum(table1$count[(table1$Arm=="Choice Arm")&(table1$Trt==1)])
    M2L<-sum(table1$count[(table1$Arm=="Choice Arm")&(table1$Trt==0)])
    ML<-M1L+M2L
    EST_0_sel=c(EST_0_sel,(sum(table1$Y[(table1$Arm=="Choice Arm")&(table1$Trt==1)])-M1L*mean(table1$Y[(table1$Arm=="Random Arm")&(table1$Trt==1)])-(sum(table1$Y[(table1$Arm=="Choice Arm")&(table1$Trt==0)])-M2L*mean(table1$Y[(table1$Arm=="Random Arm")&(table1$Trt==0)])))/(2*(PHI)*(1-PHI)*ML))
    EST_0_pre=c(EST_0_pre,(sum(table1$Y[(table1$Arm=="Choice Arm")&(table1$Trt==1)])-M1L*mean(table1$Y[(table1$Arm=="Random Arm")&(table1$Trt==1)])+(sum(table1$Y[(table1$Arm=="Choice Arm")&(table1$Trt==0)])-M2L*mean(table1$Y[(table1$Arm=="Random Arm")&(table1$Trt==0)])))/(2*(PHI)*(1-PHI)*ML))
    
    EST_tPS_trt=c(EST_tPS_trt, mean(est.stratum.trt))
    EST_tPS_sel=c(EST_tPS_sel, mean(est.stratum.sel))
    EST_tPS_pre=c(EST_tPS_pre, mean(est.stratum.pre))
    PVAL_tPS_trt=c(PVAL_tPS_trt, 2*(1-pnorm(abs(mean(est.stratum.trt))/sqrt(sum(var.stratum.trt)/strata^2))))
    PVAL_tPS_sel=c(PVAL_tPS_sel, 2*(1-pnorm(abs(mean(est.stratum.sel))/sqrt(sum(var.stratum.sel)/strata^2))))
    PVAL_tPS_pre=c(PVAL_tPS_pre, 2*(1-pnorm(abs(mean(est.stratum.pre))/sqrt(sum(var.stratum.pre)/strata^2))))
    
    EST_ePS_X1_X5_trt=c(EST_ePS_X1_X5_trt, mean(est.stratum.X1_X5.trt))
    EST_ePS_X1_X5_sel=c(EST_ePS_X1_X5_sel, mean(est.stratum.X1_X5.sel))
    EST_ePS_X1_X5_pre=c(EST_ePS_X1_X5_pre, mean(est.stratum.X1_X5.pre))
    PVAL_ePS_X1_X5_trt=c(PVAL_ePS_X1_X5_trt, 2*(1-pnorm(abs(mean(est.stratum.X1_X5.trt))/sqrt(sum(var.stratum.X1_X5.trt)/strata^2))))
    PVAL_ePS_X1_X5_sel=c(PVAL_ePS_X1_X5_sel, 2*(1-pnorm(abs(mean(est.stratum.X1_X5.sel))/sqrt(sum(var.stratum.X1_X5.sel)/strata^2))))
    PVAL_ePS_X1_X5_pre=c(PVAL_ePS_X1_X5_pre, 2*(1-pnorm(abs(mean(est.stratum.X1_X5.pre))/sqrt(sum(var.stratum.X1_X5.pre)/strata^2))))
    
    EST_ePS_X1_X4_trt=c(EST_ePS_X1_X4_trt, mean(est.stratum.X1_X4.trt))
    EST_ePS_X1_X4_sel=c(EST_ePS_X1_X4_sel, mean(est.stratum.X1_X4.sel))
    EST_ePS_X1_X4_pre=c(EST_ePS_X1_X4_pre, mean(est.stratum.X1_X4.pre))
    PVAL_ePS_X1_X4_trt=c(PVAL_ePS_X1_X4_trt, 2*(1-pnorm(abs(mean(est.stratum.X1_X4.trt))/sqrt(sum(var.stratum.X1_X4.trt)/strata^2))))
    PVAL_ePS_X1_X4_sel=c(PVAL_ePS_X1_X4_sel, 2*(1-pnorm(abs(mean(est.stratum.X1_X4.sel))/sqrt(sum(var.stratum.X1_X4.sel)/strata^2))))
    PVAL_ePS_X1_X4_pre=c(PVAL_ePS_X1_X4_pre, 2*(1-pnorm(abs(mean(est.stratum.X1_X4.pre))/sqrt(sum(var.stratum.X1_X4.pre)/strata^2))))
    
    EST_ePS_X1_X3_trt=c(EST_ePS_X1_X3_trt, mean(est.stratum.X1_X3.trt))
    EST_ePS_X1_X3_sel=c(EST_ePS_X1_X3_sel, mean(est.stratum.X1_X3.sel))
    EST_ePS_X1_X3_pre=c(EST_ePS_X1_X3_pre, mean(est.stratum.X1_X3.pre))
    PVAL_ePS_X1_X3_trt=c(PVAL_ePS_X1_X3_trt, 2*(1-pnorm(abs(mean(est.stratum.X1_X3.trt))/sqrt(sum(var.stratum.X1_X3.trt)/strata^2))))
    PVAL_ePS_X1_X3_sel=c(PVAL_ePS_X1_X3_sel, 2*(1-pnorm(abs(mean(est.stratum.X1_X3.sel))/sqrt(sum(var.stratum.X1_X3.sel)/strata^2))))
    PVAL_ePS_X1_X3_pre=c(PVAL_ePS_X1_X3_pre, 2*(1-pnorm(abs(mean(est.stratum.X1_X3.pre))/sqrt(sum(var.stratum.X1_X3.pre)/strata^2))))
    
    EST_ePS_X1_X2_trt=c(EST_ePS_X1_X2_trt, mean(est.stratum.X1_X2.trt))
    EST_ePS_X1_X2_sel=c(EST_ePS_X1_X2_sel, mean(est.stratum.X1_X2.sel))
    EST_ePS_X1_X2_pre=c(EST_ePS_X1_X2_pre, mean(est.stratum.X1_X2.pre))
    PVAL_ePS_X1_X2_trt=c(PVAL_ePS_X1_X2_trt, 2*(1-pnorm(abs(mean(est.stratum.X1_X2.trt))/sqrt(sum(var.stratum.X1_X2.trt)/strata^2))))
    PVAL_ePS_X1_X2_sel=c(PVAL_ePS_X1_X2_sel, 2*(1-pnorm(abs(mean(est.stratum.X1_X2.sel))/sqrt(sum(var.stratum.X1_X2.sel)/strata^2))))
    PVAL_ePS_X1_X2_pre=c(PVAL_ePS_X1_X2_pre, 2*(1-pnorm(abs(mean(est.stratum.X1_X2.pre))/sqrt(sum(var.stratum.X1_X2.pre)/strata^2))))
    
    EST_ePS_X1_trt=c(EST_ePS_X1_trt, mean(est.stratum.X1.trt))
    EST_ePS_X1_sel=c(EST_ePS_X1_sel, mean(est.stratum.X1.sel))
    EST_ePS_X1_pre=c(EST_ePS_X1_pre, mean(est.stratum.X1.pre))
    PVAL_ePS_X1_trt=c(PVAL_ePS_X1_trt, 2*(1-pnorm(abs(mean(est.stratum.X1.trt))/sqrt(sum(var.stratum.X1.trt)/strata^2))))
    PVAL_ePS_X1_sel=c(PVAL_ePS_X1_sel, 2*(1-pnorm(abs(mean(est.stratum.X1.sel))/sqrt(sum(var.stratum.X1.sel)/strata^2))))
    PVAL_ePS_X1_pre=c(PVAL_ePS_X1_pre, 2*(1-pnorm(abs(mean(est.stratum.X1.pre))/sqrt(sum(var.stratum.X1.pre)/strata^2))))
    
    EST_ePS_none_trt=c(EST_ePS_none_trt, mean(est.stratum.none.trt))
    EST_ePS_none_sel=c(EST_ePS_none_sel, mean(est.stratum.none.sel))
    EST_ePS_none_pre=c(EST_ePS_none_pre, mean(est.stratum.none.pre))
    PVAL_ePS_none_trt=c(PVAL_ePS_none_trt, 2*(1-pnorm(abs(mean(est.stratum.none.trt))/sqrt(sum(var.stratum.none.trt)/strata^2))))
    PVAL_ePS_none_sel=c(PVAL_ePS_none_sel, 2*(1-pnorm(abs(mean(est.stratum.none.sel))/sqrt(sum(var.stratum.none.sel)/strata^2))))
    PVAL_ePS_none_pre=c(PVAL_ePS_none_pre, 2*(1-pnorm(abs(mean(est.stratum.none.pre))/sqrt(sum(var.stratum.none.pre)/strata^2))))
    
    # loop control
    if(i%%2000==0) print(i)
    
  }
  
  # construct summary table
  # the true effect, true variance and true power for treatment, selection, preference effects
  gold.trt=c(trt,pred.var.trt,pred.power.trt)
  gold.sel=c(sel,pred.var.sel,pred.power.sel)
  gold.pre=c(pre,pred.var.pre,pred.power.pre)
  
  EST_tPS_trt[is.infinite(EST_tPS_trt)]<-NA
  EST_tPS_sel[is.infinite(EST_tPS_sel)]<-NA
  EST_tPS_pre[is.infinite(EST_tPS_pre)]<-NA
  
  EST_ePS_X1_X5_trt[is.infinite(EST_ePS_X1_X5_trt)]<-NA
  EST_ePS_X1_X5_sel[is.infinite(EST_ePS_X1_X5_sel)]<-NA
  EST_ePS_X1_X5_pre[is.infinite(EST_ePS_X1_X5_pre)]<-NA
  
  EST_ePS_X1_X4_trt[is.infinite(EST_ePS_X1_X4_trt)]<-NA
  EST_ePS_X1_X4_sel[is.infinite(EST_ePS_X1_X4_sel)]<-NA
  EST_ePS_X1_X4_pre[is.infinite(EST_ePS_X1_X4_pre)]<-NA
  
  EST_ePS_X1_X3_trt[is.infinite(EST_ePS_X1_X3_trt)]<-NA
  EST_ePS_X1_X3_sel[is.infinite(EST_ePS_X1_X3_sel)]<-NA
  EST_ePS_X1_X3_pre[is.infinite(EST_ePS_X1_X3_pre)]<-NA
  
  EST_ePS_X1_X2_trt[is.infinite(EST_ePS_X1_X2_trt)]<-NA
  EST_ePS_X1_X2_sel[is.infinite(EST_ePS_X1_X2_sel)]<-NA
  EST_ePS_X1_X2_pre[is.infinite(EST_ePS_X1_X2_pre)]<-NA
  
  EST_ePS_X1_trt[is.infinite(EST_ePS_X1_trt)]<-NA
  EST_ePS_X1_sel[is.infinite(EST_ePS_X1_sel)]<-NA
  EST_ePS_X1_pre[is.infinite(EST_ePS_X1_pre)]<-NA
  
  EST_ePS_none_trt[is.infinite(EST_ePS_none_trt)]<-NA
  EST_ePS_none_sel[is.infinite(EST_ePS_none_sel)]<-NA
  EST_ePS_none_pre[is.infinite(EST_ePS_none_pre)]<-NA
  
  
  #results for treatment effect with different stratification strategies (empirial estimated effects, empirical variance, empirical power (alpha))
  res_tPS_trt=c(mean(EST_tPS_trt,na.rm=TRUE),var(EST_tPS_trt,na.rm=TRUE),mean(PVAL_tPS_trt<alpha,na.rm=TRUE))
  res_ePS_X1_X5_trt=c(mean(EST_ePS_X1_X5_trt,na.rm=TRUE),var(EST_ePS_X1_X5_trt,na.rm=TRUE),mean(PVAL_ePS_X1_X5_trt<alpha,na.rm=TRUE))
  res_ePS_X1_X4_trt=c(mean(EST_ePS_X1_X4_trt,na.rm=TRUE),var(EST_ePS_X1_X4_trt,na.rm=TRUE),mean(PVAL_ePS_X1_X4_trt<alpha,na.rm=TRUE))
  res_ePS_X1_X3_trt=c(mean(EST_ePS_X1_X3_trt,na.rm=TRUE),var(EST_ePS_X1_X3_trt,na.rm=TRUE),mean(PVAL_ePS_X1_X3_trt<alpha,na.rm=TRUE))
  res_ePS_X1_X2_trt=c(mean(EST_ePS_X1_X2_trt,na.rm=TRUE),var(EST_ePS_X1_X2_trt,na.rm=TRUE),mean(PVAL_ePS_X1_X2_trt<alpha,na.rm=TRUE))
  res_ePS_X1_trt=c(mean(EST_ePS_X1_trt,na.rm=TRUE),var(EST_ePS_X1_trt,na.rm=TRUE),mean(PVAL_ePS_X1_trt<alpha,na.rm=TRUE))
  res_ePS_none_trt=c(mean(EST_ePS_none_trt,na.rm=TRUE),var(EST_ePS_none_trt,na.rm=TRUE),mean(PVAL_ePS_none_trt<alpha,na.rm=TRUE))
  
  res_tPS_sel=c(mean(EST_tPS_sel,na.rm=TRUE),var(EST_tPS_sel,na.rm=TRUE),mean(PVAL_tPS_sel<alpha,na.rm=TRUE))
  res_ePS_X1_X5_sel=c(mean(EST_ePS_X1_X5_sel,na.rm=TRUE),var(EST_ePS_X1_X5_sel,na.rm=TRUE),mean(PVAL_ePS_X1_X5_sel<alpha,na.rm=TRUE))
  res_ePS_X1_X4_sel=c(mean(EST_ePS_X1_X4_sel,na.rm=TRUE),var(EST_ePS_X1_X4_sel,na.rm=TRUE),mean(PVAL_ePS_X1_X4_sel<alpha,na.rm=TRUE))
  res_ePS_X1_X3_sel=c(mean(EST_ePS_X1_X3_sel,na.rm=TRUE),var(EST_ePS_X1_X3_sel,na.rm=TRUE),mean(PVAL_ePS_X1_X3_sel<alpha,na.rm=TRUE))
  res_ePS_X1_X2_sel=c(mean(EST_ePS_X1_X2_sel,na.rm=TRUE),var(EST_ePS_X1_X2_sel,na.rm=TRUE),mean(PVAL_ePS_X1_X2_sel<alpha,na.rm=TRUE))
  res_ePS_X1_sel=c(mean(EST_ePS_X1_sel,na.rm=TRUE),var(EST_ePS_X1_sel,na.rm=TRUE),mean(PVAL_ePS_X1_sel<alpha,na.rm=TRUE))
  res_ePS_none_sel=c(mean(EST_ePS_none_sel,na.rm=TRUE),var(EST_ePS_none_sel,na.rm=TRUE),mean(PVAL_ePS_none_sel<alpha,na.rm=TRUE))
  
  res_tPS_pre=c(mean(EST_tPS_pre,na.rm=TRUE),var(EST_tPS_pre,na.rm=TRUE),mean(PVAL_tPS_pre<alpha,na.rm=TRUE))
  res_ePS_X1_X5_pre=c(mean(EST_ePS_X1_X5_pre,na.rm=TRUE),var(EST_ePS_X1_X5_pre,na.rm=TRUE),mean(PVAL_ePS_X1_X5_pre<alpha,na.rm=TRUE))
  res_ePS_X1_X4_pre=c(mean(EST_ePS_X1_X4_pre,na.rm=TRUE),var(EST_ePS_X1_X4_pre,na.rm=TRUE),mean(PVAL_ePS_X1_X4_pre<alpha,na.rm=TRUE))
  res_ePS_X1_X3_pre=c(mean(EST_ePS_X1_X3_pre,na.rm=TRUE),var(EST_ePS_X1_X3_pre,na.rm=TRUE),mean(PVAL_ePS_X1_X3_pre<alpha,na.rm=TRUE))
  res_ePS_X1_X2_pre=c(mean(EST_ePS_X1_X2_pre,na.rm=TRUE),var(EST_ePS_X1_X2_pre,na.rm=TRUE),mean(PVAL_ePS_X1_X2_pre<alpha,na.rm=TRUE))
  res_ePS_X1_pre=c(mean(EST_ePS_X1_pre,na.rm=TRUE),var(EST_ePS_X1_pre,na.rm=TRUE),mean(PVAL_ePS_X1_pre<alpha,na.rm=TRUE))
  res_ePS_none_pre=c(mean(EST_ePS_none_pre,na.rm=TRUE),var(EST_ePS_none_pre,na.rm=TRUE),mean(PVAL_ePS_none_pre<alpha,na.rm=TRUE))
  
  #summary table
  sumtable=rbind(gold.trt,res_tPS_trt,res_ePS_X1_X5_trt,res_ePS_X1_X4_trt,res_ePS_X1_X3_trt,res_ePS_X1_X2_trt,res_ePS_X1_trt,res_ePS_none_trt,gold.sel,res_tPS_sel,res_ePS_X1_X5_sel,res_ePS_X1_X4_sel,res_ePS_X1_X3_sel,res_ePS_X1_X2_sel,res_ePS_X1_sel,res_ePS_none_sel,gold.pre,res_tPS_pre,res_ePS_X1_X5_pre,res_ePS_X1_X4_pre,res_ePS_X1_X3_pre,res_ePS_X1_X2_pre,res_ePS_X1_pre,res_ePS_none_pre)
  rownames(sumtable)=c("predicted","truePS","estimatedPS-X1-X5","estimatedPS-X1-X4","estimatedPS-X1-X3","estimatedPS-X1-X2","estimatedPS-X1","estimatedPS-no-covariate","predicted","truePS","estimatedPS-X1-X5","estimatedPS-X1-X4","estimatedPS-X1-X3","estimatedPS-X1-X2","estimatedPS-X1","estimatedPS-no-covariate","predicted","truePS","estimatedPS-X1-X5","estimatedPS-X1-X4","estimatedPS-X1-X3","estimatedPS-X1-X2","estimatedPS-X1","estimatedPS-no-covariate")
  colnames(sumtable)=c("dtau/dnu/dpi","var","power")
  
  # messages
  cat("true trt, sel, and pre effect verified as ", c(restrue$truth), "\n")
  cat("dtau with one strata has mean ", mean(EST_0_trt), "\n")
  cat("dnu with one strata has mean ", mean(EST_0_sel), "\n")
  cat("dpi with one strata has mean ", mean(EST_0_pre), "\n")
  return(sumtable)
}

#pred.power follows similar description above
pred.power<-function(phi,theta,strata,n,trt,sel,pre,alpha=0.05,seed){
  
  #set seed to reproduce results
  set.seed(seed)
  #compute predicted variance, power for treatment, selection, preference effects
  #first generate a dataset with 100,000 individuals in order to get true values for stratum-specfic theta, stratum-specific sigma
  #stratum-specific phi and stratum-specific selection, preference effects
  restrue=gendata(phi=phi,theta=theta,strata=strata,n=100000,trt=trt,sel=sel,pre=pre)
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
























