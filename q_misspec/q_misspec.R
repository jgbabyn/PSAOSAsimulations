args = commandArgs(trailingOnly=TRUE)
nseed = as.numeric(args[1]) ## Number
do_plot = args[2] ## True or False
switch_ord = args[3] ## ay, ya or ray or ray2
setwd("./q_misspec")


#install.packages("randtests")
#install.packages("lmtest")
#install.packages("copula")
#install.packages("EnvStats")
 
library(TMB)  
library(ggplot2)   
library(patchwork)  
library(xtable)  
library(randtests)  
library(lmtest)     
##library(copula)
library(EnvStats)  
library(reshape2)
library(ggExtra)  
library(ggpmisc)   
library(ggpubr)   
library(xtable)

library(numDeriv)  
library(Matrix)     
library(expm)      
library(matrixcalc)
library(orderstats)
require(grid)
library(orderstats) 
library(cowplot)

##Load our package to perform OTSA residuals
library(oneTimeStepPredict)

compile("../fit.cpp")

dyn.load(dynlib("../fit"))
#dyn.unload(dynlib("fit")) 


subdir.name =  '.'
if(switch_ord == "ya"){
    subdir.name2 = "./sort_year_age/"
}else if(switch_ord == "ay"){
    subdir.name2 = "./sort_age_year"
}else if(switch_ord == "ray"){
    subdir.name2 = "./sort_year_rev_age"
}else if(switch_ord == "ray2"){
    subdir.name2 = "./sort_year_rev_age2"
}



source("../plottingcode.R")



T = 40 
age=1:5
A=max(age)

tau = 1
sigma_u = 1
sigma_e = sigma_u/tau

beta=0 



nT = A*T

omat.meansd = matrix(NA,nseed,5) 
omat.sd = matrix(NA,nseed,5)       
omat.Npv = matrix(NA,nseed,5)       
omat.ARpv = array(NA,c(nseed,5,A))     
omat.outlier = matrix(NA,nseed,8)  

outputs = list()

for(seed in 1:nseed){
                                        #seed=1
    print(seed)
set.seed(seed) 

year.eff = rep(NA,T)
year.eff[1] = beta
RW_dev = rnorm(T,0,1)
for(t in 2:T){year.eff[t] = year.eff[t-1] + sigma_u*RW_dev[t-1]} 

age.eff = seq(-2,2,by=1) 

simdata=expand.grid(age=age,year=1:T,survey=1) 

simdata$year.eff = year.eff[simdata$year]  
simdata$age.eff = age.eff[simdata$age]
ind = simdata$year>30
simdata$age.eff[ind]=-simdata$age.eff[ind]  

simdata$mu =  simdata$year.eff + simdata$age.eff  

simdata$Y = simdata$mu + sigma_e*rnorm(nT,0)  

if(switch_ord == "ya"){
    simdata1=expand.grid(year=1:T,age=age,survey=1)
}else if(switch_ord == "ay"){
    simdata1=expand.grid(age=age,year=1:T,survey=1)
}else if(switch_ord == "ray"){
    simdata1=expand.grid(age=rev(age),year=1:T,survey=1)
}else if(switch_ord == "ray2"){
    simdata1=expand.grid(year=1:T,age=rev(age),survey=1)
}
         
#simdata1=simdata  
ind = simdata1$age + (simdata1$year-1)*A 
simdata = simdata[ind,];

nT = length(simdata$mu)

tmb.data = list(y=simdata$Y,iy=simdata$year-1,ia=simdata$age-1,is = simdata$survey-1,T=T)
tmb.data$wt = rep(1,nT)
tmb.data$adr_yhat=0

parameters = list(
  beta=beta,
  log_sd_err=log(rep(sigma_e,1)),
  log_sd_N=log(sigma_u),
  q=matrix(age.eff,1,A,byrow=T),
  Ny=rep(beta,T-1))

map = list(q = factor(matrix(c('q1','q2',NA,'q4','q5'),1,A,byrow=T))) 

obj <- MakeADFun(tmb.data, parameters, map=map, DLL="fit",random=c("Ny"),                
                 control = list(trace=0),silent = TRUE)

opt<-nlminb(obj$par,obj$fn,obj$gr,
     control = list(trace=0,eval.max=2000,iter.max=10000)) 

         sd.rep = sdreport(obj)              
    
    se_ran = sqrt(sd.rep$diag.cov.random) 
  
    repp = obj$report()

    simC = data.frame(y=simdata$Y,
                      iy= simdata$year-1,
                      ia = simdata$age-1,
                      is = simdata$survey-1,
                      T=T,
                      it = simdata$survey-1,
                      mu = simdata$mu)

      osa = oneStepPredict(obj,"y",method="fullGaussian",trace=FALSE)
    ##osaG = oneStepPredict(obj,"y","keep",method="oneStepGeneric")
    otsa = oneTimeStepPredict(obj,"y","keep",tmb.data$iy)

    library(rfhelpers)

    simC$otsa = otsa$residuals
    simC$osa = osa$residual
    simC$emu = repp$mu
    simC$resid = simC$y-simC$emu
    simC$std_resid = simC$resid/repp$sd_err[simC$it+1]

    ##For the delete resids
    simC$std_resid_del = 0
    simC$resid_del = 0

    for(i in 1:length(tmb.data$wt)){
        tmb.data$wt = rep(1,length(tmb.data$wt))
        tmb.data$wt[i] = 0
        obj_del = MakeADFun(tmb.data,obj$env$parList(),map=map, DLL="fit",random=c("Ny"),control=list(trace=0),silent=TRUE)
        opt_del = nlminb(obj_del$par,obj_del$fn,obj_del$gr,control=list(trace=0,eval.max=2000,iter.max=1000))
        rep_del = obj_del$report(obj_del$env$last.par)
        simC$resid_del[i] = tmb.data$y[i]-rep_del$mu[i]
        type = tmb.data$is[i]
        simC$std_resid_del[i] = simC$resid_del[i]/rep_del$sd_err[type+1]
    }

    ##For the approx delete residuals
    p = length(opt$par)
    Hess=optimHess(opt$par,obj$fn,obj$gr)
    Hess.inv = solve(Hess)
    
    tpar = obj$env$parList(par=obj$env$last.par.best)
    tpar$wt = tmb.data$wt
    tmb.data.nowt = tmb.data
    tmb.data.nowt$wt=NULL
    objwt <- MakeADFun(tmb.data.nowt, tpar, map=map, DLL="fit",random=c("Ny"),                
                       control = list(trace=0),silent = TRUE)
    ind = names(objwt$par)=="wt" 
    ind1 = names(objwt$par)!="wt"
    Hess1 = optimHess(objwt$par,objwt$fn,objwt$gr)
    Del = Hess1[ind1,ind]   
    
    tfunc = function(x){  
        objt = obj
        indx = objt$env$lrandom()
        parx= objt$env$last.par
        parx[1:p]=x  
        objt$gr(parx[!indx]) 
        repx = objt$report(obj$env$last.par)
        std.resid.x = (tmb.data$y - repx$mu)/repx$sd_err[tmb.data$is+1] 
        return(std.resid.x)     
    }
    x = opt$par
    ## this is same as using adgrad=TRUE ... 
    der.resid.theta = jacobian(tfunc,x)
    
    tfunc1 = function(x){
        objt = objwt 
        indx = objt$env$lrandom()
        parx= objt$env$last.par
        parx[1:length(tmb.data$wt)]=x
        objt$gr(parx[!indx]) 
        repx = objt$report(objt$env$last.par)
        std.resid.x = (tmb.data$y - repx$mu)/repx$sd_err[tmb.data$is+1] 
        return(std.resid.x)      
    }
    nT = length(tmb.data$wt)
    x = rep(1,nT)
    der.resid.wt = jacobian(tfunc1,x)
    der.resid.dwti = diag(der.resid.wt) 
    
    lslope = rep(NA,)
    simC$std_resid_adel = NA
    for(i in 1:length(tmb.data$wt)){
        tdeli = t(Del[,i])
        dg = der.resid.theta[i,]
        lslope[i] = der.resid.dwti[i] -tdeli%*%Hess.inv%*%dg
        simC$std_resid_adel[i] = simC$std_resid[i] - lslope[i]
    }
    
    simCs = split(simC,simC$it)

    ##doing the Shapiro Wilks test
    swres = apply(simC[,c("std_resid","std_resid_del","std_resid_adel","otsa","osa")],2,shapiro.test)

    ##Durbin-Watson test
    dwstr = aggregate(std_resid~ia,data=simC,function(x){dwtest(x~1)$p.value})
    dwdstr = aggregate(std_resid_del~ia,data=simC,function(x){dwtest(x~1)$p.value})[,-1]
    dwadstr = aggregate(std_resid_adel~ia,data=simC,function(x){dwtest(x~1)$p.value})[,-1]
    dwosa = aggregate(osa~ia,data=simC,function(x){dwtest(x~1)$p.value})[,-1]
    dwotsa = aggregate(otsa~ia,data=simC,function(x){dwtest(x~1)$p.value})[,-1]
    dwres = cbind(dwstr,dwdstr,dwadstr,dwosa,dwotsa)
    names(dwres) = c("ia","std_resid","std_resid_del","std_resid_adel","osa","otsa")
    dwres$age = dwres$ia +1

    dat = data.frame(year=2:T,YE=year.eff[2:T],YEhat=repp$Ny)
    dat$L_ran = repp$Ny - qnorm(0.975)*se_ran     
    dat$U_ran = repp$Ny + qnorm(0.975)*se_ran
    
    ret = list(simC=simC,simdata=simdata,swres=swres,dwres=dwres,dat=dat,opt=opt,report=repp,obj=obj)
    outputs[[seed]] = ret
    
}

if(!dir.exists("./sims/")){
    dir.create("./sims/")
}

onam = paste0("./sims/outputs",nseed,"-",switch_ord,".rds")
saveRDS(outputs,file=onam)
