
args = commandArgs(trailingOnly=TRUE)
nseed = as.numeric(args[1]) ## Number
doplot = args[2] ## True or False
switch_ord = args[3] ## True or FALSE
setwd("./aggregate")
library(TMB)
library(tidyverse)
library(oneTimeStepPredict)
library(numDeriv)
library(lmtest)
library(reshape2)
library(cowplot)
compile("../fitA.cpp")
dyn.load(dynlib("../fitA"))

doplot = TRUE
##Controls if survey is first or not

T = 40 
age=1:5
A=max(age)
##Number of seasons
S = 12

tau = 1
sigma_u = 1
sigma_w = 1
sigma_e = sigma_u/tau
sigma_ec = sigma_w/tau

beta=0 



nT = A*T

outputs = list()

print(switch_ord)

for(seed in 1:nseed){
    print(seed)
    set.seed(seed)

    year.eff = rep(NA,T)
    year.eff[1] = beta
    RW_dev = rnorm(T,0,1)
    for(t in 2:T){year.eff[t] = year.eff[t-1] + sigma_u*RW_dev[t-1]} 

    
    age.eff = seq(-2,2,by=1) 
    season.eff = seq(0,0,length.out = S)

    simdata=expand.grid(age=age,year=1:T,survey=1,type=0,season=S-1) 

    simdata$year.eff = year.eff[simdata$year]  
    simdata$age.eff = age.eff[simdata$age]
    simdata$season.eff = season.eff[simdata$season]

    simdata$mu =  simdata$year.eff + simdata$age.eff + simdata$season.eff  

    simdata$Y = simdata$mu + sigma_e*rnorm(A*T,0)
    nT = length(simdata$mu)

    simdataA = expand.grid(age=age,year=1:T,survey=2,type=1,season=1:S)
    simdataA$year.eff = year.eff[simdataA$year]  
    simdataA$age.eff = age.eff[simdataA$age]
    simdataA$season.eff = season.eff[simdataA$season]

    simdataA$season.eff[simdataA$season < 7 & simdataA$year >= 30] = -1.5
##    simdataA$season.eff[simdataA$season >= 7 & simdataA$year >= 31] = 5


    simdataA$mu =  simdataA$year.eff + simdataA$age.eff + simdataA$season.eff 

    simdataA$Y = simdataA$mu + sigma_ec*rnorm(nrow(simdataA),0)

    simdataAC = simdataA |>
        group_by(age,year,survey,type) |>
        summarise(year.eff=sum(year.eff),age.eff=sum(age.eff),season.eff=sum(season.eff),mu=sum(mu),Y=sum(Y))

    simC = data.frame(y=c(simdata$Y,simdataAC$Y),
                      iy=c(simdata$year,simdataAC$year)-1,
                      ia=c(simdata$age,simdataAC$age)-1,
                      is = c(simdata$survey,simdataAC$survey)-1,
                      T=T,
                      it=c(simdata$type,simdataAC$type),
                      mu = c(simdata$mu,simdataAC$mu))
    if(switch_ord == TRUE){
        simC = simC |>
            arrange(iy,desc(it))
    }else{
        simC = simC |>
            arrange(iy)
    }
    
    tmb.data = list(y=simC$y,
                    iy=simC$iy,
                    ia=simC$ia,
                    is = simC$is,
                    T=T,
                    it=simC$it,
                    S=S)


    tmb.data$wt = rep(1,length(tmb.data$y))
    tmb.data$adr_yhat=0

    parameters = list(
        beta=beta,
        log_sd_err=log(rep(sigma_e,2)),
        log_sd_N=log(c(sigma_u)),
        q=matrix(age.eff,2,A,byrow=T),
        Ny=rep(beta,T-1))

    map = list(q = factor(matrix(c('q11','q21',NA,'q41','q51',"q110","q210",NA,"q410","q510"),2,A,byrow=T))) 

    obj <- MakeADFun(tmb.data, parameters, map=map, DLL="fitA",random=c("Ny"),                
                     control = list(trace=0),silent = TRUE)

    opt<-nlminb(obj$par,obj$fn,obj$gr,
                control = list(trace=0,eval.max=2000,iter.max=10000)) 

      sd.rep = sdreport(obj)              
    
    se_ran = sqrt(sd.rep$diag.cov.random) 
  
    repp = obj$report()

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
        obj_del = MakeADFun(tmb.data,obj$env$parList(),map=map, DLL="fitA",random=c("Ny"),control=list(trace=0),silent=TRUE)
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
    objwt <- MakeADFun(tmb.data.nowt, tpar, map=map, DLL="fitA",random=c("Ny"),                
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
        std.resid.x = (tmb.data$y - repx$Nay)/repx$sd_err[tmb.data$is+1] 
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
        std.resid.x = (tmb.data$y - repx$Nay)/repx$sd_err[tmb.data$is+1] 
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
    
    ret = list(simC=simC,simdata=simdata,simdataAC=simdataAC,swres=swres,dwres=dwres,dat=dat,opt=opt,report=repp,obj=obj)
    outputs[[seed]] = ret
    
}

if(!dir.exists("./sims/")){
    dir.create("./sims/")
}

if(switch_ord == FALSE){
    onam = paste0("./sims/outputs",nseed,".rds")
    saveRDS(outputs,file=onam)
}else{
    onam = paste0("./sims/outputsSO",nseed,".rds")
    saveRDS(outputs,file=onam)
}

neo_resid_col_names = c("osa","otsa","std_resid","std_resid_del","std_resid_adel")


if(doplot == TRUE){
    ##Read in the functions to build the plots
    source("../plottingcode.R")
    
    for(i in 1:length(outputs)){
        dat = outputs[[i]]$dat
        simC = outputs[[i]]$simC
        simC$year = simC$iy + 1
        simC$age = simC$ia + 1
        simC$type = ifelse(simC$it == 0,"survey","aggregate")
        simCs = split(simC,simC$type)
        dwres = outputs[[i]]$dwres
        swres = outputs[[i]]$swres

        ##Stuff for plot 1
        p11 = Nay_p(dat)

        p120 = resid_plot(simCs[[1]],std_resid,std_resid,'black','Pearson')
        p210 = resid_plot(simCs[[1]],osa,std_resid,clr1,"OSA")
        p220 = resid_plot(simCs[[1]],otsa,std_resid,clr2,"OTSA") 
        p310 = resid_plot(simCs[[1]],std_resid_del,std_resid,clr3,"LOO")
        p320 = resid_plot(simCs[[1]],std_resid_adel,std_resid,clr4,"ALOO")
        p410 = resid_comp_p(simCs[[1]],std_resid,osa,otsa,std_resid_del,std_resid_adel,
                            c("Pearson","OSA","OTSA","LOO","ALOO"),c(clr1,clr2,clr3,clr4))
        sdat11 = fix_d(simCs[[1]],neo_resid_col_names)
        p420 = dens_p(sdat11[[1]],sdat11[[2]],c("black",clr1,clr2,clr3,clr4))

        p121 = resid_plot(simCs[[2]],std_resid,std_resid,'black','Pearson')
        p211 = resid_plot(simCs[[2]],osa,std_resid,clr1,"OSA")
        p221 = resid_plot(simCs[[2]],otsa,std_resid,clr2,"OTSA") 
        p311 = resid_plot(simCs[[2]],std_resid_del,std_resid,clr3,"LOO")
        p321 = resid_plot(simCs[[2]],std_resid_adel,std_resid,clr4,"ALOO")
        p411 = resid_comp_p(simCs[[2]],std_resid,osa,otsa,std_resid_del,std_resid_adel,
                            c("Pearson","OSA","OTSA","LOO","ALOO"),c(clr1,clr2,clr3,clr4))
        sdat22 = fix_d(simCs[[2]],neo_resid_col_names)
        p421 = dens_p(sdat22[[1]],sdat22[[2]],c("black",clr1,clr2,clr3,clr4))

        pall0 = plot_grid(p11,p120,p210,p220,p310,p320,p410,p420,ncol=2)
        pall1 = plot_grid(p11,p121,p211,p221,p311,p321,p411,p421,ncol=2)

        gname0 = paste0("./p10/ex_sim",i,".pdf")
        gname1 = paste0("./p11/ex_sim",i,".pdf")
        if(switch_ord == TRUE){
            gname0 = paste0("./p10/ex_simSO",i,".pdf")
            gname1 = paste0("./p11/ex_simSO",i,".pdf")
        }

        if(!dir.exists("./p10/")){
            dir.create("./p10/")
        }
        
        
        if(!dir.exists("./p11/")){
            dir.create("./p11/")
        }
        

        pdf(file=gname0,width=5,height=8)
        print(pall0)
        dev.off()

        pdf(file=gname1,width=5,height=8)
        print(pall1)
        dev.off()

        ##Stuff for plot 2
        library(tidyverse)



        pb110 =  bubble_plot(simCs[[1]],std_resid,std_resid,dwres,"Pearson")
        pb210 =  bubble_plot(simCs[[1]],osa,osa,dwres,"OSA")
        pb310 = bubble_plot(simCs[[1]],otsa,otsa,dwres,"OTSA")

        pb410 = bubble_plot(simCs[[1]],std_resid_adel,std_resid_adel,dwres,"ALOO")

        pb120 = qq_p(simCs[[1]],std_resid,"black",swres$std_resid$p.value)
        pb220 = qq_p(simCs[[1]],osa,clr1,swres$osa$p.value)
        pb320 = qq_p(simCs[[1]],otsa,clr2,swres$otsa$p.value)
        pb420 = qq_p(simCs[[1]],std_resid_adel,clr4,swres$std_resid_adel$p.value)

        
        pb111 =  bubble_plot(simCs[[2]],std_resid,std_resid,dwres,"Pearson")
        pb211 =  bubble_plot(simCs[[2]],osa,std_resid,dwres,"OSA")
        pb311 = bubble_plot(simCs[[2]],otsa,std_resid,dwres,"OTSA")
        pb411 = bubble_plot(simCs[[2]],std_resid_adel,std_resid,dwres,"ALOO")

        pb121 = qq_p(simCs[[2]],std_resid,"black",swres$std_resid$p.value)
        pb221 = qq_p(simCs[[2]],osa,clr1,swres$osa$p.value)
        pb321 = qq_p(simCs[[2]],otsa,clr2,swres$otsa$p.value)
        pb421 = qq_p(simCs[[2]],std_resid_adel,clr4,swres$std_resid_adel$p.value)

        pball0 = plot_grid(pb110,pb120,pb210,pb220,pb310,pb320,pb410,pb420,ncol=2,labels=c(unique(simCs[[1]]$type)))
        pball1 = plot_grid(pb111,pb121,pb211,pb221,pb311,pb321,pb411,pb421,ncol=2)
        
        
        gname0b = paste0("./p20/ex_sim",i,".pdf")
        gname1b = paste0("./p21/ex_sim",i,".pdf")
        if(switch_ord == TRUE){
            gname0b = paste0("./p20/ex_simSO",i,".pdf")
            gname1b = paste0("./p21/ex_simSO",i,".pdf")
   
        }

        if(!dir.exists("./p20/")){
            dir.create("./p20/")
        }
        
        
        if(!dir.exists("./p21/")){
            dir.create("./p21/")
        }
        

        pdf(file=gname0b,width=5,height=8)
        print(pball0)
        dev.off()

        pdf(file=gname1b,width=5,height=8)
        print(pball1)
        dev.off()

        
    }

}

