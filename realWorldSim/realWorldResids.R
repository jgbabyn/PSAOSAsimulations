args = commandArgs(trailingOnly=TRUE)
nsim = as.numeric(args[1]) ## Number
wsim = args[2]
nordp = args[3]
confp = args[4]
ssim = args[5]

library(stockassessment)
library(TMB)

##Functions for non-standard SAM residuals
sam_pearson_resid <- function(fit){
    obsSd = exp(fit$opt$par)
    obsSd = obsSd[names(obsSd) == "logSdLogObs"]
    auxy = as.data.frame(fit$data$aux)
    sdmap = fit$conf$keyVarObs[cbind(auxy$fleet,auxy$age)]+1
    raw_resid = fit$data$logobs-fit$rep$predObs
    pearson_resid = raw_resid/obsSd[sdmap]
    ret = data.frame(year=auxy$year,fleet=auxy$fleet,age=auxy$age,raw=raw_resid,pearson=pearson_resid)
    ret
}

sam_delete_resid <- function(fit,silent=TRUE,pred_error=TRUE){
    obj = fit$obj
    obs.name = "logobs"
    keep.vector = "keep"
    pred.name = "predObs"
    auxy = as.data.frame(fit$data$aux)
 
    
    args <- as.list(obj$env)[intersect(formalArgs(TMB::MakeADFun),ls(obj$env))]
    obsY = args$data[[obs.name]]

    ts.indicators = 1:length(obsY)

    unique_ts = unique(ts.indicators)

    if(length(obj$env$random) > 0){
        args$parameters = obj$env$parList(par=obj$env$last.par.best)
    }else{
        args$parameters = obj$env$parList(obj$env$last.par.best)
    }
    parm_names = names(args$parameters)
    random_names = unique(obj$env$.random)
    fixed_names = setdiff(parm_names,random_names)
    o_map = args$map
    f_map = lapply(args$parameters[fixed_names], function(x) as.factor(x*NA))
    ##o_map[fixed_names] = f_map

    time_index = data.frame(data = args$data[obs.name],ts=ts.indicators)
    time_index$index = 1:nrow(time_index)
    sdmap = fit$conf$keyVarObs[cbind(auxy$fleet,auxy$age)]+1
   
    args$parameters[[keep.vector]] = rep(1,length(obsY))
    ##args$data[obs.name] = NULL
    ##args$map = o_map
    args$random = random_names
    args$silent = TRUE

    del_resid = numeric(nrow(time_index))
    raw_del_resid = numeric(nrow(time_index))

   
    for(i in 1:length(unique_ts)){
        if(silent == FALSE){
            print(i)
        }
        
        args2 = args
        ##Turn off the obs. and map everything else
        cur_ts = i
        args2$parameters[[keep.vector]][cur_ts] = 0
        o_map[[keep.vector]] = as.factor(rep(NA,nrow(time_index)))
        args2$map = o_map

        
        n_obj = do.call(MakeADFun,args2)
        n_opt = nlminb(n_obj$par,n_obj$fn,n_obj$gr,control=list(eval.max=1000,iter.max=1000))
        rep = n_obj$report()
        mu = rep[[pred.name]]
        obsSd = exp(n_opt$par)
        obsSd = obsSd[names(obsSd) == "logSdLogObs"]
        sd_errs = obsSd[sdmap][i]
        
        if(pred_error == TRUE){
            sdr = sdreport(n_obj)
            ssdr = summary(sdr)
            whpred =  which(names(sdr$value) == pred.name)
            diagcov = diag(sdr$cov)
            cur_cov = diagcov[cur_ts]
            sd_errs = cur_cov
        }
        raw_del_resid[cur_ts] = (obsY[cur_ts]-mu[cur_ts])
        del_resid[cur_ts] = (obsY[cur_ts]-mu[cur_ts])/sd_errs
        
    }
    ret = data.frame(year=auxy$year,fleet=auxy$fleet,age=auxy$age,delete_resid=del_resid,raw_del=raw_del_resid)
    ret

    
}

sam_adel_resid <- function(fit){
    obj = fit$obj
    opt = fit$opt
    auxy = as.data.frame(fit$data$aux)
 
    pearsonR = sam_pearson_resid(fit)
    
    obs.name = "logobs"
    keep.vector = "keep"
    pred.name = "predObs"
    
    p = length(opt$par)
    Hess = fit$opt$he
    Hess.inv = solve(Hess)

    args <- as.list(obj$env)[intersect(formalArgs(TMB::MakeADFun),ls(obj$env))]
    obsY = args$data[[obs.name]]


    if(length(obj$env$random) > 0){
        args$parameters = obj$env$parList(par=obj$env$last.par.best)
    }else{
        args$parameters = obj$env$parList(obj$env$last.par.best)
    }
    parm_names = names(args$parameters)
    random_names = unique(obj$env$.random)
    fixed_names = setdiff(parm_names,random_names)
    o_map = args$map
    args$random = random_names
    f_map = lapply(args$parameters[fixed_names], function(x) as.factor(x*NA))
    ##o_map[fixed_names] = f_map

    sdmap = fit$conf$keyVarObs[cbind(auxy$fleet,auxy$age)]+1
   
    args$parameters[[keep.vector]] = rep(1,length(obsY))

    n_obj = do.call(MakeADFun,args)
    ind = names(n_obj$par) == keep.vector
    ind1 = names(n_obj$par) != keep.vector
    Hess1 = optimHess(n_obj$par,n_obj$fn,n_obj$gr)
    Del = Hess1[ind1,ind]

    tfunc <- function(x){
        objt = obj
        indx = objt$env$lrandom()
        parx = objt$env$last.par
        parx[1:p] = x
        objt$gr(parx[!indx])
        repx = objt$report(obj$env$last.par)
        sds = obj$env$last.par
        sds = exp(sds[names(sds) == "logSdLogObs"])
        std.resid.x = (obsY-repx$predObs)/sds[sdmap]
        return(std.resid.x)
    }
    x = opt$par
    der.resid.theta = numDeriv::jacobian(tfunc,x)

    tfunc1 <- function(x){
        objt = n_obj
        indx = objt$env$lrandom()
        parx = objt$env$last.par
        parx[names(parx) == keep.vector] = x
        objt$gr(parx[!indx])
        repx = objt$report(objt$env$last.par)
           sds = obj$env$last.par
        sds = exp(sds[names(sds) == "logSdLogObs"])
        std.resid.x = (obsY-repx$predObs)/sds[sdmap]
        return(std.resid.x)
    }
    x = rep(1,length(obsY))
    der.resid.wt = numDeriv::jacobian(tfunc1,x)
    der.resid.dwti = diag(der.resid.wt)

    lslope = rep(NA,length(obsY))
    std.resid.adel = rep(NA,length(obsY))
    for(i in 1:length(obsY)){
        tdeli = t(Del[,i])
        dg = der.resid.theta[i,]
        lslope[i] = der.resid.dwti[i] - tdeli%*%Hess.inv%*%dg
        std.resid.adel[i] = pearsonR$pearson[i]-lslope[i]
    }

    ret = pearsonR
    ret$pearson = NULL
    ret$raw = NULL
    ret$aloo = std.resid.adel
    ret
    
    
}






              
##Load the possible major orders of fleets in a year (see realWorldSim.R)
norders = readRDS(nordp)

aorder = lapply(norders,function(x){x$order})

##We only want these orders 1234,2134,2314,2341,4132 to reduce time on full batch
norders = norders[c(1,5,24,17,18)]
## print("Hello!")
## print(length(norders))

##Load the sims
##MSsim = readRDS("RWsims.rds")
##NMSsim = readRDS("RWNMSsims.rds")

sim = readRDS(wsim)
simname = tools::file_path_sans_ext(basename(wsim))

##The config
rwconf = readRDS(confp)

## if(wsim == "nomiss"){
##     sim = NMSsim
## }else{
##     sim = MSsim
## }

all_resids = list()
for(i in ssim:nsim){
    print(i)

    ##Fit the sim
    fit = sam.fit(sim[[i]],rwconf$conf,rwconf$par)

    ##OSA residuals all major fleet orders
    osaAO = lapply(1:length(norders),function(x){
        resid = residuals(fit,trace=FALSE,subset=norders[[x]]$norder)
        udo = norders[[x]]$undoorder
        df = data.frame(year=resid$year[udo],fleet=resid$fleet[udo],age=resid$age[udo],residual=resid$residual[udo],order=norders[[x]]$order,type="osa")
        df
    })

    pearsonR = sam_pearson_resid(fit)
    pearsonC = data.frame(year=pearsonR$year,fleet=pearsonR$fleet,age=pearsonR$age,residual=pearsonR$pearson,order=NA,type="pearson")
    
    ##LOO residuals
    deleteR = sam_delete_resid(fit)
    deleteC = data.frame(year=deleteR$year,fleet=deleteR$fleet,age=deleteR$age,residual=deleteR$delete_resid,order=NA,type="LOO")
    
    ##ALOO residuals
    adeleteR = sam_adel_resid(fit)
    adeleteC = data.frame(year=adeleteR$year,fleet=adeleteR$fleet,age=adeleteR$age,residual=adeleteR$aloo,order=NA,type="ALOO")
    

    combind = do.call(rbind,osaAO)
    combind2 = rbind(combind,pearsonC,deleteC,adeleteC)
    combind2$sim = i
    all_resids[[i]] = combind2

}

saveRDS(all_resids,file=paste0("all_resids-",simname,"-",ssim,"-",nsim,".rds"))
