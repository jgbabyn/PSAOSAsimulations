library(stockassessment)

##Fit the North Sea cod example to get an idea of parameters
fit =  sam.fit(nscodData,nscodConf,nscodParameters)

chres <- function(residuals){
    rmm = residuals
    class(rmm) = "data.frame"
    rmm
}

##write the old data
write.data.files(nscodData)

##Read them back in
cn <- read.ices("cn.dat")
cw <- read.ices("cw.dat")
dw <- read.ices("dw.dat")
lf <- read.ices("lf.dat")
lw <- read.ices("lw.dat")
mo <- read.ices("mo.dat")
nm <- read.ices("nm.dat")
pf <- read.ices("pf.dat")
pm <- read.ices("pm.dat")
sw <- read.ices("sw.dat")


make_fake_survey <- function(ages,years,time,twofirst){
    fs = matrix(1,nrow=length(years),ncol=length(ages))
    rownames(fs) = years
    colnames(fs) = ages
    attr(fs,"time") = time
    attr(fs,"twofirst") = twofirst
    fs
}

##Create our fake surveys
fakeSurveys = list()
## To simulate a changing Q over time in the first survey we create a new "survey" every year

##"Survey 1" 

##The first fixed part
fakeSurveys$sur1fixed = make_fake_survey(1:6,1983:2000,c(0.125,0.125),c(1,1))
##For the remaining 15 years a new survey made every year
for(i in 2001:2015){
    nname = paste0("sur1y",i)
    fakeSurveys[[nname]] = make_fake_survey(1:6,i,c(0.125,0.125),c(1,1))
}

##Survey 2

fakeSurveys$sur2 = make_fake_survey(1:6,1983:2015,c(0.6,0.6),c(1,1))

##Survey 3
fakeSurveys$sur3 = make_fake_survey(1:6,1983:2015,c(0.9,0.9),c(1,1))



##Create our new data setup
ndat <- setup.sam.data(surveys=fakeSurveys,
                       residual.fleets=cn,
                       prop.mature = mo,
                       stock.mean.weight = sw,
                       catch.mean.weight = cw,
                       dis.mean.weight = dw,
                       land.mean.weight = lw,
                       prop.f =pf,
                       prop.m = pm,
                       natural.mortality = nm,
                       land.frac = lf)


##Create the config file
conf <- defcon(ndat)
##Update fishing variance to be like example
## Seperate variance for age one and older
conf$keyVarF[1,2:6] = 1

##Same fbar range
conf$fbarRange = nscodConf$fbarRange
##Change the observation variance to match
## Catch has 3 variances
conf$keyVarObs[1,2] = 1
conf$keyVarObs[1,3:6] = 2
## Surveys have two, age 1 different

##Survey 1
conf$keyVarObs[2:17,1] = 3
conf$keyVarObs[2:17,2:6] = 4

##Survey 2
conf$keyVarObs[18,1] = 5
conf$keyVarObs[18,2:6] = 6
##Survey 3
conf$keyVarObs[19,1] = 7
conf$keyVarObs[19,2:6] = 8

##Setup the default parameters
pars = defpar(ndat,conf)

##Change the values of the parameters to what we are using in the sim

##Can't just copy directly these for some reason
for(i in 1:nrow(pars$logF)){
    for(j in 1:ncol(pars$logF)){
        pars$logF[i,j] = fit$pl$logF[i,j]
    }
}


for(i in 1:nrow(pars$logN)){
    for(j in 1:ncol(pars$logN)){
        pars$logN[i,j] = fit$pl$logN[i,j]
    }
}


pars$logSdLogFsta = fit$pl$logSdLogFsta
pars$logSdLogN = fit$pl$logSdLogN
pars$itrans_rho = fit$pl$itrans_rho


##Maps the SAM key to the parameters so you can see what's what
get_ptab <- function(fit,confName,parameter){
    key = fit$conf[[confName]]+1
    key[key == 0] = NA
    cf = coef(fit)
    qt <- matrix(cf[names(cf) == parameter][key],nrow=nrow(key),ncol=ncol(key))
    qt
}
get_ptab(fit,"keyVarObs","logSdLogObs")

pars$logSdLogObs[1:3] = fit$pl$logSdLogObs[1:3]
##Assume all the survey ones are just the same?
pars$logSdLogObs[c(4,6,8)] = fit$pl$logSdLogObs[4]
pars$logSdLogObs[c(5,7,9)] = fit$pl$logSdLogObs[5]

##Setup the survey Qs

qtable(fit)
conf$keyLogFpar[-c(1,18,19),1]+1
##Based off qtable(fit)

pars$logFpar[conf$keyLogFpar[c(18),]+1]
##Surveys 2 and 3
pars$logFpar[unique(conf$keyLogFpar[c(18),]+1)] = log(c(0.008,0.035,0.066,0.067,0.09))
pars$logFpar[unique(conf$keyLogFpar[c(19),]+1)] = log(c(0.008,0.035,0.066,0.067,0.09))

##First survey
pars$logFpar[unique(conf$keyLogFpar[c(2),]+1)] = log(c(0.008,0.035,0.066,0.067,0.09))
##increase Q by 5% each year
for(i in 3:17){
    pars$logFpar[unique(conf$keyLogFpar[c(i),]+1)] = log(exp(pars$logFpar[unique(conf$keyLogFpar[c(i-1),]+1)])*1.05)
}

##Create the SAM object but don't fit

sFit = sam.fit(ndat,conf,pars,run=FALSE)
##Make it be a sam class to simulate from
class(sFit) = "sam"

##Create the sims
sim = simulate(sFit,nsim=500,full.data = TRUE,set.seed=42)

##Smash survey 1 back into collected survey and put the surveys back into a list
fix_surveys <- function(dat,times,twofirsts){
    fit <- list(data = dat)
    survey1s = lapply(2:17,function(x){getFleet(fit,x)})
    survey1 = do.call(rbind,survey1s)
    attr(survey1,"time") = c(times[1],times[1])
    attr(survey1,"twofirst") = c(twofirsts[1],twofirsts[1])
    survey2 = getFleet(fit,18)
    attr(survey2,"time") = c(times[2],times[2])
    attr(survey2,"twofirst") = c(twofirsts[2],twofirsts[2])
    survey3 = getFleet(fit,19)
    attr(survey3,"time") = c(times[3],times[3])
    attr(survey3,"twofirst") = c(twofirsts[3],twofirsts[3])
    surveys = list()
    surveys$survey1 = survey1
    surveys$survey2 = survey2
    surveys$survey3 = survey3
    surveys
}



exsurv = fix_surveys(sim[[1]],c(0.125,0.6,0.9),c(1,1,1))

exdat <- setup.sam.data(surveys=exsurv,
                       residual.fleets=cn,
                       prop.mature = mo,
                       stock.mean.weight = sw,
                       catch.mean.weight = cw,
                       dis.mean.weight = dw,
                       land.mean.weight = lw,
                       prop.f =pf,
                       prop.m = pm,
                       natural.mortality = nm,
                       land.frac = lf)

##Fix all the surveys and catch
fsims = list()
for(i in 1:length(sim)){
    fitty = list(data=sim[[i]])
    cnc = getFleet(fitty,fleet=1)
    surc = fix_surveys(sim[[i]],c(0.125,0.6,0.9),c(1,1,1))
    fsims[[i]] = setup.sam.data(surveys=surc,
                                residual.fleets=cnc,
                                prop.mature=mo,
                                stock.mean.weight = sw,
                                catch.mean.weight = cw,
                                dis.mean.weight = dw,
                                land.mean.weight = lw,
                                prop.f = pf,
                                prop.m = pm,
                                natural.mortality = nm,
                                land.frac=lf)
}


##Create the config used for model fitting
exconf <- defcon(exdat)
    
##Update fishing variance to be like example
## Seperate variance for age one and older
exconf$keyVarF[1,2:6] = 1

##Same fbar range
exconf$fbarRange = nscodConf$fbarRange
##Change the observation variance to match
## Catch has 3 variances
exconf$keyVarObs[1,2] = 1
exconf$keyVarObs[1,3:6] = 2
## Surveys have two, age 1 different

##Survey 1
exconf$keyVarObs[2,1] = 3
exconf$keyVarObs[2,2:6] = 4

##Survey 2
exconf$keyVarObs[3,1] = 5
exconf$keyVarObs[3,2:6] = 6
##Survey 3
exconf$keyVarObs[4,1] = 7
exconf$keyVarObs[4,2:6] = 8

expar = defpar(exdat,exconf)

exfit = sam.fit(exdat,exconf,expar)
exresid = residuals(exfit)
exresidm = chres(exresid)

##No Misspecification sim, uses same exconf config + pars but with Q fixed 
NMSpars = pars
NMSpars$logFpar = expar$logFpar
NMSpars$logFpar[unique(exconf$keyLogFpar[c(2),]+1)] = log(c(0.008,0.035,0.066,0.067,0.09))
NMSpars$logFpar[unique(exconf$keyLogFpar[c(3),]+1)] = log(c(0.008,0.035,0.066,0.067,0.09))
NMSpars$logFpar[unique(exconf$keyLogFpar[c(4),]+1)] = log(c(0.008,0.035,0.066,0.067,0.09))

sFitNMS = sam.fit(exdat,exconf,NMSpars,run=FALSE)
##Make it be a sam class to simulate from
class(sFitNMS) = "sam"

##Create the sims
simNMS = simulate(sFitNMS,nsim=500,full.data = TRUE,set.seed=42)


##Make all the major fleet orders within a year
library(combinat)
posorders = permn(c(1,2,3,4))

norders = list()
for(i in 1:length(posorders)){
    exresid$ffleet = factor(exresid$fleet,posorders[[i]])
    ret = list()
    ret$norder = order(exresid$year,exresid$ffleet,exresid$age)
    ret$undoorder = match(1:nrow(exresidm),ret$norder)
    ret$order = paste0(posorders[[i]],collapse="")
    norders[[i]] = ret
}


saveRDS(norders,"norders.rds")

##Save the sims

saveRDS(fsims,"RWsims.rds")
saveRDS(simNMS,"RWNMSsims.rds")

rwconf = list(conf=exconf,par=expar)
saveRDS(rwconf,"rwconf.rds")


