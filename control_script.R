##These system2 calls just run Rscript with the arguments passed in the second argument to system2

##Base Case, 1000 sims, no plotting
system2("Rscript","--vanilla ./no_misspec/no_misspec.R 3 FALSE")

##High Index, 1000 sims
system2("Rscript","--vanilla ./highIndexVar/highIndexVar.R 3 FALSE")

##Ultra high index, 1000 sims
system2("Rscript","--vanilla ./UltrahighIndexVar/UltrahighIndexVar.R 3 FALSE")

##Ultra low index, 1000 sims
system2("Rscript","--vanilla ./UltraLowIndexVar/UlowIndexVar.R 3 FALSE")

##Aggregate 1000 sims, first ordering, no plots
system2("Rscript","--vanilla ./aggregate/aggregate.R 3 FALSE FALSE")
##Aggregate 1000 sims, second ordering, no plots
system2("Rscript","--vanilla ./aggregate/aggregate.R 3 FALSE TRUE")

##Equal Q, 1000, ay ordering
system2("Rscript",'--vanilla ./equalQ/equalQ.R 3 FALSE "ay"')
##Equal Q, 1000, ya ordering
system2("Rscript",'--vanilla ./equalQ/equalQ.R 3 FALSE "ya"')
##Equal Q, 1000, ray ordering
system2("Rscript",'--vanilla ./equalQ/equalQ.R 3 FALSE "ray"')

##q_misspec 1000 ay ordering
system2("Rscript",'--vanilla ./q_misspec/q_misspec.R 3 FALSE "ay"')
##q_misspec 1000 ya ordering
system2("Rscript",'--vanilla ./q_misspec/q_misspec.R 3 FALSE "ya"')
##q_misspec 1000 ray ordering
system2("Rscript",'--vanilla ./q_misspec/q_misspec.R 3 FALSE "ray"')

##surv2 1000 say ordering
system2("Rscript",'--vanilla ./surv2/surv2.R 3 FALSE "say"')
##surv2 1000 asy ordering
system2("Rscript",'--vanilla ./surv2/surv2.R 3 FALSE "asy"')
##surv2 1000 sya ordering
system2("Rscript",'--vanilla ./surv2/surv2.R 3 FALSE "sya"')
##surv2 1000 yas ordering
system2("Rscript",'--vanilla ./surv2/surv2.R 3 FALSE "yas"')
##surv2 1000 dyas ordering
system2("Rscript",'--vanilla ./surv2/surv2.R 3 FALSE "dsay"')

## noX 1000 ay ordering
system2("Rscript",'--vanilla ./noX/noX.R 3 FALSE')

##outlier 1000 ay ordering
system2("Rscript",'--vanilla ./outlier1/outlier1.R 3 FALSE')
