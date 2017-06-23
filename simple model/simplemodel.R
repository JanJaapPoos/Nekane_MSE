options(width=200)
library(RDynState5NAsigmaseason6Age)
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)

setwd("~/wur/N/Projecten/MARIFLEET/Nekane/simple\ model/")
source("functions.r")


##############################################################################
# POPULATION DYNAMICS
##############################################################################

population_dynamics <- function(pop, startyear, endyear, season, natmortality, catches, recruitment){
  #pop[age,year, season, area]                                                                                                        
  for (y in (startyear:endyear)){     # need to think this, but maybe is better change last year (y in (endy:(endy +10)))           
    for (ss in (1:length(season))){
      # move time ---------------------                                                             
      if (ss ==1){
        for(age in 2:(dim(pop)[1])){
          pop[age,as.character(y),as.character(ss),] <- pop[age-1,as.character(y-1),as.character(length(season)),];
        }
      } else {
        pop[,as.character(y),as.character(ss),] <- pop[,as.character(y),as.character(ss-1),];
      }
      
      # birth/recruitment ---------------------                                                     
      if (ss ==1){
        pop[1,as.character(y),as.character(ss),] <- recruitment;
      }
      # natural mortality  ---------------------                                                             
      pop[,as.character(y),as.character(ss),] <- pop[,as.character(y),as.character(ss),]*(1-natmortality)
      pop[pop < 1e-20 ] <- 1e-20
      
      # remove catches (dims of catches here is ages,season, area, just like in main )
      pop[,as.character(y),as.character(ss),] <- pop[,as.character(y),as.character(ss),] - catches[,,as.character(ss),]
      pop[pop < 1e-20 ] <- 1e-20 
    }
  }
  return(pop)
}



##############################################################################
# YIELD CURVE
##############################################################################

yield_curve <- function(hr,wts, natmortality, R=1, sequence = seq(0.001,2,0.001), verbose=F ){
  
  # note that definition of hr is not completely correct (should be sum over seasons, and mean over ages), but as long as consistently incorect in code it should not matter
  res <- data.frame("hr"=mean(hr) *sequence,"yield"=NA) 
  iii <- 1
  for (ii in sequence){
    respop <- yld <-  matrix(0,nrow=length(ages), ncol=length(season), dimnames=list("cat"=ages,"season"=season))  
    respop[1,1] <- R
    for(aa in ages){
      
      if (aa==1){
        # respop[aa,1] <- respop[aa,1] * (1-natmortality*1)
        respop[aa,1] <- respop[aa,1] * (1-(hr[aa,1]*ii))
        yld[aa,1]    <- R            - respop[aa,1] 
        respop[aa,1] <- respop[aa,1] * (1-natmortality*1)
      }
      if (aa > 1){
        respop[aa,1] <- respop[aa-1,max(season)] * (1-natmortality*0)
        respop[aa,1] <- respop[aa,1]             * (1-(hr[aa,1] *ii))
        yld[aa,1]    <- respop[aa-1,max(season)] - respop[aa,1]
        respop[aa,1] <- respop[aa,1]             * (1-natmortality*1)
      }
      for (ss in 2:(max(season))){
        respop[aa,ss] <- respop[aa,ss-1] * (1-natmortality*0)
        respop[aa,ss] <- respop[aa,ss]   * (1-(hr[aa,ss]*ii))
        yld[aa,ss]    <- respop[aa,ss-1] - respop[aa,ss]
        respop[aa,ss] <- respop[aa,ss]   * (1-natmortality*1)
      }
    }
    res[iii, ]$yield <- sum(yld*wts)    
    iii <- iii + 1
    if (verbose == T){ 
      print("yields (in numbers)")
      print(yld)
      print(" ")
      print("population (in numbers)")
      print(respop)
    }
  }
  return(res)
}


ages          <- 1:4
season        <- 1:6
areas         <- c("a", "b")
stab.model    <- 10
NUMRUNS       <- 10
SIMNUMBER     <- 300
SIGMA         <- 1e-20
SPP1DSCSTEPS  <- SPP2DSCSTEPS <- 0
endy          <- stab.model + NUMRUNS
Linf          <- 20
K             <- 0.3
wts           <- Linf*(1-exp(-K*ages))
q             <- 0.0005
natmortality  <- 0.0001
recs          <- c(100,0) 

pop1                <- array(0, dim=c(length(ages),endy + 1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.n.dsvm      <- array(0, dim=c(length(ages),endy + 1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.wt.dsvm     <- array(0, dim=c(length(ages),endy + 1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.wt.dsvm.tot <- array(0, dim=c(1           ,endy + 1,              1,            1), dimnames=list(cat="all", year=as.character(1:(endy+1)), season="all",                option ="all"))

#run population for 15 year
pop1 <- population_dynamics(pop=pop1, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm[,1,,, drop=F], recruitment=recs)

#calculated catches can then be used for input to DSVM (has same dims as pop (1: endyr), endyr=stabmodel+numruns)
pos_catches1 <- pop1 *q*wts

#set up dsvm
sp1<- sp2 <- sp3 <- sp4 <- sp5 <-    new("DynStateInput")
catchMean(sp2)  <- catchMean(sp3) <- catchMean(sp4) <- catchMean(sp5) <- array(0.01,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))
catchSigma(sp2) <- catchSigma(sp3)<- catchSigma(sp4)<- catchSigma(sp5)<- array(0.001,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))

#this is where our loop starts, after we set up stable population
for(yy in (stab.model):(stab.model+NUMRUNS-1)){
  
  print("====== year yy ========")
  print(yy)
  print("====== POP in year yy=")
  print(pop1[,yy,,,drop=F])
  
  catchMean(sp1)  <- array(pos_catches1[,yy,,], dim=c(length(ages), length(season),length(areas)),  dimnames=list("cat"=ages,"season"= season,"option"=areas))
  # ---No way of estimating sigma, therefore we assume that is 8% of the CPUE (note slight repetion in code for dims and dimnames of 0 catch arrays for spec 3,4,5)                                                                  
  catchSigma(sp1) <- catchMean(sp1) *0.08
  
  effort <- array(c(1), dim=c(length(areas), length(season)), dimnames=list(option =areas,season=as.character(season)))
  
  print("====== catchmean input to DSVM in year yy==")
  print(catchMean(sp1))
  
  
  sp1Price <-  sp2Price <- sp3Price <- sp4Price <- sp5Price <- array(c(4), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
  #---effort and prices used (note that now c is removed (but that if other runs, then make sure to fix/remove code that removes "c" option)                                                                                         
  control     <- DynState.control(spp1LndQuota= 800,  spp2LndQuota=500, spp1LndQuotaFine= 2000, spp2LndQuotaFine= 2000,
                                  fuelUse = 0.001, fuelPrice = 2.0, landingCosts= 0,gearMaintenance= 0, addNoFishing= TRUE, increments= 30,
                                  spp1DiscardSteps= SPP1DSCSTEPS, spp2DiscardSteps= SPP2DSCSTEPS, sigma= SIGMA, simNumber= SIMNUMBER, numThreads= 20)
  
  z <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
  
  # Extract DSVM results
  dsvm_res   <- extract_dsvm_res (z, control, ages, season)
  
  #some checks
  print("====== output catches (wt) from DSVM in weight in year yy =====")
  print(catches.wt.dsvm[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop1")) #apply(z@sim@spp1Landings + z@sim@spp1Discards,c(1,3),sum))
  catches.wt.dsvm.tot[] <- apply(catches.wt.dsvm,c(2),"sum")
  
  print("====== output catches (wt) tot from DSVM in weight in year yy =")
  print(catches.wt.dsvm.tot[,yy,,])
  
  aperm(apply(pop1[,yy,,,drop=F],1:3,sum),c(1,3,2))
  
  print("====== output catches (n) from DSVM in weight in year yy ======")
  catches.n.dsvm <- catches.wt.dsvm/wts
  print(catches.n.dsvm[,yy,,])
  
  pop1 <- population_dynamics(pop=pop1, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm[,yy,,,drop=F], recruitment=recs)
  pos_catches1 <- pop1 *q*wts
  
}



#what are the weights?
wts


hr <- apply(catches.n.dsvm,1:3,sum)/    (apply(catches.n.dsvm,1:3,sum) +    apply(pop1,1:3,sum) )

pyr <- 16 

#what happens in our yield curve for this hr?
yield_curve(hr=hr[,pyr,], wts, natmortality, R=100, sequence = 1, verbose=T)

round(catches.n.dsvm[,pyr,,],2)
round(pop1[,pyr,,],2)


#next, what happens to theoretical pop for our harvest, what is harvest?
hr[,pyr,]
mean(hr[,pyr,])



yc <- yield_curve(hr=hr[,pyr,], wts, natmortality, R=100, verbose=F)

ylim <- c(0,900)
#to check
par(mfrow=c(1,2))
plot(catches.wt.dsvm.tot, type="l", ylim=ylim)

plot(x=yc$hr, y=yc$yield, ylim=ylim)
points(mean(hr[,pyr,]),yc$yield[yc$hr>mean(hr[,pyr,])][1], col="red", pch=19)
points(mean(hr[,pyr-2,]),catches.wt.dsvm.tot[,pyr-2,,], col="blue", pch=19)
points(mean(hr[,pyr-1,]),catches.wt.dsvm.tot[,pyr-1,,], col="blue", pch=19)
points(mean(hr[,pyr,]),catches.wt.dsvm.tot[,pyr,,], col="blue", pch=19)
points(mean(hr[,pyr+1,]),catches.wt.dsvm.tot[,pyr+1,,], col="blue", pch=19)
points(mean(hr[,pyr+2,]),catches.wt.dsvm.tot[,pyr+2,,], col="blue", pch=19)



