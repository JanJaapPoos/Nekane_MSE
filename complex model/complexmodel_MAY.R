options(width=200)
library(RDynState5NAsigmaseason6Age)
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)

setwd("~/Nekane_MSE/complex model/")
source("functions.R")


##############################################################################
# POPULATION DYNAMICS
##############################################################################

population_dynamics <- function(pop, startyear, endyear, season, natmortality, catches, recruitment, migration){
  #pop[age,year, season, area]                                                                                                        
  MigToArea <- array(0, dim=dim(migration), dimnames= dimnames(migration))
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
      
      #migration
      for (age in (1:dim(pop)[1])){
        MigToArea[age,1,as.character(ss),,] <- 0
        for (toarea in (dimnames(pop)[4][[1]])){
          for (fromarea in (dimnames(pop)[4][[1]])){
            MigToArea[age,1,as.character(ss),toarea, fromarea] <-  MigToArea[age,1,as.character(ss),toarea,fromarea] + ( pop[age,as.character(y),as.character(ss),fromarea] * migration[age,1,as.character(ss),fromarea, toarea])
          }
        }
        for (toarea in (dimnames(pop)[4][[1]])){
          for (fromarea in (dimnames(pop)[4][[1]])){
            pop[age,as.character(y),as.character(ss),toarea]  <-  pop[age,as.character(y),as.character(ss),toarea] +  MigToArea[age,1,as.character(ss),toarea, fromarea]
          }
        }
      }
    }
  }
  return(pop)
}



##############################################################################
# YIELD CURVE
##############################################################################

yield_curve <- function(hr,lratio, wts, natmortality, R=1, sequence = seq(0.001,2,0.001), verbose=F ){
  
  # note that definition of hr is not completely correct (should be sum over seasons, and mean over ages), but as long as consistently incorect in code it should not matter
  res <- data.frame("hr"=mean(hr) *sequence,"catch"=NA, "landings"=NA) 
  iii <- 1
  sumR <- sum(R)
  if (verbose == T){ 
    print("total Recruitment")
    print(R)
  }
  for (ii in sequence){
    respop <- yld <-  matrix(0,nrow=length(ages), ncol=length(season), dimnames=list("cat"=ages,"season"=season))  
    respop[1,1] <- sumR
    for(aa in ages){
      if (aa==1){
        # respop[aa,1] <- respop[aa,1] * (1-natmortality*1)
        respop[aa,1] <- respop[aa,1] * (1-(hr[aa,1]*ii))
        yld[aa,1]    <- sumR   - respop[aa,1] 
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
    res[iii, ]$catch    <- sum(yld*wts[,1,,1])
    res[iii, ]$landings <- sum(yld*lratio*wts[,1,,1])
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

#SIGMAS        <- c(10, 40, 70, 100)   #SIGMA         <- 100 #comes from 2// chanheg from 300 to 200
#NVESSELS      <- c(600, 700, 800)     #SIMNUMBER     <- 600

#for (pos in NVESSELS){
#  for (sig in SIGMAS){

ages          <- 1:4
season        <- 1:6
areas         <- c("a", "b")
stab.model    <- 10
NUMRUNS       <- 40
MPstart       <- 20
MPstartLO     <- 30
SIMNUMBER     <- 700 #pos
SIGMA         <- 20 #sig 
SPP1DSCSTEPS  <- 1
SPP2DSCSTEPS  <- 1
endy          <- stab.model + NUMRUNS
Linf          <- 50
K             <- 0.4
alpha         <- 0.00005
beta          <- 3
sages         <- array(seq(min(ages)+((1/max(season))/2), max(ages+1),1/max(season)), dim=c(length(season),1,length(ages),1), dimnames=list(season=as.character(season),   year="all", cat=ages, option ="all"))
lens          <- Linf*(1-exp(-K*(sages)))
wts           <- alpha * lens ^ beta
wts           <- aperm(wts, c(3,2,1,4))

q             <- 0.0005
natmortality  <- 0.0001

migconstant   <- 0.2
sp1price      <- sp2price      <- 150
slope1price <- 100
slope2price <- 100 # 0.50*150

# scenario I: discarding is not allowed, YPR based in C (C=L)
# scenario II: discarding is allowed, YPR based in L, hr wanted based in catches
# scenario III: discarding ocurred but not perceived, YPR based in L, hr wanted based in landings

recs1          <- c(400,0) 
mig1     <- array(0, dim=c(length(ages),1,length(season),length(areas), length(areas)), dimnames=list(cat=ages,year="all",season=as.character(season), from =areas, to=areas)) 
mig1[,,,"a","a"] <- -migconstant
mig1[,,,"b","b"] <- -migconstant
mig1[,,,"a","b"] <- migconstant
mig1[,,,"b","a"] <- migconstant
aperm( mig1,c(1,3,2,4,5))

recs2          <- c(0,400) 
mig2     <- array(0, dim=c(length(ages),1,length(season),length(areas), length(areas)), dimnames=list(cat=ages,year="all",season=as.character(season), from =areas, to=areas)) 
mig2[,,,"a","a"] <- -migconstant
mig2[,,,"b","b"] <- -migconstant
mig2[,,,"a","b"] <- migconstant
mig2[,,,"b","a"] <- migconstant
aperm( mig2,c(1,3,2,4,5))

effort <- array(c(1), dim=c(length(areas), length(season)), dimnames=list(option =areas,season=as.character(season)))


pop1  <- pop2         <-array(0, dim=c(length(ages),endy+1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.n.dsvm1       <- catches.n.dsvm2       <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
landings.n.dsvm1      <- landings.n.dsvm2      <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
catches.wt.dsvm1      <- catches.wt.dsvm2      <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
landings.wt.dsvm1     <- landings.wt.dsvm2     <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
catches.wt.dsvm.tot1  <- catches.wt.dsvm.tot2  <- array(0, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))
landings.wt.dsvm.tot1 <- landings.wt.dsvm.tot2 <- array(0, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))
quota1                <- quota2                <- array(1.6, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))

pos_catches1<- pos_catches2 <- pop1

#run population for 15 year
pop1 <- population_dynamics(pop=pop1, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs1, migration=mig1)
pop2 <- population_dynamics(pop=pop2, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs2, migration=mig2)

#calculated catches can then be used for input to DSVM (has same dims as pop (1: endyr), endyr=stabmodel+numruns)
for (ii in 1:endy){
  for(jj in areas){
    pos_catches1[,as.character(ii),,as.character(jj)] <- pop1[,as.character(ii),,as.character(jj)] *q*wts[,1,,1]
    pos_catches2[,as.character(ii),,as.character(jj)]<- pop2[,as.character(ii),,as.character(jj)] *q*wts[,1,,1]
  }
}

#set up dsvm
sp1<- sp2 <- sp3 <- sp4 <- sp5 <-    new("DynStateInput")
catchMean(sp3)  <- catchMean(sp4) <- catchMean(sp5) <- array(0.01,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))
catchSigma(sp3) <- catchSigma(sp4)<- catchSigma(sp5)<- array(0.0000001,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))

#SIZE DEPENDENT PRICING, following Zimmermann et al. (2011)
sp1Price <- array(c(sp1price + slope1price*(((wts-mean(wts))/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
sp2Price <- array(c(sp2price + slope2price*(((wts-mean(wts))/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
sp3Price <- sp4Price <- sp5Price <- array(c(0), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
#---effort and prices used (note that now c is removed (but that if other runs, then make sure to fix/remove code that removes "c" option)                                                                                         
control     <- DynState.control(spp1LndQuota= 200,  spp2LndQuota=200, spp1LndQuotaFine= 3e6, spp2LndQuotaFine= 3e6, fuelUse = 1, fuelPrice = 150.0, landingCosts= 0,gearMaintenance= 0, addNoFishing= TRUE, increments= 25, spp1DiscardSteps= SPP1DSCSTEPS, spp2DiscardSteps= SPP2DSCSTEPS, sigma= SIGMA, simNumber= SIMNUMBER, numThreads= 68)

#this is where our loop starts, after we set up stable population
for(yy in (stab.model):(stab.model+NUMRUNS)){
  
  print("====== year yy ========")
  print(yy)

  catchMean(sp1)  <- array(apply(pos_catches1[,(yy-2):yy,,,drop=F],c(1,3,4),mean), dim=c(length(ages), length(season),length(areas)),  dimnames=list("cat"=ages,"season"= season,"option"=areas))
  catchMean(sp2)  <- array(apply(pos_catches2[,(yy-2):yy,,,drop=F],c(1,3,4),mean), dim=c(length(ages), length(season),length(areas)),  dimnames=list("cat"=ages,"season"= season,"option"=areas))
  
  # ---No way of estimating sigma, therefore we assume that is 8% of the CPUE (note slight repetion in code for dims and dimnames of 0 catch arrays for spec 3,4,5)                                                                  
  catchSigma(sp1) <- catchMean(sp1) *0.08
  catchSigma(sp2) <- catchMean(sp2) *0.08
  
  # if we are in MP period, then set quota based on last year
  if (yy > MPstart){
    control@spp1LndQuota <-  quota1[,yy,,]
    control@spp2LndQuota <-  quota2[,yy,,]
  }

  # if we are in MPLO period, then set discardsteps =0 
  if (yy > MPstartLO){
    control@spp1DiscardSteps <-  0
    control@spp2DiscardSteps <-  0
  }

  #run DSVM (wiht quota constraining if in MP time)   
  z <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
  
  #extract DSVM results
  dsvm_res <-  extract_dsvm_res (z, control, ages, season)
  #Net revenuw from DSVM
  economics_res             <- as.data.frame(as.matrix(netRev(z)))
  names(economics_res )     <- "NetRev"
  economics_res$Grossrev    <- as.data.frame(as.matrix(grossRev(z)))$V1 
  economics_res$Fuelcosts   <- apply(effort(sim(z)) * control(z)@fuelUse * control(z)@fuelPrice, 2, sum)
  economics_res$Annualfine  <- as.data.frame(as.matrix(annualFine(z)))$V1
  
  if (yy == stab.model){ 
    dsvm_res_allyrs       <- cbind("year"= yy,dsvm_res)
    economics_res_allyrs  <- cbind("year"= yy,economics_res)
  } else {
    dsvm_res_allyrs <- rbind(dsvm_res_allyrs, (cbind("year"= yy,dsvm_res)))
    economics_res_allyrs <- rbind(economics_res_allyrs, (cbind("year"= yy,economics_res)))
  } 
  
  #get catches in wts from DSVM 
  catches.wt.dsvm1[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop1", catch_option="catch.wt")
  catches.wt.dsvm2[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop2", catch_option="catch.wt") 
  
  landings.wt.dsvm1[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop1", catch_option="landings.wt")
  landings.wt.dsvm2[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop2", catch_option="landings.wt") 
  
  #calculate total catches (by summing over seasons and ages)
  catches.wt.dsvm.tot1[] <- apply(catches.wt.dsvm1,c(2),"sum")
  catches.wt.dsvm.tot2[] <- apply(catches.wt.dsvm2,c(2),"sum")
  
  landings.wt.dsvm.tot1[] <- apply(landings.wt.dsvm1,c(2),"sum")
  landings.wt.dsvm.tot2[] <- apply(landings.wt.dsvm2,c(2),"sum")
  
  #calculate numbers caught from weight caught 
  for(jj in areas){
    catches.n.dsvm1 [,yy,,as.character(jj)]  <- catches.wt.dsvm1[,yy,,as.character(jj)]%/%wts[,1,,1]
    catches.n.dsvm2 [,yy,,as.character(jj)]  <- catches.wt.dsvm2[,yy,,as.character(jj)]%/%wts[,1,,1]
    
    landings.n.dsvm1 [,yy,,as.character(jj)] <- landings.wt.dsvm1[,yy,,as.character(jj)]%/%wts[,1,,1]
    landings.n.dsvm2 [,yy,,as.character(jj)] <- landings.wt.dsvm2[,yy,,as.character(jj)]%/%wts[,1,,1]
  }
 
  # calculatae what happens to population based on catches
  pop1 <- population_dynamics(pop=pop1, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,yy,,,drop=F], recruitment=recs1, migration=mig1)
  pop2 <- population_dynamics(pop=pop2, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm2[,yy,,,drop=F], recruitment=recs2, migration=mig2)
  
  #calculate the catches that can be input into DSVM based on updated pop
  for(jj in areas){
      pos_catches1[,as.character(yy),,as.character(jj)] <- pop1[,as.character(yy),,as.character(jj)] *q*wts[,1,,1]
      pos_catches2[,as.character(yy),,as.character(jj)] <- pop2[,as.character(yy),,as.character(jj)] *q*wts[,1,,1]
    }

  #------------------------------------------
  #MANAGEMENT PROCEDURE
  #------------------------------------------
  
  hr1 <- (apply(catches.n.dsvm1,1:3,sum)+1e-20)/    (apply(catches.n.dsvm1,1:3,sum) +  apply(pop1[,1:endy,,],1:3,sum))
  hr2 <- (apply(catches.n.dsvm2,1:3,sum)+1e-20)/    (apply(catches.n.dsvm2,1:3,sum) +  apply(pop2[,1:endy,,],1:3,sum))
  
  landings.ratio1 <- (apply(landings.n.dsvm1,1:3,sum)+1e-20)/ (apply(catches.n.dsvm1,1:3,sum)+1e-20)
  landings.ratio2 <- (apply(landings.n.dsvm2,1:3,sum)+1e-20)/ (apply(catches.n.dsvm2,1:3,sum)+1e-20)
  
  if (yy == MPstartLO){
    landings.ratio1[,yy,]<- 1
    landings.ratio2[,yy,]<- 1
  }
  
  yc1 <- yield_curve(hr=hr1[,yy,], landings.ratio1[,yy,], wts, natmortality, R=recs1, verbose=F)
  yc2 <- yield_curve(hr=hr2[,yy,], landings.ratio2[,yy,], wts, natmortality, R=recs2, verbose=F)
      
  hr1wanted <- yc1[yc1$landings==max(yc1$landings),]$hr
  hr2wanted <- yc2[yc2$landings==max(yc2$landings),]$hr
  
  #We take the first value of hr1wanted vector to guarantee that we always get a single value in case landingratio is 0
  if (yy > (MPstart-1) & yy < endy){ #if (yy > (MPstart-1))
    prel.quota1 <- sum((hr1wanted[1]/mean(hr1[,yy,]))* hr1[,yy,]*landings.ratio1[,yy,]*apply(pop1[,yy+1,,],c(1,2), sum)*wts[,1,,1])/SIMNUMBER 
      #sum(sweep((hr1wanted[1]/mean(hr1[,yy,]))* hr1[,yy,]*landings.ratio1[,yy,]*apply(pop1[,yy+1,,],c(1,2), sum) ,1,wts[,1,,1],"*"))/SIMNUMBER
    prel.quota2 <- sum((hr2wanted[1]/mean(hr2[,yy,]))* hr2[,yy,]*landings.ratio2[,yy,]*apply(pop2[,yy+1,,],c(1,2), sum)*wts[,1,,1])/SIMNUMBER
      #sum(sweep((hr2wanted[1]/mean(hr2[,yy,]))* hr2[,yy,]*landings.ratio2[,yy,]*apply(pop2[,yy+1,,],c(1,2), sum) ,1,wts,"*"))/SIMNUMBER
    
    # We constrained +- 15% TAC change, until the LO implementation, where we allow the HR wanted just for the transition
    # In the transition we assume landing ratio equal 1
    if (yy == MPstartLO){ # landing.ratio are equal to 1 for all years before estimating them
      quota1[,yy+1,,] <- prel.quota1
      quota2[,yy+1,,] <- prel.quota2
    }else{
      tac.constrained1 <- c(quota1[,yy,,]*0.85, quota1[,yy,,]*1.15)
      tac.constrained2 <- c(quota2[,yy,,]*0.85, quota2[,yy,,]*1.15)
      quota1[,yy+1,,] <- max(min(prel.quota1, tac.constrained1[2]), tac.constrained1[1]) 
      quota2[,yy+1,,] <- max(min(prel.quota2, tac.constrained2[2]), tac.constrained2[1]) 
    }
       
  }
}


#what are the weights?
wts

pyrnoMP   <- MPstart- 4
pyrMP     <- MPstartLO - 4
pyrMPLO   <- endy   - 4


# Scenarios II: full landings selectivity, when discarding is allowed
hr1 <- (apply(catches.n.dsvm1,1:3,sum)+1e-20)/    ((apply(catches.n.dsvm1,1:3,sum)+1e-20) +    apply(pop1[,1:endy,,],1:3,sum))
hr2 <- (apply(catches.n.dsvm2,1:3,sum)+1e-20)/    ((apply(catches.n.dsvm2,1:3,sum)+1e-20) +    apply(pop2[,1:endy,,],1:3,sum))
landings.ratio1 <- (apply(landings.n.dsvm1,1:3,sum)+1e-20)/ (apply(catches.n.dsvm1,1:3,sum)+1e-20)
landings.ratio2 <- (apply(landings.n.dsvm2,1:3,sum)+1e-20)/ (apply(catches.n.dsvm2,1:3,sum)+1e-20)

yc1noMP <- yield_curve(hr=hr1[,pyrnoMP,], landings.ratio1[,pyrnoMP,], wts, natmortality, R=recs1, verbose=F)
yc2noMP <- yield_curve(hr=hr2[,pyrnoMP,], landings.ratio2[,pyrnoMP,], wts, natmortality, R=recs2, verbose=F)
#what happens in our yield curve for this hr?
#yield_curve(hr=hr1[,pyrnoMP,], landings.ratio1[,pyrnoMP,], wts, natmortality, R=recs1, sequence = 1, verbose=T)
#yield_curve(hr=hr2[,pyrnoMP,], landings.ratio2[,pyrnoMP,], wts, natmortality, R=recs2, sequence = 1, verbose=T)

yc1MP <- yield_curve(hr=hr1[,pyrMP,], landings.ratio1[,pyrMP,], wts, natmortality, R=recs1, verbose=F)
yc2MP <- yield_curve(hr=hr2[,pyrMP,], landings.ratio2[,pyrMP,], wts, natmortality, R=recs2, verbose=F)

yc1MPLO <- yield_curve(hr=hr1[,pyrMPLO,], landings.ratio1[,pyrMPLO,], wts, natmortality, R=recs1, verbose=F)
yc2MPLO <- yield_curve(hr=hr2[,pyrMPLO,], landings.ratio2[,pyrMPLO,], wts, natmortality, R=recs2, verbose=F)
 
Fmsy1noMP <- yc1noMP[yc1noMP$landings==max(yc1noMP$landings),]$hr
Fmsy2noMP <- yc2noMP[yc2noMP$landings==max(yc2noMP$landings),]$hr

Fmsy1MP <- yc1MP[yc1MP$landings==max(yc1MP$landings),]$hr
Fmsy2MP <- yc2MP[yc2MP$landings==max(yc2MP$landings),]$hr

Fmsy1MPLO <- yc1MPLO[yc1MPLO$landings==max(yc1MPLO$landings),]$hr
Fmsy2MPLO <- yc2MPLO[yc2MPLO$landings==max(yc2MPLO$landings),]$hr

#round(catches.n.dsvm1[,pyrnoMP,,],2)
#round(pop1[,pyrnoMP,,],2)

#next, what happens to theoretical pop for our harvest, what is harvest?
#hr1[,pyrnoMP,]
#mean(hr1[,pyrnoMP,])

# Effort pattern and economics
#--------------------------------------------------------------------------------------
#effort_plot_dsvm(SIMNUMBER, dsvm_res_allyrs, stab.model, economics_res_allyrs)

#effort
effort      <- subset(dsvm_res_allyrs,(spp %in% "sp1"))
effort      <- subset(effort,(cat %in% 1))
effort$year <- factor(effort$year, levels= 1:endy)
effort      <- aggregate(cbind(landings.wt, discards.wt, catch.wt, effort,trip)~ spp+cat+option+year, FUN=sum, data=effort)

trip      <- with(effort,data.frame(year, trip, option))
trip      <- dcast(trip ,option~year, value.var="trip")
trip[is.na(trip)]<- 0
rownames(trip)  <- trip$option
trip            <- trip[,-1]
trip_percentage <- matrix(apply(trip, 2, function(x){x*100/sum(x,na.rm=T)}), nrow=nrow(trip))
rownames(trip_percentage) <- rownames(trip)
colnames(trip_percentage) <- colnames(trip)

days      <- with(effort,data.frame(year, effort, option))
days<- days[!days$option=="Stay in port",]
days      <- dcast(days ,option~year, value.var="effort")
days[is.na(days)]<- 0
rownames(days)  <- days$option
days            <- days[,-1]

#economics
netrev      <-  melt(economics_res_allyrs, id.vars = "year", measure.vars = c("NetRev"))
grossrev    <-  melt(economics_res_allyrs, id.vars = "year", measure.vars = c("Grossrev"))
annualfine  <-  melt(economics_res_allyrs, id.vars = "year", measure.vars = c("Annualfine"))
fuelcosts   <-  melt(economics_res_allyrs, id.vars = "year", measure.vars = c("Fuelcosts"))
