options(width=200)
library(RDynState5NAsigmaseason6Age)
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)

setwd("~/Dropbox/BoB/MSE/Git/Nekane/complex model/")
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
    res[iii, ]$catch    <- sum(yld*wts)
    res[iii, ]$landings <- sum(yld*lratio*wts)
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
NUMRUNS       <- 80
MPstart       <- 40
MPstartLO     <- 65
SIMNUMBER     <- 700 #pos
SIGMA         <- 40 #sig 
SPP1DSCSTEPS  <- 1
SPP2DSCSTEPS  <- 1
endy          <- stab.model + NUMRUNS
Linf          <- 20
K             <- 0.3
wts           <- Linf*(1-exp(-K*ages))
q             <- 0.0005
natmortality  <- 0.0001
migconstant   <- 0.2
sp1price      <- sp2price      <- 150
slope1price <- 150
slope2price <- 0.50*150
# scenario I: discarding is not allowed, YPR based in C (C=L)
# scenario II: discarding is allowed, YPR based in L, hr wanted based in catches
# scenario III: discarding ocurred but not perceived, YPR based in L, hr wanted based in landings

recs1          <- c(100,0) 
mig1     <- array(0, dim=c(length(ages),1,length(season),length(areas), length(areas)), dimnames=list(cat=ages,year="all",season=as.character(season), from =areas, to=areas)) 
mig1[,,,"a","a"] <- -migconstant
mig1[,,,"b","b"] <- -migconstant
mig1[,,,"a","b"] <- migconstant
mig1[,,,"b","a"] <- migconstant
aperm( mig1,c(1,3,2,4,5))

recs2          <- c(0,100) 
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
quota1                <- quota2                <- array(1.2, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))


#run population for 15 year
pop1 <- population_dynamics(pop=pop1, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs1, migration=mig1)
pop2 <- population_dynamics(pop=pop2, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs2, migration=mig2)

#calculated catches can then be used for input to DSVM (has same dims as pop (1: endyr), endyr=stabmodel+numruns)
pos_catches1 <- pop1 *q*wts
pos_catches2 <- pop2 *q*wts

#set up dsvm
sp1<- sp2 <- sp3 <- sp4 <- sp5 <-    new("DynStateInput")
catchMean(sp3)  <- catchMean(sp4) <- catchMean(sp5) <- array(0.01,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))
catchSigma(sp3) <- catchSigma(sp4)<- catchSigma(sp5)<- array(0.0000001,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))

#SIZE DEPENDENT PRICING, following Zimmermann et al. (2011)
sp1Price <- array(c(sp1price+ slope1price*((wts-mean(wts)/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
sp2Price <- array(c(sp2price+ slope2price*((wts-mean(wts)/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
sp3Price <- sp4Price <- sp5Price <- array(c(0), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
#---effort and prices used (note that now c is removed (but that if other runs, then make sure to fix/remove code that removes "c" option)                                                                                         
control     <- DynState.control(spp1LndQuota= 200,  spp2LndQuota=200, spp1LndQuotaFine= 3e6, spp2LndQuotaFine= 3e6, fuelUse = 1, fuelPrice = 150.0, landingCosts= 0,gearMaintenance= 0, addNoFishing= TRUE, increments= 25, spp1DiscardSteps= SPP1DSCSTEPS, spp2DiscardSteps= SPP2DSCSTEPS, sigma= SIGMA, simNumber= SIMNUMBER, numThreads= 20)

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
  catches.n.dsvm1 <- catches.wt.dsvm1/wts
  catches.n.dsvm2 <- catches.wt.dsvm2/wts
  
  landings.n.dsvm1 <- landings.wt.dsvm1/wts
  landings.n.dsvm2 <- landings.wt.dsvm2/wts
  
  # calculatae what happens to population based on catches
  pop1 <- population_dynamics(pop=pop1, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,yy,,,drop=F], recruitment=recs1, migration=mig1)
  pop2 <- population_dynamics(pop=pop2, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm2[,yy,,,drop=F], recruitment=recs2, migration=mig2)
  
  #calculate the catches that can be input into DSVM based on updated pop
  pos_catches1 <- pop1 *q*wts
  pos_catches2 <- pop2 *q*wts
  
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
    prel.quota1 <-  sum(sweep((hr1wanted[1]/mean(hr1[,yy,]))* hr1[,yy,]*landings.ratio1[,yy,]*apply(pop1[,yy+1,,],c(1,2), sum) ,1,wts,"*"))/SIMNUMBER
    prel.quota2 <-  sum(sweep((hr2wanted[1]/mean(hr2[,yy,]))* hr2[,yy,]*landings.ratio2[,yy,]*apply(pop2[,yy+1,,],c(1,2), sum) ,1,wts,"*"))/SIMNUMBER
    
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



ylim <- c(0,1000)
xlimYPR <- c(0,0.2)

setwd("~/Dropbox/BoB/MSE/Git/Nekane/doc/Figures")
png(filename=paste("MIXEDMP2_CATCH_", paste0("SIGMA ", SIGMA, "; INCREMENTS ",control@increments, "; SIMNUMBER ",SIMNUMBER, "; DISCARDSTEPS ",SPP1DSCSTEPS, "; MIGRATION ",migconstant, "; REC1 ",paste(recs1, collapse = " "), "; REC2 ",paste(recs2, collapse = " "), "; SP1PRICE ",paste0(round(sp1Price[,1]), sep = ', ', collapse = ''),";SP2PRICE ",paste0(round(sp2Price[,1]), sep = ', ', collapse = ''), "; FUELPRICE ",control@fuelPrice,".png"), sep=""),width=20, height=8, units="cm", res=500, pointsize=6.5)
#to check

par(mfrow=c(2,5),oma = c(3,0,0,0) + 0.1, mar = c(4,4,1,1) + 0.1)
#,oma = c(4,3,1,1) + 0.1, mar = c(4,4,1,1) + 0.1)

plot(catches.wt.dsvm.tot1, type="l",  xlim=c(10,90), ylim=ylim, xaxs='i', yaxs='i', xlab= "Year", ylab = "Total catches (weight)", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="ivory3")
polygon(x=c(pyrMPLO-2,pyrMPLO+2,pyrMPLO+2,pyrMPLO-2), border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1, at = c(10, 30, 40, 50, 65, 70, 90),
     labels = c(10, 30, "MP", 50, "LO" ,70,90))
axis(1, at = c(40, 65),labels = c("MP",  "LO" ))
abline(v=c(30, 50, 70), lty="dotted", col = "ivory2")
abline(v=MPstart, lty=2)
abline(v=MPstartLO, lty=2)
lines((quota1* SIMNUMBER), col="red")
lines(catches.wt.dsvm.tot1,  type="l", ylim=ylim)
lines(landings.wt.dsvm.tot1, type="l", lty=2, ylim=ylim)
legend("bottomright", inset=.05, legend=c("Catches","Landings","TAC"), pch=c(1,46,1), col=c("black","black","red"), bty='n', cex=0.8)

plot(apply(hr1,c(1,2),mean)[1,], type="p",  xlim=c(10,90), ylim=c(0,0.2), xaxs='i', yaxs='i', xlab= "Year", ylab = "Harvest rates", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="ivory3")
polygon(x=c(pyrMPLO-2,pyrMPLO+2,pyrMPLO+2,pyrMPLO-2), border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1, at = c(10, 30, 40, 50, 65, 70, 90),
     labels = c(10, 30, "MP", 50, "LO" ,70,90))
axis(1, at = c(40, 65),labels = c("MP",  "LO" ))
abline(v=c(30, 50, 70), lty="dotted", col = "ivory2")
abline(h=c(0.05, 0.1, 0.15), lty="dotted", col = "ivory2")
abline(v=MPstart, lty=2)
abline(v=MPstartLO, lty=2)
points(apply(hr1,c(1,2),mean)[1,], col="black", pch=19)
points(apply(hr1,c(1,2),mean)[2,], col="red", pch=19)
points(apply(hr1,c(1,2),mean)[3,], col="black", pch=21)
points(apply(hr1,c(1,2),mean)[4,], col="red", pch=21)
legend("topright", inset=.05, legend=c("Age 1","Age 2","Age 3","Age 4"), pch=c(19,19,21,21), col=c("black","red", "black","red"), bty='n', cex=0.8)

plot(x=yc1noMP$hr, y=yc1noMP$landings, type="l", xlim=xlimYPR, ylim=ylim,xaxs='i', yaxs='i',  xlab="Harvest rate", ylab = "Yield per recruit", panel.first=grid(col = "ivory2"))
text(xlimYPR[2]*0.8, yc1noMP$landings[length(yc1noMP$hr)]+5, "Unconstrained")
abline(v=Fmsy1noMP, col="ivory4")
#text(xlimYPR[2]*0.8, ylim[2]*0.9, paste0("SIMNUMBER ",SIMNUMBER))
points(mean(hr1[,pyrnoMP,]),yc1noMP$landings[yc1noMP$hr>mean(hr1[,pyrnoMP,])][1], col="red", pch=19)
points(mean(hr1[,pyrnoMP-2,]),landings.wt.dsvm.tot1[,pyrnoMP-2,,],col="ivory4", pch=19)
points(mean(hr1[,pyrnoMP-1,]),landings.wt.dsvm.tot1[,pyrnoMP-1,,], col="ivory4", pch=19)
points(mean(hr1[,pyrnoMP,]),landings.wt.dsvm.tot1[,pyrnoMP,,], col="ivory4", pch=19)
points(mean(hr1[,pyrnoMP+1,]),landings.wt.dsvm.tot1[,pyrnoMP+1,,], col="ivory4", pch=19)
points(mean(hr1[,pyrnoMP+2,]),landings.wt.dsvm.tot1[,pyrnoMP+2,,], col="ivory4", pch=19)
lines(x=yc1MP$hr, y=yc1MP$landings, ylim=ylim, col="ivory3")
text(0.15, yc1MP$landings[length(yc1MP$hr)]+80, "Constrained \n 15% TAC change")
abline(v=Fmsy1MP, col="ivory3")
points(mean(hr1[,pyrMP,]),yc1MP$landings[yc1MP$hr>mean(hr1[,pyrMP,])][1], col="red", pch=21, bg="white")
points(mean(hr1[,pyrMP-2,]),landings.wt.dsvm.tot1[,pyrMP-2,,], col="ivory3", pch=21, bg="white")
points(mean(hr1[,pyrMP-1,]),landings.wt.dsvm.tot1[,pyrMP-1,,], col="ivory3", pch=21, bg="white")
points(mean(hr1[,pyrMP,]),landings.wt.dsvm.tot1[,pyrMP,,], col="ivory3", pch=21, bg="white")
points(mean(hr1[,pyrMP+1,]),landings.wt.dsvm.tot1[,pyrMP+1,,], col="ivory3", pch=21, bg="white")
points(mean(hr1[,pyrMP+2,]),landings.wt.dsvm.tot1[,pyrMP+2,,], col="ivory3", pch=21, bg="white")
lines(x=yc1MPLO$hr, y=yc1MPLO$landings, ylim=ylim, col="darkgreen")
text(yc1MPLO$hr[length(yc1MPLO$hr)]-0.02, yc1MPLO$landings[length(yc1MPLO$hr)]+80, "Constrained LO (MSY)")
abline(v=Fmsy1MPLO, col="darkgreen")
points(mean(hr1[,pyrMPLO,]),yc1MPLO$landings[yc1MPLO$hr>mean(hr1[,pyrMPLO,])][1], col="red", pch=21, bg="white")
points(mean(hr1[,pyrMPLO-2,]),landings.wt.dsvm.tot1[,pyrMPLO-2,,], col="darkgreen", pch=21, bg="white")
points(mean(hr1[,pyrMPLO-1,]),landings.wt.dsvm.tot1[,pyrMPLO-1,,], col="darkgreen", pch=21, bg="white")
points(mean(hr1[,pyrMPLO,]),landings.wt.dsvm.tot1[,pyrMPLO,,], col="darkgreen", pch=21, bg="white")
points(mean(hr1[,pyrMPLO+1,]),landings.wt.dsvm.tot1[,pyrMPLO+1,,], col="darkgreen", pch=21, bg="white")
points(mean(hr1[,pyrMPLO+2,]),landings.wt.dsvm.tot1[,pyrMPLO+2,,], col="darkgreen", pch=21, bg="white")

plot(rowMeans(hr1[,pyrnoMP,]), type="b", ylim=c(0,.2),  xlab="Age", ylab = "Selectivity", panel.first=grid(col = "ivory2"), xaxt="n")
text(1,rowMeans(hr1[,pyrnoMP,])[1]+0.01, "Unconstrained", pos=4)
lines(rowMeans(hr1[,pyrMP,]), type="b", ylim=c(0,.2), col="grey")
text(1,rowMeans(hr1[,pyrMP,])[1]+0.01, "Constrained 15% TAC change", pos=4)
lines(rowMeans(hr1[,pyrMPLO,]), type="b", ylim=c(0,.2), col="darkgreen")
text(1,rowMeans(hr1[,pyrMPLO,])[1]-0.01, "Constrained LO (MSY)", pos=4)
axis(1, at = seq(1, 4, by = 1))

plot(apply(catches.wt.dsvm1,c(2,4),sum)[,1], col="blue", type="l",  ylim=ylim,xaxs='i', yaxs='i', xlim=c(10,90),xlab = "Year", ylab = "Catches by area (weight)",xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="ivory3")
polygon(x=c(pyrMPLO-2,pyrMPLO+2,pyrMPLO+2,pyrMPLO-2), border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1, at = c(10, 30, 40, 50, 65, 70, 90),
     labels = c(10, 30, "MP", 50, "LO" ,70,90))
axis(1, at = c(40, 65),labels = c("MP",  "LO" ))
abline(v=c(30, 50, 70), lty="dotted", col = "ivory2")
abline(v=MPstart, lty=2)
abline(v=MPstartLO, lty=2)
lines(apply(catches.wt.dsvm1,c(2,4),sum)[,1], col="blue", type="l")
lines(apply(catches.wt.dsvm1,c(2,4),sum)[,2], col="red")
lines(apply(landings.wt.dsvm1,c(2,4),sum)[,1], col="blue", lty=2)
lines(apply(landings.wt.dsvm1,c(2,4),sum)[,2], col="red", lty=2)
#lines(apply(catches.wt.dsvm1,c(2,4),sum)[,3], col="black")
legend("topright", inset=.05, c("a","b"),  pch=c(1,1), col=c("blue","red"), bty='n', cex=0.8)
abline(v=MPstart, lty=2)


#round(pop1,0)

plot(catches.wt.dsvm.tot2, type="l",  xlim=c(10,90), ylim=ylim, xaxs='i', yaxs='i', xlab= "Year", ylab = "Total catches (weight)", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="ivory3")
polygon(x=c(pyrMPLO-2,pyrMPLO+2,pyrMPLO+2,pyrMPLO-2), border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1, at = c(10, 30, 40, 50, 65, 70, 90),
     labels = c(10, 30, "MP", 50, "LO" ,70,90))
axis(1, at = c(40, 65),labels = c("MP",  "LO" ))
abline(v=c(30, 50, 70), lty="dotted", col = "ivory2")
abline(v=MPstart, lty=2)
abline(v=MPstartLO, lty=2)
lines((quota2* SIMNUMBER), col="red")
lines(catches.wt.dsvm.tot2,  type="l", ylim=ylim)
lines(landings.wt.dsvm.tot2, type="l", ylim=ylim, lty= 2)

plot(apply(hr2,c(1,2),mean)[1,], type="p",  xlim=c(10,90), ylim=c(0,0.2), xaxs='i', yaxs='i', xlab= "Year", ylab = "Harvest rates", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="ivory3")
polygon(x=c(pyrMPLO-2,pyrMPLO+2,pyrMPLO+2,pyrMPLO-2), border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1, at = c(10, 30, 40, 50, 65, 70, 90),
     labels = c(10, 30, "MP", 50, "LO" ,70,90))
axis(1, at = c(40, 65),labels = c("MP",  "LO" ))
abline(v=c(30, 50, 70), lty="dotted", col = "ivory2")
abline(v=MPstart, lty=2)
abline(v=MPstartLO, lty=2)
points(apply(hr2,c(1,2),mean)[1,], col="black", pch=19)
points(apply(hr2,c(1,2),mean)[2,], col="red", pch=19)
points(apply(hr2,c(1,2),mean)[3,], col="black", pch=21)
points(apply(hr2,c(1,2),mean)[4,], col="red", pch=21)

plot(x=yc2noMP$hr, y=yc2noMP$landings, type="l", xlim=xlimYPR, ylim=ylim,xaxs='i', yaxs='i', xlab = "Harvest rate", ylab = "Yield per recruit", panel.first=grid(col = "ivory2"))
text(xlimYPR[2]*0.8, yc2noMP$landings[length(yc2noMP$hr)]+5, "Unconstrained")
abline(v=Fmsy2noMP)
points(mean(hr2[,pyrnoMP,]),yc2noMP$landings[yc2noMP$hr>mean(hr2[,pyrnoMP,])][1], col="red", pch=19)
points(mean(hr2[,pyrnoMP-2,]),landings.wt.dsvm.tot2[,pyrnoMP-2,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP-1,]),landings.wt.dsvm.tot2[,pyrnoMP-1,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP,]),landings.wt.dsvm.tot2[,pyrnoMP,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP+1,]),landings.wt.dsvm.tot2[,pyrnoMP+1,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP+2,]),landings.wt.dsvm.tot2[,pyrnoMP+2,,], col="blue", pch=19)
lines(x=yc2MP$hr, y=yc2MP$landings, ylim=ylim, col="grey")
text(yc2MP$hr[length(yc2MP$hr)]-0.02, yc2MP$landings[length(yc2MP$hr)]+80, "Constrained 15% TAC change")
abline(v=Fmsy2MP, col="grey")
points(mean(hr2[,pyrMP,]),yc2MP$landings[yc2MP$hr>mean(hr2[,pyrMP,])][1], col="red", pch=21, bg="white")
points(mean(hr2[,pyrMP-2,]),landings.wt.dsvm.tot2[,pyrMP-2,,], col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP-1,]),landings.wt.dsvm.tot2[,pyrMP-1,,], col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP,]),landings.wt.dsvm.tot2[,pyrMP,,], col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP+1,]),landings.wt.dsvm.tot2[,pyrMP+1,,], col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP+2,]),landings.wt.dsvm.tot2[,pyrMP+2,,], col="blue", pch=21, bg="white")
lines(x=yc2MPLO$hr, y=yc2MPLO$landings, ylim=ylim, col="darkgreen")
text(yc2MPLO$hr[length(yc2MPLO$hr)]-0.02, yc2MPLO$landings[length(yc2MPLO$hr)]+80, "Constrained LO (MSY)")
abline(v=Fmsy2MP, col="darkgreen")
points(mean(hr2[,pyrMPLO,]),yc2MPLO$landings[yc2MPLO$hr>mean(hr2[,pyrMPLO,])][1], col="red", pch=21, bg="white")
points(mean(hr2[,pyrMPLO-2,]),landings.wt.dsvm.tot2[,pyrMPLO-2,,], col="darkgreen", pch=21, bg="white")
points(mean(hr2[,pyrMPLO-1,]),landings.wt.dsvm.tot2[,pyrMPLO-1,,], col="darkgreen", pch=21, bg="white")
points(mean(hr2[,pyrMPLO,]),landings.wt.dsvm.tot2[,pyrMPLO,,], col="darkgreen", pch=21, bg="white")
points(mean(hr2[,pyrMPLO+1,]),landings.wt.dsvm.tot2[,pyrMPLO+1,,], col="darkgreen", pch=21, bg="white")
points(mean(hr2[,pyrMPLO+2,]),landings.wt.dsvm.tot2[,pyrMPLO+2,,], col="darkgreen", pch=21, bg="white")


plot(rowMeans(hr2[,pyrnoMP,]), type="b", ylim=c(0,.2), xlab = "Age", ylab = "Selectivity", panel.first=grid(col = "ivory2"), xaxt="n")
text(1,rowMeans(hr2[,pyrnoMP,])[1]+0.01, "Unconstrained",pos=4)
lines(rowMeans(hr2[,pyrMP,]), type="b", ylim=c(0,.2), col="grey")
text(1,rowMeans(hr2[,pyrMP,])[1]+0.01, "Constrained 15% TAC change",pos=4)
lines(rowMeans(hr2[,pyrMPLO,]), type="b", ylim=c(0,.2), col="darkgreen")
text(1,rowMeans(hr2[,pyrMPLO,])[1]-0.03, "Constrained LO (MSY)",pos=4)
axis(1, at = seq(1, 4, by = 1))

plot(apply(catches.wt.dsvm2,c(2,4),sum)[,1], col="blue", type="l",  ylim=ylim,xaxs='i', yaxs='i', xlim=c(10,90),xlab = "Year", ylab = "Catches by area (weight)",xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="ivory3")
polygon(x=c(pyrMPLO-2,pyrMPLO+2,pyrMPLO+2,pyrMPLO-2), border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1, at = c(10, 30, 40, 50, 65, 70, 90),
     labels = c(10, 30, "MP", 50, "LO" ,70,90))
axis(1, at = c(40, 65),labels = c("MP",  "LO" ))

abline(v=MPstart, lty=2)
abline(v=MPstartLO, lty=2)
lines(apply(catches.wt.dsvm2,c(2,4),sum)[,1], col="blue", type="l")
lines(apply(catches.wt.dsvm2,c(2,4),sum)[,2], col="red")
lines(apply(landings.wt.dsvm2,c(2,4),sum)[,1], col="blue", lty=2)
lines(apply(landings.wt.dsvm2,c(2,4),sum)[,2], col="red", lty=2)
#lines(apply(catches.wt.dsvm2,c(2,4),sum)[,3], col="black")

add_legend("bottomright", legend=paste0("SIGMA ", SIGMA, "; INCREMENTS ",control@increments, "; SIMNUMBER ",SIMNUMBER, "; DISCARDSTEPS1 ",SPP1DSCSTEPS,"; DISCARDSTEPS2 ",SPP2DSCSTEPS, "; MIGRATION ",migconstant, "; REC1 ",paste(recs1, collapse = " "), "; REC2 ",paste(recs2, collapse = " "), "; SP1PRICE ",paste0(round(sp1Price[,1]), sep = ', ', collapse = ''),";SP2PRICE ",paste0(round(sp2Price[,1]), sep = ', ', collapse = ''), "; FUELPRICE ",control@fuelPrice), col="black", horiz=TRUE, bty='n', cex=1)
#round(pop2,0)
#
dev.off()

# Effort pattern and economics
#effort_plot_dsvm(SIMNUMBER, dsvm_res_allyrs, stab.model, economics_res_allyrs)

#effort
effort      <- subset(dsvm_res_allyrs,(spp %in% "sp1"))
effort      <- subset(effort,(cat %in% 1))
effort$year <- factor(effort$year, levels= 1:90)
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


png(filename=paste("MIXEDMP2_EFFORT_", paste0("SIGMA ", SIGMA, "; INCREMENTS ",control@increments,"; SIMNUMBER ",SIMNUMBER, "; DISCARDSTEPS ",SPP1DSCSTEPS, "; MIGRATION ",migconstant, "; REC1 ",paste(recs1, collapse = " "), "; REC2 ",paste(recs2, collapse = " "), "; SP1PRICE ",paste0(round(sp1Price[,1]), sep = ', ', collapse = ''),";SP2PRICE ",paste0(round(sp2Price[,1]), sep = ', ', collapse = ''), "; FUELPRICE ",control@fuelPrice,".png"), sep=""),width=20, height=8, units="cm", res=500, pointsize=6.5)
#to check
par(mfrow=c(2,5),oma = c(3,0,0,0) + 0.1, mar = c(4,4,1,1) + 0.1)

#layout(matrix(c(1,1,2,2, 3,4,5, 6), 2, 4, byrow = TRUE))

mypalette <- c("#808080","#CCCCCC","#D55E00")
names(mypalette) <- c("a", "b", "Stay in port")
barplot(trip_percentage, col= mypalette, border=NA, xlim = c(1,81), xlab = "Year", ylab = "Effort pattern (%)", xaxt = "n",space = 0)
axis(1, at = c(0.5, 20.5, 30.5, 40.5, 60.5, 80.5),labels = c(10, 30,"MP",50,70,90))
legend("topright", inset=.05, legend=c("a", "b", "Stay in port"), fill=mypalette, cex=0.6)
abline(v=MPstart-9.5, lty=2)

boxplot(value ~ year, netrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean net revenue", 
        ylim=c(0,4000), xlim=c(1, 81), xaxt = "n")
axis(1, at = c(1, 21, 31, 41, 61, 81),labels = c(10, 30,"MP", 50,70,90))
grid(NA, NULL, col = "ivory2")
abline(v=c(1, 21, 41, 61, 81), lty="dotted", col = "ivory2")
abline(v=MPstart-9, lty=2)
boxplot(value ~ year, netrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean net revenue", 
        ylim=c(0,4000), xlim=c(1, 81), xaxt = "n", add=TRUE)

boxplot(value ~ year, grossrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean gross revenue", 
        ylim=c(0,4000), xlim=c(1,81), xaxt = "n")
axis(1, at = c(1, 21, 31, 41, 61, 81),labels = c(10, 30,"MP", 50,70,90))
grid(NA, NULL, col = "ivory2")
abline(v=c(1, 21, 41, 61, 81), lty="dotted", col = "ivory2")
abline(v=MPstart-9, lty=2)
boxplot(value ~ year, grossrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean gross revenue", 
        ylim=c(0,4000), xlim=c(1,81), xaxt = "n", add=TRUE)

boxplot(value ~ year, fuelcosts, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean fuel costs", 
        ylim=c(0,4000), xlim=c(1,81), xaxt = "n")
axis(1, at = c(1, 21, 31, 41, 61, 81),labels = c(10, 30,"MP", 50,70,90))
grid(NA, NULL, col = "ivory2")
abline(v=c(1, 21, 41, 61, 81), lty="dotted", col = "ivory2")
abline(v=MPstart-9, lty=2)
boxplot(value ~ year, fuelcosts, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean fuel costs", 
        ylim=c(0,4000), xlim=c(1,81), xaxt = "n", add=TRUE)

boxplot(value ~ year, annualfine, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean annual fine", 
        ylim=c(0,4000), xlim=c(1,81),  xaxt = "n")
axis(1, at = c(1, 21, 31, 41, 61, 81),labels = c(10, 30,"MP", 50,70,90))
grid(NA, NULL, col = "ivory2")
abline(v=c(1, 21, 41, 61, 81), lty="dotted", col = "ivory2")
abline(v=MPstart-9, lty=2)
boxplot(value ~ year, annualfine, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean annual fine", 
        ylim=c(0,4000), xlim=c(1,81), xaxt = "n", add=TRUE)

barplot(as.matrix(days), col=c("#808080","#CCCCCC"),border=NA, xlim = c(1,81), xlab = "Year", ylab = "Total days at sea", xaxt = "n",space = 0)
axis(1, at = c(0.5, 20.5, 30.5, 40.5, 60.5, 80.5),labels = c(10, 30,"MP",50,70,90))
abline(v=MPstart-9.5, lty=2)

plot(aggregate(value ~ year, FUN=sum, data=netrev)$year,aggregate(value ~ year, FUN=sum, data=netrev)$value, type="l", ylim=c(0,2091351), xlim=c(10,90), xlab = "Year", ylab = "Total net revenue", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
axis(1, at = c(10, 30, 40, 50, 70, 90),labels = c(10, 30, "MP", 50,70,90))
abline(v=c(10, 30, 50, 70, 90), lty="dotted", col = "ivory2")
abline(v=MPstart, lty=2)
lines(aggregate(value ~ year, FUN=sum, data=netrev)$year,aggregate(value ~ year, FUN=sum, data=netrev)$value, type="l")

plot(aggregate(value ~ year, FUN=sum, data=grossrev)$year, aggregate(value ~ year, FUN=sum, data=grossrev)$value, type="l", ylim=c(0,2091351), xlim=c(10,90), xlab = "Year", ylab = "Total gross revenue", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
axis(1, at = c(10, 30, 40, 50, 70, 90),labels = c(10, 30, "MP", 50,70,90))
abline(v=c(10, 30, 50, 70, 90), lty="dotted", col = "ivory2")
abline(v=MPstart, lty=2)
lines(aggregate(value ~ year, FUN=sum, data=grossrev)$year, aggregate(value ~ year, FUN=sum, data=grossrev)$value, type="l")

plot(aggregate(value ~ year, FUN=sum, data=fuelcosts)$year, aggregate(value ~ year, FUN=sum, data=fuelcosts)$value, type="l", ylim=c(0,2091351), xlim=c(10,90), xlab = "Year", ylab = "Total fuel costs", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
axis(1, at = c(10, 30, 40, 50, 70, 90),labels = c(10, 30, "MP", 50,70,90))
abline(v=c(10, 30, 50, 70, 90), lty="dotted", col = "ivory2")
abline(v=MPstart, lty=2)
lines(aggregate(value ~ year, FUN=sum, data=fuelcosts)$year, aggregate(value ~ year, FUN=sum, data=fuelcosts)$value, type="l")

plot(aggregate(value ~ year, FUN=sum, data=annualfine)$year, aggregate(value ~ year, FUN=sum, data=annualfine)$value, type="l", ylim=c(0,209130), xlim=c(10,90),xlab = "Year", ylab = "Total annual fine", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
axis(1, at = c(10, 30, 40, 50, 70, 90),labels = c(10, 30, "MP", 50,70,90))
abline(v=c(10, 30, 50, 70, 90), lty="dotted", col = "ivory2")
abline(v=MPstart, lty=2)
lines(aggregate(value ~ year, FUN=sum, data=annualfine)$year, aggregate(value ~ year, FUN=sum, data=annualfine)$value, type="l")
add_legend("bottomright", legend=paste0("SIGMA ", SIGMA, "; SIMNUMBER ",SIMNUMBER, "; DISCARDSTEPS ",SPP1DSCSTEPS, "; MIGRATION ",migconstant, "; REC1 ",paste(recs1, collapse = " "), "; REC2 ",paste(recs2, collapse = " "), "; SP1PRICE ",paste0(round(sp1Price[,1]), sep = ', ', collapse = ''),";SP2PRICE ",paste0(round(sp2Price[,1]), sep = ', ', collapse = ''), "; FUELPRICE ",control@fuelPrice), col="black", horiz=TRUE, bty='n', cex=1)

dev.off()

#  }
#}
