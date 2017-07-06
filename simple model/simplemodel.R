options(width=200)
library(RDynState5NAsigmaseason6Age)
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)

setwd("~/wur/N/Projecten/MARIFLEET/Nekane/simple\ model/")
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

yield_curve <- function(hr,wts, natmortality, R=1, sequence = seq(0.001,2,0.001), verbose=F ){
  
  # note that definition of hr is not completely correct (should be sum over seasons, and mean over ages), but as long as consistently incorect in code it should not matter
  res <- data.frame("hr"=mean(hr) *sequence,"yield"=NA) 
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
NUMRUNS       <- 80
MPstart       <- 40
SIMNUMBER     <- 850
SIGMA         <- 300 #comes from 2
SPP1DSCSTEPS  <- SPP2DSCSTEPS <- 0
endy          <- stab.model + NUMRUNS
Linf          <- 20
K             <- 0.3
wts           <- Linf*(1-exp(-K*ages))
q             <- 0.0005
natmortality  <- 0.0001

recs1          <- c(100,0) 
mig1     <- array(0, dim=c(length(ages),1,length(season),length(areas), length(areas)), dimnames=list(cat=ages,year="all",season=as.character(season), from =areas, to=areas)) 
mig1[,,,"a","a"] <- -0.2
mig1[,,,"b","b"] <- -0.2
mig1[,,,"a","b"] <- 0.2
mig1[,,,"b","a"] <- 0.2
aperm( mig1,c(1,3,2,4,5))

recs2          <- c(100,0) 
mig2     <- array(0, dim=c(length(ages),1,length(season),length(areas), length(areas)), dimnames=list(cat=ages,year="all",season=as.character(season), from =areas, to=areas)) 
mig2[,,,"a","a"] <- -0.2
mig2[,,,"b","b"] <- -0.2
mig2[,,,"a","b"] <- 0.2
mig2[,,,"b","a"] <- 0.2
aperm( mig2,c(1,3,2,4,5))

effort <- array(c(1), dim=c(length(areas), length(season)), dimnames=list(option =areas,season=as.character(season)))


pop1  <- pop2   <-array(0, dim=c(length(ages),endy + 1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.n.dsvm1      <- catches.n.dsvm2      <- array(0, dim=c(length(ages),endy + 1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.wt.dsvm1     <- catches.wt.dsvm2     <- array(0, dim=c(length(ages),endy + 1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.wt.dsvm.tot1 <- catches.wt.dsvm.tot2 <- array(0, dim=c(1           ,endy + 1,              1,            1), dimnames=list(cat="all", year=as.character(1:(endy+1)), season="all",                option ="all"))
quota1               <- quota2               <- array(NA, dim=c(1           ,endy + 1,              1,            1), dimnames=list(cat="all", year=as.character(1:(endy+1)), season="all",                option ="all"))


#run population for 15 year
pop1 <- population_dynamics(pop=pop1, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs1, migration=mig1)
pop2 <- population_dynamics(pop=pop2, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs2, migration=mig2)

#calculated catches can then be used for input to DSVM (has same dims as pop (1: endyr), endyr=stabmodel+numruns)
pos_catches1 <- pop1 *q*wts
pos_catches2 <- pop2 *q*wts

#set up dsvm
sp1<- sp2 <- sp3 <- sp4 <- sp5 <-    new("DynStateInput")
catchMean(sp3)  <- catchMean(sp4) <- catchMean(sp5) <- array(0.01,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))
catchSigma(sp3) <- catchSigma(sp4)<- catchSigma(sp5)<- array(0.001,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))

control     <- DynState.control(spp1LndQuota= 800,  spp2LndQuota=800, spp1LndQuotaFine= 20000, spp2LndQuotaFine= 20000,
                                fuelUse = 0.001, fuelPrice = 1.0, landingCosts= 0,gearMaintenance= 0, addNoFishing= TRUE, increments= 25,
                                spp1DiscardSteps= SPP1DSCSTEPS, spp2DiscardSteps= SPP2DSCSTEPS, sigma= SIGMA, simNumber= SIMNUMBER, numThreads= 20)

#this is where our loop starts, after we set up stable population
for(yy in (stab.model):(stab.model+NUMRUNS-1)){
  
  print("====== year yy ========")
  print(yy)

  catchMean(sp1)  <- array(apply(pos_catches1[,(yy-2):yy,,,drop=F],c(1,3,4),mean), dim=c(length(ages), length(season),length(areas)),  dimnames=list("cat"=ages,"season"= season,"option"=areas))
  catchMean(sp2)  <- array(apply(pos_catches1[,(yy-2):yy,,,drop=F],c(1,3,4),mean), dim=c(length(ages), length(season),length(areas)),  dimnames=list("cat"=ages,"season"= season,"option"=areas))
  
  # ---No way of estimating sigma, therefore we assume that is 8% of the CPUE (note slight repetion in code for dims and dimnames of 0 catch arrays for spec 3,4,5)                                                                  
  catchSigma(sp1) <- catchMean(sp1) *0.08
  catchSigma(sp2) <- catchMean(sp2) *0.08
  
  sp1Price <-  sp2Price <- sp3Price <- sp4Price <- sp5Price <- array(c(1000), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
  #---effort and prices used (note that now c is removed (but that if other runs, then make sure to fix/remove code that removes "c" option)                                                                                         

  # if we are in MP period, then set quota based on last year
  if (yy > (stab.model + MPstart))
    control@spp1LndQuota <-  quota1[,yy,,]
  
  #run DSVM (wiht quota constraining if in MP time)   
  z <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
  
  #extract DSVM results
  dsvm_res <-  extract_dsvm_res (z, control, ages, season)
  
  if (yy == stab.model){ 
   dsvm_res_allyrs  <- cbind("year"= yy,dsvm_res)
  } else {
    dsvm_res_allyrs <- rbind(dsvm_res_allyrs, (cbind("year"= yy,dsvm_res)))
  } 

    
  #get catches in wts from DSVM 
  catches.wt.dsvm1[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop1") 
  catches.wt.dsvm2[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop2") 
  
  #calculate total catches (by summing over seasons and ages)
  catches.wt.dsvm.tot1[] <- apply(catches.wt.dsvm1,c(2),"sum")
  catches.wt.dsvm.tot2[] <- apply(catches.wt.dsvm2,c(2),"sum")
  
  #calculate numbers caught from weight caught 
  catches.n.dsvm1 <- catches.wt.dsvm1/wts
  catches.n.dsvm2 <- catches.wt.dsvm2/wts

  # calculatae what happens to population based on catches
  pop1 <- population_dynamics(pop=pop1, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,yy,,,drop=F], recruitment=recs1, migration=mig1)
  pop2 <- population_dynamics(pop=pop2, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm2[,yy,,,drop=F], recruitment=recs2, migration=mig2)

  #calculate the catches that can be input into DSVM based on updated pop
  pos_catches1 <- pop1 *q*wts
  pos_catches2 <- pop2 *q*wts
  
  #MANAGEMENT PROCEDURE
  hr1 <- apply(catches.n.dsvm1,1:3,sum)/    (apply(catches.n.dsvm1,1:3,sum) +    apply(pop1,1:3,sum) )
  hr2 <- apply(catches.n.dsvm2,1:3,sum)/    (apply(catches.n.dsvm2,1:3,sum) +    apply(pop2,1:3,sum) )

  yc1 <- yield_curve(hr=hr1[,yy,], wts, natmortality, R=recs1, verbose=F)
  yc2 <- yield_curve(hr=hr2[,yy,], wts, natmortality, R=recs2, verbose=F)
  
  hr1wanted <- yc1[yc1$yield==max(yc1$yield),]$hr
  hr2wanted <- yc2[yc2$yield==max(yc2$yield),]$hr
  
  if (yy > (stab.model + MPstart-1))
    quota1[,yy+1,,] <- sum(sweep((hr1wanted/mean(hr1[,yy,])) * hr1[,yy,] * apply(pop1[,yy,,],c(1,2), sum) ,1,wts,"*"))/SIMNUMBER 
  
}



#what are the weights?
wts

hr1 <- apply(catches.n.dsvm1,1:3,sum)/    (apply(catches.n.dsvm1,1:3,sum) +    apply(pop1,1:3,sum) )
hr2 <- apply(catches.n.dsvm2,1:3,sum)/    (apply(catches.n.dsvm2,1:3,sum) +    apply(pop2,1:3,sum) )

pyrnoMP <- 46 
pyrMP   <- 86 

#what happens in our yield curve for this hr?
yield_curve(hr=hr1[,pyr,], wts, natmortality, R=recs1, sequence = 1, verbose=T)
yield_curve(hr=hr2[,pyr,], wts, natmortality, R=recs2, sequence = 1, verbose=T)

round(catches.n.dsvm1[,pyrnoMP,,],2)
round(pop1[,pyrnoMP,,],2)

#next, what happens to theoretical pop for our harvest, what is harvest?
hr1[,pyrnoMP,]
mean(hr1[,pyrnoMP,])

yc1noMP <- yield_curve(hr=hr1[,pyrnoMP,], wts, natmortality, R=recs1, verbose=F)
yc2noMP <- yield_curve(hr=hr2[,pyrnoMP,], wts, natmortality, R=recs2, verbose=F)

Fmsy1noMP <- yc1noMP[yc1noMP$yield==max(yc1noMP$yield),]$hr

yc1MP <- yield_curve(hr=hr1[,pyrMP,], wts, natmortality, R=recs1, verbose=F)
yc2MP <- yield_curve(hr=hr2[,pyrMP,], wts, natmortality, R=recs2, verbose=F)

Fmsy1MP <- yc1MP[yc1MP$yield==max(yc1MP$yield),]$hr

ylim    <- c(0,900)
xlimYPR <- c(0,0.2)

#to check
par(mfrow=c(2,3))
plot(catches.wt.dsvm.tot1, type="l", ylim=ylim,xaxs='i', yaxs='i')
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="grey")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="grey")

lines((quota1*SIMNUMBER), col="red" )
abline(v=stab.model + MPstart, lty=2)
lines(catches.wt.dsvm.tot1, type="l", ylim=ylim)

plot(x=yc1noMP$hr, y=yc1noMP$yield, type="l", xlim=xlimYPR, ylim=ylim)
abline(v=Fmsy1noMP)
text(xlimYPR[2]*0.9, ylim[2]*0.9, paste0("SIMNUMBER ",SIMNUMBER))
points(mean(hr1[,pyrnoMP,]),yc1noMP$yield[yc1noMP$hr>mean(hr1[,pyrnoMP,])][1], col="red", pch=19)
points(mean(hr1[,pyrnoMP-2,]),catches.wt.dsvm.tot1[,pyrnoMP-2,,], col="blue", pch=19)
points(mean(hr1[,pyrnoMP-1,]),catches.wt.dsvm.tot1[,pyrnoMP-1,,], col="blue", pch=19)
points(mean(hr1[,pyrnoMP,]),catches.wt.dsvm.tot1[,pyrnoMP,,], col="blue", pch=19)
points(mean(hr1[,pyrnoMP+1,]),catches.wt.dsvm.tot1[,pyrnoMP+1,,], col="blue", pch=19)
points(mean(hr1[,pyrnoMP+2,]),catches.wt.dsvm.tot1[,pyrnoMP+2,,], col="blue", pch=19)
lines(x=yc1MP$hr, y=yc1MP$yield, ylim=ylim, col="grey")
abline(v=Fmsy1MP, col="grey")
points(mean(hr1[,pyrMP,]),yc1MP$yield[yc1MP$hr>mean(hr1[,pyrMP,])][1], col="red", pch=19)

plot(apply(catches.wt.dsvm1,c(2,4),sum)[,1], col="blue", type="l",  ylim=ylim, xaxs='i', yaxs='i')
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="grey")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="grey")
lines(apply(catches.wt.dsvm1,c(2,4),sum)[,2], col="red")
lines(apply(catches.wt.dsvm1,c(2,4),sum)[,3], col="black")
legend("topleft",c("a","b"), col=c("blue","red"), lty=c(1,1))
abline(v=stab.model + MPstart, lty=2)

dsvm_res_allyrs[dsvm_res_allyrs$year %in% ((pyr-1):(pyr+1))  & dsvm_res_allyrs$spp == "sp1",]

round(pop1,0)

#to check
plot(catches.wt.dsvm.tot2, type="l", ylim=ylim)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="grey")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="grey")
abline(v=stab.model + MPstart, lty=2)
lines(catches.wt.dsvm.tot2, type="l", ylim=ylim)

plot(x=yc2noMP$hr, y=yc2noMP$yield, type="l", xlim=xlimYPR, ylim=ylim)
points(mean(hr2[,pyrnoMP,]),yc2$yield[yc2$hr>mean(hr2[,pyrnoMP,])][1], col="red", pch=19)
points(mean(hr2[,pyrnoMP-2,]),catches.wt.dsvm.tot2[,pyrnoMP-2,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP-1,]),catches.wt.dsvm.tot2[,pyrnoMP-1,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP,]),catches.wt.dsvm.tot2[,pyrnoMP,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP+1,]),catches.wt.dsvm.tot2[,pyrnoMP+1,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP+2,]),catches.wt.dsvm.tot2[,pyrnoMP+2,,], col="blue", pch=19)
lines(x=yc2MP$hr, y=yc2MP$yield, ylim=ylim, col="grey")
points(mean(hr2[,pyrMP,]),yc2MP$yield[yc2MP$hr>mean(hr1[,pyrMP,])][1], col="red", pch=19)

plot(apply(catches.wt.dsvm2,c(2,4),sum)[,1], col="blue", type="l",  ylim=ylim)
lines(apply(catches.wt.dsvm2,c(2,4),sum)[,2], col="red")
lines(apply(catches.wt.dsvm2,c(2,4),sum)[,3], col="black")
abline(v=stab.model + MPstart, lty=2)

dsvm_res_allyrs[dsvm_res_allyrs$year %in% ((pyr-1):(pyr+1))  & dsvm_res_allyrs$spp == "sp2",]

round(pop2,0)
