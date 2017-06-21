options(width=200)
library(RDynState5NAsigmaseason6Age)
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)


###############################################################################
# Extract the results from DSVM
##############################################################################

extract_dsvm_res <- function(z, control, ages, season){
  
  simNumber <-control@simNumber 
  sp        <- c("sp1","sp2","sp3","sp4","sp5")
  
  dsvm_res             <- as.data.frame(rbind(as.matrix(spp1Landings(sim(z))),
                                              as.matrix(spp2Landings(sim(z))),
                                              as.matrix(spp3Landings(sim(z))),
                                              as.matrix(spp4Landings(sim(z))),
                                              as.matrix(spp5Landings(sim(z)))))      
  names(dsvm_res)      <- "landings.wt"
  dsvm_res$discards.wt <- c(rbind(as.matrix(spp1Discards(sim(z))), 
                                  as.matrix(spp2Discards(sim(z))), 
                                  as.matrix(spp3Discards(sim(z))), 
                                  as.matrix(spp4Discards(sim(z))), 
                                  as.matrix(spp5Discards(sim(z)))))
  
  dsvm_res$catch.wt    <- dsvm_res$ landings + dsvm_res$discards
  dsvm_res$effort      <- rep(rep(as.matrix(effort(sim(z))),each=length(ages)),length(sp))
  dsvm_res$option      <- rep(rep(as.matrix(choice(sim(z))),each=length(ages)),length(sp))
  dsvm_res$spp         <- as.factor(c(rep(sp, each=(simNumber*length(ages)*length(season)))))
  dsvm_res$cat         <- ages
  dsvm_res$season      <- c(rep(season, each=simNumber*length(ages)))
  dsvm_res$vessel      <- rep(1:simNumber,each=length(ages))  
  dsvm_res$option[is.na(dsvm_res$option)] <- "Stay in port"
  
  dsvm_res[c(1:4)]     <- lapply(dsvm_res[c(1:4)], function(x) as.numeric(as.character(x)))
  dsvm_res[c(5:9)]     <- lapply(dsvm_res[c(5:9)], function(x) as.factor(x))
  
  is.num               <- sapply(dsvm_res, is.numeric)
  dsvm_res[is.num]     <- lapply(dsvm_res[is.num], round, 6)
  
  dsvm_res             <- subset(dsvm_res,(spp %in% c("sp1", "sp2")))
  
  trip                 <- count(dsvm_res,c("spp","cat","season","option"))
  names(trip)[5]       <- "trip"
  
  dsvm_res             <- aggregate(cbind(landings.wt, discards.wt, catch.wt, effort)~ spp+cat+season+option, FUN=sum, data=dsvm_res)  
  
  dsvm_res             <- merge(dsvm_res, trip, by=c("spp","cat", "season","option"),all.x=TRUE)
  
  return(dsvm_res)
}

##############################################################################
# Data frame to array
##############################################################################

catch_dataframe_to_array <- function(dsvm_result, ages, season, areas, stock){
  
  # need to revisit this function if we want to subset landings or discards separetly
  dsvm_catch <- with(dsvm_result,data.frame(cat, season, option, spp, catch.wt)) 
  dsvm_catch <- melt(dsvm_catch, id=c("cat", "season", "option","spp"))
  levels(dsvm_catch$spp) <- c(levels(dsvm_catch$spp), "pop1", "pop2")
  dsvm_catch$spp[dsvm_catch$spp=="sp1"]<-"pop1"
  dsvm_catch$spp[dsvm_catch$spp=="sp2"]<-"pop2"
  
  #subset the desired stock population
  dsvm_catch_wt  <-subset(dsvm_catch,(spp %in% stock))[-c(4,5)]
  
  #trick to avoid deciding what to do with staying at port, it is just to sum catches
  dsvm_catch_wt            <- dsvm_catch_wt[dsvm_catch_wt$option!="Stay in port",] 
  
  dim          <- c(length(ages),length(season),length(areas)) 
  dimnames     <- list(cat=ages,season=as.character(season),option =areas)
  empty.df     <-melt(array(c(0), dim=dim, dimnames=dimnames))
  
  catch        <- merge(empty.df,dsvm_catch_wt,by=c("cat", "season","option"),all.x=TRUE)
  
  catch        <- catch[-4]
  catch[is.na(catch)] <- 0
  names(catch) <- c("cat", "season","option","data")
  
  catch        <- tapply(catch[,"data"], list(factor(x = catch[,"cat"], levels = unique(catch[,"cat"])), 
                                              factor(x = catch[,"season"], levels = unique(catch[,"season"])), 
                                              factor(x = catch[,"option"], levels = unique(areas))), sum)
  return(catch)
}

##############################################################################
# POPULATION DYNAMICS
##############################################################################

population_dynamics <- function(pop, recruitarea, startyear, endyear, season, natmortality, catches, recruitment){
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
        pop[1,as.character(y),as.character(ss),recruitarea] <- recruitment;
      }
      # natural mortality  ---------------------                                                             
      pop[,as.character(y),as.character(ss),] <- pop[,as.character(y),as.character(ss),]*(1-natmortality)
      pop[pop < 1e-20 ] <- 1e-20
      
      # maturation ---------------------                                                            
      nummatures <- pop[,as.character(y),as.character(ss),recruitarea] * c(0,0.1,0.3,0.5) # last is vector of maturity              
      
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
stab.model    <- 20
NUMRUNS       <- 50
SIMNUMBER     <- 410
SIGMA         <- 1e-20
SPP1DSCSTEPS  <- SPP2DSCSTEPS <- 0
endy          <- stab.model + NUMRUNS
Linf          <- 20
K             <- 0.3
wts           <- Linf*(1-exp(-K*ages))
q             <- 0.0005
natmortality  <- 0.0001

pop1                <- array(0, dim=c(length(ages),endy + 1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.n.dsvm      <- array(0, dim=c(length(ages),endy + 1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.wt.dsvm     <- array(0, dim=c(length(ages),endy + 1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.wt.dsvm.tot <- array(0, dim=c(1           ,endy + 1,              1,            1), dimnames=list(cat="all", year=as.character(1:(endy+1)), season="all",                option ="all"))

#run population for 15 year
pop1 <- population_dynamics(pop=pop1, recruitarea="a", startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm[,1,,, drop=F], recruitment=100)

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
                                  fuelUse = 0.001, fuelPrice = 2.0, landingCosts= 0,gearMaintenance= 0, addNoFishing= TRUE, increments= 40,
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
  
  pop1 <- population_dynamics(pop=pop1, recruitarea = "a", startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm[,yy,,,drop=F], recruitment=100)
  pos_catches1 <- pop1 *q*wts
  
}



#what are the weights?
wts


hr <- apply(catches.n.dsvm,1:3,sum)/    (apply(catches.n.dsvm,1:3,sum) +    apply(pop1,1:3,sum) )


#what happens in our yield curve for this hr?
yield_curve(hr=hr[,64,], wts, natmortality, R=100, sequence = 1, verbose=T)

round(catches.n.dsvm[,64,,],2)
round(pop1[,64,,],2)


#next, what happens to theoretical pop for our harvest, what is harvest?
hr[,64,]
mean(hr[,64,])



yc <- yield_curve(hr=hr[,64,], wts, natmortality, R=100, verbose=F)

ylim <- c(0,900)
#to check
par(mfrow=c(1,2))
plot(catches.wt.dsvm.tot, type="l", ylim=ylim)

plot(x=yc$hr, y=yc$yield, ylim=ylim)
points(mean(hr[,64,]),yc$yield[yc$hr>mean(hr[,64,])][1], col="red", pch=19)
points(mean(hr[,62,]),catches.wt.dsvm.tot[,62,,], col="blue", pch=19)
points(mean(hr[,63,]),catches.wt.dsvm.tot[,63,,], col="blue", pch=19)
points(mean(hr[,64,]),catches.wt.dsvm.tot[,64,,], col="blue", pch=19)
points(mean(hr[,65,]),catches.wt.dsvm.tot[,65,,], col="blue", pch=19)
points(mean(hr[,66,]),catches.wt.dsvm.tot[,66,,], col="blue", pch=19)



