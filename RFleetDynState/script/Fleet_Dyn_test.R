library(RFleetDynState)
hke <-  mac <-  meg <-  hom <- mon <- new("DynStateInput")

catchMean(hke)  <- array(c(rep(10,2), rep(20,2),rep(10,2), rep(5,2)), dim=c(6,4,2), dimnames=list(cat=1:6,season=as.character(1:4),option=c("a","b")))
catchMean(hke)[,,"b"] <- catchMean(hke)[,,"b"]/2

catchSigma(hke) <- array(2                    ,dim=c(6,4,2), dimnames=dimnames(catchMean(hke)))

catchMean(mac)  <- array(0                      ,dim=c(6,4,2), dimnames=dimnames(catchMean(hke)))
catchSigma(mac) <- array(0.1                    ,dim=c(6,4,2), dimnames=dimnames(catchMean(hke)))

catchMean(mon) <- catchMean(hom)  <- catchMean(meg)  <- catchMean(mac)  
catchSigma(mon)  <- catchSigma(hom) <- catchSigma(meg) <- catchSigma(mac) 


effort <- array(c(rep(8,4)), dim=c(2,4), dimnames=list(option=c("a","b"), seasons=as.character(1:4)))

hkePrice <- array(c(8,10,12,16,20,22),dim=c(6,4), dimnames=list(cat=1:6, seasons=as.character(1:4)))
macPrice <- array(c(0),dim=c(6,4), dimnames=dimnames(hkePrice))
monPrice <- homPrice <- megPrice <- macPrice

##############################################################################
# test 1 unconstrained hake quotum, high hake price
###############################################################################

control <- DynState.control(spp1LndQuota= 600,  spp2LndQuota= 50,#tn
                            spp1LndQuotaFine= 1E8, spp2LndQuotaFine= 1E8, # tn
                            fuelUse = 1, 
                            spp1DiscardSteps=1,
                            fuelPrice = 12, # euros/day
                            gearMaintenance= 2, # 287 euros/day 
                            landingCosts= 1, # euros/tn
                            sigma=0.000001,
                            addNoFishing= TRUE, increments= 10, 
		            simNumber= 20, numThreads= 60)


z <- DynState(hke, mac, meg, hom, mon, hkePrice, macPrice, megPrice, homPrice, monPrice, effort, control) 

##############################################################################
# Check results
###############################################################################
layout (matrix(c(1, 4, 1, 4, 1, 5, 1, 5, 2, 6, 2, 6, 2, 7, 3, 7, 3, 8, 3, 8), 10, 2, byrow = TRUE))
plotChoice(hke,z)
hist(netRev(z),20,  main="Net Revenue",xlab = "", xlim= c(0,max(netRev(z))))
plot(apply(effort(sim(z))[,,],1,cumsum)[,1], type="s", ylim=c(0,max((apply(effort(sim(z))[,,],1,cumsum)[,1]))), main="Effort",xlab = "season",ylab = "")
for (ii in 2:dim(effort(sim(z)))[2]) lines(apply(effort(sim(z))[,,],1,cumsum)[,ii], type="s", ylim=c(0,50))
hist(apply(spp1Landings(sim(z)),2,sum),20, col="grey", main="Landings HKE",xlab = "", xlim= c(0, max(c(control@spp1LndQuota),max(apply(spp1Landings(sim(z)),2,sum)))))
abline(v=control@spp1LndQuota, lty=2)
hist(apply(spp2Landings(sim(z)),2,sum),20, col="grey", main="Landings MAC",xlab = "", xlim= c(0,60))
abline(v=control@spp2LndQuota, lty=2)
hist(apply(spp3Landings(sim(z)),2,sum),20, col="grey", main="Landings MEG",xlab = "", xlim= c(0,60))
hist(apply(spp4Landings(sim(z)),2,sum),20, col="grey", main="Landings HOM",xlab = "", xlim= c(0,60))
hist(apply(spp5Landings(sim(z)),2,sum),20, col="grey", main="Landings MON",xlab = "", xlim= c(0,60))

z@Sim@Derivative[,,]
