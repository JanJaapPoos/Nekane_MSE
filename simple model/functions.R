
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
  
  # Just focus on sp1 and sp2
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
  
  #subset the desired stock population catch in weight
  dsvm_catch_wt  <-subset(dsvm_catch,(spp %in% stock))[-c(4,5)]
  
  #trick to avoid deciding what to do with staying at port, it is just to sum catches
  dsvm_catch_wt  <- dsvm_catch_wt[dsvm_catch_wt$option!="Stay in port",] 
  
  dim            <- c(length(ages),length(season),length(areas)) 
  dimnames       <- list(cat=ages,season=as.character(season),option =areas)
  empty.df       <-melt(array(c(0), dim=dim, dimnames=dimnames))
  
  catch          <- merge(empty.df,dsvm_catch_wt,by=c("cat", "season","option"),all.x=TRUE)
  
  catch          <- catch[-4]
  catch[is.na(catch)] <- 0
  names(catch)   <- c("cat", "season","option","data")
  
  catch          <- tapply(catch[,"data"], list(factor(x = catch[,"cat"], levels = unique(catch[,"cat"])), 
                                                factor(x = catch[,"season"], levels = unique(catch[,"season"])), 
                                                factor(x = catch[,"option"], levels = unique(areas))), sum)
  return(catch)
}
