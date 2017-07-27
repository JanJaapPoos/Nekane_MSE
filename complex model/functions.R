
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

catch_dataframe_to_array <- function(dsvm_result, ages, season, areas, stock, catch_option){
  
  # need to revisit this function if we want to subset landings or discards separetly
  dsvm_catch <- with(dsvm_result,data.frame(cat, season, option, spp, get(catch_option)))
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

##############################################################################
# EFFORT PLOT
##############################################################################

effort_plot_dsvm <- function(NVESSELS, effort_dsvm_res_allyrs, stab.model, economics_res_allyrs){
  
  # create a list to store the results
  efforts<- list()
  netrev<-  melt(economics_res_allyrs,, id.vars = "year", measure.vars = c("NetRev"))
  grossrev<-  melt(economics_res_allyrs,, id.vars = "year", measure.vars = c("Grossrev"))
  annualfine<-  melt(economics_res_allyrs,, id.vars = "year", measure.vars = c("Annualfine"))
  
  mypalette <-c("#808080","#CCCCCC","#D55E00")
  names(mypalette) <- c("a", "b", "Stay in port")
  
  for (n in (1: length(NVESSELS))){
    
    a<- subset(effort_dsvm_res_allyrs,(spp %in% "sp1"))
    a<- subset(a,(cat %in% 1))
    a$year<- factor(a$year, levels= 1:90)
    #a<- subset(a,(nvessels %in% NVESSELS[n]))
    #a<- within(a,year <- year- stab.model)
    a<- aggregate(cbind(landings.wt, discards.wt, catch.wt, effort,trip)~ spp+cat+option+year, FUN=sum, data=a)
    
    efforts[[n]]<-a
    
  }    
  
  efforts <- ldply(efforts, rbind)
  
  p <- ggplot(a, aes(x=as.factor(year), y=trip, fill=option)) + 
    geom_bar(stat="identity", position = "fill", colour="black")+
    #annotate("text", label =paste0("SIMNUMBER ",SIMNUMBER), x =84,y = - .04)+
    scale_y_continuous( labels = percent)+
    scale_x_discrete(breaks=seq(0,90, 2), labels = c(seq(0,90, 2)), drop=FALSE)+
    geom_vline(xintercept = MPstart, linetype=2, color = "black", size=2)+
    annotate("text", label ="MP", x = MPstart+2, y = - .04)+
    scale_fill_manual(name= "option",values=mypalette,breaks=c("a", "b", "Stay in port"))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(),legend.title=element_blank(),
          legend.position="bottom",axis.text=element_text(size=8),
          strip.text = element_text(size = 8),
          text = element_text(size=10))+
    guides(fill= guide_legend(nrow=1, byrow=TRUE))+
    xlab("model time") +
    ylab("Effort pattern")
  
  q<- ggplot()+
    stat_summary(data= netrev, aes(year, value), geom="ribbon", fun.ymin = function(x) quantile(x, 0.025), fun.ymax = function(x) quantile(x, 0.975), alpha=0.2)+
    stat_summary(data= netrev, aes(year, value), geom="line", size= 0.25, fun.y=mean)+
    stat_summary(data= netrev, aes(year, value), geom="point", size= 1, , fun.y=mean)+
    scale_x_discrete(breaks=seq(0,90, 2), labels = c(seq(0,90, 2)), drop=FALSE)+
    ylim(0, 4000)+
    geom_vline(xintercept = MPstart, linetype=2, color = "black", size=2)+
    annotate("text", label ="MP", x = MPstart+2, y = - .04)+
     theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(),legend.title=element_blank(),
          legend.position="bottom",axis.text=element_text(size=8),
          strip.text = element_text(size = 8),
          text = element_text(size=10))+
    guides(fill= guide_legend(nrow=1, byrow=TRUE))+
    ylab("Net Revenue")
    
    r<- ggplot()+
    stat_summary(data= grossrev, aes(year, value), geom="ribbon", fun.ymin = function(x) quantile(x, 0.025), fun.ymax = function(x) quantile(x, 0.975), alpha=0.2)+
    stat_summary(data=grossrev, aes(year, value), geom="line", size= 0.25, fun.y=mean)+
    stat_summary(data= grossrev, aes(year, value), geom="point", size= 1, , fun.y=mean)+
    scale_x_discrete(breaks=seq(0,90, 2), labels = c(seq(0,90, 2)), drop=FALSE)+
    ylim(0, 4000)+
    geom_vline(xintercept = MPstart, linetype=2, color = "black", size=2)+
    annotate("text", label ="MP", x = MPstart+2, y = - .04)+
     theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(),legend.title=element_blank(),
          legend.position="bottom",axis.text=element_text(size=8),
          strip.text = element_text(size = 8),
          text = element_text(size=10))+
    guides(fill= guide_legend(nrow=1, byrow=TRUE))+
    ylab("Gross Revenue")
    
    s<- ggplot()+
    stat_summary(data= annualfine, aes(year, value), geom="ribbon", fun.ymin = function(x) quantile(x, 0.025), fun.ymax = function(x) quantile(x, 0.975), alpha=0.2)+
    stat_summary(data= annualfine, aes(year, value), geom="line", size= 0.25, fun.y=mean)+
    stat_summary(data= annualfine, aes(year, value), geom="point", size= 1, , fun.y=mean)+
    scale_x_discrete(breaks=seq(0,90, 2), labels = c(seq(0,90, 2)), drop=FALSE)+
    ylim(0, 4000)+
    geom_vline(xintercept = MPstart, linetype=2, color = "black", size=2)+
    annotate("text", label ="MP", x = MPstart+2, y = - .04)+
     theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(),legend.title=element_blank(),
          legend.position="bottom",axis.text=element_text(size=8),
          strip.text = element_text(size = 8),
          text = element_text(size=10))+
    guides(fill= guide_legend(nrow=1, byrow=TRUE))+
    ylab("Annual Fine")
    
    
    
  
  return(grid.arrange(p, q, r, s, ncol=1))
  
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
