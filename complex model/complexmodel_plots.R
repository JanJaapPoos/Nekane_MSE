###############################################################################
# Figures from results
##############################################################################
library(fields)
#We need to check if legend is correct (i.e. high and low abundances are correctly coloured for the two species.
#setwd("~/Dropbox/BoB/MSE/Git/Nekane/doc/Figures")
setwd("~/Nekane_MSE/doc/Figures")

#DISTRIBUTIONS
##############################################################################
postscript("distributions.eps")
par(oma=c(0,0,0,0), mar=c(4.1, 4.1, 3.1, 1.1))
split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
split.screen(c(2,1), screen=1)-> ind
screen( ind[1])

image(matrix(round(aperm(pop1[,5,,],c(2,1,3))), ncol=2), col=gray((0:64)/64), axes=F, main= "Species 1", ylab="Area" )
box()
axis(2,at= seq(0.,1,1),labels= c("S","N"), las=1)
axis(1,at= seq(0.,0.833,0.166),labels= seq(1:max(ages)))

screen( ind[2])

image(matrix(round(aperm(pop2[,5,,],c(2,1,3))), ncol=2), col=gray((0:64)/64), axes=F, main="Species 2", xlab="Age (years)", ylab="Area")
box()
axis(2,at= seq(0.,1,1),labels= c("S","N"), las=1)
axis(1,at= seq(0.,0.833,0.166),labels= seq(1:max(ages)))

#image.plot(matrix(round(aperm(pop1[,5,,],c(2,1,3))), ncol=2), col=gray((0:64)/64),  main= "Species 1" )
screen(2)
image.plot(legend.only=T,zlim=c(0,400),  col=gray((0:64)/64),  smallplot=c(.2,.4, .3,.7))
close.screen( all=TRUE)
dev.off()


#PRICES
##############################################################################
postscript("prices.eps")
par(mar = c(5,5,2,5))
plot(x=seq(min(ages)+((1/max(season))/2), max(ages+1),1/max(season)),y=sort(c(sp1Price)),
     xlab="Age (year)", ylab= "Price (euro per ton)", ylim=c(0,max(sp1Price)), type="l", lty=1, las=1)
par(new = T)
plot(x=seq(min(ages)+((1/max(season))/2), max(ages+1),1/max(season)),
y=c(lens),ylim=c(0,max(lens)), axes=F, xlab=NA,ylab=NA, type="l",lty=2)
axis(side = 4)
mtext(side = 4, line = 3, 'Fish length (cm)')
legend("topleft",
legend=c("Price", "Length"),
lty=c(1,2), pch=c(NA, NA), col=c("black", "black"))
dev.off()


#POPULATION AND CATCH COHORTS
#by sizing points according the value of a variable
##############################################################################
postscript("Popcatch.eps")

data           <- melt(pop1)
levels(data$option) <- c(levels(data$option),"North", "South")
data$option[data$option=="a"] <- c("North")
data$option[data$option=="b"] <- c("South")
data$time      <- paste0(data$year, sep=".", data$season)
time           <-levels(as.factor(as.numeric(data$time)))
data$years     <- match(data$time,time)
data$ variable <- "Population numbers"

data1           <- melt(sweep( catches.wt.dsvm1 ,c(1,3),  wts , FUN="/"))
levels(data1$option) <- c(levels(data1$option),"North", "South")
data1$option[data1$option=="a"] <- c("North")
data1$option[data1$option=="b"] <- c("South")
data1$time      <- paste0(data1$year, sep=".", data1$season)
data1$years     <- match(data1$time,time)
data1$ variable <- "Catch numbers"

data<- rbind(data, data1)
rm(data1)
#label for x axis
lab <-c(seq(from = 0, to = stab.model+NUMRUNS+1, by = 5))
lab[1] <-1
lab[which(lab %in% stab.model)]<- "Stab.model"
lab[which(lab %in% MPstart)]<- "MP" 

brk<-c(seq(from = 1, to = length(time)+max(season), by = max(season)*5))
brk[-1]<- (brk-max(season))[-1]

ggplot(data, aes(years, as.factor(cat), size = value)) +
     geom_point(shape = 21, aes(colour = factor(sign(value)),fill = factor(sign(value)))) + scale_size(range = c(1, 15)) + scale_colour_manual(values = c("black","white")) + 
     scale_fill_manual(values = c("lightgray", "black"))+ 
     facet_grid(option~ variable)+
     geom_vline(xintercept=(MPstart*max(season))+1,linetype=2, color="red")+ 
     geom_vline(xintercept=(stab.model*max(season))-(max(season)-1),linetype=2, color="red")+
     theme(legend.position = "none", panel.background = element_blank()) + ylab("age")+
     scale_x_continuous(name="year", limits=c(60, length(time)),breaks=brk, labels=lab)
      
dev.off()


#CATCHES AND YIELD
##############################################################################
postscript(file=paste("MIXEDMP2_CATCH_", paste0("SIGMA ", SIGMA, "; INCREMENTS ",control@increments, "; SIMNUMBER ",SIMNUMBER, "; DISCARDSTEPS ",SPP1DSCSTEPS, "; MIGRATION ",migconstant, "; REC1 ",paste(recs1, collapse = " "), "; REC2 ",paste(recs2, collapse = " "), "; SP1PRICE ",paste0(round(sp1Price[,1]), sep = ', ', collapse = ''),";SP2PRICE ",paste0(round(sp2Price[,1]), sep = ', ', collapse = ''), "; FUELPRICE ",control@fuelPrice,".eps"), sep=""),width=20, height=8, units="cm", res=500, pointsize=6.5)

par(mfrow=c(2,5),oma = c(3,0,0,0) + 0.1, mar = c(4,4,1,1) + 0.1)

ylim <- c(0, ceiling(max(sort(catches.wt.dsvm.tot1,partial=length(catches.wt.dsvm.tot1)-1)[length(catches.wt.dsvm.tot1)-1])/100)*100)
xlimYPR <- c(0,0.3)

plot(catches.wt.dsvm.tot1, type="l",  xlim=c(stab.model,endy), ylim=ylim, xaxs='i', yaxs='i', xlab= "Year", ylab = "Total catches (weight)", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines((quota1* SIMNUMBER), col="red")
lines(catches.wt.dsvm.tot1,  type="l", ylim=ylim)
legend("bottomright", inset=.05, legend=c("Catches","TAC"), pch=c(1,46,1), col=c("black","red"), bty='n', cex=0.8)

#harvest rate plots
plot(apply(hr1,c(1,2),mean)[1,], type="p",  xlim=c(stab.model,endy), ylim=xlimYPR, xaxs='i', yaxs='i', xlab= "Year", ylab = "Harvest rates", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
abline(h=c(0.05, 0.1, 0.15), lty="dotted", col = "ivory2")
points(apply(hr1,c(1,2),mean)[1,], col="black", pch=19)
points(apply(hr1,c(1,2),mean)[2,], col="red", pch=19)
points(apply(hr1,c(1,2),mean)[3,], col="black", pch=21)
points(apply(hr1,c(1,2),mean)[4,], col="red", pch=21)
points(apply(hr1,c(1,2),mean)[5,], col="black", pch=23)
points(apply(hr1,c(1,2),mean)[6,], col="red", pch=23)
lines(hr1wanted[1,,1,1], col="red")    #wanted hr
lines(apply(hr1,c(2),mean), col="black") #observed hr

legend("topright", inset=.05, legend=c("Age 1","Age 2","Age 3","Age 4",,"Age 5","Age 6"), pch=c(19,19,21,21,23,23), col=c("black","red", "black","red", "black","red"), bty='n', cex=0.8)

#YPR plots
plot(x=yc1noMP$hr, y=yc1noMP$landings, type="l", xlim=xlimYPR, ylim=ylim,xaxs='i', yaxs='i',  xlab="Harvest rate", ylab = "Yield per recruit", panel.first=grid(col = "ivory2"))
abline(v=Fmsy1noMP, col="ivory4")
#text(xlimYPR[2]*0.8, ylim[2]*0.9, paste0("SIMNUMBER ",SIMNUMBER))
points(mean(hr1[,pyrnoMP,]),yc1noMP$landings[yc1noMP$hr>mean(hr1[,pyrnoMP,])][1], col="red", pch=19) # current hr on yield curve 
points(mean(hr1[,pyrnoMP-2,]),landings.wt.dsvm.tot1[,pyrnoMP-2,,],col="ivory4", pch=19)
points(mean(hr1[,pyrnoMP-1,]),landings.wt.dsvm.tot1[,pyrnoMP-1,,], col="ivory4", pch=19)
points(mean(hr1[,pyrnoMP,]),landings.wt.dsvm.tot1[,pyrnoMP,,], col="ivory4", pch=19)
points(mean(hr1[,pyrnoMP+1,]),landings.wt.dsvm.tot1[,pyrnoMP+1,,], col="ivory4", pch=19)
points(mean(hr1[,pyrnoMP+2,]),landings.wt.dsvm.tot1[,pyrnoMP+2,,], col="ivory4", pch=19)
lines(x=yc1MP$hr, y=yc1MP$landings, ylim=ylim, col="ivory3")
text(xlimYPR[2]*0.8, yc1MP$landings[length(yc1noMP$hr)/3], "Unconstrained")
text(xlimYPR[2]*0.8, yc1MP$landings[length(yc1MP$hr)/3], "Constrained")
text(1,max(rowMeans(hr1[,pyrMP,]))+0.01, "Constrained", pos=4)
abline(v=Fmsy1MP, col="ivory3")
points(mean(hr1[,pyrMP,]),yc1MP$landings[yc1MP$hr>mean(hr1[,pyrMP,])][1], col="red", pch=21, bg="white")
points(mean(hr1[,pyrMP-2,]),landings.wt.dsvm.tot1[,pyrMP-2,,], col="ivory3", pch=21, bg="white")
points(mean(hr1[,pyrMP-1,]),landings.wt.dsvm.tot1[,pyrMP-1,,], col="ivory3", pch=21, bg="white")
points(mean(hr1[,pyrMP,]),landings.wt.dsvm.tot1[,pyrMP,,], col="ivory3", pch=21, bg="white")
points(mean(hr1[,pyrMP+1,]),landings.wt.dsvm.tot1[,pyrMP+1,,], col="ivory3", pch=21, bg="white")
points(mean(hr1[,pyrMP+2,]),landings.wt.dsvm.tot1[,pyrMP+2,,], col="ivory3", pch=21, bg="white")

plot(rowMeans(hr1[,pyrnoMP,]), type="b", ylim= xlimYPR,  xlab="Age", ylab = "Selectivity", panel.first=grid(col = "ivory2"), xaxt="n")
text(1,max(rowMeans(hr1[,pyrnoMP,]))+0.01, "Unconstrained", pos=4)
lines(rowMeans(hr1[,pyrMP,]), type="b", ylim= xlimYPR, col="grey")
text(1,max(rowMeans(hr1[,pyrMP,]))+0.01, "Constrained", pos=4)
axis(1, at = seq(1, 6, by = 1))

plot(apply(catches.wt.dsvm1,c(2,4),sum)[,1], col="blue", type="l",  ylim=ylim,xaxs='i', yaxs='i', xlim=c(stab.model,endy),xlab = "Year", ylab = "Catches by area (weight)",xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(apply(catches.wt.dsvm1,c(2,4),sum)[,1], col="blue", type="l")
lines(apply(catches.wt.dsvm1,c(2,4),sum)[,2], col="red")
lines(apply(landings.wt.dsvm1,c(2,4),sum)[,1], col="blue", lty=2)
lines(apply(landings.wt.dsvm1,c(2,4),sum)[,2], col="red", lty=2)
#lines(apply(catches.wt.dsvm1,c(2,4),sum)[,3], col="black")
legend("topright", inset=.05, c("North","South"),  pch=c(1,1), col=c("blue","red"), bty='n', cex=0.8)
abline(v=MPstart, lty=2)


#round(pop1,0)

plot(catches.wt.dsvm.tot2, type="l",  xlim=c(stab.model,endy), ylim=ylim, xaxs='i', yaxs='i', xlab= "Year", ylab = "Total catches (weight)", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines((quota2* SIMNUMBER), col="red")
lines(catches.wt.dsvm.tot2,  type="l", ylim=ylim)
lines(landings.wt.dsvm.tot2, type="l", ylim=ylim, lty= 2)

plot(apply(hr2,c(1,2),mean)[1,], type="p",  xlim=c(stab.model,endy), ylim=xlimYPR, xaxs='i', yaxs='i', xlab= "Year", ylab = "Harvest rates", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(hr2wanted[1,,1,1], col="red")    #wanted hr
points(apply(hr2,c(1,2),mean)[1,], col="black", pch=19)
points(apply(hr2,c(1,2),mean)[2,], col="red", pch=19)
points(apply(hr2,c(1,2),mean)[3,], col="black", pch=21)
points(apply(hr2,c(1,2),mean)[4,], col="red", pch=21)
points(apply(hr2,c(1,2),mean)[5,], col="black", pch=23)
points(apply(hr2,c(1,2),mean)[6,], col="red", pch=23)


#YPR
plot(x=yc2noMP$hr, y=yc2noMP$landings, type="l", xlim=xlimYPR, ylim=ylim,xaxs='i', yaxs='i', xlab = "Harvest rate", ylab = "Yield per recruit", panel.first=grid(col = "ivory2"))
abline(v=Fmsy2noMP)
points(mean(hr2[,pyrnoMP,]),yc2noMP$landings[yc2noMP$hr>mean(hr2[,pyrnoMP,])][1], col="red", pch=19)
points(mean(hr2[,pyrnoMP-2,]),landings.wt.dsvm.tot2[,pyrnoMP-2,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP-1,]),landings.wt.dsvm.tot2[,pyrnoMP-1,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP,]),landings.wt.dsvm.tot2[,pyrnoMP,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP+1,]),landings.wt.dsvm.tot2[,pyrnoMP+1,,], col="blue", pch=19)
points(mean(hr2[,pyrnoMP+2,]),landings.wt.dsvm.tot2[,pyrnoMP+2,,], col="blue", pch=19)
lines(x=yc2MP$hr, y=yc2MP$landings, ylim=ylim, col="grey")
text(xlimYPR[2]*0.8, yc2MP$landings[length(yc2noMP$hr)/3], "Unconstrained")
text(xlimYPR[2]*0.8, yc2MP$landings[length(yc2MP$hr)/3], "Constrained")
abline(v=Fmsy2MP, col="grey")
points(mean(hr2[,pyrMP,]),yc2MP$landings[yc2MP$hr>mean(hr2[,pyrMP,])][1], col="red", pch=21, bg="white")
points(mean(hr2[,pyrMP-2,]),landings.wt.dsvm.tot2[,pyrMP-2,,], col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP-1,]),landings.wt.dsvm.tot2[,pyrMP-1,,], col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP,]),landings.wt.dsvm.tot2[,pyrMP,,], col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP+1,]),landings.wt.dsvm.tot2[,pyrMP+1,,], col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP+2,]),landings.wt.dsvm.tot2[,pyrMP+2,,], col="blue", pch=21, bg="white")

plot(rowMeans(hr2[,pyrnoMP,]), type="b", ylim= xlimYPR,  xlab="Age", ylab = "Selectivity", panel.first=grid(col = "ivory2"), xaxt="n")
text(1,max(rowMeans(hr2[,pyrnoMP,]))+0.01, "Unconstrained", pos=4)
lines(rowMeans(hr2[,pyrMP,]), type="b", ylim=  xlimYPR, col="grey")
text(1,max(rowMeans(hr2[,pyrMP,]))+0.01, "Constrained", pos=4)
axis(1, at = seq(1, 6, by = 1))

plot(apply(catches.wt.dsvm2,c(2,4),sum)[,1], col="blue", type="l",  ylim=ylim,xaxs='i', yaxs='i', xlim=c(stab.model,endy),xlab = "Year", ylab = "Catches by area (weight)",xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="ivory4")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkgreen")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(apply(catches.wt.dsvm2,c(2,4),sum)[,1], col="blue", type="l")
lines(apply(catches.wt.dsvm2,c(2,4),sum)[,2], col="red")
lines(apply(landings.wt.dsvm2,c(2,4),sum)[,1], col="blue", lty=2)
lines(apply(landings.wt.dsvm2,c(2,4),sum)[,2], col="red", lty=2)


add_legend("bottomright", legend=paste0("SIGMA ", SIGMA, "; INCREMENTS ",control@increments, "; SIMNUMBER ",SIMNUMBER, "; DISCARDSTEPS1 ",SPP1DSCSTEPS,"; DISCARDSTEPS2 ",SPP2DSCSTEPS, "; MIGRATION ",migconstant, "; REC1 ",paste(recs1, collapse = " "), "; REC2 ",paste(recs2, collapse = " "), "; SP1PRICE ",paste0(round(sp1Price[,1]), sep = ', ', collapse = ''),";SP2PRICE ",paste0(round(sp2Price[,1]), sep = ', ', collapse = ''), "; FUELPRICE ",control@fuelPrice), col="black", horiz=TRUE, bty='n', cex=1)

dev.off()


#EFFORT AND ECONOMICS
##############################################################################
postscript(file=paste("MIXEDMP2_EFFORT_", paste0("SIGMA ", SIGMA, "; INCREMENTS ",control@increments,"; SIMNUMBER ",SIMNUMBER, "; DISCARDSTEPS ",SPP1DSCSTEPS, "; MIGRATION ",migconstant, "; REC1 ",paste(recs1, collapse = " "), "; REC2 ",paste(recs2, collapse = " "), "; SP1PRICE ",paste0(round(sp1Price[,1]), sep = ', ', collapse = ''),";SP2PRICE ",paste0(round(sp2Price[,1]), sep = ', ', collapse = ''), "; FUELPRICE ",control@fuelPrice,".eps"), sep=""),width=20, height=8, units="cm", res=500, pointsize=6.5)

par(mfrow=c(2,5),oma = c(3,0,0,0) + 0.1, mar = c(4,4,1,1) + 0.1)

mypalette <- c("#808080","#CCCCCC","#D55E00")
names(mypalette) <- c("North", "South", "Stay in port")
ylim <- c(0, ceiling(max(sort(netrev$value,partial=length(netrev$value)-1)[length(netrev$value)-1])/100)*100)


barplot(trip_percentage, col= mypalette, border=NA, xlim = c(1,NUMRUNS+1), xlab = "Year", ylab = "Effort pattern (%)", xaxt = "n",space = 0)
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
    labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
legend("topright", inset=.05, legend=c("North", "South", "Stay in port"), fill=mypalette, cex=0.6)
abline(v=MPstart-9.5, lty=2)

boxplot(value ~ year, netrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean net revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n")
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart-9.5, lty=2)
boxplot(value ~ year, netrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean net revenue", xlim=c(1,  NUMRUNS+1), ylim=ylim, xaxt = "n", add=TRUE)

boxplot(value ~ year, grossrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean gross revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n")
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart-9.5, lty=2)
boxplot(value ~ year, grossrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean gross revenue", xlim=c(1,NUMRUNS+1), ylim=ylim, xaxt = "n", add=TRUE)

boxplot(value ~ year, fuelcosts, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean fuel costs",  xlim=c(1,NUMRUNS+1), ylim=ylim, xaxt = "n")
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart-9.5, lty=2)
boxplot(value ~ year, fuelcosts, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean fuel costs", 
        ylim=ylim, xlim=c(1,NUMRUNS+1), xaxt = "n", add=TRUE)

boxplot(value ~ year, annualfine, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean annual fine", 
        ylim=ylim, xlim=c(1,NUMRUNS+1),  xaxt = "n")
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart-9.5, lty=2)
boxplot(value ~ year, annualfine, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean annual fine", 
        ylim=ylim, xlim=c(1,NUMRUNS+1), xaxt = "n", add=TRUE)

barplot(as.matrix(days), col=c("#808080","#CCCCCC"),border=NA, xlim = c(1,NUMRUNS+1), xlab = "Year", ylab = "Total days at sea", xaxt = "n",space = 0)
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart-9.5, lty=2)

ylim= c(0, ceiling(max(aggregate(value ~ year, FUN=sum, data=netrev)[-1,])/100)*100)
plot(aggregate(value ~ year, FUN=sum, data=netrev), type="l", ylim=ylim, xlab = "Year", ylab = "Total net revenue", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(aggregate(value ~ year, FUN=sum, data=netrev), type="l")

plot(aggregate(value ~ year, FUN=sum, data=grossrev), type="l", ylim=ylim, xlab = "Year", ylab = "Total gross revenue", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(aggregate(value ~ year, FUN=sum, data=grossrev), type="l")

plot(aggregate(value ~ year, FUN=sum, data=fuelcosts), type="l", ylim=ylim, xlab = "Year", ylab = "Total fuel costs", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(aggregate(value ~ year, FUN=sum, data=fuelcosts), type="l")

plot(aggregate(value ~ year, FUN=sum, data=annualfine),type="l", ylim=ylim, xlab = "Year", ylab = "Total annual fine", xaxt = "n", panel.first=grid(NA, NULL,col = "ivory2"))
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(aggregate(value ~ year, FUN=sum, data=annualfine), type="l")

add_legend("bottomright", legend=paste0("SIGMA ", SIGMA, "; SIMNUMBER ",SIMNUMBER, "; DISCARDSTEPS ",SPP1DSCSTEPS, "; MIGRATION ",migconstant, "; REC1 ",paste(recs1, collapse = " "), "; REC2 ",paste(recs2, collapse = " "), "; SP1PRICE ",paste0(round(sp1Price[,1]), sep = ', ', collapse = ''),";SP2PRICE ",paste0(round(sp2Price[,1]), sep = ', ', collapse = ''), "; FUELPRICE ",control@fuelPrice), col="black", horiz=TRUE, bty='n', cex=1)

dev.off()


