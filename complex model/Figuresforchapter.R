###############################################################################
# Figures from results
##############################################################################
library(fields)
require(graphics)

setwd("~/Dropbox/BoB/MSE/Git/Nekane_MSE/doc/Figures")

#DISTRIBUTIONS
##############################################################################
postscript(file="Distributions.eps", onefile=FALSE, horizontal=FALSE,width=8.5,height=3.25)
load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/1quota12seasons.RData")
par(oma=c(0,0,0,0), mar=c(4.1, 4.5, 1.1, 1.1))
split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))

split.screen(c(1,2), screen=1)-> ind
screen( ind[1])
image(matrix(round(aperm(pop1[,(stab.model-1),,],c(2,1,3))), ncol=2)[ , c(2,1)], col=gray((0:72)/72), axes=F, main= "", ylab="" )
box()
mtext("(a)", side=1, line = -1.3, adj = 0.96, font=2, cex = 1)
mtext(text="Area"  ,side=2,line=3.5    , adj = 0.50, font=2)

axis(2,at= seq(0.,1,1),labels= c("South","North"), las=1)
axis(1,at= seq(0.,1,0.1666),labels= seq(1:max(ages+1)))

screen( ind[2])
image(matrix(round(aperm(pop2[,(stab.model-1),,],c(2,1,3))), ncol=2)[ , c(2,1)], col=gray((0:72)/72), axes=F, main="", xlab="", ylab="")
box()
mtext("(b)", side=1, line = -1.3, adj = 0.96, font=2, cex = 1)
mtext(text="Age (year)"  ,side=1,line=2, adj = -0.60, font=2)
axis(2,at= seq(0.,1,1),labels= c("South","North"), las=1)
axis(1,at= seq(0.,1,0.1666),labels= seq(1:max(ages+1)))

#image.plot(matrix(round(aperm(pop1[,5,,],c(2,1,3))), ncol=2), col=gray((0:64)/64),  main= "Species 1" )
screen(2)
image.plot(legend.only=T,zlim=c(0,500),  col=gray((0:72)/72), smallplot=c(.15,.35, .4,.7),
           legend.args=list(text='  Abundance \n (numbers)', side=3, font=2, line=1, adj = 0.2, cex=0.8))
close.screen( all=TRUE)
dev.off()

#PRICES and FISH LENGTHS
##############################################################################
postscript(file="weights.eps", onefile=FALSE, horizontal=FALSE,width=4.25,height=3.25)
par(mar = c(3,5,2,4))
plot(x=seq(min(ages)+((1/max(season))/2), max(ages+1),1/max(season)),y=sort(c(sp1Price/1000)), xlab="", ylab= "", ylim=c(0,8), type="l", lty=1, las=1)
mtext(side = 2, line = 3.5, 'Price (euro per kg)', font=2, cex = 1)
mtext(side = 1, line = 2, 'Age (year)', font=2, cex = 1)
par(new = T)
plot(x=seq(min(ages)+((1/max(season))/2), max(ages+1),1/max(season)),
     y=c(lens),ylim=c(0,max(lens)), axes=F, xlab=NA,ylab=NA, type="l",lty=2)
axis(side = 4, las=1)
mtext(side = 4, line = 2.5, 'Fish length (cm)', font=2, cex = 1)
legend("bottomright", inset=.05, legend=c("Price", "Length"), pch=c(NA, NA), lty=c(1,2), lwd=c(1,1),col=c("black","black"), bty='n',text.font =1, cex=0.7)
dev.off()

#FISH WEIGHTS
##############################################################################
postscript(file="Prices.eps", onefile=FALSE, horizontal=FALSE,width=4.25,height=3.25)
par(mar = c(3,5,2,4))
plot(x=seq(min(ages)+((1/max(season))/2), max(ages+1),1/max(season)),y=sort(wts[,,,]), xlab="", ylab= "", ylim=c(0,3), type="l", lty=1, las=1)
mtext(side = 2, line = 3.5, 'Fish weight (kg)', font=2, cex = 1)
mtext(side = 1, line = 2, 'Age (year)', font=2, cex = 1)
dev.off()

#CATCHES
##############################################################################
postscript(file="Catches.eps", onefile=FALSE, horizontal=FALSE,width=8.5,height=6.5)
ylim=c(0,8)
par(par(no.readonly=TRUE))
par(oma=c(3,3,0,3),mar=c(3,3,2,3),mfrow=c(2,2))

load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/1quota12seasons.RData")
plot(catches.wt.dsvm.tot1/1000, type="l",  xlim=c(stab.model,endy), ylim=ylim, xaxs='i', yaxs='i', xlab="", ylab="", xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart, lty=2)
q1<- quota1
q1[1:25]<-NA
lines((q1* SIMNUMBER)/1000, col="red")
lines(catches.wt.dsvm.tot1/1000,  type="l", ylim=ylim)
mtext("(a)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1)
par(new = T)
plot(apply(hr1,c(2),mean), type="p",  xlim=c(stab.model,endy), ylim=c(0,0.08), axes=F, xlab=NA,ylab=NA, col="black", pch=20)
points(hr1wanted, col="blue", pch="*") #observed hr
axis(side = 4, las=1)

plot(catches.wt.dsvm.tot2/1000, type="l",  xlim=c(stab.model,endy), ylim=ylim, xaxs="i", yaxs="i", xlab="", ylab="", xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart, lty=2)
lines(catches.wt.dsvm.tot2/1000,  type="l", ylim=ylim)
mtext("(b)",  side=1, line = -1.3, adj = 0.04, font=2, cex = 1)
par(new = T)
plot(apply(hr2,c(2),mean), type="p",  xlim=c(stab.model,endy), ylim=c(0,0.08), axes=F, xlab=NA,ylab=NA, col="black", pch=20)
points(hr2wanted, col="blue", pch="*") #observed hr
axis(side = 4, las=1)
e <- bquote(expression(H[.("max")]))
legend("bottomright", inset=c(0,0), legend=c("Catches","TAC", "observed H", eval(e) ),lty = c(1,1,NA,NA), lwd = c(2,2,NA,NA),pch = c(NA,NA,20,8),
       pt.bg = c(NA,NA,"black","blue"),col=c("black","red","black","blue"), bty='n',text.font =1, cex=0.7)

load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/2quotas12seasons.RData")

plot(catches.wt.dsvm.tot1/1000, type="l",  xlim=c(stab.model,endy), ylim=ylim, xaxs='i', yaxs='i', xlab="", ylab="", xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart, lty=2)
q1<- quota1
q1[1:25]<-NA
lines((q1* SIMNUMBER)/1000, col="red")
lines(catches.wt.dsvm.tot1/1000,  type="l", ylim=ylim)
mtext("(c)",  side=1, line = -1.3, adj = 0.04, font=2, cex = 1)
par(new = T)
plot(apply(hr1,c(2),mean), type="p",  xlim=c(stab.model,endy), ylim=c(0,0.08), axes=F, xlab=NA,ylab=NA, col="black", pch=20)
points(hr1wanted, col="blue", pch="*") #observed hr
axis(side = 4, las=1)

plot(catches.wt.dsvm.tot2/1000, type="l",  xlim=c(stab.model,endy), ylim=ylim, xaxs="i", yaxs="i", xlab="", ylab="", xaxt = "n", panel.first=grid(col = "ivory3"))
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart, lty=2)
q2<- quota2
q2[1:25]<-NA
lines((q2* SIMNUMBER)/1000, col="red")
lines(catches.wt.dsvm.tot2/1000,  type="l", ylim=ylim)
mtext("(d)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1)
par(new = T)
plot(apply(hr2,c(2),mean), type="p",  xlim=c(stab.model,endy), ylim=c(0,0.08), axes=F, xlab=NA,ylab=NA, col="black", pch=20)
points(hr2wanted, col="blue", pch="*") #observed hr
axis(side = 4, las=1)

mtext(text="Year",side=1,line=0,font=2,outer=TRUE)
mtext(text="Total catches (tonnes)",side=2,line=0,font=2,outer=TRUE)
mtext(text="Mean harvest rate",side = 4, line=0,font=2,outer=TRUE)
dev.off()

#YIELDS
##############################################################################
postscript(file="Yields.eps", onefile=FALSE, horizontal=FALSE,width=8.5,height=6.5)
load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/1quota12seasons.RData")
ylim=c(0,8)
xlimYPR <- c(0,0.08)
par(par(no.readonly=TRUE))
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(2,2))

plot(x=yc1noMP$hr, y=yc1noMP$landings/1000, type="l", lwd=2, xlim=xlimYPR, ylim=ylim,xaxs='i', yaxs='i',  xlab="Harvest rate", ylab = "Yield per recruit", panel.first=grid(col = "ivory3"), las=1)
lines(x=yc1MP$hr, y=yc1MP$landings/1000, lwd=2, ylim=ylim, col="#808080")
abline(v=Fmsy1noMP, col="black",lwd=1, lty=2)
points(mean(hr1[,pyrnoMP-2,]),landings.wt.dsvm.tot1[,pyrnoMP-2,,]/1000,col="blue", pch=19)
points(mean(hr1[,pyrnoMP-1,]),landings.wt.dsvm.tot1[,pyrnoMP-1,,]/1000, col="blue", pch=19)
points(mean(hr1[,pyrnoMP,]),landings.wt.dsvm.tot1[,pyrnoMP,,]/1000, col="blue", pch=19)
points(mean(hr1[,pyrnoMP+1,]),landings.wt.dsvm.tot1[,pyrnoMP+1,,]/1000, col="blue", pch=19)
points(mean(hr1[,pyrnoMP+2,]),landings.wt.dsvm.tot1[,pyrnoMP+2,,]/1000, col="blue", pch=19)
points(mean(hr1[,pyrnoMP,]),yc1noMP$landings[yc1noMP$hr>mean(hr1[,pyrnoMP,])][1]/1000, col="red", pch=19) # current hr on yield curve 
text(xlimYPR[2]*0.8, (yc1noMP$landings[round(yc1noMP$hr,2)==0.06][1]/1000)-1.2, "Unconstrained")
text(xlimYPR[2]*0.8, (yc1MP$landings[round(yc1MP$hr,2)==0.06][1]/1000)+0.9, "Constrained")
#text(1,max(rowMeans(hr1[,pyrMP,]))+0.01, "Constrained", pos=4)
abline(v=Fmsy1MP, col="ivory3",lwd=1, lty=2)
points(mean(hr1[,pyrMP-2,]),landings.wt.dsvm.tot1[,pyrMP-2,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr1[,pyrMP-1,]),landings.wt.dsvm.tot1[,pyrMP-1,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr1[,pyrMP,]),landings.wt.dsvm.tot1[,pyrMP,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr1[,pyrMP+1,]),landings.wt.dsvm.tot1[,pyrMP+1,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr1[,pyrMP+2,]),landings.wt.dsvm.tot1[,pyrMP+2,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr1[,pyrMP,]),yc1MP$landings[yc1MP$hr>mean(hr1[,pyrMP,])][1]/1000, col="red", pch=21, bg="white")
mtext("(a)", side=3, line = -1.3, adj = 0.04, font=2, cex = 1)

plot(x=yc2noMP$hr, y=yc2noMP$landings/1000, type="l", lwd=2, xlim=xlimYPR, ylim=ylim,xaxs='i', yaxs='i', xlab = "Harvest rate", ylab = "Yield per recruit", panel.first=grid(col = "ivory3"), las=1)
lines(x=yc2MP$hr, y=yc2MP$landings/1000, lwd=2, ylim=ylim, col="#808080")
abline(v=Fmsy2noMP, col="black",lwd=1, lty=2)
points(mean(hr2[,pyrnoMP-2,]),landings.wt.dsvm.tot2[,pyrnoMP-2,,]/1000, col="blue", pch=19)
points(mean(hr2[,pyrnoMP-1,]),landings.wt.dsvm.tot2[,pyrnoMP-1,,]/1000, col="blue", pch=19)
points(mean(hr2[,pyrnoMP,]),landings.wt.dsvm.tot2[,pyrnoMP,,]/1000, col="blue", pch=19)
points(mean(hr2[,pyrnoMP+1,]),landings.wt.dsvm.tot2[,pyrnoMP+1,,]/1000, col="blue", pch=19)
points(mean(hr2[,pyrnoMP+2,]),landings.wt.dsvm.tot2[,pyrnoMP+2,,]/1000, col="blue", pch=19)
points(mean(hr2[,pyrnoMP,]),yc2noMP$landings[yc2noMP$hr>mean(hr2[,pyrnoMP,])][1]/1000, col="red", pch=19)
text(xlimYPR[2]*0.8, (yc2noMP$landings[round(yc2noMP$hr,2)==0.06][1]/1000)+0.9, "Unconstrained")
text(xlimYPR[2]*0.8, (yc2MP$landings[round(yc2MP$hr,2)==0.06][1]/1000)-1.2, "Constrained")
abline(v=Fmsy2MP, col="ivory3",lwd=1, lty=2)
points(mean(hr2[,pyrMP-2,]),landings.wt.dsvm.tot2[,pyrMP-2,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP-1,]),landings.wt.dsvm.tot2[,pyrMP-1,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP,]),landings.wt.dsvm.tot2[,pyrMP,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP+1,]),landings.wt.dsvm.tot2[,pyrMP+1,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP+2,]),landings.wt.dsvm.tot2[,pyrMP+2,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP,]),yc2MP$landings[yc2MP$hr>mean(hr2[,pyrMP,])][1]/1000, col="red", pch=21, bg="white")
mtext("(b)",  side=3, line = -1.3, adj = 0.04, font=2, cex = 1)

load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/2quotas12seasons.RData")

plot(x=yc1noMP$hr, y=yc1noMP$landings/1000, type="l", lwd=2, xlim=xlimYPR, ylim=ylim,xaxs='i', yaxs='i',  xlab="Harvest rate", ylab = "Yield per recruit", panel.first=grid(col = "ivory3"), las=1)
lines(x=yc1MP$hr, y=yc1MP$landings/1000, lwd=2, ylim=ylim, col="#808080")
abline(v=Fmsy1noMP, col="black",lwd=1, lty=2)
points(mean(hr1[,pyrnoMP-2,]),landings.wt.dsvm.tot1[,pyrnoMP-2,,]/1000,col="blue", pch=19)
points(mean(hr1[,pyrnoMP-1,]),landings.wt.dsvm.tot1[,pyrnoMP-1,,]/1000, col="blue", pch=19)
points(mean(hr1[,pyrnoMP,]),landings.wt.dsvm.tot1[,pyrnoMP,,]/1000, col="blue", pch=19)
points(mean(hr1[,pyrnoMP+1,]),landings.wt.dsvm.tot1[,pyrnoMP+1,,]/1000, col="blue", pch=19)
points(mean(hr1[,pyrnoMP+2,]),landings.wt.dsvm.tot1[,pyrnoMP+2,,]/1000, col="blue", pch=19)
points(mean(hr1[,pyrnoMP,]),yc1noMP$landings[yc1noMP$hr>mean(hr1[,pyrnoMP,])][1]/1000, col="red", pch=19) # current hr on yield curve 
text(xlimYPR[2]*0.8, (yc1noMP$landings[round(yc1noMP$hr,2)==0.06][1]/1000)+0.9, "Unconstrained")
text(xlimYPR[2]*0.8, (yc1MP$landings[round(yc1MP$hr,2)==0.06][1]/1000)-1.2, "Constrained")
#text(1,max(rowMeans(hr1[,pyrMP,]))+0.01, "Constrained", pos=4)
abline(v=Fmsy1MP, col="ivory3",lwd=1, lty=2)
points(mean(hr1[,pyrMP-2,]),landings.wt.dsvm.tot1[,pyrMP-2,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr1[,pyrMP-1,]),landings.wt.dsvm.tot1[,pyrMP-1,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr1[,pyrMP,]),landings.wt.dsvm.tAnnual net revenue by vessel (thousand euroot1[,pyrMP,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr1[,pyrMP+1,]),landings.wt.dsvm.tot1[,pyrMP+1,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr1[,pyrMP+2,]),landings.wt.dsvm.tot1[,pyrMP+2,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr1[,pyrMP,]),yc1MP$landings[yc1MP$hr>mean(hr1[,pyrMP,])][1]/1000, col="red", pch=21, bg="white")
mtext("(c)", side=3, line = -1.3, adj = 0.04, font=2, cex = 1)

plot(x=yc2noMP$hr, y=yc2noMP$landings/1000, type="l", lwd=2, xlim=xlimYPR, ylim=ylim,xaxs='i', yaxs='i', xlab = "Harvest rate", ylab = "Yield per recruit", panel.first=grid(col = "ivory3"), las=1)
lines(x=yc2MP$hr, y=yc2MP$landings/1000, lwd=2, ylim=ylim, col="#808080")
abline(v=Fmsy2noMP, col="black",lwd=1, lty=2)
points(mean(hr2[,pyrnoMP-2,]),landings.wt.dsvm.tot2[,pyrnoMP-2,,]/1000, col="blue", pch=19)
points(mean(hr2[,pyrnoMP-1,]),landings.wt.dsvm.tot2[,pyrnoMP-1,,]/1000, col="blue", pch=19)
points(mean(hr2[,pyrnoMP,]),landings.wt.dsvm.tot2[,pyrnoMP,,]/1000, col="blue", pch=19)
points(mean(hr2[,pyrnoMP+1,]),landings.wt.dsvm.tot2[,pyrnoMP+1,,]/1000, col="blue", pch=19)
points(mean(hr2[,pyrnoMP+2,]),landings.wt.dsvm.tot2[,pyrnoMP+2,,]/1000, col="blue", pch=19)
points(mean(hr2[,pyrnoMP,]),yc2noMP$landings[yc2noMP$hr>mean(hr2[,pyrnoMP,])][1]/1000, col="red", pch=19)
text(xlimYPR[2]*0.8, (yc2noMP$landings[round(yc2noMP$hr,2)==0.06][1]/1000)+0.9, "Unconstrained")
text(xlimYPR[2]*0.8, (yc2MP$landings[round(yc2MP$hr,2)==0.06][1]/1000)-1.2, "Constrained")
abline(v=Fmsy2MP, col="ivory3",lwd=1, lty=2)
points(mean(hr2[,pyrMP-2,]),landings.wt.dsvm.tot2[,pyrMP-2,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP-1,]),landings.wt.dsvm.tot2[,pyrMP-1,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP,]),landings.wt.dsvm.tot2[,pyrMP,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP+1,]),landings.wt.dsvm.tot2[,pyrMP+1,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP+2,]),landings.wt.dsvm.tot2[,pyrMP+2,,]/1000, col="blue", pch=21, bg="white")
points(mean(hr2[,pyrMP,]),yc2MP$landings[yc2MP$hr>mean(hr2[,pyrMP,])][1]/1000, col="red", pch=21, bg="white")
mtext("(d)",  side=3, line = -1.3, adj = 0.04, font=2, cex = 1)

mtext(text="Mean harvest rate (ages)",side=1,line=0,font=2,outer=TRUE)
mtext(text="Expected yield (tonnes)",side=2,line=0,font=2,outer=TRUE)
dev.off()

#SELECTIVITY
##############################################################################
postscript(file="Selectivity.eps", onefile=FALSE, horizontal=FALSE,width=8.5,height=6.5)
load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/1quota12seasons.RData")
ylim=c(0,8000)
xlimYPR <- c(0,0.08)
par(par(no.readonly=TRUE))
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(2,2))

plot(rowMeans(hr1[,pyrnoMP,]), type="b", ylim= xlimYPR,  xlab="Age", ylab = "Selectivity", panel.first=grid(col = "ivory3"), xaxt="n",las=1)
#rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="#CCCCCC")
#par(new=TRUE)
text(1,0.078, "Unconstrained", pos=4)
#lines(rowMeans(hr1[,pyrnoMP,]), type="b", ylim= xlimYPR)
lines(rowMeans(hr1[,pyrMP,]), type="b", ylim= xlimYPR, col="#808080")
text(1,max(rowMeans(hr1[,pyrMP,]))+0.01, "Constrained", pos=4)
axis(1, at = seq(1, 6, by = 1))
mtext("(a)", side=3, line = -1.3, adj = 0.96, font=2, cex = 1)

plot(rowMeans(hr2[,pyrnoMP,]), type="b", ylim= xlimYPR,  xlab="Age", ylab = "Selectivity", panel.first=grid(col = "ivory3"), xaxt="n",las=1)
#rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="darkseagreen3")
#par(new=TRUE)
text(1,0.078, "Unconstrained", pos=4)
#lines(rowMeans(hr2[,pyrnoMP,]), type="b", ylim= xlimYPR)
lines(rowMeans(hr2[,pyrMP,]), type="b", ylim=  xlimYPR, col="#808080")
text(1,max(rowMeans(hr2[,pyrMP,]))-0.02, "Constrained", pos=4)
axis(1, at = seq(1, 6, by = 1))
mtext("(b)",  side=3, line = -1.3, adj = 0.96, font=2, cex = 1)

load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/2quotas12seasons.RData")

plot(rowMeans(hr1[,pyrnoMP,]), type="b", ylim= xlimYPR,  xlab="Age", ylab = "Selectivity", panel.first=grid(col = "ivory3"), xaxt="n",las=1)
text(1,0.078, "Unconstrained", pos=4)
lines(rowMeans(hr1[,pyrMP,]), type="b", ylim= xlimYPR, col="#808080")
text(1,max(rowMeans(hr1[,pyrMP,]))+0.01, "Constrained", pos=4)
axis(1, at = seq(1, 6, by = 1))
mtext("(c)", side=3, line = -1.3, adj = 0.96, font=2, cex = 1)

plot(rowMeans(hr2[,pyrnoMP,]), type="b", ylim= xlimYPR,  xlab="Age", ylab = "Selectivity", panel.first=grid(col = "ivory3"), xaxt="n",las=1)
text(1,0.078, "Unconstrained", pos=4)
lines(rowMeans(hr2[,pyrMP,]), type="b", ylim=  xlimYPR, col="#808080")
text(1,max(rowMeans(hr2[,pyrMP,]))+0.01, "Constrained", pos=4)
axis(1, at = seq(1, 6, by = 1))
mtext("(d)",  side=3, line = -1.3, adj = 0.96, font=2, cex = 1)

mtext(text="Age (years)",side=1,line=0,font=2,outer=TRUE)
mtext(text="Mean harvest rate (month)",side=2,line=0,font=2,outer=TRUE)
dev.off()

#Effort- Economics and Mean Effort
##############################################################################
postscript(file="Efforteconomics.eps", onefile=FALSE, horizontal=FALSE,width=8.5,height=6.5)

load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/1quota12seasons.RData")
mypalette <- c("#808080","#CCCCCC","#D55E00")
names(mypalette) <- c("North", "South", "Stay in port")

#nf=layout(matrix(c(1,2,1,3,4,5,4,6), 4, 2, byrow = TRUE))
#par(oma = c(3,0,0,0) + 0.1, mar = c(4,5,1,1) + 0.1)
par(par(no.readonly=TRUE))
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(2,2))

barplot(trip_percentage, col= mypalette, border=NA, xlim = c(1,NUMRUNS+1), xlab = "", ylab = "", xaxt = "n",space = 0,las=1)
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(2, at =seq(0,100,20), label=seq(0,100,20),las=1)
legend("topright", inset=.05, legend=c("North", "South", "Stay in port"), fill=mypalette, cex=0.6)
abline(v=MPstart-9.5, lty=2)
box()
mtext("(a)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1) 
mtext(text="Effort pattern (%)",side=2, line=1,adj = 0.5,font=2,outer=TRUE)

a<- aggregate(value ~ year, FUN=sum, data=netrev)
a[2]<-a[2]/1000000
b<- aggregate(value ~ year, FUN=sum, data=grossrev)
b[2]<-b[2]/1000000
c<- aggregate(value ~ year, FUN=sum, data=fuelcosts)
c[2]<-c[2]/1000000
d<- aggregate(value ~ year, FUN=sum, data=annualfine)
d[2]<-d[2]/1000000
ylim= c(0, 350)
plot(a, type="l",  xlim=c(stab.model,endy), ylim=ylim, xaxs='i', yaxs='i', xlab="", ylab="", xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
lines(a, type="l")
lines(b, type="l", col="#808080")
lines(c, type="l",col="blue")
lines(d, type="l",col="red")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart, lty=2)
mtext("(b)",side=1, line = -1.3, adj = 0.04, font=2, cex = 1) 
mtext(text="Economics (thousand euro)",side=2, line=3,adj =0.5,at=-75,font=2)#side=2, line=3,adj =45,font=2,outer=TRUE)

load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/2quotas12seasons.RData")

barplot(trip_percentage, col= mypalette, border=NA, xlim = c(1,NUMRUNS+1), xlab = "", ylab = "", xaxt = "n",space = 0,las=1)
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(2, at =seq(0,100,20), label=seq(0,100,20),las=1)
legend("topright", inset=.05, legend=c("North", "South", "Stay in port"), fill=mypalette, cex=0.6)
abline(v=MPstart-9.5, lty=2)
box()
mtext("(c)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1)

a<- aggregate(value ~ year, FUN=sum, data=netrev)
b<- aggregate(value ~ year, FUN=sum, data=grossrev)
c<- aggregate(value ~ year, FUN=sum, data=fuelcosts)
d<- aggregate(value ~ year, FUN=sum, data=annualfine)
# maxval_mp<- b$value[16]
# a[2]<-a[2]/maxval_mp*100
# b[2]<-b[2]/maxval_mp*100
# c[2]<-c[2]/maxval_mp*100
# d[2]<-d[2]/maxval_mp*100
# ylim= c(0, 150)
# a[2]<-a[2]/1000000
a[2]<-a[2]/1000000
b[2]<-b[2]/1000000
c[2]<-c[2]/1000000
d[2]<-d[2]/1000000
ylim= c(0, 350)
plot(a, type="l",  xlim=c(stab.model,endy), ylim=ylim, xaxs='i', yaxs='i', xlab="", ylab="", xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
lines(a, type="l")
lines(b, type="l", col="#808080")
lines(c, type="l",col="blue")
lines(d, type="l",col="red")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart, lty=2)
mtext("(d)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1) 
legend("bottomright", inset=.01, legend=c("Net revenue","Gross revenue", "Fuel cost", "Annual fine"), pch=c(1,1,1,1), lwd=c(2,2,2,2),col=c("black","#808080","blue","red"), bty='n',text.font =2, cex=0.75)

mtext(text="Year",side=1,line=0,font=2,outer=TRUE)

dev.off()


postscript(file="Meaneffort.eps", onefile=FALSE, horizontal=FALSE,width=8.5,height=6.5)

load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/1quota12seasons.RData")
mypalette <- c("#808080","#CCCCCC")
names(mypalette) <- c("North", "South")

#nf=layout(matrix(c(1,2,1,3,4,5,4,6), 4, 2, byrow = TRUE))
#par(oma = c(3,0,0,0) + 0.1, mar = c(4,5,1,1) + 0.1)
par(par(no.readonly=TRUE))
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(2,2))
barplot(as.matrix(daysmen[1:2,])/1000, col=c("#808080","#CCCCCC"),xlim = c(0,max(season)), ylim=c(0, 7),xlab = "", ylab = "", xaxt = "n",space = 0, las=1)
axis( 1, at =seq(0.5,max(season)-0.5,1), labels = c(month.abb))
box()
mtext("(a)", side=3, line = -1.5, adj = 0.96, font=2, cex = 1)

barplot(as.matrix(daysmen[3:4,])/1000, col=c("#808080","#CCCCCC"), xlim = c(0,max(season)), ylim=c(0, 7),xlab = " ", ylab = "", xaxt = "n",space = 0, las=1)
axis( 1, at =seq(0.5,max(season)-0.5,1), labels = c(month.abb))
box()
mtext("(b)", side=3, line = -1.5, adj = 0.96, font=2, cex = 1)
legend(10,6, legend=c("North", "South"), fill=mypalette, cex=0.6)

load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/2quotas12seasons.RData")

barplot(as.matrix(daysmen[1:2,])/1000, col=c("#808080","#CCCCCC"),xlim = c(0,max(season)), ylim=c(0, 7),xlab = "", ylab = "", xaxt = "n",space = 0, las=1)
axis( 1, at =seq(0.5,max(season)-0.5,1), labels = c(month.abb))
box()
mtext("(c)", side=3, line = -1.5, adj = 0.96, font=2, cex = 1)

barplot(as.matrix(daysmen[3:4,])/1000, col=c("#808080","#CCCCCC"), xlim = c(0,max(season)), ylim=c(0, 7),xlab = "", ylab = "", xaxt = "n",space = 0, las=1)
axis( 1, at =seq(0.5,max(season)-0.5,1), labels = c(month.abb))
box()
mtext("(d)", side=3, line = -1.5, adj = 0.96, font=2, cex = 1)


mtext(text="Month",side=1,line=0,font=2,outer=TRUE)
mtext(text="Monthly fleet effort (thousand days)",side=2,line=0,font=2,outer=TRUE)

dev.off()

#CATCHES BY AREA
##############################################################################
postscript(file="Catchesbyarea.eps", onefile=FALSE, horizontal=FALSE,width=8.5,height=6.5)
load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/1quota12seasons.RData")
ylim=c(0,8)
par(par(no.readonly=TRUE))
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(2,2))

plot(apply(catches.wt.dsvm1/1000,c(2,4),sum)[,1], col="black", type="l",  ylim=ylim,xaxs='i', yaxs='i', xlim=c(stab.model,endy),xlab = "Year", ylab = "Catches by area (weight)",xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(apply(catches.wt.dsvm1/1000,c(2,4),sum)[,1], col="black", type="l")
lines(apply(catches.wt.dsvm1/1000,c(2,4),sum)[,2], col="#808080")
lines(apply(landings.wt.dsvm1/1000,c(2,4),sum)[,1], col="black", lty=2)
lines(apply(landings.wt.dsvm1/1000,c(2,4),sum)[,2], col="#808080", lty=2)
mtext("(a)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1)

plot(apply(catches.wt.dsvm2/1000,c(2,4),sum)[,1], col="black", type="l",  ylim=ylim,xaxs='i', yaxs='i', xlim=c(stab.model,endy),xlab = "Year", ylab = "Catches by area (weight)",xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(apply(catches.wt.dsvm2/1000,c(2,4),sum)[,1], col="black", type="l")
lines(apply(catches.wt.dsvm2/1000,c(2,4),sum)[,2], col="#808080")
lines(apply(landings.wt.dsvm2/1000,c(2,4),sum)[,1], col="black", lty=2)
lines(apply(landings.wt.dsvm2/1000,c(2,4),sum)[,2], col="#808080", lty=2)
mtext("(b)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1)
legend("topright", inset=.05, legend=c("North","South"),  pch=c(1,1), lwd=c(2.5,2.5), col=c("black","#808080"), bg="white",text.font =2, cex=0.6)

load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/2quotas12seasons.RData")

plot(apply(catches.wt.dsvm1/1000,c(2,4),sum)[,1], col="black", type="l",  ylim=ylim,xaxs='i', yaxs='i', xlim=c(stab.model,endy),xlab = "Year", ylab = "Catches by area (weight)",xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(apply(catches.wt.dsvm1/1000,c(2,4),sum)[,1], col="black", type="l")
lines(apply(catches.wt.dsvm1/1000,c(2,4),sum)[,2], col="#808080")
lines(apply(landings.wt.dsvm1/1000,c(2,4),sum)[,1], col="black", lty=2)
lines(apply(landings.wt.dsvm1/1000,c(2,4),sum)[,2], col="#808080", lty=2)
mtext("(c)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1)

plot(apply(catches.wt.dsvm2/1000,c(2,4),sum)[,1], col="black", type="l",  ylim=ylim,xaxs='i', yaxs='i', xlim=c(stab.model,endy),xlab = "Year", ylab = "Catches by area (weight)",xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(apply(catches.wt.dsvm2/1000,c(2,4),sum)[,1], col="black", type="l")
lines(apply(catches.wt.dsvm2/1000,c(2,4),sum)[,2], col="#808080")
lines(apply(landings.wt.dsvm2/1000,c(2,4),sum)[,1], col="black", lty=2)
lines(apply(landings.wt.dsvm2/1000,c(2,4),sum)[,2], col="#808080", lty=2)
mtext("(d)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1)

mtext(text="Year",side=1,line=0,font=2,outer=TRUE)
mtext(text="Catches by area (tonnes)",side=2,line=0,font=2,outer=TRUE)
dev.off()

#CPUE BY AREA
##############################################################################
postscript(file="CPUEbyarea.eps", onefile=FALSE, horizontal=FALSE,width=8.5,height=6.5)
load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/1quota12seasons.RData")
ylim=c(0,0.35)
par(par(no.readonly=TRUE))
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(2,2))

#Need to sum the catches for all categories but keep the effort as a mean
#xx<- aggregate(cbind(catch.wt, effort)~ spp+option+year+cat, FUN=sum, data=dsvm_res_allyrs)
xx<- aggregate(catch.wt~ spp+option+year, FUN=sum, data=dsvm_res_allyrs)
yy<- aggregate(effort~ spp+option+year+season, FUN=mean, data=dsvm_res_allyrs)
yy<- aggregate(effort~ spp+option+year, FUN=sum, data=yy)

# xx<- aggregate(catch.wt~ spp+option+year+season, FUN=sum, data=dsvm_res_allyrs)
# xx$effort<- aggregate(effort~ spp+option+year+season, FUN=mean, data=dsvm_res_allyrs)$effort
# xx$cpue<- xx$catch.wt/xx$effort
# yy<- aggregate(cpue~ spp+option+year, FUN=mean, data=xx)

xx$effort<- yy$effort
xx<- xx[!xx$option=="Stay in port",]
xx$cpue<- xx$catch.wt/xx$effort
# sp 1
xx1<- subset(xx,(spp %in% "sp1"))
# sp 2
xx2<- subset(xx,(spp %in% "sp2"))

plot(xx1$year[xx1$option=="a"],xx1$cpue[xx1$option=="a"],  col="black", type="l",  ylim=ylim, xaxs='i', yaxs='i', xlim=c(stab.model,endy),xlab = "Year", ylab = "CPUE by area",xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(xx1$year[xx1$option=="a"],xx1$cpue[xx1$option=="a"], col="black", type="l")
lines(xx1$year[xx1$option=="b"],xx1$cpue[xx1$option=="b"], col="#808080")
lines(xx1$year[xx1$option=="a"],xx1$cpue[xx1$option=="a"], col="black", lty=2)
lines(xx1$year[xx1$option=="b"],xx1$cpue[xx1$option=="b"], col="#808080", lty=2)
mtext("(a)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1)

plot(xx2$year[xx2$option=="a"],xx2$cpue[xx2$option=="a"],  col="black", type="l",  ylim=ylim, xaxs='i', yaxs='i', xlim=c(stab.model,endy),xlab = "Year", ylab = "CPUE by area",xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(xx2$year[xx2$option=="a"],xx2$cpue[xx2$option=="a"], col="black", type="l")
lines(xx2$year[xx2$option=="b"],xx2$cpue[xx2$option=="b"], col="#808080")
lines(xx2$year[xx2$option=="a"],xx2$cpue[xx2$option=="a"], col="black", lty=2)
lines(xx2$year[xx2$option=="b"],xx2$cpue[xx2$option=="b"], col="#808080", lty=2)
mtext("(b)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1)
legend("topright", inset=.05, legend=c("North","South"),  pch=c(1,1), lwd=c(2.5,2.5), col=c("black","#808080"), bg="white",text.font =2, cex=0.6)

load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/2quotas12seasons.RData")
#all categories together
xx<- aggregate(catch.wt~ spp+option+year, FUN=sum, data=dsvm_res_allyrs)
yy<- aggregate(effort~ spp+option+year+season, FUN=mean, data=dsvm_res_allyrs)
yy<- aggregate(effort~ spp+option+year, FUN=sum, data=yy)
xx$effort<- yy$effort
xx<- xx[!xx$option=="Stay in port",]
xx$cpue<- xx$catch.wt/xx$effort
# sp 1
xx1<- subset(xx,(spp %in% "sp1"))
# sp 2
xx2<- subset(xx,(spp %in% "sp2"))

plot(xx1$year[xx1$option=="a"],xx1$cpue[xx1$option=="a"],  col="black", type="l",  ylim=ylim, xaxs='i', yaxs='i', xlim=c(stab.model,endy),xlab = "Year", ylab = "CPUE by area",xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(xx1$year[xx1$option=="a"],xx1$cpue[xx1$option=="a"], col="black", type="l")
lines(xx1$year[xx1$option=="b"],xx1$cpue[xx1$option=="b"], col="#808080")
lines(xx1$year[xx1$option=="a"],xx1$cpue[xx1$option=="a"], col="black", lty=2)
lines(xx1$year[xx1$option=="b"],xx1$cpue[xx1$option=="b"], col="#808080", lty=2)
mtext("(c)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1)

plot(xx2$year[xx2$option=="a"],xx2$cpue[xx2$option=="a"],  col="black", type="l",  ylim=ylim, xaxs='i', yaxs='i', xlim=c(stab.model,endy),xlab = "Year", ylab = "CPUE by area",xaxt = "n", panel.first=grid(col = "ivory3"), las=1)
polygon(x=c(pyrnoMP-2,pyrnoMP+2,pyrnoMP+2,pyrnoMP-2), border=NA, y=c(rep(ylim,each=2)), col="#CCCCCC")
polygon(x=c(pyrMP-2,pyrMP+2,pyrMP+2,pyrMP-2)        , border=NA, y=c(rep(ylim,each=2)), col="darkseagreen3")
axis(1,  at = c(stab.model, stab.model+((MPstart-stab.model)/2), MPstart, MPstart+((endy-MPstart)/2), endy),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
axis(1, at = c(MPstart), labels = c("MP"))
abline(v=MPstart, lty=2)
lines(xx2$year[xx2$option=="a"],xx2$cpue[xx2$option=="a"], col="black", type="l")
lines(xx2$year[xx2$option=="b"],xx2$cpue[xx2$option=="b"], col="#808080")
lines(xx2$year[xx2$option=="a"],xx2$cpue[xx2$option=="a"], col="black", lty=2)
lines(xx2$year[xx2$option=="b"],xx2$cpue[xx2$option=="b"], col="#808080", lty=2)
mtext("(d)", side=1, line = -1.3, adj = 0.04, font=2, cex = 1)

mtext(text="Year",side=1,line=0,font=2,outer=TRUE)
mtext(text="CPUE by area (kg/vessel)",side=2,line=0,font=2,outer=TRUE)
dev.off()

#ECONOMICS
##############################################################################
postscript(file="Mean_Economics.eps", onefile=FALSE, horizontal=FALSE,width=8.5,height=6.5)
load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/1quota12seasons.RData")
#par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(2,3))
par(par(no.readonly=TRUE))
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(2,2))
ylim <- c(0, 50)
boxplot((value/1000)~ year, data=netrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean net revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n",outline=FALSE, las=1)
grid(nx=NA, ny=NULL,col = "ivory3") #grid over boxplot
par(new=TRUE)
boxplot(value/1000~ year, netrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean net revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n",outline=FALSE, las=1)
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart-9.5, lty=2)
abline(h=seq(0,60,10),col="ivory3", lty="dotted",lwd = 0.4)
mtext("(a)", side=3, line = -1.3, adj = 0.04, font=2, cex = 1)
mtext(text="Annual net revenue by vessel (euro)",side=2, line=0,adj = 0.5,font=2,outer=TRUE)

boxplot(value/1000~ year, grossrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean gross revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n",outline=FALSE, las=1)
grid(nx=NA, ny=NULL,col = "ivory3") #grid over boxplot
par(new=TRUE)
boxplot(value/1000~ year, grossrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean gross revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n",outline=FALSE, las=1)
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart-9.5, lty=2)
abline(h=seq(0,60,10),col="ivory3", lty="dotted",lwd = 0.4)
mtext("(b)", side=3, line = -1.3, adj = 0.04, font=2, cex = 1)
mtext(text="Annual gross revenue by vessel (euro)",side=2, line=3,adj =0.5,at=-15,font=2)

load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/2quotas12seasons.RData")

boxplot(value/1000~ year, netrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean net revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n",outline=FALSE, las=1)
grid(nx=NA, ny=NULL,col = "ivory3") #grid over boxplot
par(new=TRUE)
boxplot(value/1000~ year, netrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean net revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n",outline=FALSE, las=1)
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart-9.5, lty=2)
abline(h=seq(0,60,10),col="ivory3", lty="dotted",lwd = 0.4)
mtext("(c)", side=3, line = -1.3, adj = 0.04, font=2, cex = 1)

boxplot(value/1000 ~ year, grossrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean gross revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n",outline=FALSE, las=1)
grid(nx=NA, ny=NULL,col = "ivory3") #grid over boxplot
par(new=TRUE)
boxplot(value/1000 ~ year, grossrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean gross revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n",outline=FALSE, las=1)
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart-9.5, lty=2)
abline(h=seq(0,60,10),col="ivory3", lty="dotted",lwd = 0.4)
mtext("(d)", side=3, line = -1.3, adj = 0.04, font=2, cex = 1)

mtext(text="Year",side=1,line=0,font=2,outer=TRUE)
#mtext(text="Economics by vessel (thousand euro)",side=2,line=0,font=2,outer=TRUE)
dev.off()



load("~/Dropbox/BoB/MSE/Git/Nekane_MSE/complex model/2quotas12seasons.RData")
ylim <- c(0, 1.5)
boxplot(value/mean(value) ~ year, netrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean net revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n",outline=FALSE, las=1)
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart-9.5, lty=2)
abline(h=seq(0,60,10),col="ivory3", lty="dotted",lwd = 0.4)
mtext("(c)", side=3, line = -1.3, adj = 0.04, font=2, cex = 1)

boxplot(value/mean(value) ~ year, grossrev, cex = 0.6, col="grey",boxwex=1, xlab = "Year", ylab = "Mean gross revenue", xlim=c(1, NUMRUNS+1), ylim=ylim, xaxt = "n",outline=FALSE, las=1)
axis(1, at =c(0.5,(MPstart-stab.model)/2+0.5, MPstart-stab.model+0.5, MPstart+((endy-MPstart)/2)-stab.model+0.5, NUMRUNS+0.5),
     labels = c(stab.model, stab.model+((MPstart-stab.model)/2),    "MP", MPstart+((endy-MPstart)/2), endy))
abline(v=MPstart-9.5, lty=2)
abline(h=seq(0,60,10),col="ivory3", lty="dotted",lwd = 0.4)
mtext("(d)", side=3, line = -1.3, adj = 0.04, font=2, cex = 1)

mtext(text="Year",side=1,line=0,font=2,outer=TRUE)
mtext(text="Economics by vessel (thousand euro)",side=2,line=0,font=2,outer=TRUE)
dev.off()


