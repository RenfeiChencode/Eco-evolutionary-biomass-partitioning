basedir <- "C:\\Users\\Lenovo\\Desktop\\EBToutput\\"

datadir <- paste0(basedir, "fig4equilibrium\\")
modeldir <- paste0(basedir, "fig4equilibrium\\")
fname <- paste0(basedir, "fig4equilibrium\\fig4.pdf")
setwd(basedir)

#setwd("C:\\Users\\Lenovo\\Desktop\\EBToutput\\fig4equilibrium")
library(ggplot2);library(gridExtra);library(lmodel2);library(scales);library(grid)
library(tseries)# test the stability of time series
library(sqldf);library(dplyr);library(patchwork);library(ggsignif)
library(PSPManalysis)
library(latex2exp)
library(shape)
library(igraph)
library(DescTools)
####################################################################
ToPdf <- T
EBT_PSPM_cmp <- F
#1 nursery; 2 harvested area; 3 protected area
EBTtcol         <-  1
EBTres1col      <-  2
EBTres2col      <-  3
EBTres3col      <-  4
EBTtotnumcol    <-  5
EBTtotbiocol    <-  6
EBTh1numcol     <-  7
EBTh1biocol     <-  8
EBTh2numcol     <-  9
EBTh2biocol     <- 10
EBTh3numcol     <- 11
EBTh3biocol     <- 12
EBTh2juvnumcol  <- 13
EBTh2juvbiocol  <- 14
EBTh3juvnumcol  <- 15
EBTh3juvbiocol  <- 16
EBTh2adunumcol  <- 17
EBTh2adubiocol  <- 18
EBTh3adunumcol  <- 19
EBTh3adubiocol  <- 20
EBTsmoltsizecol <- 21   
EBTtotreprocol  <- 22
EBTnumcohcol    <- 23
EBTyieldjuv     <- 24
EBTyieldadu     <- 25
EBTyieldtot     <- 26
EBTbifparcol    <- 27
EBTperiodcol    <- 28

PSPMbifparcol    <-  1
PSPMres1col      <-  2
PSPMres2col      <-  3
PSPMres3col      <-  4
PSPMh1biocol     <-  9
PSPMh2juvbiocol  <- 10
PSPMh3juvbiocol  <- 11
PSPMh2adubiocol  <- 12
PSPMh3adubiocol  <- 13
PSPMh1numcol     <- 14
PSPMh2numcol     <- 15
PSPMh3numcol     <- 16

iArrows <- igraph:::igraph.Arrows
DefaultPars <- c(Rho1 = 0.5, Rho2 =    0.5, Delta = 1.0, 
                 Eta1 = 1.0, Eta2 =    0.8, Eta3  = 0.8, 
                 ETSJ = 1.4, ETSA =    1.4, WS    = 0.35, 
                 Q    = 1.0, Beta = 2000.0, SR    = 0.0)

if (ToPdf) pdf(file = fname, width = (44.0 / 2.54), height = 10.0)
layout(matrix(1:6, nrow = 2, ncol = 3), widths = c(1.3, rep(1.0, 6)), heights = rep(c(0.98, 1.3), 6))
xlimc <- c(-0.0, 1.0)
xlimw <- c(0.05, 0.45)
ylim1 <- c(-0.05, 1.05)
ylim2 <- c(-0.02, 0.35)

cexlab <- 2.0
cexaxs <- 2.0
cexleg <- 1.8
cexttl <- 2.2
cexpnt <- 2.0
cexbt  <- 1.6

axislwd <- 1

par(tcl = 0.6, mgp = c(3, 1, 0))
########## Panel A
lwd=2
pars <- DefaultPars
init=c(0.99, 0.0473386150	,	0.4993470305	,	0.4885173488	, 1)
if (!exists("Equi_SR_Ws035")) {
  Equi_SR_Ws035 <- PSPMequi(modelname = paste0(modeldir, "HarvestWithReserve"), 
                            biftype = "EQ", startpoint = init, 
                            stepsize = -0.01, parbnds = c(11, 0.01, 1.2), 
                            parameters = pars, clean=TRUE, 
                            options = c("popEVO", "0", "parEVO", "8"))
}
dt          <- Equi_SR_Ws035$curvepoints
EBTdt       <- read.table(paste0(datadir, "fig4noevo14dwn.minmax2.out"))
EBTmax      <- EBTdt[2 * (1:(nrow(EBTdt) / 2.0)),]
EBTmin      <- EBTdt[2 * (1:(nrow(EBTdt) / 2.0)) - 1,]

stable      <- EBTmax[, EBTperiodcol] > 1.0E-4
indx2       <- nrow(EBTmax)
indx1       <- max((1:indx2)[stable]) + 1
EBTdwnM     <- rbind(EBTmax[(1:indx1),], EBTmin[(indx1:1),])
stable      <- dt[, PSPMbifparcol] < min(EBTdwnM[, EBTbifparcol])

par(mar = c(0, 10, 4.0, 2.5))
plot(NULL, NULL, xlim = xlimc, ylim = ylim1, 
     xlab = "", ylab = "", xaxs = "i", yaxs = "i",
     xaxt = "n", yaxt = "n")
axis(2, label = T, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext("Number \nin stage 2 ", 2, cex = cexlab, line = 4.8)
mtext("A: no evolution", 3, cex = cexttl, line = 1.0)
lines(dt[ stable,PSPMbifparcol], dt[ stable,PSPMh3numcol], lwd = 2 * lwd, col = "blue")
lines(dt[!stable,PSPMbifparcol], dt[!stable,PSPMh3numcol], lwd = lwd, col = "blue", lty = 2)
lines(EBTdwnM[,EBTbifparcol], EBTdwnM[,EBTh3numcol], lwd = lwd, col = "blue")
lines(dt[ stable,PSPMbifparcol], dt[ stable,PSPMh2numcol], lwd = 2 * lwd, col = "red")
lines(dt[!stable,PSPMbifparcol], dt[!stable,PSPMh2numcol], lwd = lwd, col = "red", lty = 2)
lines(EBTdwnM[,EBTbifparcol], EBTdwnM[,EBTh2numcol], lwd = lwd, col = "red")



par(mar = c(10, 10, 0, 2.5))
plot(NULL, NULL, xlim = xlimc, ylim = ylim2, 
     xlab = "", ylab = "", xaxs = "i", yaxs = "i",
     xaxt = "n", yaxt = "n")
axis(1, label = T, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd)
axis(2, at = (0:3) * 0.1, label = c("0.0", "0.1", "0.2", "0.3"), 
     cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext("Biomass\nin stage 2", 2, cex = cexlab, line = 4.8)
mtext("Protected fraction", 1, cex = cexlab, line = 4.5)

lines(dt[stable,PSPMbifparcol], dt[stable,PSPMh2juvbiocol] + dt[stable,PSPMh2adubiocol], lwd = 2 * lwd, col = "red")
lines(dt[stable,PSPMbifparcol], dt[stable,PSPMh3juvbiocol] + dt[stable,PSPMh3adubiocol], lwd = 2 * lwd, col = "blue")
lines(dt[!stable,PSPMbifparcol], dt[!stable,PSPMh2juvbiocol] + dt[!stable,PSPMh2adubiocol], lwd = lwd, lty = 2, col = "red")
lines(dt[!stable,PSPMbifparcol], dt[!stable,PSPMh3juvbiocol] + dt[!stable,PSPMh3adubiocol], lwd = lwd, lty = 2, col = "blue")
lines(EBTdwnM[,EBTbifparcol], EBTdwnM[,EBTh2biocol], lwd = lwd, col = "red")
lines(EBTdwnM[,EBTbifparcol], EBTdwnM[,EBTh3biocol], lwd = lwd, col = "blue")

########## Panel B
init=c(0.5, 0.1363040134,	0.4840562681,	0.5000000000,1)#######
pars <- DefaultPars
if (!exists("Equi_SR00_Ws")) {
  Equi_SR00_Ws <- PSPMequi(modelname = paste0(modeldir, "HarvestWithReserve"), 
                           biftype = "EQ", startpoint = init, 
                           stepsize = -0.01, parbnds = c(8, 0.05, 0.5), 
                           parameters = pars, clean=TRUE, 
                           options = c("popEVO", "0"))
}

dt          <- Equi_SR00_Ws$curvepoints
dtbp        <- Equi_SR00_Ws$bifpoints
dtbt        <- Equi_SR00_Ws$biftypes

EBTdt       <- read.table(paste0(datadir, "fig4nopro14dwn.minmax2.out"))
EBTmax      <- EBTdt[2 * (1:(nrow(EBTdt) / 2.0)),]
EBTmin      <- EBTdt[2 * (1:(nrow(EBTdt) / 2.0)) - 1,]

stable      <- EBTmax[, EBTperiodcol] >1.0E-4
indx2       <- nrow(EBTmax)
indx1       <- max((1:indx2)[stable]) + 1
EBTmax_upS  <- EBTmax[(1:indx1),]
EBTmin_upS  <- EBTmin[(1:indx1),]
EBTupM      <- rbind(EBTmax[(indx2:indx1),], EBTmin[(indx1:indx2),])
stable      <- dt[, PSPMbifparcol] < max(EBTupM[, EBTbifparcol])

par(mar = c(0, 1, 4.0, 2.5))
plot(NULL, NULL, xlim = xlimw, ylim = ylim1, 
     xlab = "", ylab = "", xaxs = "i", yaxs = "i",
     xaxt = "n", yaxt = "n")
mtext("B: no protection", 3, cex = cexttl, line = 1.0)
lines(dt[ stable,PSPMbifparcol], dt[ stable,PSPMh2numcol], lwd = 2 * lwd, col = "red")
lines(dt[!stable,PSPMbifparcol], dt[!stable,PSPMh2numcol], lwd = lwd, col = "red", lty = 2)
#lines(EBTupM[,EBTbifparcol], EBTupM[,EBTh2numcol], lwd = lwd, col = "green")
lines(EBTmax_upS[,EBTbifparcol],  EBTmax_upS[,EBTh2numcol],  lwd = lwd, col = "red")
lines(EBTmin_upS [,EBTbifparcol], EBTmin_upS [,EBTh2numcol], lwd = lwd, col = "red")

points(dtbp[,PSPMbifparcol], dtbp[,PSPMh2numcol], pch = 19, col = "yellow", cex = cexpnt)
text(dtbp[,PSPMbifparcol], dtbp[,PSPMh2numcol], "ESS", cex = cexbt, adj = c(1.3, 0.45))

Arrows(0.09, 0.07, 0.16, 0.07, arr.type="triangle", lwd = 5, arr.width = 0.2, col = "black")
Arrows(0.44, 0.85, 0.39, 0.85, arr.type="triangle", lwd = 5, arr.width = 0.2, col = "black")
legend("topleft", legend = c("Stable", "Unstable","Minimum / maximum"), col = "black", lty = c(1,2,1), cex = cexbt,lwd = c(4,2,2))

par(mar = c(10, 1, 0, 2.5))
plot(NULL, NULL, xlim = xlimw, ylim = ylim2, 
     xlab = "", ylab = "", xaxs = "i", yaxs = "i",
     xaxt = "n", yaxt = "n")
axis(1, at = 0.05 + (0:4) * 0.1, label = c("0.05", "0.15", "0.25", "0.35","0.45"), 
     cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd)
mtext("Minimum reproductive \n threshold size", 1, cex = cexlab, line = 7)

lines(dt[stable,PSPMbifparcol], dt[stable,PSPMh2juvbiocol] + dt[stable,PSPMh2adubiocol], lwd = 2 * lwd, col = "red")
lines(dt[!stable,PSPMbifparcol], dt[!stable,PSPMh2juvbiocol] + dt[!stable,PSPMh2adubiocol], lwd = lwd, lty = 2, col = "red")
#lines(EBTupM[,EBTbifparcol], EBTupM[,EBTh2biocol], lwd = lwd, col = "red")
lines(EBTmax_upS[,EBTbifparcol],  EBTmax_upS[,EBTh2biocol],  lwd = lwd, col = "red")
lines(EBTmin_upS [,EBTbifparcol], EBTmin_upS [,EBTh2biocol], lwd = lwd, col = "red")

points(dtbp[,PSPMbifparcol], dtbp[,PSPMh2juvbiocol] + dtbp[,PSPMh2adubiocol], pch = 19, col = "yellow", cex = cexpnt)

########## Panel C
init=c(0.05, 0.9142577824,	0.2843496419,	0.1620576624, 1)#######
pars <- DefaultPars
pars["SR"] <- 0.25
if (!exists("Equi_SR03_Ws")) {
  Equi_SR03_Ws <- PSPMequi(modelname = paste0(modeldir, "HarvestWithReserve"), 
                           biftype = "EQ", startpoint = init, 
                           stepsize = 0.01, parbnds = c(8, 0.05, 0.5), 
                           parameters = pars, clean=TRUE, 
                           options = c("popEVO", "0"))
}

dt          <- Equi_SR03_Ws$curvepoints
dtbp        <- Equi_SR03_Ws$bifpoints
dtbt        <- Equi_SR03_Ws$biftypes

EBTdt       <- read.table(paste0(datadir, "fig4c025evo14dwn.minmax2.out"))
EBTmax      <- EBTdt[2 * (1:(nrow(EBTdt) / 2.0)),]
EBTmin      <- EBTdt[2 * (1:(nrow(EBTdt) / 2.0)) - 1,]

stable      <- EBTmax[, EBTperiodcol] > 1.0E-4
indx2       <- nrow(EBTmax)
indx1       <- max((1:indx2)[stable]) + 1
EBTdwnM     <- rbind(EBTmax[(1:indx1),], EBTmin[(indx1:1),])
stable      <- dt[, PSPMbifparcol] < min(EBTdwnM[, EBTbifparcol])
par(mar = c(0, 1, 4.0, 2.5))
plot(NULL, NULL, xlim = xlimw, ylim = ylim1, 
     xlab = "", ylab = "", xaxs = "i", yaxs = "i",
     xaxt = "n", yaxt = "n")
mtext("C: 25% protected", 3, cex = cexttl, line = 1.0)
# text(xlim1[2], ylimb[2], "A", cex = 4.0, xpd = T, adj = c(0.5, 0))


indx1 <- which.min(abs(dt[,PSPMbifparcol] - dtbp[1,PSPMbifparcol]))
indx2 <- which.min(abs(dt[,PSPMbifparcol] - dtbp[3,PSPMbifparcol]))
indx3 <- max((indx2:nrow(dt))[dt[(indx2:nrow(dt)),PSPMbifparcol] < min(EBTdwnM[,EBTbifparcol])])

lines(dt[(1:indx1),PSPMbifparcol], dt[(1:indx1),PSPMh3numcol], lwd = 2 * lwd, col = "blue")
lines(dt[(indx1:indx2),PSPMbifparcol], dt[(indx1:indx2),PSPMh3numcol], lwd = lwd, lty = 3, col = "blue")
lines(dt[(indx2:indx3),PSPMbifparcol], dt[(indx2:indx3),PSPMh3numcol], lwd = 2 * lwd, col = "blue")
lines(dt[(indx3:nrow(dt)),PSPMbifparcol], dt[(indx3:nrow(dt)),PSPMh3numcol], lwd = lwd, col = "blue", lty = 2)
lines(EBTdwnM[,EBTbifparcol], EBTdwnM[,EBTh3numcol], lwd = lwd, col = "blue")

lines(dt[(1:indx1),PSPMbifparcol], dt[(1:indx1),PSPMh2numcol], lwd = 2 * lwd, col = "red")
lines(dt[(indx1:indx2),PSPMbifparcol], dt[(indx1:indx2),PSPMh2numcol], lwd = lwd, lty = 3, col = "red")
lines(dt[(indx2:indx3),PSPMbifparcol], dt[(indx2:indx3),PSPMh2numcol], lwd = 2 * lwd, col = "red")
lines(dt[(indx3:nrow(dt)),PSPMbifparcol], dt[(indx3:nrow(dt)),PSPMh2numcol], lwd = lwd, col = "red", lty = 2)
lines(EBTdwnM[,EBTbifparcol], EBTdwnM[,EBTh2numcol], lwd = lwd, col = "red")


points(dtbp[2,PSPMbifparcol], dtbp[2,PSPMh3numcol], pch = 19, col = "yellow", cex = cexpnt)
points(dtbp[1,PSPMbifparcol], dtbp[1,PSPMh3numcol], pch = 19, col = "blue", cex = cexpnt)
points(dtbp[3,PSPMbifparcol], dtbp[3,PSPMh3numcol], pch = 19, col = "blue", cex = cexpnt)
text(dtbp[1,PSPMbifparcol]+0.04, dtbp[1,PSPMh3numcol], "LP", cex = cexbt, adj = c(1.3, -0.25))

points(dtbp[2,PSPMbifparcol], dtbp[2,PSPMh2numcol], pch = 19, col = "yellow", cex = cexpnt)
points(dtbp[1,PSPMbifparcol], dtbp[1,PSPMh2numcol], pch = 19, col = "red", cex = cexpnt)
points(dtbp[3,PSPMbifparcol], dtbp[3,PSPMh2numcol], pch = 19, col = "red", cex = cexpnt)
text(dtbp[3,PSPMbifparcol], dtbp[3,PSPMh2numcol], "LP", cex = cexbt, adj = c(1.3, -0.25))

Arrows(0.09, 0.09, 0.16, 0.09, arr.type="triangle", lwd = 5, arr.width = 0.2, col = "black")
Arrows(0.44, 0.85, 0.37, 0.85, arr.type="triangle", lwd = 5, arr.width = 0.2, col = "black")
DrawEllipse(dtbp[1,PSPMbifparcol] + 0.012, 0.195, radius.x = 0.008, radius.y = 0.168, col = NA, lwd = 1.5, border = "red")
DrawEllipse(dtbp[1,PSPMbifparcol] + 0.052, 0.199, radius.x = 0.008, radius.y = 0.1805, col = NA, lwd = 1.5, border = "blue")
Arrows(dtbp[1,PSPMbifparcol] + 0.0047, 0.15, 
       dtbp[1,PSPMbifparcol] + 0.0047, 0.14, 
       arr.type="triangle", lwd = 3, arr.width = 0.1, col = "red")
Arrows(dtbp[1,PSPMbifparcol] + 0.0192, 0.27, 
       dtbp[1,PSPMbifparcol] + 0.0192, 0.28, 
       arr.type="triangle", lwd = 3, arr.width = 0.1, col = "red")
Arrows(dtbp[1,PSPMbifparcol] + 0.0447, 0.15, 
       dtbp[1,PSPMbifparcol] + 0.0447, 0.14, 
       arr.type="triangle", lwd = 3, arr.width = 0.1, col = "blue")
Arrows(dtbp[1,PSPMbifparcol] + 0.0592, 0.27, 
       dtbp[1,PSPMbifparcol] + 0.0592, 0.28, 
       arr.type="triangle", lwd = 3, arr.width = 0.1, col = "blue")
legend("topleft", legend = c("Harvested fraction", "Protected fraction"), col = c("red","blue"), lty =1, cex = cexbt,lwd = 4)

par(mar = c(10, 1, 0, 2.5))
plot(NULL, NULL, xlim = xlimw, ylim = ylim2, 
     xlab = "", ylab = "", xaxs = "i", yaxs = "i",
     xaxt = "n", yaxt = "n")
axis(1, at = 0.05 + (0:4) * 0.1, label = c("0.05", "0.15", "0.25", "0.35","0.45"), 
     cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd)
mtext("Minimum reproductive \n threshold size", 1, cex = cexlab, line = 7)

lines(dt[(1:indx1),PSPMbifparcol], dt[(1:indx1),PSPMh2juvbiocol] + dt[(1:indx1),PSPMh2adubiocol], lwd = 2 * lwd, col = "red")
lines(dt[(indx1:indx2),PSPMbifparcol], dt[(indx1:indx2),PSPMh2juvbiocol] + dt[(indx1:indx2),PSPMh2adubiocol], lwd = lwd, lty = 3, col = "red")
lines(dt[(indx2:indx3),PSPMbifparcol], dt[(indx2:indx3),PSPMh2juvbiocol] + dt[(indx2:indx3),PSPMh2adubiocol], lwd = 2 * lwd, col = "red")
lines(dt[(indx3:nrow(dt)),PSPMbifparcol], dt[(indx3:nrow(dt)),PSPMh2juvbiocol] + dt[(indx3:nrow(dt)),PSPMh2adubiocol], lwd = lwd, col = "red", lty = 2)

lines(dt[(1:indx1),PSPMbifparcol], dt[(1:indx1),PSPMh3juvbiocol] + dt[(1:indx1),PSPMh3adubiocol], lwd = 2 * lwd, col = "blue")
lines(dt[(indx1:indx2),PSPMbifparcol], dt[(indx1:indx2),PSPMh3juvbiocol] + dt[(indx1:indx2),PSPMh3adubiocol], lwd = lwd, lty = 3, col = "blue")
lines(dt[(indx2:indx3),PSPMbifparcol], dt[(indx2:indx3),PSPMh3juvbiocol] + dt[(indx2:indx3),PSPMh3adubiocol], lwd = 2 * lwd, col = "blue")
lines(dt[(indx3:nrow(dt)),PSPMbifparcol], dt[(indx3:nrow(dt)),PSPMh3juvbiocol] + dt[(indx3:nrow(dt)),PSPMh3adubiocol], lwd = lwd, col = "blue", lty = 2)

lines(EBTdwnM[,EBTbifparcol], EBTdwnM[,EBTh2biocol], lwd = lwd, col = "red")
lines(EBTdwnM[,EBTbifparcol], EBTdwnM[,EBTh3biocol], lwd = lwd, col = "blue")

points(dtbp[2,PSPMbifparcol], dtbp[2,PSPMh2juvbiocol] + dtbp[2,PSPMh2adubiocol], pch = 19, col = "yellow", cex = cexpnt)
points(dtbp[2,PSPMbifparcol], dtbp[2,PSPMh3juvbiocol] + dtbp[2,PSPMh3adubiocol], pch = 19, col = "yellow", cex = cexpnt)
points(dtbp[1,PSPMbifparcol],dtbp[1,PSPMh2juvbiocol] + dtbp[1,PSPMh2adubiocol], pch = 19, col = "red", cex = cexpnt)
points(dtbp[3,PSPMbifparcol], dtbp[3,PSPMh2juvbiocol] + dtbp[3,PSPMh2adubiocol], pch = 19, col = "red", cex = cexpnt)
points(dtbp[1,PSPMbifparcol],dtbp[1,PSPMh3juvbiocol] + dtbp[1,PSPMh3adubiocol], pch = 19, col = "blue", cex = cexpnt)
points(dtbp[3,PSPMbifparcol], dtbp[3,PSPMh3juvbiocol] + dtbp[3,PSPMh3adubiocol], pch = 19, col = "blue", cex = cexpnt)

if (ToPdf) {
  dev.off()
  system(paste("open ", fname))
}


setwd("C:\\Users\\Lenovo\\Desktop\\EBToutput\\fig2dynamic")
library(ggplot2);library(gridExtra);library(lmodel2);library(scales);library(grid)
library(tseries)# test the stability of time series
library(sqldf);library(dplyr);library(patchwork);library(ggsignif)
####################################################################
c025noevo=read.table("fig2c025noevo14.out")
c005noevo=read.table("fig2c005noevo14.out")
c025evo=read.table("fig2c025evo14.out")
c005evo=read.table("fig2c005evo14.out")
tcol         <-  1
res1col      <-  2
res2col      <-  3
res3col      <-  4
totnumcol    <-  5
totbiocol    <-  6
h1numcol     <-  7
h1biocol     <-  8
h2numcol     <-  9
h2biocol     <- 10
h3numcol     <- 11
h3biocol     <- 12
h2juvnumcol  <- 13
h2juvbiocol  <- 14
h3juvnumcol  <- 15
h3juvbiocol  <- 16
h2adunumcol  <- 17
h2adubiocol  <- 18
h3adunumcol  <- 19
h3adubiocol  <- 20
smoltsizecol <- 21   
totreprocol  <- 22
numcohcol    <- 23
yieldjuvcol  <- 24
yieldaducol  <- 25
yieldtotcol  <- 26

numnoevo=c(c005noevo[,h2numcol],c025noevo[,h2numcol])
bionoevo=c(c005noevo[,h2biocol],c025noevo[,h2biocol])
numevo=c(c005evo[,h2numcol],c025evo[,h2numcol])
bioevo=c(c005evo[,h2biocol],c025evo[,h2biocol])
NN=length(c005noevo[,tcol])
class=c(rep("5% protected",NN),rep("25% protected",NN))
time=rep((c005noevo[,tcol]),2)
dfnoevo=data.frame(numnoevo=numnoevo,bionoevo=bionoevo,class=class,time=time)

NN=length(c005evo[,tcol])
class=c(rep("5% protected",NN),rep("25% protected",NN))
time=rep((c005evo[,tcol]),2)
dfevo=data.frame(numevo=numevo,bioevo=bioevo,class=class,time=time)

fig2A=ggplot(data = dfnoevo,aes(x=time,y=bionoevo,group=class,color=class))+geom_path(size=1)+
  geom_vline(xintercept=500,colour="grey",lty="dashed",size=2)+
  annotate("text",x=460, y=0.2, label="Protection start",colour="grey",angle = 90,size=10,family="serif")+
  annotate("text",x=1, y=0.37, label="A",angle = 0,size=10,family="serif")+
  labs(title="No evolution",colour="",x="Time",y="Biomass")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.9), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(0,1000,by=250),limits=c(0,1000))+
  scale_y_continuous(breaks=seq(0,0.4,by=0.1),limits=c(0,0.4))

fig2B=ggplot(data = dfevo,aes(x=time,y=bioevo,group=class,color=class))+geom_path(size=1)+
  geom_vline(xintercept=500,colour="grey",lty="dashed",size=2)+
  annotate("text",x=390, y=0.2, label="Protection start",colour="grey",angle = 90,size=10 ,family="serif")+
  annotate("text",x=1, y=0.38, label="B",angle = 0,size=10,family="serif")+
  labs(title="Evolution",colour="",x="Time",y="Biomass")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = "none", 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(0,3000,by=1000),limits=c(0,3000))+
  scale_y_continuous(breaks=seq(0,0.4,by=0.1),limits=c(0,0.4))


fig2C=ggplot(data = dfnoevo,aes(x=time,y=numnoevo,group=class,color=class))+geom_path(size=1)+
  geom_vline(xintercept=500,colour="grey",lty="dashed",size=2)+
  annotate("text",x=460, y=0.5, label="Protection start",colour="grey",angle = 90,size=10,family="serif" )+
  annotate("text",x=1, y=0.98, label="C",angle = 0,size=10,family="serif")+
  labs(title="No evolution",colour="",x="Time",y="Number")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = "none", 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(0,1000,by=250),limits=c(0,1000))+
  scale_y_continuous(breaks=seq(0,1,by=0.25),limits=c(0,1))


fig2D=ggplot(data = dfevo,aes(x=time,y=numevo,group=class,color=class))+geom_path(size=1)+geom_vline(xintercept=500,colour="grey",lty="dashed",size=2)+
  annotate("text",x=390, y=0.5, label="Protection start",colour="grey",angle = 90,size=10,family="serif")+
  annotate("text",x=1, y=0.98, label="D",angle = 0,size=10,family="serif")+
  labs(title="Evolution",colour="",x="Time",y="Number")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = "none", 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(0,3000,by=1000),limits=c(0,3000))+
  scale_y_continuous(breaks=seq(0,1,by=0.25),limits=c(0,1))


####################################
#figure 5
setwd("C:\\Users\\Lenovo\\Desktop\\EBToutput\\fig5yield")
c025noevo=read.table("fig5c025noevofishing14.out")
c005noevo=read.table("fig5c005noevofishing14.out")
c025evo=read.table("fig5c025evofishing14.out")
c005evo=read.table("fig5c005evofishing14.out")

bef025evo=mean(c025evo[1:500,yieldjuvcol])
aft025evo=c025evo[501:3001,yieldjuvcol]
ReYjuvevo=aft025evo/bef025evo
bef025evo=mean(c025evo[1:500,yieldaducol])
aft025evo=c025evo[501:3001,yieldaducol]
ReYaduevo=aft025evo/bef025evo

bef025noevo=mean(c025noevo[1:500,yieldjuvcol])
aft025noevo=c025noevo[501:3001,yieldjuvcol]
ReYjuvnoevo=aft025noevo/bef025noevo

bef025noevo=mean(c025noevo[1:500,yieldaducol])
aft025noevo=c025noevo[501:3001,yieldaducol]
ReYadunoevo=aft025noevo/bef025noevo
ReY025=c(ReYjuvevo,ReYaduevo,ReYjuvnoevo,ReYadunoevo)
NNevo=length(ReYjuvevo);NNnoevo=length(ReYjuvnoevo)
group1=c(rep("1Moderate",NNevo),rep("2Large",NNevo),rep("4Moderate",NNnoevo),rep("5Large",NNnoevo))
group2=c(rep("Evolution",2*NNevo),rep("No evolution",2*NNnoevo))
meanReY025=c(mean(ReYjuvevo),mean(ReYaduevo),mean(ReYjuvnoevo),mean(ReYadunoevo))
sdReY025=c(sd(ReYjuvevo),sd(ReYaduevo),sd(ReYjuvnoevo),sd(ReYadunoevo))
yerrormin=meanReY025-sdReY025;yerrormax=meanReY025+sdReY025
ymin=c(rep(yerrormin[1],NNevo),rep(yerrormin[2],NNevo),
       rep(yerrormin[3],NNnoevo),rep(yerrormin[4],NNnoevo))
ymax=c(rep(yerrormax[1],NNevo),rep(yerrormax[2],NNevo),
       rep(yerrormax[3],NNnoevo),rep(yerrormax[4],NNnoevo))
df=data.frame(yield=ReY025,group1=group1,group2=group2,ymin=ymin,ymax=ymax)
fig5A=ggplot(df,aes(group1,yield))+ geom_bar(aes(fill=group2),stat="summary",fun=mean,position="dodge")+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.2,size=1, colour="orange",position = position_dodge(c(0.9)),data=df)+
  geom_signif(comparisons=list(c("1Moderate","2Large"), c("2Large","5Large"),c("4Moderate","5Large")),
              map_signif_level=T,##是否使用*显示显著性
              tip_length=c(0,0,0,0),y_position=c(2200,1300,850),size=1,textsize=7,test="t.test")+
  #scale_fill_manual(values=c("#037ef3","#f85a40"))+#自定义颜色
  geom_hline(yintercept=1,colour="grey",lty="dashed",size=2)+
  annotate("text",x=0.75, y=3000, label="A",angle = 0,size=10,family="serif")+
  scale_fill_discrete(name="")+#添加图例名字
  labs(title="25% protected",colour="",x="Plant body size",y= expression("Relative yield"))+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0,size=20),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.9), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm")) +
  scale_y_continuous(breaks=seq(0,3000,by=1000),limits=c(0,3000))
################################################
bef005evo=mean(c005evo[1:500,yieldjuvcol])
aft005evo=c005evo[501:3001,yieldjuvcol]
ReYjuvevo=aft005evo/bef005evo
bef005evo=mean(c005evo[1:500,yieldaducol])
aft005evo=c005evo[501:3001,yieldaducol]
ReYaduevo=aft005evo/bef005evo
bef005noevo=mean(c005noevo[1:500,yieldjuvcol])
aft005noevo=c005noevo[501:3001,yieldjuvcol]
ReYjuvnoevo=aft005noevo/bef005noevo
bef005noevo=mean(c005noevo[1:500,yieldaducol])
aft005noevo=c005noevo[501:3001,yieldaducol]
ReYadunoevo=aft005noevo/bef005noevo
ReY005=c(ReYjuvevo,ReYaduevo,ReYjuvnoevo,ReYadunoevo)
NNevo=length(ReYjuvevo);NNnoevo=length(ReYjuvnoevo)
group1=c(rep("1Moderate",NNevo),rep("2Large",NNevo),rep("4Moderate",NNnoevo),rep("5Large",NNnoevo))
group2=c(rep("Evolution",2*NNevo),rep("No evolution",2*NNnoevo))
meanReY005=c(mean(ReYjuvevo),mean(ReYaduevo),mean(ReYjuvnoevo),mean(ReYadunoevo))
sdReY005=c(sd(ReYjuvevo),sd(ReYaduevo),sd(ReYjuvnoevo),sd(ReYadunoevo))
yerrormin=meanReY005-sdReY005;yerrormax=meanReY005+sdReY005
ymin=c(rep(yerrormin[1],NNevo),rep(yerrormin[2],NNevo),
       rep(yerrormin[3],NNnoevo),rep(yerrormin[4],NNnoevo))
ymax=c(rep(yerrormax[1],NNevo),rep(yerrormax[2],NNevo),
       rep(yerrormax[3],NNnoevo),rep(yerrormax[4],NNnoevo))
df=data.frame(yield=ReY005,group1=group1,group2=group2,ymin=ymin,ymax=ymax)
fig5B=ggplot(df,aes(group1,yield))+ geom_bar(aes(fill=group2),stat="summary",fun=mean,position="dodge")+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.2,size=1, colour="orange",position = position_dodge(c(0.9)),data=df)+
  geom_signif(comparisons=list(c("1Moderate","2Large"), c("2Large","5Large"),c("4Moderate","5Large")),
              map_signif_level=T,##是否使用*显示显著性
              tip_length=c(0,0,0,0),y_position=c(2500,1800,1500),size=1,textsize=7,test="t.test")+
  #scale_fill_manual(values=c("#037ef3","#f85a40"))+#自定义颜色
  geom_hline(yintercept=1,colour="grey",lty="dashed",size=2)+
  annotate("text",x=0.75, y=3000, label="B",angle = 0,size=10,family="serif")+
  scale_fill_discrete(name="")+#添加图例名字
  labs(title="5% protected",colour="",x="Plant body size",y= expression("Relative yield"))+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0,size=20),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.9), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm")) +
  scale_y_continuous(breaks=seq(0,3000,by=1000),limits=c(0,3000))
#ggsave("fig2 harvested area.pdf",grid.arrange(fig2A,fig2B,fig2C,fig2D,ncol=2,nrow=2),width = 40, height = 30, units = "cm", dpi = 300) 
#ggsave("fig5.pdf",grid.arrange(fig5A,fig5B,ncol=2,nrow=1),width = 40, height = 20, units = "cm", dpi = 300) 

