# create figure 4--partitioning of coexistence mechanisms

library(tidyverse)

# source the file that paritions coexistence
source("invader_resident_comparison.R")

avena_bootstrap <- read.csv("avena_bootstrap_results.csv")
erodium_bootstrap <- read.csv("erodium_bootstrap_results.csv")

# 95% confidence interval includes sorted values [25:975]
min <- (dim(avena_bootstrap)[1] - .95*dim(avena_bootstrap)[1])/2
max <- dim(avena_bootstrap)[1] - min

# extract 95% CI
avena_ldgr <- sort(avena_bootstrap[,2])
avena_ldgr_min <- avena_ldgr[min]
avena_ldgr_max <- avena_ldgr[max]
avena_ldgr_sd <- sd(avena_ldgr)

erodium_ldgr <- sort(erodium_bootstrap[,2])
erodium_ldgr_min <- erodium_ldgr[min]
erodium_ldgr_max <- erodium_ldgr[max]
erodium_ldgr_sd <- sd(erodium_ldgr)

avena_0 <- sort(avena_bootstrap[,3])
avena_0_min <- avena_0[min]
avena_0_max <- avena_0[max]
avena_0_sd <- sd(avena_0)

erodium_0 <- sort(erodium_bootstrap[,3])
erodium_0_min <- erodium_0[min]
erodium_0_max <- erodium_0[max]
erodium_0_sd <- sd(erodium_0)

avena_a <- sort(avena_bootstrap[,4])
avena_a_min <- avena_a[min]
avena_a_max <- avena_a[max]
avena_a_sd <- sd(avena_a)

erodium_a <- sort(erodium_bootstrap[,4])
erodium_a_min <- erodium_a[min]
erodium_a_max <- erodium_a[max]
erodium_a_sd <- sd(erodium_a)

avena_l <- sort(avena_bootstrap[,5])
avena_l_min <- avena_l[min]
avena_l_max <- avena_l[max]
avena_l_sd <- sd(avena_l)

erodium_l <- sort(erodium_bootstrap[,5])
erodium_l_min <- erodium_l[min]
erodium_l_max <- erodium_l[max]
erodium_l_sd <- sd(erodium_l)

avena_al <- sort(avena_bootstrap[,6])
avena_al_min <- avena_al[min]
avena_al_max <- avena_al[max]
avena_al_sd <- sd(avena_al)

erodium_al <- sort(erodium_bootstrap[,6])
erodium_al_min <- erodium_al[min]
erodium_al_max <- erodium_al[max]
erodium_al_sd <- sd(erodium_al)

avena_lower <- c(avena_ldgr_min, avena_0_min, avena_a_min, avena_l_min, avena_al_min)
avena_upper <- c(avena_ldgr_max, avena_0_max, avena_a_max, avena_l_max, avena_al_max)
avena_sd <- c(avena_ldgr_sd, avena_0_sd, avena_a_sd, avena_l_sd, avena_al_sd)
avena_mean <- c(mean(avena_ldgr), mean(avena_0), mean(avena_a), mean(avena_l), mean(avena_al))

erodium_lower <- c(erodium_ldgr_min, erodium_0_min, erodium_a_min, erodium_l_min, erodium_al_min)
erodium_upper <- c(erodium_ldgr_max, erodium_0_max, erodium_a_max, erodium_l_max, erodium_al_max)
erodium_sd <- c(erodium_ldgr_sd, erodium_0_sd, erodium_a_sd, erodium_l_sd, erodium_al_sd)
erodium_mean <- c(mean(erodium_ldgr), mean(erodium_0), mean(erodium_a), mean(erodium_l), mean(erodium_al))

params <- c("a", "b", "c", "d", "e")

dat <- as.data.frame(t(rbind(params, ir_avena_results_weighted, ir_erodium_results_weighted, 
                             avena_lower, avena_upper, erodium_lower, erodium_upper)))
names(dat) = c("params", "Avena", "Erodium", "Avena_lower", "Avena_upper", "Erodium_lower", "Erodium_upper")

 # -------------------------------------------------------------------------------------------------------------

quartz(width=6, height=4)
par(mfrow=c(1,2), oma=c(4,2, 1.5, 1), mar=c(0,1,0,0))
x <- barplot(ir_avena_results_weighted, ylim=c(-3, 3.5), xlab="", ylab=c("Growth Rate When Rare"),
             col=c("grey40", "grey70", "grey70", "grey70", "grey70"))

abline(h=0)
axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
                                                 "b" = expression(bar(Delta)[i]^0),
                                                 "c" = expression(bar(Delta)[i]^alpha),
                                                 "d" = expression(bar(Delta)[i]^lambda),
                                                 "e" = expression(bar(Delta)[i]^{alpha*lambda})))


box(which = "plot", lty = "solid")
mtext(expression("Avena"), side=3, outer=FALSE, adj=0.5)
text(x=5.9, y=1.85, "A)")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=avena_lower, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=avena_upper, length=.05,
       angle=90, col=c("black"), code=3)

x <- barplot(ir_erodium_results_weighted, ylim=c(-3, 3.5), xlab="", ylab=c("Growth Rate When Rare"),
             col=c("grey40", "grey70", "grey70", "grey70", "grey70"), yaxt="n")

abline(h=0)
axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
                                                 "b" = expression(bar(Delta)[i]^0),
                                                 "c" = expression(bar(Delta)[i]^alpha),
                                                 "d" = expression(bar(Delta)[i]^lambda),
                                                 "e" = expression(bar(Delta)[i]^{alpha*lambda})))


box(which = "plot", lty = "solid")
mtext(expression("Erodium"), side=3, outer=FALSE, adj=0.5)
text(x=5.9, y=1.85, "B)")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=erodium_lower, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=erodium_upper, length=.05,
       angle=90, col=c("black"), code=3)

mtext("Mechanistic Partitioning", side=1, outer=TRUE, adj=0.5, line=2.25)
mtext("Growth Rate When Rare", side=2, outer=TRUE, adj=0.5, line=.75)


# -------------------------------------------------------------------------------------------------------------

quartz(width=6, height=4)
par(mfrow=c(1,2), oma=c(4,2, 1.5, 1), mar=c(0,1,0,0))
x <- barplot(ir_avena_results_weighted, ylim=c(-2, 2.75), xlab="", ylab=c("Growth Rate When Rare"),
             col=c("grey40", "grey70", "grey70", "grey70", "grey70"))

abline(h=0)
axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
                                                 "b" = expression(bar(Delta)[i]^0),
                                                 "c" = expression(bar(Delta)[i]^alpha),
                                                 "d" = expression(bar(Delta)[i]^lambda),
                                                 "e" = expression(bar(Delta)[i]^{alpha*lambda})))


box(which = "plot", lty = "solid")
mtext(expression("Avena"), side=3, outer=FALSE, adj=0.5)
text(x=.3, y=2.5, "(a)")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=ir_avena_results_weighted-avena_sd, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=ir_avena_results_weighted+avena_sd, length=.05,
       angle=90, col=c("black"), code=3)

x <- barplot(ir_erodium_results_weighted, ylim=c(-2, 2.75), xlab="", ylab=c("Growth Rate When Rare"),
             col=c("grey40", "grey70", "grey70", "grey70", "grey70"), yaxt="n")

abline(h=0)
axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
                                                 "b" = expression(bar(Delta)[i]^0),
                                                 "c" = expression(bar(Delta)[i]^alpha),
                                                 "d" = expression(bar(Delta)[i]^lambda),
                                                 "e" = expression(bar(Delta)[i]^{alpha*lambda})))


box(which = "plot", lty = "solid")
mtext(expression("Erodium"), side=3, outer=FALSE, adj=0.5)
text(x=.3, y=2.5, "(b)")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=ir_erodium_results_weighted-erodium_sd, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=ir_erodium_results_weighted+erodium_sd, length=.05,
       angle=90, col=c("black"), code=3)

mtext("Mechanistic Partitioning", side=1, outer=TRUE, adj=0.5, line=2.25)
mtext("Growth Rate When Rare", side=2, outer=TRUE, adj=0.5, line=.95)
