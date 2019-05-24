# create figure 4--partitioning of coexistence mechanisms

library(tidyverse)

# source the file that paritions coexistence
source("invader_resident_comparison.R")

avena_bootstrap <- read.csv("avena_leave_one_out.csv")
erodium_bootstrap <- read.csv("erodium_leave_one_out.csv")

avena_lower <- rep(NA,5)
avena_upper <- erodium_lower <- erodium_upper <- avena_lower

for(temp in 1:5) {
  avena_lower[temp] <- min(avena_bootstrap[,temp+1])
  avena_upper[temp] <- max(avena_bootstrap[,temp+1])
  erodium_lower[temp] <- min(erodium_bootstrap[,temp+1])
  erodium_upper[temp] <- max(erodium_bootstrap[,temp+1])
}

params <- c("a", "b", "c", "d", "e")

dat <- as.data.frame(t(rbind(params, ir_avena_results_weighted, ir_erodium_results_weighted, 
                             avena_lower, avena_upper, erodium_lower, erodium_upper)))
names(dat) = c("params", "Avena", "Erodium", "Avena_lower", "Avena_upper", "Erodium_lower", "Erodium_upper")

 # -------------------------------------------------------------------------------------------------------------

quartz(width=6, height=4)
par(mfrow=c(1,2), oma=c(4,2, 1.5, 1), mar=c(0,1,0,0))
x <- barplot(ir_avena_results_weighted, ylim=c(-2, 3), xlab="", ylab=c("Growth Rate When Rare"),
             col=c("grey40", "grey70", "grey70", "grey70", "grey70"))

abline(h=0)
axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
                                                 "b" = expression(bar(Delta)[i]^0),
                                                 "c" = expression(bar(Delta)[i]^alpha),
                                                 "d" = expression(bar(Delta)[i]^lambda),
                                                 "e" = expression(bar(Delta)[i]^{alpha*lambda})))


box(which = "plot", lty = "solid")
mtext(expression("Avena"), side=3, outer=FALSE, adj=0.5)
text(x=5.9, y=2.8, "A)")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=avena_lower, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=avena_upper, length=.05,
       angle=90, col=c("black"), code=3)

x <- barplot(ir_erodium_results_weighted, ylim=c(-2, 3), xlab="", ylab=c("Growth Rate When Rare"),
             col=c("grey40", "grey70", "grey70", "grey70", "grey70"), yaxt="n")

abline(h=0)
axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
                                                 "b" = expression(bar(Delta)[i]^0),
                                                 "c" = expression(bar(Delta)[i]^alpha),
                                                 "d" = expression(bar(Delta)[i]^lambda),
                                                 "e" = expression(bar(Delta)[i]^{alpha*lambda})))


box(which = "plot", lty = "solid")
mtext(expression("Erodium"), side=3, outer=FALSE, adj=0.5)
text(x=5.9, y=2.8, "B)")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=erodium_lower, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=erodium_upper, length=.05,
       angle=90, col=c("black"), code=3)

mtext("Mechanistic Partitioning", side=1, outer=TRUE, adj=0.5, line=2.25)
mtext("Growth Rate When Rare", side=2, outer=TRUE, adj=0.5, line=.9)


# -------------------------------------------------------------------------------------------------------------
# 
# quartz(width=6, height=4)
# par(mfrow=c(1,2), oma=c(4,2, 1.5, 1), mar=c(0,1,0,0))
# x <- barplot(ir_avena_results_weighted, ylim=c(-2, 2.75), xlab="", ylab=c("Growth Rate When Rare"),
#              col=c("grey40", "grey70", "grey70", "grey70", "grey70"))
# 
# abline(h=0)
# axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
#                                                  "b" = expression(bar(Delta)[i]^0),
#                                                  "c" = expression(bar(Delta)[i]^alpha),
#                                                  "d" = expression(bar(Delta)[i]^lambda),
#                                                  "e" = expression(bar(Delta)[i]^{alpha*lambda})))
# 
# 
# box(which = "plot", lty = "solid")
# mtext(expression("Avena"), side=3, outer=FALSE, adj=0.5)
# text(x=.3, y=2.5, "(a)")
# arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=ir_avena_results_weighted-avena_sd, 
#        x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=ir_avena_results_weighted+avena_sd, length=.05,
#        angle=90, col=c("black"), code=3)
# 
# x <- barplot(ir_erodium_results_weighted, ylim=c(-2, 2.75), xlab="", ylab=c("Growth Rate When Rare"),
#              col=c("grey40", "grey70", "grey70", "grey70", "grey70"), yaxt="n")
# 
# abline(h=0)
# axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
#                                                  "b" = expression(bar(Delta)[i]^0),
#                                                  "c" = expression(bar(Delta)[i]^alpha),
#                                                  "d" = expression(bar(Delta)[i]^lambda),
#                                                  "e" = expression(bar(Delta)[i]^{alpha*lambda})))
# 
# 
# box(which = "plot", lty = "solid")
# mtext(expression("Erodium"), side=3, outer=FALSE, adj=0.5)
# text(x=.3, y=2.5, "(b)")
# arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=ir_erodium_results_weighted-erodium_sd, 
#        x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=ir_erodium_results_weighted+erodium_sd, length=.05,
#        angle=90, col=c("black"), code=3)
# 
# mtext("Mechanistic Partitioning", side=1, outer=TRUE, adj=0.5, line=2.25)
# mtext("Growth Rate When Rare", side=2, outer=TRUE, adj=0.5, line=.95)
