##########################################################################################
#                       Replicating Franks and Deneubourg                               #
##########################################################################################


library(readr)
stone_n_results <- read_csv("results/stone_n_results.csv")
stone_n_results$model <- as.factor(stone_n_results$model)
View(stone_n_results)

lists <- read_csv("results/Rbar_and_CoV_values.csv")

# CoV plot
FD1000 <- lists$SO1000_CoV
FD3000 <- lists$SO3000_CoV
FD5000 <- lists$SO5000_CoV

T1000 <- lists$T1000_CoV
T3000 <- lists$T3000_CoV
T5000 <- lists$T5000_CoV

# Rbar
FD1000R <- lists$SO1000_Rbar
FD3000R <- lists$SO3000_Rbar
FD5000R <- lists$SO5000_Rbar

T1000R <- lists$T1000_Rbar
T3000R <- lists$T3000_Rbar
T5000R <- lists$T5000_Rbar

opar = par()
par(family = "serif", mfrow = c(2,1), oma = c(2,0,3,0) + 0.1, mar = c(0,5,1,1) + 0.1) # default serif family is Times New Roman

boxplot(FD1000,T1000, FD3000,T3000, FD5000,T5000, col = c("orange","red"), xaxt = 'n', ylab = "Distance Spread (r)")
#axis(1, at=c(1.5,3.5,5.5), labels = c("stones = 1000","stones = 3000", "stones = 5000"), tick = FALSE)
abline(v = 2.5, lty = 2)
abline(v = 4.5, lty = 2)
legend("bottomright", legend = c("Template + Feedback", "Template Only"), col = c("orange","red"), pch = 15, bty = 'n')

boxplot(FD1000R,T1000R, FD3000R,T3000R, FD5000R,T5000R, ylim = c(-0.002,0.40), col = c("orange","red"), xaxt = 'n', ylab = "Circular Spread (theta)")
axis(1, at=c(1.5,3.5,5.5), labels = c("Stones = 1000","Stones = 3000", "Stones = 5000"), tick = FALSE)
abline(v = 2.5, lty = 2)
abline(v = 4.5, lty = 2)
#legend("topright", legend = c("template + feedback", "template only"), col = c("orange","red"), pch = 15, bty = 'n')

title(main = "Final Stone Distribution - Model Comparison", outer = TRUE, line = 1)

#############################
# delayed pop size increase #
#############################

library(readr)
dlists <- read_csv("results/delayed/Rbar_and_CoV_values_FD_delayed.csv")
r6000 <- dlists$SO_T6000
r10000 <- dlists$SO_T10000
r75000 <- dlists$SO_T75000

par = opar
par(family = "serif", mfrow = c(1,1), mar = c(3,5,3,1), oma = c(0,0,0,0))
boxplot(r6000, r10000, r75000, xaxt = 'n', ylab = "Mean Distance r of Deposition Locations from Cluster", main = "Shift to New Nest Size Over Time", col = c("orange"))
axis(1, at=c(1,2,3), labels = c("t = 6000","t = 10000", "t = 75000"), tick = FALSE)


movedStones.dt <- read_csv("results/stone_movement.csv")
yearlyMean <- rowMeans(movedStones.dt)

par(family = "serif", mfrow = c(1,1), mar = c(6.1, 4.1, 4.1, 2.1))
plot(x = c(1:4999), y = yearlyMean, type = 'l', main = "Average Colony Deposition Rate Across Simulation Time", ylab = "Average Depositions/Round", xlab = "Time (Rounds)")

######
# SA #
######
library(ggplot2)
library(cowplot)
library(readr)

SA1.dt <- read_csv("SA/Franks and Deneubourg/FD - pars and statistics.csv")
# CoVCoL <- heat.colors(10)[as.numeric(cut(SA.dt$CoV, 10))] 

# CoV
P.CoV <- ggplot(SA1.dt, aes(x = Pmax, y = CoV, color = CoV)) +
  geom_point() +
  labs(x = "Pmax", y = "Distance spread (r)")
P.CoV + theme_bw()

F.CoV <- ggplot(SA1.dt, aes(x = Fmax, y = CoV, color = CoV)) +
  geom_point() +
  labs(x = "Fmax", y= "Distance spread (r)", color = "Distance spread")
F.CoV + scale_color_continuous() + theme_bw()

PF.CoV <- ggplot(SA1.dt, aes(x = Pmax, y = Fmax, color = CoV)) +
  geom_point() +
  labs(x = "Pmax", y = "Fmax", color = "Distance spread")
PF.CoV + scale_color_continuous() + theme_bw()

#Rbar
Gm.Rbar <- ggplot(SA1.dt, aes(x = Gmin, y = Rbar, color = Rbar)) +
  geom_point() +
  labs(x = "Gmin", y = "Circular spread\n(theta)", color = "Circular\nspread")
Gm.Rbar + scale_color_gradient(low = "yellow", high = "red") + theme_bw()

PM.Rbar <- ggplot(SA1.dt, aes(x = Pmax, y = Rbar, color = Rbar)) +
  geom_point() +
  labs(x = "Pmax", y = "Circular spread (theta)", color = "Circular spread")
PM.Rbar + scale_color_gradient(low = "yellow", high = "red") + theme_bw()

FM.Rbar <- ggplot(SA1.dt, aes(x = Fmax, y = Rbar, color = Rbar)) +
  geom_point() +
  labs(x = "Fmax", y = "Circular spread (theta)", color = "Circular spread")
FM.Rbar + scale_color_gradient(low = "yellow", high = "red") + theme_bw()

GM.Rbar <- ggplot(SA1.dt, aes(x = Gmax, y = Rbar, color = Rbar)) +
  geom_point() +
  labs(x = "Gmax", y = "Circular spread (theta)", color = "Circular spread")
GM.Rbar + scale_color_gradient(low = "yellow", high = "red") + theme_bw()

tau.Rbar <- ggplot(SA1.dt, aes(x = tau, y = Rbar, color = Rbar)) +
  geom_point() +
  labs(x = "tau", y = "Circular spread (theta)", color = "Circular spread")
tau.Rbar + scale_color_gradient(low = "yellow", high = "red") + theme_bw()



PMGm.Rbar <- ggplot(SA1.dt, aes(x = Pmax, y = Gmin, color = Rbar)) +
  geom_point() +
  labs(x = "Pmax", y = "Gmin", color = "Circular spread")
PMGm.Rbar + scale_color_continuous() + theme_bw()

FMGm.Rbar <- ggplot(SA1.dt, aes(x = Fmax, y = Gmin, color = Rbar)) +
  geom_point() +
  labs(x = "Fmax", y = "Gmin", color = "Circular spread")
FMGm.Rbar + scale_color_continuous() + theme_bw()

FmGm.Rbar <- ggplot(SA1.dt, aes(x = Fmin, y = Gmin, color = Rbar)) +
  geom_point() +
  labs(x = "Fmin", y = "Gmin", color = "Circular spread")
FmGm.Rbar + scale_color_continuous() + theme_bw()

tauGm.Rbar <- ggplot(SA1.dt, aes(x = tau, y = Gmin, color = Rbar)) +
  geom_point() +
  labs(x = "tau", y = "Gmin", color = "Circular spread")
tauGm.Rbar + scale_color_continuous() + theme_bw()

title <- ggdraw() + draw_label("Change in circular spread across parameter space", fontface = 'bold')
bottom_rows <- plot_grid(PMGm.Rbar + scale_color_gradient(low = "yellow", high = "red", guide = F) + theme_bw(), FMGm.Rbar + scale_color_gradient(low = "yellow", high = "red", guide = F) + theme_bw(), FmGm.Rbar + scale_color_gradient(low = "yellow", high = "red", guide = F) + theme_bw(), tauGm.Rbar + scale_color_gradient(low = "yellow", high = "red", guide = F) + theme_bw(), ncol = 2, labels = c("","","",""), scale = 0.9)
plot_grid(title, Gm.Rbar + scale_color_gradient(low = "yellow", high = "red") + theme_bw(), bottom_rows, nrow = 3, labels = c("","",""), rel_heights = c(0.2,1,2), scale = 0.9)


########################################################################################
#                                  rugatulus model                                     #
########################################################################################

library(readr)
stone_n_results <- read_csv("results/stone_n_results_with_rug.csv")
stone_n_results$model <- as.factor(stone_n_results$model)
View(stone_n_results)

lists <- read_csv("results/Rbar_and_CoV_values_with_rug.csv")

# CoV plot
FD1000 <- lists$SO1000_CoV
FD3000 <- lists$SO3000_CoV

T1000 <- lists$T1000_CoV
T3000 <- lists$T3000_CoV

r1000 <- lists$r1000_CoV
r3000 <- lists$r3000_CoV

# Rbar
FD1000R <- lists$SO1000_Rbar
FD3000R <- lists$SO3000_Rbar

T1000R <- lists$T1000_Rbar
T3000R <- lists$T3000_Rbar

r1000R <- lists$r1000_Rbar
r3000R <- lists$r3000_Rbar

opar = par()
par(mfrow = c(2,1), oma = c(2,0,3,0) + 0.1, mar = c(0,5,1,1) + 0.1)

boxplot(FD1000,T1000,r1000, FD3000,T3000,r3000, col = c("orange","red","dodgerblue"), xaxt = 'n', ylab = "distance spread (r)")
#axis(1, at=c(1.5,3.5,5.5), labels = c("stones = 1000","stones = 3000", "stones = 5000"), tick = FALSE)
abline(v = 3.5, lty = 2)

boxplot(FD1000R,T1000R,r1000R, FD3000R,T3000R,r3000R, ylim = c(-0.002,0.40), col = c("orange","red","dodgerblue"), xaxt = 'n', ylab = "circular spread (theta)")
axis(1, at=c(2,5), labels = c("stones = 1000","stones = 3000"), tick = FALSE)
abline(v = 3.5, lty = 2)
legend("topright", legend = c("Franks & Deneubourg\n(template + feedback)\n", "template only", "gradual model"), col = c("orange","red", "dodgerblue"), pch = 15, bty = 'n')

title(main = "Final stone distribution - model comparison", outer = TRUE, line = 1)

#############################
# delayed pop size increase #
#############################

library(readr)
dlists <- read_csv("results/delayed/Rbar_and_CoV_values_rug_delayed.csv")
r6000 <- dlists$SO_T6000
r10000 <- dlists$SO_T10000
r75000 <- dlists$SO_T75000

par = opar
par(mfrow = c(1,1), mar = c(3,5,3,1), oma = c(0,0,0,0))
boxplot(r6000, r10000, r75000, xaxt = 'n', ylab = "mean distance r of deposition locations from cluster (mm)", main = "Shift to new nest size over time", col = c("orange"))
axis(1, at=c(1,2,3), labels = c("t = 6000","t = 10000", "t = 75000"), tick = FALSE)

######
# SA #
######
library(ggplot2)
library(cowplot)
library(readr)

SA2.dt <- read_csv("SA/rug - pars and statistics.csv")
dev.off()
#mean and variance
boxplot(SA2.dt$CoV, SA2.dt$Rbar, xaxt = 'n', col = c('deepskyblue3', 'darkgoldenrod1'), ylab = "distribution across parameter space", cex.lab = 1.5)
abline(v = 1.5, lty = 2)
axis(1, at=c(1,2), labels = c("distance spread", "circular spread"), tick = F, cex.axis = 1.5)

# CoV
PM.CoV <- ggplot(SA2.dt, aes(x = Pmax, y = CoV, color = CoV)) +
  geom_point() +
  labs(x = "Pmax", y = "Distance spread (r)", color = "Distance\nspread")
PM.CoV + scale_color_continuous() + theme_bw()

FM.CoV <- ggplot(SA2.dt, aes(x = Fmax, y = CoV, color = CoV)) +
  geom_point() +
  labs(x = "Fmax", y = "Distance spread (r)", color = "Distance\nspread")
FM.CoV + scale_color_continuous() + theme_bw()

PMFM.CoV <- ggplot(SA2.dt, aes(x = Pmax, y = Fmax, color = CoV)) +
  geom_point() +
  labs(x = "Pmax", y = "Fmax", color = "Distance\nspread")
PMFM.CoV + scale_color_continuous() + theme_bw()
  

## for higher order interactions
# library(scatterplot3d) 
# colors <- heat.colors(10)[as.numeric(cut(SA.dt$Rbar, 10))]
# 
# # below the code for a double plot including the custon-made legend
# layout(matrix(1:2,nrow=1),widths=c(0.8,0.2))
# par(mar=c(5.1,4.1,4.1,2.1))
# scatterplot3d(x = SA.dt$Pmax, y = SA.dt$Dmax, z = SA.dt$tau, pch = 16, color = colors, xlab = "Pmax", ylab = "Dmax", zlab = "tau")
# 
# xl <- 1
# yb <- 1
# xr <- 1.5
# yt <- 2
# 
# par(mar=c(5.1,0.5,4.1,0.5))
# plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
# rect(
#   xl,
#   head(seq(yb,yt,(yt-yb)/10),-1),
#   xr,
#   tail(seq(yb,yt,(yt-yb)/10),-1),
#   col=heat.colors(10)
# )
# 
# mtext(sapply(seq(min(SA.dt$Rbar), max(SA.dt$Rbar)+0.01, length.out = 10), round, 3), side=2,at=tail(seq(yb,yt,(yt-yb)/10),-1)-0.05,las=2,cex=0.7)
# dev.off()


#######################################################################################
                            # Within-population variation #
#######################################################################################

########################
# Mutation in germline #
########################

library(readr)
library(ggplot2)

fn <- "C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/mutation in germline/Rq_sweep_3000stones.csv"
mut.dt <- read_csv(fn)

T.dt <- subset(mut.dt, model == "T")
View(T.dt)

T.CoV <- ggplot(T.dt, aes(x = q, y = R2, color = CoV_ave)) +
  geom_point() +
  labs(title = "Template model", x = "proportion of mutant workers", y = "R2 value (mm)", color = "distance spread")
T.CoV + scale_color_continuous() + theme_bw()

T.Rbar <- ggplot(T.dt, aes(x = q, y = R2, color = Rbar_ave)) +
  geom_point() +
  labs(title = "Template model",x = "proportion of mutant workers", y = "R2 value (mm)", color = "circular spread")
T.Rbar + scale_color_gradient(low = "yellow", high = "red") + theme_bw()


FD.dt <- subset(mut.dt, model == "FD")

FD.CoV <- ggplot(FD.dt, aes(x = q, y = R2, color = CoV_ave)) +
  geom_point() +
  labs(title = "FD model", x = "proportion of mutant workers", y = "R2 value (mm)", color = "distance spread")
FD.CoV + scale_color_continuous() + theme_bw()

FD.Rbar <- ggplot(FD.dt, aes(x = q, y = R2, color = Rbar_ave)) +
  geom_point() +
  labs(title = "FD model",x = "proportion of mutant workers", y = "R2 value (mm)", color = "circular spread")
FD.Rbar + scale_color_gradient(low = "yellow", high = "red") + theme_bw()


rug.dt <- subset(mut.dt, model == 'rug')
View(rug.dt)
max(rug.dt$CoV_SD)
max(rug.dt$Rbar_SD)

rug.CoV <- ggplot(rug.dt, aes(x = q, y = R2, color = CoV_ave)) +
  geom_point() +
  labs(title = "gradual model", x = "proportion of mutant workers", y = "R2 value (mm)", color = "distance spread")
rug.CoV + scale_color_continuous() + theme_bw()

rug.Rbar <- ggplot(rug.dt, aes(x = q, y = R2, color = Rbar_ave)) +
  geom_point() +
  labs(title = "rug model",x = "proportion of mutant workers", y = "R2 value (mm)", color = "circular spread")
rug.Rbar + scale_color_gradient(low = "yellow", high = "red") + theme_bw()


## visualise selected q cases
library(dplyr)

T.q10.dt <- subset(T.dt, q == 0.1)
T.q30.dt <- subset(T.dt, q == 0.3)
T.q50.dt <- subset(T.dt, q == 0.5)
FD.q10.dt <- subset(FD.dt, q == 0.1)
FD.q30.dt <- subset(FD.dt, q == 0.3)
FD.q50.dt <- subset(FD.dt, q == 0.5)
rug.q10.dt <- subset(rug.dt, q == 0.1)
rug.q30.dt <- subset(rug.dt, q == 0.3)
rug.q50.dt <- subset(rug.dt, q == 0.5)

Rdiff.v <- c(-10,-8,-6,-4,-2,0,+2,+4,+6)
Trug.subset <- rbind(T.q10.dt,T.q30.dt,T.q50.dt,rug.q10.dt,rug.q30.dt,rug.q50.dt)
Trug.subset$Rdiff <- rep(Rdiff.v,6)
Trug.subset$CoV_norm <- ifelse(Trug.subset$Rdiff == 0, 0, Trug.subset$CoV_ave/((Trug.subset$Rdiff)^2))
#zero.subset <- subset(Trug.subset, Rdiff ==0)

CoV.p <- ggplot(Trug.subset, aes(x = Rdiff, y = CoV_norm, color = model)) +
  geom_line()
CoV.p

p <- ggplot(Trug.subset, aes(x = Rdiff, color = model)) +
  geom_line(aes(y = CoV_ave), linetype = 'solid') +
  geom_ribbon(aes(ymin = CoV_ave - CoV_SD, ymax = CoV_ave + CoV_SD), linetype=0, alpha=0.1) +
  geom_line(aes(y = Rbar_ave/0.6), linetype = 'dashed') +
  geom_ribbon(aes(ymin = Rbar_ave/0.6 - Rbar_SD/0.6, ymax = Rbar_ave/0.6 + Rbar_SD/0.6), linetype=0, alpha=0.1) +
  scale_color_manual(name = "model", labels = c("gradual", "template-only"), values = c("#F8766D", "#00BFC4")) +
  #geom_point(data = zero.subset, aes(x = Rdiff, y = CoV_ave), color = 'gray', size = 3) +
  #geom_point(data = zero.subset, aes(x = Rdiff, y = Rbar_ave/0.6), color = 'gray', size = 3) +
  geom_vline(xintercept = 0,linetype = 5, size = 0.2, color = "Gray") +
  ylim(0,0.5) +
  scale_y_continuous(name = "distance dispersal", sec.axis = sec_axis(~.*0.6, name = "circular dispersal")) +
  labs(x = "difference between resident and mutant distance rule (mm)") +
  theme_bw() +
  theme(axis.title.y = element_text(hjust = 0.95, margin = unit(c(0, 5, 0, 0), "mm")), axis.title.y.right = element_text(margin = unit(c(0, 0, 0, 5), 'mm')), axis.title.x = element_text(margin = unit(c(5,0,0,0), 'mm')))
p + facet_grid(cols = vars(q), labeller = label_both)


########################
# Individual variation #
########################

library(readr)
library(ggplot2)

fn <- "C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/individual variation/ind_var_mean18_3000stones.csv"
var.dt <- read_csv(fn)

v.p <- ggplot(var.dt, aes(x = as.factor(sd), color = model)) +
  geom_point(aes(y = CoV_ave), shape = 16, size = 2.5) +
  geom_point(aes(y = Rbar_ave/0.8), shape = 15, size = 2.5) +
  scale_color_manual(name = "model", labels = c("gradual", "template-only"), values = c("#F8766D", "#00BFC4")) +
  scale_y_continuous(name = "distance dispersal", sec.axis = sec_axis(~.*0.8, name = "circular dispersal")) +
  labs(x = "standard deviation (mm) of worker population around mean distance rule") +
  theme_bw() +
  theme(axis.title.y = element_text(hjust = 0.95, margin = unit(c(0, 5, 0, 0), "mm")), axis.title.y.right = element_text(margin = unit(c(0, 0, 0, 5), 'mm')), axis.title.x = element_text(margin = unit(c(5,0,0,0), 'mm')))
v.p

ggplot(var.dt, aes(x = as.factor(sd), fill = model)) +
  geom_boxplot(aes(upper = CoV_ave + CoV_SD, lower = CoV_ave-CoV_SD, middle = CoV_ave)) +
  theme_bw()


# boxplot
library(gridExtra)
library(ggpubr)

sd50.dt <- read_csv("C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/individual variation/ind_var_list-sd50-mean18-3000stones.csv")
sd100.dt <- read_csv("C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/individual variation/ind_var_list-sd100-mean18-3000stones.csv")
sd150.dt <- read_csv("C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/individual variation/ind_var_list-sd150-mean18-3000stones.csv")
sd200.dt <- read_csv("C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/individual variation/ind_var_list-sd200-mean18-3000stones.csv")
sd500.dt <- read_csv("C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/individual variation/ind_var_list-sd500-mean18-3000stones.csv")
var2.dt <- rbind(sd50.dt,sd100.dt,sd150.dt,sd200.dt,sd500.dt)
var2.dt$SD <- as.factor(var2.dt$SD)

CoV.bp <- ggplot(var2.dt, aes(x = SD, y = CoV, fill = model)) +
  geom_boxplot() +
  ylim(0,0.45) +
  scale_fill_discrete(name = "model", labels = c("gradual", "template-only")) +
  labs(x = "", y = "distance dispersal") +
  theme_bw() +
  theme(legend.position = "none")
CoV.bp

Rbar.bp <- ggplot(var2.dt, aes(x = SD, y = Rbar, fill = model)) +
  geom_boxplot() +
  ylim(0,0.45) +
  scale_fill_discrete(name = "model", labels = c("gradual", "template-only")) +
  labs(x = "", y = "circular dispersal") +
  theme_bw() +
  theme(plot.margin = margin(l = 30, t = 6, r = 6, b = 5))
Rbar.bp

bottom <- text_grob("standard deviation (mm) of worker population around mean distance rule                        ", size = 11)
grid.arrange(CoV.bp, Rbar.bp, widths=c(0.38,0.62), bottom = bottom, ncol = 2)

######################################
# genetic variation among patrilines #
######################################

library(readr)
library(ggplot2)

fn <- "C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/multiple patrilines/four-patrilines_mean18_3000stones.csv"
pat.dt <- read_csv(fn)
pat.dt$sd <- as.factor(pat.dt$sd)
view(pat.dt)

f25f25.dt <- subset(pat.dt, f2 == 0.25)
f25.p <- ggplot(f25f25.dt, aes(x = as.factor(sd), color = model)) +
  geom_point(aes(y = CoV_ave), shape = 16, size = 2.5) +
  geom_point(aes(y = Rbar_ave), shape = 15, size = 2.5) +
  scale_color_manual(name = "model", labels = c("gradual", "template-only"), values = c("#F8766D", "#00BFC4")) +
  scale_y_continuous(name = "distance dispersal", sec.axis = sec_axis(~., name = "circular dispersal")) +
  labs(x = "standard deviation of the building rule Ropt distribution (mm)") +
  theme_bw() +
  theme(axis.title.y = element_text(hjust = 0.95, margin = unit(c(0, 5, 0, 0), "mm")), axis.title.y.right = element_text(margin = unit(c(0, 0, 0, 5), 'mm')), axis.title.x = element_text(margin = unit(c(5,0,0,0), 'mm')))
f25.p



f50f40.dt <- subset(pat.dt, f2 == 0.40)
f40.p <- ggplot(f50f40.dt, aes(x = as.factor(sd), color = model)) +
  geom_point(aes(y = CoV_ave), shape = 16, size = 2.5) +
  geom_point(aes(y = Rbar_ave), shape = 15, size = 2.5) +
  scale_color_manual(name = "model", labels = c("gradual", "template-only"), values = c("#F8766D", "#00BFC4")) +
  scale_y_continuous(name = "distance dispersal", sec.axis = sec_axis(~., name = "circular dispersal")) +
  labs(x = "standard deviation of the building rule Ropt distribution (mm)") +
  theme_bw() +
  theme(axis.title.y = element_text(hjust = 0.95, margin = unit(c(0, 5, 0, 0), "mm")), axis.title.y.right = element_text(margin = unit(c(0, 0, 0, 5), 'mm')), axis.title.x = element_text(margin = unit(c(5,0,0,0), 'mm')))
f40.p

f50f40.sd2 <- subset(f50f40.dt, sd = "2.0")
Rplot <- ggplot(f50f40.sd2, aes())



f50f20.dt <- subset(pat.dt, f2 == 0.20)
f20.p <- ggplot(f50f20.dt, aes(x = as.factor(sd), color = model)) +
  geom_point(aes(y = CoV_ave), shape = 16, size = 2.5) +
  geom_point(aes(y = Rbar_ave), shape = 15, size = 2.5) +
  scale_color_manual(name = "model", labels = c("gradual", "template-only"), values = c("#F8766D", "#00BFC4")) +
  scale_y_continuous(name = "distance dispersal", sec.axis = sec_axis(~., name = "circular dispersal")) +
  labs(x = "standard deviation of the building rule Ropt distribution (mm)") +
  theme_bw() +
  theme(axis.title.y = element_text(hjust = 0.95, margin = unit(c(0, 5, 0, 0), "mm")), axis.title.y.right = element_text(margin = unit(c(0, 0, 0, 5), 'mm')), axis.title.x = element_text(margin = unit(c(5,0,0,0), 'mm')))
f20.p



fn200 <- "C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/multiple patrilines/sdsimvals-sd200.csv"
pat200.dt <- read_csv(fn200)
pat200.dt$f2 <- as.factor(pat200.dt$f2)
View(pat200.dt)

fn500 <- "C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/multiple patrilines/sdsimvals-sd500.csv"
pat500.dt <- read_csv(fn500)
pat500.dt$f2 <- as.factor(pat500.dt$f2)
View(pat500.dt)

b200 <- ggplot(pat200.dt, aes(x = model, y = CoV, fill = f2)) +
  geom_boxplot() +
  ylim(0,0.5) +
  scale_x_discrete(labels = c("gradual model", "template-only model")) +
  labs(x = "", y = "distance dispersal") +
  scale_fill_discrete(name = "patriline frequencies\n(f1, f2, f3, f4)", labels = c("0.5, 0.2, 0.2, 0.1", "0.25, 0.25, 0.25, 0.25", "0.5, 0.4, 0.05, 0.05")) +
  theme_bw()
b200

b500 <- ggplot(pat500.dt, aes(x = model, y = CoV, fill = f2)) +
  geom_boxplot() +
  ylim(0,0.5) +
  scale_x_discrete(labels = c("gradual model", "template-only model")) +
  labs(x = "", y = "distance dispersal") +
  scale_fill_discrete(name = "patriline frequencies\n(f1, f2, f3, f4)", labels = c("0.5, 0.2, 0.2, 0.1", "0.25, 0.25, 0.25, 0.25", "0.5, 0.4, 0.05, 0.05")) +
  theme_bw()
b500


pat200500.dt <- rbind(pat200.dt,pat500.dt)
pat200500.dt$SD <- as.factor(pat200500.dt$SD)
b.CoV <- ggplot(pat200500.dt, aes(x = model, y = CoV, fill = f2)) +
  geom_boxplot() +
  ylim(0,0.5) +
  facet_grid(~SD, labeller = label_both) +
  scale_x_discrete(labels = c("gradual model", "template-only model")) +
  labs(x = "", y = "distance dispersal") +
  scale_fill_discrete(name = "patriline frequencies\n[f1, f2, f3, f4]", labels = c("[0.5, 0.2, 0.2, 0.1]", "[0.25, 0.25, 0.25, 0.25]", "[0.5, 0.4, 0.05, 0.05]")) +
  theme_bw()
b.CoV

b.Rbar <- ggplot(pat200500.dt, aes(x = model, y = Rbar, fill = f2)) +
  geom_boxplot() +
  ylim(0,0.5) +
  facet_grid(~SD, labeller = label_both) +
  scale_x_discrete(labels = c("gradual model", "template-only model")) +
  labs(x = "", y = "circular dispersal") +
  scale_fill_discrete(name = "patriline frequencies\n[f1, f2, f3, f4]", labels = c("[0.5, 0.2, 0.2, 0.1]", "[0.25, 0.25, 0.25, 0.25]", "[0.5, 0.4, 0.05, 0.05]")) +
  theme_bw()
b.Rbar


#####################################
# spatially specialised patrilines #
#####################################

library(readr)
library(ggplot2)
library(dplyr)

fn200 <- "C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/spatially specialised patrilines/spsp-sdsimvals-sd200.csv"
sppat200.dt <- read_csv(fn200)
sppat200.dt$f2 <- as.factor(sppat200.dt$f2)
sppat200.dt$seed <- rep(c(1:10),6)-1
View(sppat200.dt)

fn500 <- "C:/Users/FastRun/OneDrive - University of St Andrews/Projects/T_albipennis_ABM/results/spatially specialised patrilines/spsp-sdsimvals-sd500.csv"
sppat500.dt <- read_csv(fn500)
sppat500.dt$f2 <- as.factor(sppat500.dt$f2)
sppat500.dt$seed <- rep(c(1:10),6)-1
View(sppat500.dt)

b200.sp <- ggplot(sppat200.dt, aes(x = model, y = CoV, fill = f2)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("gradual model", "template-only model")) +
  labs(x = "", y = "distance dispersal") +
  scale_fill_discrete(name = "patriline frequencies\n(f1, f2, f3, f4)", labels = c("0.5, 0.2, 0.2, 0.1", "0.25, 0.25, 0.25, 0.25", "0.5, 0.4, 0.05, 0.05")) +
  theme_bw()
b200.sp

b500 <- ggplot(pat500.dt, aes(x = model, y = CoV, fill = f2)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("gradual model", "template-only model")) +
  labs(x = "", y = "distance dispersal") +
  scale_fill_discrete(name = "patriline frequencies\n(f1, f2, f3, f4)", labels = c("0.5, 0.2, 0.2, 0.1", "0.25, 0.25, 0.25, 0.25", "0.5, 0.4, 0.05, 0.05")) +
  theme_bw()
b500


sppat200500.dt <- rbind(sppat200.dt,sppat500.dt)
sppat200500.dt$SD <- as.factor(sppat200500.dt$SD)
View(sppat200500.dt)

un1 <- unique(sppat200500.dt[which(sppat200500.dt$CoV>0.35 & sppat200500.dt$model == "rug"), "seed"])
un2 <- unique(sppat200500.dt[which(sppat200500.dt$Rbar>0.2 & sppat200500.dt$model == "rug"), "seed"])
un3 <- unique(sppat200.dt[which(sppat200.dt$Rbar>0.07 & sppat200.dt$model == "rug" & sppat200.dt$f2 != 0.4), "seed"])

sppat500f40.dt <- subset(sppat500.dt, f2 == 0.4 & model == 'rug') 
median(sppat500f40.dt$CoV)
sppat500f25.dt <- subset(sppat500.dt, f2 == 0.25 & model == 'rug') 
median(sppat500f25.dt$CoV)
sppat500f20.dt <- subset(sppat500.dt, f2 == 0.2 & model == 'rug') 
median(sppat500f20.dt$CoV)

seedList <- c(2,8,5,4) # first two are outliers, third is close to distance disp median and fourth is hgih quality wall
points.dt <- subset(sppat200500.dt, seed %in% seedList )
points.dt$label <- c()
points.dt[which(points.dt$seed == 2), "label"] <- 1
points.dt[which(points.dt$seed == 8), "label"] <- 2
points.dt[which(points.dt$seed == 5), "label"] <- 3
points.dt[which(points.dt$seed == 4), "label"] <- 4


b.sp.CoV <- ggplot(sppat200500.dt, aes(x = model, y = CoV, fill = f2)) +
  geom_boxplot() +
  geom_point(data = points.dt, aes(x = model, y = CoV, fill = f2), colour = "black", position = position_dodge(width = 0.75), pch = 21) +
  geom_text(data = points.dt, aes(x = model, y = CoV, label = label), position = position_jitterdodge(jitter.width = 0.6, jitter.height = 0.005, seed = 0), size = 3) +
  ylim(0,0.5) +
  facet_grid(~SD, labeller = label_both) +
  scale_x_discrete(labels = c("gradual model", "template-only model")) +
  labs(x = "", y = "distance dispersal") +
  scale_fill_discrete(name = "patriline frequencies\n[f1, f2, f3, f4]", labels = c("[0.5, 0.2, 0.2, 0.1]", "[0.25, 0.25, 0.25, 0.25]", "[0.5, 0.4, 0.05, 0.05]")) +
  theme_bw()
b.sp.CoV

b.sp.Rbar <- ggplot(sppat200500.dt, aes(x = model, y = Rbar, fill = f2)) +
  geom_boxplot() +
  geom_point(data = points.dt, aes(x = model, y = Rbar, fill = f2), colour = "black", position = position_dodge(width = 0.75), pch = 21) +
  geom_text(data = points.dt, aes(x = model, y = Rbar, label = label), position = position_jitterdodge(jitter.width = 0.6, jitter.height = 0.005), size = 3.5) +
  ylim(0,0.5) +
  facet_grid(~SD, labeller = label_both) +
  scale_x_discrete(labels = c("gradual model", "template-only model")) +
  labs(x = "", y = "circular dispersal") +
  scale_fill_discrete(name = "patriline frequencies\n[f1, f2, f3, f4]", labels = c("[0.5, 0.2, 0.2, 0.1]", "[0.25, 0.25, 0.25, 0.25]", "[0.5, 0.4, 0.05, 0.05]")) +
  theme_bw()
b.sp.Rbar
