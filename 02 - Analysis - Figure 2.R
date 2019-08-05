# R code from the manuscript

# Title: "A multi-locus demo-genetic simulator to model local adaptation in heterogeneous landscapes"
# Authors: Marco Andrello, Christelle Noirot, Florence Débarre and Stéphanie Manel
# Submitted to: Methods in Ecology and Evolution
# Date: 05 Aug 2019

# "02 - Analysis - Figure 2.R"
# Analyse the results of the demo-genetic simulations on Mullus surmuletus
# Global patterns presented in Figure 2
# Marco Andrello
# 05/08/2019


rm(list=ls())

library(maps)
library(mapdata)
library(RColorBrewer)
library(rgdal)

# Load only useful data
load("Simulation_parameters.RData")
rm(al, m, Proba, kappa0, delta, List_gene, mu, n, phi_F, phi_M, r, sigma_M, T_max,z)

load("Results_base.RData")
sigma_ij <- t(surv.trans) # NB : sigma_ij is the survival of individuals originating in i (row) transplanted in j (column)
rm(fst_A, fst_B, fst_C, fst_D, Globfst, Nt, Nt_A, Nt_B, Nt_C, Nt_D, ro, surv.trans)

deme <- readOGR(dsn="Demes.shp")

load("Distances.RData")

load("Salinity_data.RData")


# CalculatING suitability S_i
S_i <- apply(sigma_F[,,1,1],2,max) 


# Calculating useful metrics
v_ik <- AL * (1/rowSums(AL))                            # 1 - Phenotype frequency
u_mean <- as.vector(v_ik %*% c(0:6))                    # 2 - Mean number of adaptive alleles per individual (in each deme)
u_var <- as.vector(v_ik %*% c(0:6)^2) - u_mean^2        # Variance of the number of adaptive alleles per individual (in each deme)
u_CV <- sqrt(u_var)/u_mean                              # 3 - CV of the number of adaptive alleles per individual (in each deme)
sigma_i <- diag(sigma_ij)                               # 4 - Mean survival in each deme
A_i <- sigma_i/S_i                                      # 5 - Adaptation for each deme i



###############################################################################
# Spatial patterns in number of "+" alleles
# Figure 2
###############################################################################

# Regression of u_mean on Salinity, in demes with salinity > 36
id <- which(Salinity_per_deme[3,]>36)
id <- which(deme@data$Salinity>36)
mod<-lm(u_mean[id]~deme@data$Salinity[id])
summary(mod)

# Plot Figure 2a
png("Figure 2a.png",width=56,height=55,units="mm",res=300)
par(mar=c(3,4,1,1))
plot(deme@data$Salinity,u_mean,xlab="",ylab="",axes=F,pch=21,bg="gray",cex=0.75,col="black")
axis(1,at=c(33,35,37,39),cex.axis=0.70,lwd=0.5)
mtext("Salinity (PSU)",1, line=2, cex=0.75)
axis(2,las=2,at=c(0:6), cex.axis=0.75,lwd=0.5)
mtext("Mean number\nof \"+\" alleles",2, line=2.5, cex=0.75)
abline(coef(mod),lty=1,col="red")
dev.off()

# Plot Figure 2b
png("Figure 2b.png",width=56,height=55,units="mm",res=300)
par(mar=c(3,4,1,1))
plot(deme@data$Salinity,u_CV,xlab="",ylab="",axes=F,pch=21,bg="gray",cex=0.75,col="black")
axis(1,at=c(33,35,37,39),cex.axis=0.70,lwd=0.5)
mtext("Salinity (PSU)",1, line=2, cex=0.75)
axis(2,las=2,at=seq(0,3.5,0.5), cex.axis=0.75,lwd=0.5);
mtext("CV of number\nof \"+\" alleles",2, line=2.5, cex=0.75)
dev.off()

# Plot Figure 2c
col.level <- cut(A_i,seq(0,1,0.25))
colors.map <- brewer.pal(4,"Reds")
png("Figure 2c.png",width=112,height=55,units="mm",res=600,pointsize=6)
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47), col="gray90", fill=TRUE, mar=c(1,1,1,1), lwd=0.25)
plot(deme,col=colors.map[col.level],add=T,lwd=0.25)
legend(-8.2,33,legend=levels(col.level),fill=colors.map, cex=1,xpd=TRUE,
       title="Adaptation", box.lwd=0.25, ncol=2)
dev.off()




