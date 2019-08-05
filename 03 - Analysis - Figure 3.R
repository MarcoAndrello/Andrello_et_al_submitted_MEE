# R code from the manuscript

# Title: "A multi-locus demo-genetic simulator to model local adaptation in heterogeneous landscapes"
# Authors: Marco Andrello, Christelle Noirot, Florence Débarre and Stéphanie Manel
# Submitted to: Methods in Ecology and Evolution
# Date: 05 Aug 2019

# "03 - Analysis - Figure 3.R"
# Analyse the patterns of adaptation in the Eastern Gulf of Lion and Dardanelles
# of the demo-genetic simulations on Mullus surmuletus
# Patterns EGL and DAR presented in Figure 3 of the manuscript
# Marco Andrello
# 05/08/2019


rm(list=ls())
library(rgdal)
library(rgeos)
library(maps)
library(mapdata)
library(RColorBrewer)

# Function to jitter points in the good direction
spJitter <- function(coord.d,coord.s) {
    x <- rep(coord.d,length(coord.s))
    diff <- coord.s-x
    x + 0.1*diff
}


# Load data
deme <- readOGR("Demes.shp")
load("Dispersal_SC.RData")
d <- Dispersal_SC; rm(Dispersal_SC)



###############################################################################

# To produce the figure of Eastern Gulf of Lion (Figure 3a and 3b)
# execute the "Description of Eastern Gulf of Lion" and then jump to PRODUCE FIGURES

# To produce the figure of Dardanelles (Figure 3c and 3d)
# execute the "Description of Dardanelles" and then jump to PRODUCE FIGURES

###############################################################################


###############################################################################
# Description of Eastern Gulf of Lion
# Figure 3a and 3b
id.d <- which(deme$Name=="Eastern Gulf of Lion")
xlim.Map <- c(0,15)
ylim.Map <- c(40,46)
xlim.Leg<- -0.15
ylim.Leg <- 41.1

###############################################################################


###############################################################################
# Description of Dardanelles
# Figure 3c and 3d
id.d <- which(deme$Name=="Dardanelles")
xlim.Map <- c(22,28)
ylim.Map <- c(38,41)
xlim.Leg<- 21.9
ylim.Leg <- 38.5

###############################################################################

###############################################################################
# PRODUCE FIGURES
###############################################################################

# 
# Figure Maps
# Map: flows: Wheight dispersal probabilities by number of individuals in source demes

# Calculating flow: dispersal probability times reproductive output
load("Results_base.RData")
f <- array(0,c(100,100))
for (j in 1 : 100) {
    f[,j] <- d[,j] * ro[j]
}
which(f[id.d,]>0)
id.s <- which(f[id.d,]>0)
id.s.noLocRet <- id.s[-which(id.s==id.d)]

# Extract centroids to draw arrows
# Exclude the destination deme to improve plotting
p <- gCentroid(deme,byid=T)
lon.s <- p@coords[id.s.noLocRet,1]
lat.s <- p@coords[id.s.noLocRet,2]
lon.d <- p@coords[id.d,1]
lat.d <- p@coords[id.d,2]
 
# Stength of flows control the line type of arrows
d_strength <- f[id.d,which(f[id.d,]>0)]
bins <- c(1,20,200,2000,20000)
bins <- c(0.2,20,200,2000,20000)
cuts  <- cut(d_strength, breaks=bins,include.lowest=T)
LineType <- c(3,4,2,1)[cuts]

# Demes are colored according to their identity (one color per deme)
cols <- brewer.pal(length(id.s),"Set3") 
col.d <- rep("gray95",length(deme))
col.d[id.s] <- cols

# Plot map
png(paste0("Figure 3_Map_",id.d,".png"),width=112,height=55,units="mm",res=600,pointsize=6)
map("worldHires", xlim=xlim.Map, ylim=ylim.Map, col="gray90", fill=T, mar=c(1,1,1,1), lwd=0.25)
plot(deme,col=col.d,add=T,lwd=0.25)
lon.d.j <- spJitter(lon.d,lon.s)
lat.d.j <- spJitter(lat.d,lat.s)
arrows(lon.s, lat.s, lon.d.j, lat.d.j, lty=LineType, col="black", lwd=1, length=0.05, angle=30 )
legend(xlim.Leg, ylim.Leg, legend=levels(cuts), lty=c(3,4,2,1), cex=1, ncol=2, title="Larval flow", box.lwd=0.25)
dev.off()

names(d_strength) <- which(f[id.d,]>0)


########
# Figure Histograms
# Distribution of genotype classes in these demes, and survival probabilities of the receiving deme 
pk <- AL * (1/rowSums(AL))

load("Salinity_data.RData")
xi <- seq(36,39,0.5) # Optima i.e. phenotype
x <- Salinity_per_deme[3,]
y <- array(NA,c(100,7))
for (i in 1 : 100) {
    for (k in 1 : 7) {
        y[i,k] <- 0.8 * exp( -(x[i] - xi[k])^2 / 2 )
    }
}

# Stack barplots
height <- array(NA,c(length(d_strength),7)) # 7 is the number of phenotypes
for (i in 1 : length(d_strength)) {
    height[i,] <- pk[id.s[i],] * d_strength[i]
}
height <- height / sum(height)


png(paste0("Figure 3_Hist_",id.d,".png"),width=56,height=55,units="mm",res=600,pointsize=6)
par(mar=c(5, 4, 4, 6) + 0.1)
barplot(height,main=deme$Name[id.d],names.arg = xi,xlab="Phenotype (optimum salinity)",
    ylim=c(0,0.85),space=0,axes=F,col=cols)
axis(2,las=2)
mtext("Frequency",2,line=3)
if (id.d == 48) { # For DAR a second axis is necessary
    points(c(1:7)-0.5,y[id.d,],pch=16,cex=1)
} else {
    points(c(1:7)-0.5,y[69,]*500,pch=16,cex=1)
    axis(4,las=2,at=seq(0,0.8,0.2),labels=c("0","0.0004","0.0008","0.0012","0.0016"))
    mtext("Survival",4,line=4)
}
dev.off()


png(paste0("Figure 3_Leg_",id.d,".png"),width=56,height=55,units="mm",res=600,pointsize=6)
plot(1,xlim=c(0,1),ylim=c(0,1),type="n",axes=F,xlab="",ylab="")
legend(0,1,deme$Name[id.s],pch=15,col=cols,cex=1,bty="n",title="Sources")
dev.off()



