 ####################################################################
#### Code accompanying:
#
# Gomez-Gras et al. 2021. Climate change transforms the functional identity in Mediterranean coralligenous assemblages. Ecology Letters
#
# Code to calculate: 
# 1) Abundance distributions of traits in the multidimensional trait-space (Adapted from: Teixido et al.(2018). Functional biodiversity loss along natural CO2 gradients. Nature Communications) 
# 2) Temporal trends in functional identity (FI) in each site. 
# 3) Permutational analysis of Variance (PERMANOVA) to explore temporal changes in FI in each site  
#
# Script written by: Daniel Gomez-Gras & Viviana Brambilla 

# Functions required and written by Sebastien Villeguer: 
# quality_funct_space.R: This function was obtained in the Supplementary data of Teixido et al. 2018
#
#  
 
####################################################################

#This script needs to be run after 1.Functional space & Frich 

######################################

# set working directory
setwd("XXXX/")

####### Load libraries
library('FD')
library('tripack')
library('geometry')
library('matrixStats')
library('vegan')


###### Load data

ab <- read.csv2("./Raw_files/Abundances.csv", sep=";", dec=",", row.names=1)
fes_traits <- read.csv2("./Raw_files/FE_ordered.csv", sep=";", dec=",", row.names=1)
spe_fes <- read.csv2("./Raw_files/Species_FE.csv", sep=";", dec=",")
sites <- read.table("./Raw_files/Sites.txt", sep="\t", header=T, row.names=1)


############################### 1) Abundance distributions of traits in the functional space 

condition <- c("Pzzu_cor_2003", "Pzzu_cor_2011" , "Pzzu_cor_2018", "Pzzu_par_2006",
               "Pzzu_par_2011", "Pzzu_par_2018", "Gabin_par_1999", "Gabin_par_2007", "Gabin_par_2009",
               "Pzzinu_par_2006","Pzzinu_par_2011","Pzzinu_par_2016", "Passe_cor_2006","Passe_cor_2011","Passe_cor_2018")

##We load the coordinates generated previously when we created the functional space 

fd.coord<- read.csv2 ("FE_4D_coord.csv", sep=",", dec=",", row.names=1)
fd.coord<- as.matrix(fd.coord)

######## Data manipulation and arrangements

ab.conditions <- lapply(condition, function(x) {
  
  quad <- rownames(sites[sites$Year == x,])
  
  colSums(ab[rownames(ab) %in% quad,])
  
})#eo lapply


ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition

ab.conditions <- ab.conditions/2400*100 #number of quadrats 24 per condition and expressed as %

######### compute abundance of FEs for the different conditions 

fes <- levels(spe_fes$FE)

ab.fe.conditions <- lapply(condition, function (z) {
  
                       abund.fes <-  sapply(fes, function (x) {
  
                                            spec <- as.character(spe_fes[which(spe_fes$FE == x),]$Species)
                         
                                            sum(ab.conditions[z,spec])
                                                 
  
                                     })#eo sapply
                       
                       abund.fes

})#eo lapply

names(ab.fe.conditions) = condition

ab.fe.conditions <- do.call(rbind, ab.fe.conditions)

###################### 2) FUNCTIONAL IDENTITY (FI) TRENDS

###Obtaining FE_4D_SP_Coord from FE_4D_Coord that will be useful for further analysis

FE_4D_coord <- read.csv2 ("FE_4D_coord.csv", sep=",", dec=".", header =TRUE)
Species_FE <- read.csv2 ("./Raw_files/Species_FE.csv", sep=";", dec=".", header = TRUE)

new.df <- merge(FE_4D_coord,Species_FE, by.x = "X", by.y = "FE")

new.df <- new.df[,2:6] 
row.names(new.df) <- new.df$Species
new.df1 <-new.df[,1:4] 
new.df2 <- new.df1[ order(row.names(new.df1)), ]

coord<-new.df2
FE_4D_SP_coord <-coord
write.csv(coord, file="FE_4D_SP_coord.csv") #to use it for further analyses

##### Import Weights file to obtain weighted centroids of the assemblages in trait space (Functional identity) 

WEIGHTS. <- read.csv2 ("./Raw_files/Weights.csv", sep=";", dec=",", header = TRUE)

######Obtaining weighted centroids for each assemblage and Time point 

###################################Gabin_par

###Gabin_par_1999

weightedxPC1999 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Gabin_par_1999
xa<-sum(weightedxPC1999)
xPC1999 <- xa/100

weightedyPC1999 <-FE_4D_SP_coord$PC2 * WEIGHTS.$Gabin_par_1999
ya<-sum(weightedyPC1999)
yPC1999 <- ya/100

###Gabin_par_2007

weightedxPC2007 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Gabin_par_2007 
xb <- sum(weightedxPC2007)
xPC2007 <- xb/100

weightedyPC2007 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Gabin_par_2007
yb<- sum(weightedyPC2007)
yPC2007 <- yb/100

####Gabin_par_2009 
weightedxPC2009 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Gabin_par_2009 
xc <- sum(weightedxPC2009)
xPC2009 <- xc/100

weightedyPC2009 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Gabin_par_2009 
yc<- sum(weightedyPC2009)
yPC2009 <- yc/100


######################### Pzzu_cor 

###Pzzu_cor_2003 
weightedxPzzu_cor2003 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Pzzu_cor_2003
xd<-sum(weightedxPzzu_cor2003)
xPzzu_cor2003 <- xd/100

weightedyPzzu_cor2003 <-FE_4D_SP_coord$PC2 * WEIGHTS.$Pzzu_cor_2003
yd<-sum(weightedyPzzu_cor2003)
yPzzu_cor2003 <- yd/100

###Pzzu_cor_2011

weightedxPzzu_cor2011 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Pzzu_cor_2011
xe<-sum(weightedxPzzu_cor2011)
xPzzu_cor2011 <- xe/100

weightedyPzzu_cor2011 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Pzzu_cor_2011
ye<-sum(weightedyPzzu_cor2011)
yPzzu_cor2011 <- ye/100

###Pzzu_cor_2018

weightedxPzzu_cor2018 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Pzzu_cor_2018
xf<-sum(weightedxPzzu_cor2018)
xPzzu_cor2018 <- xf/100

weightedyPzzu_cor2018 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Pzzu_cor_2018
yf<-sum(weightedyPzzu_cor2018)
yPzzu_cor2018 <- yf/100


################# Pzzuu_par

###Pzzu_par_2006
weightedxPzzu_par2006 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Pzzu_par_2006
xg<-sum(weightedxPzzu_par2006)
xPzzu_par2006 <- xg/100

weightedyPzzu_par2006 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Pzzu_par_2006
yg<-sum(weightedyPzzu_par2006)
yPzzu_par2006 <- yg/100

###Pzzu_par_2011

weightedxPzzu_par2011 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Pzzu_par_2011
xh<-sum(weightedxPzzu_par2011)
xPzzu_par2011 <- xh/100

weightedyPzzu_par2011 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Pzzu_par_2011
yh<-sum(weightedyPzzu_par2011)
yPzzu_par2011 <- yh/100

###Pzzu_par_2018


weightedxPzzu_par2018 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Pzzu_par_2018
xi<-sum(weightedxPzzu_par2018)
xPzzu_par2018 <- xi/100

weightedyPzzu_par2018 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Pzzu_par_2018
yi<-sum(weightedyPzzu_par2018)
yPzzu_par2018 <- yi/100


############################# Pzzinu_par

###Pzzinu_par_2006 

xPalazzinu2006 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Pzzinu_par_2006
xj <- sum(xPalazzinu2006)
xPzzinu_par_2006 <- xj/100

yPalazzinu2006 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Pzzinu_par_2006
yj <- sum(yPalazzinu2006)
yPzzinu_par_2006 <- yj/100


####Pzzinu_par_2016

xPalazzinu2016 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Pzzinu_par_2016
xk <- sum(xPalazzinu2016)
xPzzinu_par_2016 <- xk/100

yPalazzinu2016 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Pzzinu_par_2016
yk <- sum(yPalazzinu2016)
yPzzinu_par_2016 <- yk/100


####Pzzinu_par_2011

xPalazzinu2011 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Pzzinu_par_2011
xl <- sum(xPalazzinu2011)
xPzzinu_par_2011 <- xl/100

yPalazzinu2011 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Pzzinu_par_2011
yl <- sum(yPalazzinu2011)
yPzzinu_par_2011 <- yl/100


######################## Passe_cor

###Passe_cor_2006 

xPassepalazzu2006 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Passe_cor_2006
xm <- sum(xPassepalazzu2006)
xPasse_cor_2006 <- xm/100

yPassepalazzu2006 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Passe_cor_2006
ym <- sum(yPassepalazzu2006)
yPasse_cor_2006 <- ym/100


####Passe_cor_2011

xPassepalazzu2011 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Passe_cor_2011
xn <- sum(xPassepalazzu2011)
xPasse_cor_2011 <- xn/100

yPassepalazzu2011 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Passe_cor_2011
yn <- sum(yPassepalazzu2011)
yPasse_cor_2011 <- yn/100


####Passe_cor_2018

xPassepalazzu2018 <- FE_4D_SP_coord$PC1 * WEIGHTS.$Passe_cor_2018
xo <- sum(xPassepalazzu2018)
xPasse_cor_2018 <- xo/100

yPassepalazzu2018 <- FE_4D_SP_coord$PC2 * WEIGHTS.$Passe_cor_2018
yo <- sum(yPassepalazzu2018)
yPasse_cor_2018 <- yo/100


# Figure S6. Relative abundance distribution of FE across the trait space with FI values 


#define number of axes

n_axes = 4

tiff(filename="Fig_S6.tif", height=27, width=15, units="cm", compression = c("lzw"), res=300, pointsize=8)

par(mfrow=c(5,3))


###Potential functional space (background polygon)

m2 <- fd.coord[rownames(fd.coord),]
tr2 <-tri.mesh(m2[,1],m2[,2])
ch2 <- convex.hull(tr2)   

##Pzzu_cor
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzu_cor_2003", col.main="#CD2626")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#CD2626" , bg="#CD262670", cex=sqrt(ab.fe.conditions["Pzzu_cor_2003",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPzzu_cor2003,yPzzu_cor2003, pch="+", col="red", cex=4)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzu_cor_2011", col.main="#F9D71C")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#F9D71C" , bg="#F9D71C70", cex=sqrt(ab.fe.conditions["Pzzu_cor_2011",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPzzu_cor2011,yPzzu_cor2011,pch="+", col="yellow", cex=4)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzu_cor_2018", col.main="#3A5FCD")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#3A5FCD" , bg="#3A5FCD70", cex=sqrt(ab.fe.conditions["Pzzu_cor_2018",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPzzu_cor2018,yPzzu_cor2018,pch="+", col="blue", cex=4)


###Pzzu_par

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzu_par_2006", col.main="#CD2626")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#CD2626" , bg="#CD262670", cex=sqrt(ab.fe.conditions["Pzzu_par_2006",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPzzu_par2006,yPzzu_par2006, pch="+", col="red", cex=4)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzu_par_2011", col.main="#F9D71C")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#F9D71C" , bg="#F9D71C70", cex=sqrt(ab.fe.conditions["Pzzu_par_2011",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPzzu_par2011,yPzzu_par2011,pch="+", col="yellow", cex=4)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzu_par_2018", col.main="#3A5FCD")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#3A5FCD" , bg="#3A5FCD70", cex=sqrt(ab.fe.conditions["Pzzu_par_2018",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPzzu_par2018,yPzzu_par2018,pch="+", col="blue", cex=4)


##Gabin_par 
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Gabin_par_1999", col.main="#CD2626")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#CD2626" , bg="#CD262670", cex=sqrt(ab.fe.conditions["Gabin_par_1999",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPC1999,yPC1999, pch="+", col="red", cex=4)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Gabin_par_2007", col.main="#F9D71C")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#F9D71C" , bg="#F9D71C70", cex=sqrt(ab.fe.conditions["Gabin_par_2007",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPC2007,yPC2007,pch="+", col="yellow", cex=4)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Gabin_par_2009", col.main="#3A5FCD")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#3A5FCD" , bg="#3A5FCD70", cex=sqrt(ab.fe.conditions["Gabin_par_2009",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPC2009,yPC2009,pch="+", col="blue", cex=4)


##Pzzinu_par 
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzinu_par_2006", col.main="#CD2626")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#CD2626" , bg="#CD262670", cex=sqrt(ab.fe.conditions["Pzzinu_par_2006",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPzzinu_par_2006,yPzzinu_par_2006, pch="+", col="red", cex=4)


plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzinu_par_2011", col.main="#F9D71C")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#F9D71C" , bg="#F9D71C70", cex=sqrt(ab.fe.conditions["Pzzinu_par_2011",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPzzinu_par_2011,yPzzinu_par_2011, pch="+", col="yellow", cex=4)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzinu_par_2016", col.main="#3A5FCD")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#3A5FCD" , bg="#3A5FCD70", cex=sqrt(ab.fe.conditions["Pzzinu_par_2016",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPzzinu_par_2016,yPzzinu_par_2016,pch="+", col="blue", cex=4)


### Passe_cor

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Passe_cor_2006", col.main="#CD2626")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#CD2626" , bg="#CD262670", cex=sqrt(ab.fe.conditions["Passe_cor_2006",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPasse_cor_2006,yPasse_cor_2006, pch="+", col="red", cex=4)


plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Passe_cor_2011", col.main="#F9D71C")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#F9D71C" , bg="#F9D71C70", cex=sqrt(ab.fe.conditions["Passe_cor_2011",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPasse_cor_2011,yPasse_cor_2011,pch="+", col="yellow", cex=4)


plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Passe_cor_2016", col.main="#3A5FCD")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#3A5FCD" , bg="#3A5FCD70", cex=sqrt(ab.fe.conditions["Passe_cor_2018",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPasse_cor_2018,yPasse_cor_2018, pch="+", col="blue", cex=4)

dev.off()


#####START-END plot; Trait abundance distributions of traits with FI values (Figure_3a_g)

#Potential functional space
m2 <- fd.coord[rownames(fd.coord),]
tr2 <-tri.mesh(m2[,1],m2[,2])
ch2 <- convex.hull(tr2)   

n_axes = 4

tiff(filename="Figure_3a_g.tif", height=12, width=16, units="cm", compression = c("lzw"), res=300, pointsize=8)

par(mfrow=c(2,3))

##Pzzu_cor
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzu_cor_2003_2018", col.main="#CD2626")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#CD2626" , bg="#CD262670", cex=sqrt(ab.fe.conditions["Pzzu_cor_2003",]/100*30)*3)
points(fd.coord[,1], fd.coord[,2], pch=21, col="#3A5FCD" , bg="#3A5FCD70", cex=sqrt(ab.fe.conditions["Pzzu_cor_2018",]/100*30)*3)#text(fd.coord[,1], fd.coord[,2],labels = rownames(fd.coord), cex = 0.5)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPzzu_cor2003,yPzzu_cor2003, pch="+", col="red", cex=3.5)
points(xPzzu_cor2018,yPzzu_cor2018,pch="+", col="blue", cex=3.5)
points(xPzzu_cor2011,yPzzu_cor2011,pch="+", col="yellow", cex=3.5)


###Pzzu_par

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzu_par_2006-2018", col.main="#CD2626")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#CD2626" , bg="#CD262670", cex=sqrt(ab.fe.conditions["Pzzu_par_2006",]/100*30)*3)
points(fd.coord[,1], fd.coord[,2], pch=21, col="#3A5FCD" , bg="#3A5FCD70", cex=sqrt(ab.fe.conditions["Pzzu_par_2018",]/100*30)*3)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPzzu_par2006,yPzzu_par2006, pch="+", col="red", cex=3.5)
points(xPzzu_par2018,yPzzu_par2018,pch="+", col="blue", cex=3.5)
points(xPzzu_par2011,yPzzu_par2011,pch="+", col="yellow", cex=3.5)


##Gabib_par
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Gabin_par_1999-2009", col.main="#CD2626")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#CD2626" , bg="#CD262670", cex=sqrt(ab.fe.conditions["Gabin_par_1999",]/100*30)*3)
points(fd.coord[,1], fd.coord[,2], pch=21, col="#3A5FCD" , bg="#3A5FCD70", cex=sqrt(ab.fe.conditions["Gabin_par_2009",]/100*30)*3)#text(fd.coord[,1], fd.coord[,2],labels = rownames(fd.coord), cex = 0.5)
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(xPC1999,yPC1999, pch="+", col="red", cex=3.5)
points(xPC2009,yPC2009,pch="+", col="blue", cex=3.5)
points(xPC2007,yPC2007,pch="+", col="yellow", cex=3.5)

####Passe_cor
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Passe_cor_2006-2018", col.main="#CD2626")
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(fd.coord[,1], fd.coord[,2], pch=21, col="#CD2626" , bg="#CD262670", cex=sqrt(ab.fe.conditions["Passe_cor_2006",]/100*30)*3)
points(fd.coord[,1], fd.coord[,2], pch=21, col="#3A5FCD" , bg="#3A5FCD70", cex=sqrt(ab.fe.conditions["Passe_cor_2018",]/100*30)*3)
points(xPasse_cor_2006,yPasse_cor_2006, pch="+", col="red", cex=3.5)
points(xPasse_cor_2011,yPasse_cor_2011,pch="+", col="yellow", cex=3.5)
points(xPasse_cor_2018,yPasse_cor_2018, pch="+", col="blue", cex=3.5)


###Pzzinu_par 
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzinu_par_2006-2016", col.main="#CD2626")
polygon(ch2, col="#CCCCCC30", border=FALSE)
points(fd.coord[,1], fd.coord[,2], pch=21, col="#CD2626" , bg="#CD262670", cex=sqrt(ab.fe.conditions["Pzzinu_par_2006",]/100*30)*3)
points(fd.coord[,1], fd.coord[,2], pch=21, col="#3A5FCD" , bg="#3A5FCD70", cex=sqrt(ab.fe.conditions["Pzzinu_par_2016",]/100*30)*3)
points(xPzzinu_par_2006,yPzzinu_par_2006, pch="+", col="red", cex=3.5)
points(xPzzinu_par_2016,yPzzinu_par_2016,pch="+", col="blue", cex=3.5)
points(xPzzinu_par_2011,yPzzinu_par_2011, pch="+", col="yellow", cex=3.5)

dev.off()

######### 3) PERMUTATIONAL ANALYSIS OF VARIANZA (PERMANOVA) WITH FUNCTIONAL IDENTITY TRENDS

##We select the first two axes (PCOA1 and PCOA2) from FE_4D_SP_coord to obtain x and y coordinates 

FE_2D_SP_coord = FE_4D_SP_coord[,1:2]

##We obtain the weighted centroids (x,y) of the community for each quadrat.

ab_P <- t(ab)
weighted_x = FE_2D_SP_coord$PC1 * ab_P
xsum <- colSums(weighted_x)
weighted_x_coordinates <- xsum/100

weighted_y = FE_2D_SP_coord$PC2 * ab_P
ysum <- colSums(weighted_y)
weighted_y_coordinates <- ysum/100


### We creat a data frame that contains the obtained weighted "x" "y" coordinates. 

weighted_centroids <- cbind(weighted_x_coordinates,weighted_y_coordinates)
weighted_centroids = as.data.frame(weighted_centroids)

###We complete the data frame with more useful info for the analysis.

##Site
weighted_centroids$Site = c(1:120)
weighted_centroids$Site[1:72] = "Pzzu_cor"
weighted_centroids$Site[73:144] = "Pzzu_par"
weighted_centroids$Site[145:216] = "Gabin_par"
weighted_centroids$Site[217:288] = "Pzzinu_par"
weighted_centroids$Site[289:360] = "Passe_cor"

##Time
weighted_centroids$Time = c(1:120)
weighted_centroids$Time[c(1:24, 73:96,144:168,216:240,288:312)] = "Time_1"
weighted_centroids$Time[c(25:48, 97:120, 169:192, 241:264,313:336)] = "Time_2"
weighted_centroids$Time[c(49:72, 121:144, 193:216, 265:288, 337:360)] = "Time_3"

##MHW-Impacted vs Not-impacted

weighted_centroids$Impact = c(1:120)
weighted_centroids$Impact[c(1:216)] = "MHW-Impacted"
weighted_centroids$Impact[c(217:360)] = "Not-Impacted"

### We use this new data frame to conduct the PERMANOVA ANALYSIS 


data <-weighted_centroids

###################Not-impacted sites 

####Passe_cor 

Passe_cor = data[data$Site =="Passe_cor",]
per_Passe_cor <- adonis2(cbind(Passe_cor$weighted_x_coordinates, Passe_cor$weighted_y_coordinates) ~ Passe_cor$Time, method="canberra")
per_Passe_cor

####Pzzinu_par
Pzzinu_par = data[data$Site=="Pzzinu_par",]
per_Pzzinu_par <- adonis2(cbind(Pzzinu_par$weighted_x_coordinates,Pzzinu_par$weighted_y_coordinates) ~ Pzzinu_par$Time, method="canberra")
per_Pzzinu_par

###################MHW-impacted sites 

###Pzzu_cor
Pzzu_cor = data[data$Site =="Pzzu_cor",]
per_Pzzu_cor <- adonis2(cbind(Pzzu_cor$weighted_x_coordinates,Pzzu_cor$weighted_y_coordinates) ~ Pzzu_cor$Time, method="canberra")
per_Pzzu_cor

###Pzzu_par
Pzzu_par  = data[data$Site=="Pzzu_par",]
per_Pzzu_par <- adonis2(cbind(Pzzu_par$weighted_x_coordinates,Pzzu_par$weighted_y_coordinates) ~ Pzzu_par$Time, method="canberra")
per_Pzzu_par

###Gabin_par
Gabin_par = data[data$Site=="Gabin_par",]
per_Gabin_par <- adonis2(cbind(Gabin_par$weighted_x_coordinates, Gabin_par$weighted_y_coordinates) ~ Gabin_par$Time, method="canberra")
per_Gabin_par




