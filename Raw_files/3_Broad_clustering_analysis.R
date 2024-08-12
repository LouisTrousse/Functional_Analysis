####################################################################
#### Code accompanying:
#
# Gomez-Gras et al. 2021. Climate change transforms the functional identity in Mediterranean coralligenous assemblages. Ecology Letters

# Code to compute: 
  #1) Broad functional classification into functional clusters (with the overall pool). Figure 4a
  #2) Broad functional classification into functional clusters (per site). Figure 4b-f
  #3) Number of species per functional claster in each site and Time point. Figure 4g-k
  #4) Temporal abundance changes in functional clusters per site. Figure 4l-p
#
# Script written by: Daniel Gomez-Gras & Viviana Brambilla 

####################################################################

# This script needs to be run after: 

#1. Functional space & Frich
#2. Abundance distribution of traits & FI

################################# Load libraries 

library(cluster) # For the Gower distance matrix and PAM
library(ggplot2) # For plotting the data
library(matrixStats) #For some data manipulation
library(dplyr) # For some data arrangements 
library(tidyr) # For some data arrangements 
library(stringr) # For some data arrangements 


######## 1) BROAD CLASSIFICATION IN FUNCTIONAL CLUSTERS

# Convert the distance list into a matrix
gower_mat = as.matrix(gower)

# How many clusters? 
# We will iterate through 2-51 clusters (52 FE - 1) and see which shows the maximum separation between centres
# The "Silhoutte Width" parameter is used by default for PAM clusters

sil_width = c(NA)
for(i in 2:51)
{
  pam_fit = pam(gower, diss = TRUE, k=i)
  sil_width[i] = pam_fit$silinfo$avg.width
}

# We plot the silhoutte widths and see how many clusters is optimal (maximum value of sil_width): Figure_S6

tiff(filename ="Figure_S5.tif", height=8, width=9, units="cm", compression = c("lzw"), res=300, pointsize=8)
plot(1:51,sil_width,xlab="Number of Clusters",ylab="Silhouette width")
lines(1:51, sil_width)
#After 2 clusters, 8 is the optimum number of broad clusters in our case because it has the maximum value of sil_width.
###We discard values higher than 20 because we are interested in a broad classification 
abline(h = 0, v = 8, col = "red")
dev.off()

# Create the clusters using PAM for 8 clusters

pam_fit = pam(gower, diss=TRUE, k=8)


##########Create the plot. Figure 4.a

# We add the "cluster N" column to our data with coordinates. This will allow us to know to which cluster has been assigned every FE. 
###We will use these data later to calculate the N of sp per cluster (functional redundancy) and if there are temporal changes in the
##relative abundance of functional clusters in the different assemblages. 


#We need fd.coord as a list again. We can directly create again the fd.coord object as in script 1. Functional space and Frich
fd.coord <- qfs$details_funct_space$mat_coord[,1:4]

pam_fit$clustering = as.data.frame(pam_fit$clustering)
binded <- cbind(pam_fit$clustering, fd.coord)
colnames(binded)[which(names(binded) == "pam_fit$clustering")] <- "groups"

#####plotting with colors 

#first we create the convex hulls of the clusters 

m2 <- fd.coord[rownames(fd.coord),]
tr2 <-tri.mesh(m2[,1],m2[,2])
ch2 <- convex.hull(tr2)   

m3 <- fd.coord[rownames(fd.coord),]
m3 <- cbind(m3,binded)
m3 <-m3[,1:5]
m3 <- m3[rownames(m3),]
m3 <- m3[m3$groups=="1",]
tr3 <-tri.mesh(m3[,1],m3[,2])
ch3 <- convex.hull(tr3)   

m4 <- fd.coord[rownames(fd.coord),]
m4 <- cbind(m4,binded)
m4 <-m4[,1:5]
m4 <- m4[rownames(m4),]
m4 <- m4[m4$groups=="2",]
tr4 <-tri.mesh(m4[,1],m4[,2])
ch4 <- convex.hull(tr4) 

m5 <- fd.coord[rownames(fd.coord),]
m5 <- cbind(m5,binded)
m5 <-m5[,1:5]
m5 <- m5[rownames(m5),]
m5 <- m5[m5$groups=="3",]
tr5 <-tri.mesh(m5[,1],m5[,2])
ch5 <- convex.hull(tr5)   

m6 <- fd.coord[rownames(fd.coord),]
m6 <- cbind(m6,binded)
m6 <-m6[,1:5]
m6 <- m6[rownames(m6),]
m6 <- m6[m6$groups=="4",]
tr6 <-tri.mesh(m6[,1],m6[,2])
ch6 <- convex.hull(tr6)   

m7 <- fd.coord[rownames(fd.coord),]
m7 <- cbind(m7,binded)
m7 <-m7[,1:5]
m7 <- m7[rownames(m7),]
m7 <- m7[m7$groups=="5",]
tr7 <-tri.mesh(m7[,1],m7[,2])
ch7 <- convex.hull(tr7)   

m8 <- fd.coord[rownames(fd.coord),]
m8 <- cbind(m8,binded)
m8 <-m8[,1:5]
m8 <- m8[rownames(m8),]
m8 <- m8[m8$groups=="6",]
tr8 <-tri.mesh(m8[,1],m8[,2])
ch8 <- convex.hull(tr8)   

m9 <- fd.coord[rownames(fd.coord),]
m9 <- cbind(m9,binded)
m9 <-m9[,1:5]
m9 <- m9[rownames(m9),]
m9 <- m9[m9$groups=="7",]
tr9 <-tri.mesh(m9[,1],m9[,2])
ch9 <- convex.hull(tr9)   


#We create the color palette 
index <- c("1", "2","3", "4", "5", "6", "7", "8")
valuesC <- c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple")
binded$color <-valuesC[match(binded$groups,index)]


####### We create Figure 4a

tiff(filename ="Figure_4a.tif", height=8, width=9, units="cm", compression = c("lzw"), res=300, pointsize=8)

par(xpd = T, mar = par()$mar + c(0,0,0,5))
plot(binded[,2], binded[,3], xlab="PCoA1", ylab="PCoA2", type="n", main="Clusters", col.main="black")
points(binded[,2], binded[,3], pch=21, col=binded$color , bg=binded$color, cex=1)
legend(0.55, 0.3, 
       c("Cluster 1", "Cluster 2","Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8"),
       pch=21, col = valuesC, pt.bg=valuesC, 
       cex = 0.8)

#We add the polygons to the figure
polygon(ch2, col="#CCCCCC10", border=FALSE)
polygon(ch3, col="#20202060", border=FALSE)
polygon(ch4, col="#FF666680", border=FALSE)
polygon(ch5, col="#B2FF6680", border=FALSE)
polygon(ch6, col="#FFFF9980", border=FALSE)
polygon(ch7, col="#66FFFF80", border=FALSE)
polygon(ch8, col="#00666680", border=FALSE)
polygon(ch9, col="#FF800080", border=FALSE)

par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()


################ 2) PLOTTING FUNCTIONAL CLUSTERS PER SITE.

######Data arragments and subsetting every site. 

binded <- binded[,1:5]
fe.conditions = t(ab.fe.conditions)
fe_binded<-cbind(binded,fe.conditions)

#Pzzu_cor 
f_b_Pzzu_cor = fe_binded[,1:8]
f_b_Pzzu_cor$Sum <- f_b_Pzzu_cor$Pzzu_cor_2003+f_b_Pzzu_cor$Pzzu_cor_2011 + f_b_Pzzu_cor$Pzzu_cor_2018
f_b_Pzzu_cor <- f_b_Pzzu_cor %>% filter(Sum > 0)
subset_Pzzu_cor = f_b_Pzzu_cor[,1:5]

#Pzzu_par

f_b_Pzzu_par_par = fe_binded[,c(1,2,3,4,5,9,10,11)]
f_b_Pzzu_par_par$Sum <- f_b_Pzzu_par_par$Pzzu_par_2006+f_b_Pzzu_par_par$Pzzu_par_2011 + f_b_Pzzu_par_par$Pzzu_par_2018
f_b_Pzzu_par_par <- f_b_Pzzu_par_par %>% filter(Sum > 0)
subset_Pzzu_par = f_b_Pzzu_par_par[,1:5]

#Gabin_par
f_b_Gabin_par = fe_binded[,c(1,2,3,4,5,12,13,14)]
f_b_Gabin_par$Sum <- f_b_Gabin_par$Gabin_par_1999+f_b_Gabin_par$Gabin_par_2007 + f_b_Gabin_par$Gabin_par_2009
f_b_Gabin_par <- f_b_Gabin_par %>% filter(Sum > 0)
subset_Gabin_par = f_b_Gabin_par[,1:5]


#Pzzinu_par
f_b_Pzzinu_par = fe_binded[,c(1,2,3,4,5,15,16,17)]
f_b_Pzzinu_par$Sum <- f_b_Pzzinu_par$Pzzinu_par_2006+f_b_Pzzinu_par$Pzzinu_par_2011 + 
  f_b_Pzzinu_par$Pzzinu_par_2016
f_b_Pzzinu_par <- f_b_Pzzinu_par %>% filter(Sum > 0)
subset_Pzzinu_par = f_b_Pzzinu_par[,1:5]

#Passe_cor
f_b_Passe_cor = fe_binded[,c(1,2,3,4,5,18,19,20)]
f_b_Passe_cor$Sum <- f_b_Passe_cor$Passe_cor_2006 + f_b_Passe_cor$Passe_cor_2011+ f_b_Passe_cor$Passe_cor_2018
f_b_Passe_cor <- f_b_Passe_cor %>% filter(Sum > 0)
subset_Passe_cor = f_b_Passe_cor[,1:5]


#### Now we make the plots (Figure 4b-f)

###Passe_cor (Fig.4b)

m2 <- fd.coord[rownames(fd.coord),]
tr2 <-tri.mesh(m2[,1],m2[,2])
ch2 <- convex.hull(tr2)   
index <- c("1", "2","3", "4", "5", "6", "7", "8")
valuesC <- c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple")
subset_Passe_cor$color <-valuesC[match(subset_Passe_cor$groups,index)]

tiff(filename ="Figure_4b.tif", height=7, width=7, units="cm", compression = c("lzw"), res=300, pointsize=8)

plot(binded[,2], binded[,3], xlab="PCoA1", ylab="PCoA2", type="n", main="Passe_cor", col.main="black") 

points(subset_Passe_cor[,2], subset_Passe_cor[,3], pch=21, col=subset_Passe_cor$color , bg=subset_Passe_cor$color, cex=1)

#We add the polygons to the figure

m3 <- subset_Passe_cor[rownames(subset_Passe_cor),]
m3 <-m3[,1:5]
m3 <- m3[rownames(m3),]
m3 <- m3[m3$groups=="1",]
tr3 <-tri.mesh(m3[,2],m3[,3])
ch3 <- convex.hull(tr3)   

m4 <- subset_Passe_cor[rownames(subset_Passe_cor),]
m4 <-m4[,1:5]
m4 <- m4[rownames(m4),]
m4 <- m4[m4$groups=="2",]
tr4 <-tri.mesh(m4[,2],m4[,3])
ch4 <- convex.hull(tr4) 

m6 <- subset_Passe_cor[rownames(subset_Passe_cor),]
m6 <-m6[,1:5]
m6 <- m6[rownames(m6),]
m6 <- m6[m6$groups=="4",]
tr6 <-tri.mesh(m6[,2],m6[,3])
ch6 <- convex.hull(tr6)   

m9 <- subset_Passe_cor[rownames(subset_Passe_cor),]
m9 <-m9[,1:5]
m9 <- m9[rownames(m9),]
m9 <- m9[m9$groups=="7",]
tr9 <-tri.mesh(m9[,2],m9[,3])
ch9 <- convex.hull(tr9)   

polygon(ch2, col="#CCCCCC10", border=FALSE)
polygon(ch3, col="#20202060", border=FALSE)
polygon(ch4, col="#FF666680", border=FALSE)
polygon(ch6, col="#FFFF9980", border=FALSE)
polygon(ch9, col="#FF800080", border=FALSE)

dev.off()



##Pzzu_cor. Fig. 4c

 
m2 <- fd.coord[rownames(fd.coord),]
tr2 <-tri.mesh(m2[,1],m2[,2])
ch2 <- convex.hull(tr2)   
index <- c("1", "2","3", "4", "5", "6", "7", "8")
valuesC <- c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple")
subset_Pzzu_cor$color <-valuesC[match(subset_Pzzu_cor$groups,index)]

tiff(filename ="Figure_4c.tif", height=7, width=7, units="cm", compression = c("lzw"), res=300, pointsize=8)

plot(binded[,2], binded[,3], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzu_cor", col.main="black")
points(subset_Pzzu_cor[,2], subset_Pzzu_cor[,3], pch=21, col=subset_Pzzu_cor$color , bg=subset_Pzzu_cor$color, cex=1)

#We add the polygons to the figure

m3 <- subset_Pzzu_cor[rownames(subset_Pzzu_cor),]
m3 <-m3[,1:5]
m3 <- m3[rownames(m3),]
m3 <- m3[m3$groups=="1",]
tr3 <-tri.mesh(m3[,2],m3[,3])
ch3 <- convex.hull(tr3)   

m4 <- subset_Pzzu_cor[rownames(subset_Pzzu_cor),]
m4 <-m4[,1:5]
m4 <- m4[rownames(m4),]
m4 <- m4[m4$groups=="2",]
tr4 <-tri.mesh(m4[,2],m4[,3])
ch4 <- convex.hull(tr4) 

m6 <- subset_Pzzu_cor[rownames(subset_Pzzu_cor),]
m6 <-m6[,1:5]
m6 <- m6[rownames(m6),]
m6 <- m6[m6$groups=="4",]
tr6 <-tri.mesh(m6[,2],m6[,3])
ch6 <- convex.hull(tr6)   

m9 <- subset_Pzzu_cor[rownames(subset_Pzzu_cor),]
m9 <-m9[,1:5]
m9 <- m9[rownames(m9),]
m9 <- m9[m9$groups=="7",]
tr9 <-tri.mesh(m9[,2],m9[,3])
ch9 <- convex.hull(tr9)   

polygon(ch2, col="#CCCCCC10", border=FALSE)
polygon(ch3, col="#20202060", border=FALSE)
polygon(ch4, col="#FF666680", border=FALSE)
polygon(ch6, col="#FFFF9980", border=FALSE)
polygon(ch9, col="#FF800080", border=FALSE)

dev.off()


##Pzzinu_par   Fig. 4d

m2 <- fd.coord[rownames(fd.coord),]
tr2 <-tri.mesh(m2[,1],m2[,2])
ch2 <- convex.hull(tr2)   
index <- c("1", "2","3", "4", "5", "6", "7", "8")
valuesC <- c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple")
subset_Pzzinu_par$color <-valuesC[match(subset_Pzzinu_par$groups,index)]

tiff(filename ="Figure_4d.tif", height=7, width=7, units="cm", compression = c("lzw"), res=300, pointsize=8)

plot(binded[,2], binded[,3], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzinu_par", col.main="black") 
points(subset_Pzzinu_par[,2], subset_Pzzinu_par[,3], pch=21, col=subset_Pzzinu_par$color , bg=subset_Pzzinu_par$color, cex=1)

#We add the polygons to the figure

m3 <- subset_Pzzinu_par[rownames(subset_Pzzinu_par),]
m3 <-m3[,1:5]
m3 <- m3[rownames(m3),]
m3 <- m3[m3$groups=="1",]
tr3 <-tri.mesh(m3[,2],m3[,3])
ch3 <- convex.hull(tr3)   

m4 <- subset_Pzzinu_par[rownames(subset_Pzzinu_par),]
m4 <-m4[,1:5]
m4 <- m4[rownames(m4),]
m4 <- m4[m4$groups=="2",]
tr4 <-tri.mesh(m4[,2],m4[,3])
ch4 <- convex.hull(tr4) 

m5 <- subset_Pzzinu_par[rownames(subset_Pzzinu_par),]
m5 <-m5[,1:5]
m5 <- m5[rownames(m5),]
m5 <- m5[m5$groups=="3",]
tr5 <-tri.mesh(m5[,2],m5[,3])
ch5 <- convex.hull(tr5)   

m6 <- subset_Pzzinu_par[rownames(subset_Pzzinu_par),]
m6 <-m6[,1:5]
m6 <- m6[rownames(m6),]
m6 <- m6[m6$groups=="4",]
tr6 <-tri.mesh(m6[,2],m6[,3])
ch6 <- convex.hull(tr6)   

m7 <- subset_Pzzinu_par[rownames(subset_Pzzinu_par),]
m7 <-m7[,1:5]
m7 <- m7[rownames(m7),]
m7 <- m7[m7$groups=="5",]
tr7 <-tri.mesh(m7[,2],m7[,3])
ch7 <- convex.hull(tr7)   

m9 <- subset_Pzzinu_par[rownames(subset_Pzzinu_par),]
m9 <-m9[,1:5]
m9 <- m9[rownames(m9),]
m9 <- m9[m9$groups=="7",]
tr9 <-tri.mesh(m9[,2],m9[,3])
ch9 <- convex.hull(tr9)   

polygon(ch2, col="#CCCCCC10", border=FALSE)
polygon(ch3, col="#20202060", border=FALSE)
polygon(ch4, col="#FF666680", border=FALSE)
polygon(ch5, col="#B2FF6680", border=FALSE)
polygon(ch6, col="#FFFF9980", border=FALSE)
polygon(ch7, col="#66FFFF80", border=FALSE)

dev.off()



####Pzzu_par    Fig. 4e

m2 <- fd.coord[rownames(fd.coord),]
tr2 <-tri.mesh(m2[,1],m2[,2])
ch2 <- convex.hull(tr2)   
index <- c("1", "2","3", "4", "5", "6", "7", "8")
valuesC <- c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple")
subset_Pzzu_par$color <-valuesC[match(subset_Pzzu_par$groups,index)]

tiff(filename ="Figure_4e.tif", height=7, width=7, units="cm", compression = c("lzw"), res=300, pointsize=8)

plot(binded[,2], binded[,3], xlab="PCoA1", ylab="PCoA2", type="n", main="Pzzu_par", col.main="black")
oints(subset_Pzzu_par[,2], subset_Pzzu_par[,3], pch=21, col=subset_Pzzu_par$color , bg=subset_Pzzu_par$color, cex=1)

#We add the polygons to the figure

m3 <- subset_Pzzu_par[rownames(subset_Pzzu_par),]
m3 <-m3[,1:5]
m3 <- m3[rownames(m3),]
m3 <- m3[m3$groups=="1",]
tr3 <-tri.mesh(m3[,2],m3[,3])
ch3 <- convex.hull(tr3)   

m4 <- subset_Pzzu_par[rownames(subset_Pzzu_par),]
m4 <-m4[,1:5]
m4 <- m4[rownames(m4),]
m4 <- m4[m4$groups=="2",]
tr4 <-tri.mesh(m4[,2],m4[,3])
ch4 <- convex.hull(tr4) 

m5 <- subset_Pzzu_par[rownames(subset_Pzzu_par),]
m5 <-m5[,1:5]
m5 <- m5[rownames(m5),]
m5 <- m5[m5$groups=="3",]
tr5 <-tri.mesh(m5[,2],m5[,3])
ch5 <- convex.hull(tr5)   

m6 <- subset_Pzzu_par[rownames(subset_Pzzu_par),]
m6 <-m6[,1:5]
m6 <- m6[rownames(m6),]
m6 <- m6[m6$groups=="4",]
tr6 <-tri.mesh(m6[,2],m6[,3])
ch6 <- convex.hull(tr6)   

m7 <- subset_Pzzu_par[rownames(subset_Pzzu_par),]
m7 <-m7[,1:5]
m7 <- m7[rownames(m7),]
m7 <- m7[m7$groups=="5",]
tr7 <-tri.mesh(m7[,2],m7[,3])
ch7 <- convex.hull(tr7)   

m9 <- subset_Pzzu_par[rownames(subset_Pzzu_par),]
m9 <-m9[,1:5]
m9 <- m9[rownames(m9),]
m9 <- m9[m9$groups=="7",]
tr9 <-tri.mesh(m9[,2],m9[,3])
ch9 <- convex.hull(tr9)   

polygon(ch2, col="#CCCCCC10", border=FALSE)
polygon(ch3, col="#20202060", border=FALSE)
polygon(ch4, col="#FF666680", border=FALSE)
polygon(ch5, col="#B2FF6680", border=FALSE)
polygon(ch6, col="#FFFF9980", border=FALSE)
polygon(ch7, col="#66FFFF80", border=FALSE)
polygon(ch9, col="#FF800080", border=FALSE)

dev.off()



##Gabin_par    Fig. 4f

m2 <- fd.coord[rownames(fd.coord),]
tr2 <-tri.mesh(m2[,1],m2[,2])
ch2 <- convex.hull(tr2)   
index <- c("1", "2","3", "4", "5", "6", "7", "8")
valuesC <- c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple")
subset_Gabin_par$color <-valuesC[match(subset_Gabin_par$groups,index)]



tiff(filename ="Figure_4f.tif", height=7, width=7, units="cm", compression = c("lzw"), res=300, pointsize=8)

plot(binded[,2], binded[,3], xlab="PCoA1", ylab="PCoA2", type="n", main="Gabin_par", col.main="black") 
points(subset_Gabin_par[,2], subset_Gabin_par[,3], pch=21, col=subset_Gabin_par$color , bg=subset_Gabin_par$color, cex=1)

#We add the polygons to the figure

m3 <- subset_Gabin_par[rownames(subset_Gabin_par),]
m3 <-m3[,1:5]
m3 <- m3[rownames(m3),]
m3 <- m3[m3$groups=="1",]
tr3 <-tri.mesh(m3[,2],m3[,3])
ch3 <- convex.hull(tr3)   

m4 <- subset_Gabin_par[rownames(subset_Gabin_par),]
m4 <-m4[,1:5]
m4 <- m4[rownames(m4),]
m4 <- m4[m4$groups=="2",]
tr4 <-tri.mesh(m4[,2],m4[,3])
ch4 <- convex.hull(tr4) 

m5 <- subset_Gabin_par[rownames(subset_Gabin_par),]
m5 <-m5[,1:5]
m5 <- m5[rownames(m5),]
m5 <- m5[m5$groups=="3",]
tr5 <-tri.mesh(m5[,2],m5[,3])
ch5 <- convex.hull(tr5)   

m6 <- subset_Gabin_par[rownames(subset_Gabin_par),]
m6 <-m6[,1:5]
m6 <- m6[rownames(m6),]
m6 <- m6[m6$groups=="4",]
tr6 <-tri.mesh(m6[,2],m6[,3])
ch6 <- convex.hull(tr6)   

m7 <- subset_Gabin_par[rownames(subset_Gabin_par),]
m7 <-m7[,1:5]
m7 <- m7[rownames(m7),]
m7 <- m7[m7$groups=="5",]
tr7 <-tri.mesh(m7[,2],m7[,3])
ch7 <- convex.hull(tr7)   

polygon(ch2, col="#CCCCCC10", border=FALSE)
polygon(ch3, col="#20202060", border=FALSE)
polygon(ch4, col="#FF666680", border=FALSE)
polygon(ch5, col="#B2FF6680", border=FALSE)
polygon(ch6, col="#FFFF9980", border=FALSE)
polygon(ch7, col="#66FFFF80", border=FALSE)

dev.off()



######
#################3) PLOTTING THE NUMBER OF SPECIES PER FUNCTIONAL CLUSTER (Figure_4g-k)

###We load the file with the N of species per cluster previously calculated and obtained in STEP 1. 

Nsp_cluster <- read.csv2("Nsp_per_cluster.csv", sep=";", dec=",")

#####Passe_cor    Fig 4g

Nsp_cluster_Passe_cor <- Nsp_cluster[13:15,]
Nsp_cluster_Passe_cor$Time = c("1", "2", "3")
Nsp_cluster_Passe_cor$Time = as.factor(Nsp_cluster_Passe_cor$Time)

Nsp_cluster_plot_1 <- Nsp_cluster_Passe_cor %>% 
  tidyr::gather(key="Cluster",value="Nsp",-Time)

Nsp_cluster_plot_1 = na.omit(Nsp_cluster_plot_1)

plot_1 = Nsp_cluster_plot_1 %>% 
  ggplot(aes(Time,Nsp)) +
  geom_point(aes(color=Cluster,shape=Cluster), position=position_jitter(w=0.1, h=0, seed=NULL),size=3) + 
  scale_color_manual(values = c("black", "red", "green", "yellow", "orange", "purple"))+
  scale_shape_manual(values = c(16,16,16,16,16,16))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24),limits = c(0, 24))+
  ggtitle("Passe_cor") + ylab("N of sp within each cluster") +
  geom_hline(yintercept=2, linetype="dashed", 
             color = "red", size=0.5)+
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=15,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.y=element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(colour="black", size=12))+
  theme(legend.position="none")+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())

tiff(filename ="Figure_4g.tif", height=7, width=8, units="cm", compression = c("lzw"), res=300, pointsize=8)

plot_1
dev.off()


###Pzzu_cor    Fig. 4h

Nsp_cluster_Pzzu_cor <- Nsp_cluster[1:3,]
Nsp_cluster_Pzzu_cor$Time = c("1", "2", "3")
Nsp_cluster_Pzzu_cor$Time = as.factor(Nsp_cluster_Pzzu_cor$Time)

Nsp_cluster_plot_2 <- Nsp_cluster_Pzzu_cor %>% 
  tidyr::gather(key="Cluster",value="Nsp",-Time)

Nsp_cluster_plot_2 = na.omit(Nsp_cluster_plot_2)

plot_2 = Nsp_cluster_plot_2 %>% 
  ggplot(aes(Time,Nsp)) +
  geom_point(aes(color=Cluster,shape=Cluster), position=position_jitter(w=0.1, h=0, seed=NULL),size=3) + 
     scale_color_manual(values = c("black", "red", "green", "yellow", "orange", "purple"))+
  scale_shape_manual(values = c(16,16,16,16,16,16))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24),limits = c(0, 24))+
  ggtitle("Pzzu_cor") + ylab("N of sp within each cluster") +
  geom_hline(yintercept=2, linetype="dashed", 
             color = "red", size=0.5)+
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=15,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.y=element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(colour="black", size=12))+
  theme(legend.position="none")+
    theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())

tiff(filename ="Figure_4h.tif", height=7, width=8, units="cm", compression = c("lzw"), res=300, pointsize=8)

plot_2
dev.off()



###Pzzinu_par    Fig. 4i

Nsp_cluster_Pzzinu_par <- Nsp_cluster[10:12,]
Nsp_cluster_Pzzinu_par$Time = c("1", "2", "3")
Nsp_cluster_Pzzinu_par$Time = as.factor(Nsp_cluster_Pzzinu_par$Time)

Nsp_cluster_plot_3 <- Nsp_cluster_Pzzinu_par %>% 
  tidyr::gather(key="Cluster",value="Nsp",-Time)

Nsp_cluster_plot_3 = na.omit(Nsp_cluster_plot_3)

plot_3 = Nsp_cluster_plot_3 %>% 
  ggplot(aes(Time,Nsp)) +
  geom_point(aes(color=Cluster,shape=Cluster), position=position_jitter(w=0.1, h=0, seed=NULL),size=3) + 
  scale_color_manual(values = c("black", "red", "green", "yellow", "#66FFFF",  "orange", "purple"))+
  scale_shape_manual(values = c(16,16,16,16,16,16,16))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24),limits = c(0, 24))+
  ggtitle("Pzzinu_par") + ylab("N of sp within each cluster") +
  geom_hline(yintercept=2, linetype="dashed", 
             color = "red", size=0.5)+
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=15,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.y=element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(colour="black", size=12))+
  theme(legend.position="none")+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())

tiff(filename ="Figure_4i.tif", height=7, width=8, units="cm", compression = c("lzw"), res=300, pointsize=8)

plot_3
dev.off()


######Pzzu_par   Fig. 4j

Nsp_cluster_Pzzu_par <- Nsp_cluster[4:6,]
Nsp_cluster_Pzzu_par$Time = c("1", "2", "3")
Nsp_cluster_Pzzu_par$Time = as.factor(Nsp_cluster_Pzzu_par$Time)

Nsp_cluster_plot_4 <- Nsp_cluster_Pzzu_par %>% 
  tidyr::gather(key="Cluster",value="Nsp",-Time)

Nsp_cluster_plot_4 = na.omit(Nsp_cluster_plot_4)

plot_4 = Nsp_cluster_plot_4 %>% 
  ggplot(aes(Time,Nsp)) +
  geom_point(aes(color=Cluster,shape=Cluster), position=position_jitter(w=0.1, h=0, seed=NULL),size=3) + 
  scale_color_manual(values = c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple"))+
  scale_shape_manual(values = c(16,16,16,16,16,16,16,16))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24),limits = c(0, 24))+
  ggtitle("Pzzu_par") + ylab("N of sp within each cluster") +
  geom_hline(yintercept=2, linetype="dashed", 
             color = "red", size=0.5)+
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=15,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.y=element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(colour="black", size=12))+
  theme(legend.position="none")+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())

tiff(filename ="Figure_4j.tif", height=7, width=8, units="cm", compression = c("lzw"), res=300, pointsize=8)

plot_4
dev.off()


######Gabin_par   Fig. 4k

Nsp_cluster_Gabin_par <- Nsp_cluster[7:9,]
Nsp_cluster_Gabin_par$Time = c("1", "2", "3")
Nsp_cluster_Gabin_par$Time = as.factor(Nsp_cluster_Gabin_par$Time)

Nsp_cluster_plot_5 <- Nsp_cluster_Gabin_par %>% 
  tidyr::gather(key="Cluster",value="Nsp",-Time)

Nsp_cluster_plot_5 = na.omit(Nsp_cluster_plot_5)

plot_5 = Nsp_cluster_plot_5 %>% 
  ggplot(aes(Time,Nsp)) +
  geom_point(aes(color=Cluster,shape=Cluster), position=position_jitter(w=0.1, h=0, seed=NULL),size=3) + 
  scale_color_manual(values = c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "purple"))+
  scale_shape_manual(values = c(16,16,16,16,16,16,16))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24),limits = c(0, 24))+
  ggtitle("Gabin_par") + ylab("N of sp within each cluster") +
  geom_hline(yintercept=2, linetype="dashed", 
             color = "red", size=0.5)+
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=15,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.y=element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(colour="black", size=12))+
  theme(legend.position="none")+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())

tiff(filename ="Figure_4k.tif", height=7, width=8, units="cm", compression = c("lzw"), res=300, pointsize=8)

plot_5
dev.off()



############### 4) PLOTTING CLUSTER ABUNDANCE CHANGES PER SITE. Figure_4l-p

###We calculate the abundance of each cluster in each site and Time point. 

clusters <- binded$groups

ab_clus<- rbind(ab.fe.conditions,clusters) 

colnames(ab_clus) = ab_clus[16,]

ab_clus = ab_clus[1:15,]
ab_clus = as.data.frame((ab_clus))

ab_clusters<-t(rowsum(t(ab_clus), group = colnames(ab_clus), na.rm = T))

colnames(ab_clusters) <- paste("Cluster", colnames(ab_clusters), sep = "_")

ab_clusters = as.data.frame((ab_clusters))

Time = c("T1","T2","T3","T1","T2","T3","T1","T2","T3","T1","T2","T3","T1","T2","T3")

ab_clusters<-cbind(ab_clusters,Time)

ab_clusters$Time = as.factor(ab_clusters$Time)

###We make the abundance plots. Figure_4l-o

###Passe_cor Fig. 4l

ab_clus_Passe_cor <-ab_clusters[13:15,]

Passe_cor_ab_cluster <- ab_clus_Passe_cor %>% 
  tidyr::gather(key="Cluster",value="ab",-Time) 

cc13 <-c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple")

plot1 <-ggplot(data=Passe_cor_ab_cluster, aes(fill=Cluster, y=ab, x=Time ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  ggtitle("Passe_cor") + ylab("Abundance (%)") +
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=15,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.y=element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(colour="black", size=12))+
  theme(legend.position="none")+
  scale_fill_manual(values = cc13)+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())


tiff(filename ="Figure_4l.tif", height=7, width=7, units="cm", compression = c("lzw"), res=300, pointsize=8)
plot1
dev.off()


####Pzzu_cor   Fig. 4m

ab_clus_Pzzu_cor <-ab_clusters[1:3,]

Pzzu_cor_ab_cluster <- ab_clus_Pzzu_cor %>% 
  tidyr::gather(key="Cluster",value="ab",-Time) 

cc13 <-c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple")

plot2 <-ggplot(data=Pzzu_cor_ab_cluster, aes(fill=Cluster, y=ab, x=Time ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  ggtitle("Pzzu_cor") + ylab("Abundance (%)") +
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=15,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.y=element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(colour="black", size=12))+
  theme(legend.position="none")+
  scale_fill_manual(values = cc13)+
   theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())


tiff(filename ="Figure_4m.tif", height=7, width=7, units="cm", compression = c("lzw"), res=300, pointsize=8)
plot2
dev.off()


####Pzzinu_par   Fig. 4n

ab_clus_Pzzinu_par <-ab_clusters[10:12,]

Pzzinu_par_ab_cluster <- ab_clus_Pzzinu_par %>% 
  tidyr::gather(key="Cluster",value="ab",-Time) 

cc13 <-c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple")

plot3 <-ggplot(data=Pzzinu_par_ab_cluster, aes(fill=Cluster, y=ab, x=Time ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  ggtitle("Pzzinu_par") + ylab("Abundance (%)") +
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=15,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.y=element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(colour="black", size=12))+
  theme(legend.position="none")+
  scale_fill_manual(values = cc13)+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())


tiff(filename ="Figure_4n.tif", height=7, width=7, units="cm", compression = c("lzw"), res=300, pointsize=8)
plot3
dev.off()



####Pzzu_par  Fig. 4o

ab_clus_Pzzu_par <-ab_clusters[4:6,]

Pzzu_par_ab_cluster <- ab_clus_Pzzu_par %>% 
  tidyr::gather(key="Cluster",value="ab",-Time) 

cc13 <-c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple")

plot4 <-ggplot(data=Pzzu_par_ab_cluster, aes(fill=Cluster, y=ab, x=Time ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  ggtitle("Pzzu_par") + ylab("Abundance (%)") +
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=15,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.y=element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(colour="black", size=12))+
  theme(legend.position="none")+
  scale_fill_manual(values = cc13)+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())


tiff(filename ="Figure_4o.tif", height=7, width=7, units="cm", compression = c("lzw"), res=300, pointsize=8)
plot4
dev.off()

####Gabin_par    Fig. 4p

ab_clus_Gabin_par <-ab_clusters[7:9,]

Gabin_par_ab_cluster <- ab_clus_Gabin_par %>% 
  tidyr::gather(key="Cluster",value="ab",-Time) 

cc13 <-c("black", "red", "green", "yellow", "#66FFFF", "cadetblue", "orange", "purple")

plot5 <-ggplot(data=Gabin_par_ab_cluster, aes(fill=Cluster, y=ab, x=Time ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  ggtitle("Gabin_par") + ylab("Abundance (%)") +
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=15,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.y=element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(colour="black", size=12))+
  theme(legend.position="none")+
  scale_fill_manual(values = cc13)+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())


tiff(filename ="Figure_4p.tif", height=7, width=7, units="cm", compression = c("lzw"), res=300, pointsize=8)
plot5
dev.off()




