####################################################################
#### Code accompanying:
#
# Gomez-Gras et al. 2021. Climate change transforms the functional identity in Mediterranean coralligenous assemblages. Ecology Letters

# Code to compute: 
  #1) The functional trait space 
  #2) Trends of functional richness (Frich): Volume filled by each assemblage in the trait space (Fig. 2a-o)
  #3) Null models of functional richness (Fig. 2p-t) 
  #4) Ordination direction of trait vectors and position of categorized traits values in trait space (Fig. 3h-l) 
#
# Script written by: Daniel Gomez-Gras & Viviana Brambilla 
# Adapted from: Teixido et al.(2018). Functional biodiversity loss along natural CO2 gradients. Nature Communications
#
# Functions required and written by Sebastien Villeguer: 
# quality_funct_space.R: This function was obtained in the Supplementary data of Teixido et al. 2018
#
# 
#
####################################################################

# set working directory

setwd("~/xxx")


#### Load required packages
# install.packages("FD") # FD is a package for functional diversity measures
library('FD')
# install.packages("tripack") # Triangulation of irregularly spaced data
library('tripack')
# install.packages("geometry") # Geometric operations on point sets
library('geometry')
# install.packages("matrixStats") # Functions that operate on rows and columns of matrices
library('matrixStats')
# install.packages("ggplot2")
library('ggplot2')


####Load DATA

# Load FEs data
fes <- read.csv2("FE_ordered.csv", sep=";", dec=",", row.names=1)

#Load Species and FEs data
spe_fes <- read.csv2("Species_FE.csv", sep=";", dec=",", row.names=1)

#Load Abundance data
ab <- read.csv2("Abundances.csv", sep=";", dec=",", row.names=1)

#Load sites and quadrats 
sites <- read.table("Sites.txt", sep="\t", header=T, row.names=1)

# Defining sites and Time points. 
condition <- c("Pzzu_cor_2003", "Pzzu_cor_2011" , "Pzzu_cor_2018", "Pzzu_par_2006",
               "Pzzu_par_2011", "Pzzu_par_2018", "Gabin_par_1999", "Gabin_par_2007", 
               "Gabin_par_2009", "Pzzinu_par_2006","Pzzinu_par_2011","Pzzinu_par_2016", 
               "Passe_cor_2006","Passe_cor_2011","Passe_cor_2018")


############################### 1) CREATE THE FUNCTIONAL SPACE

# Creating  multidimensional functional spaces (2 to 11 D) for the 12 traits 
#load additional functions

source("quality_funct_space.R")

qfs <- quality_funct_space(fes, traits_weights=NULL, nbdim=11, metric= "Gower", dendro=FALSE, plot="quality_funct_space") 

# quality of spaces (low meanSD = high quality)
round( qfs$meanSD , 4)

# keeping coordinates on the 4 dimensions, meanSD<0.004
fd.coord <- qfs$details_funct_space$mat_coord[,1:4]

write.csv(fd.coord, file="FE_4D_coord.csv") #to use it for further analyses

#see variance explained by the PCoA axes
gower<-qfs$details_funct_space$mat_dissim

fit <- cmdscale(gower,eig=TRUE, k=4) # PCoA

# variance explained by the axes
cumsum(fit$eig[fit$eig>=0]) / sum(fit$eig[fit$eig>0])


#################################  2) TRENDS IN FUNCTIONAL RICHNESS

# Data manipulation and arrangements 

ab.conditions <- lapply(condition, function(x) {
  
            quad <- rownames(sites[sites$Year == x,])
            
            colSums(ab[rownames(ab) %in% quad,])
  
})#eo lapply

ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition
  
#### Calculate convex hull (111 taxonomic units in total for calculating the relative richness)

Fric <- lapply(condition, function (x) {
  
        species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
        fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
        
        m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
        
        ch <- convhulln(m, options = "FA")
        
        chg <- convhulln(fd.coord, options = "FA")
        
        c(length(species), length(species)/111*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol)
        
})#eo lapply

names(Fric) = condition

# Fric contains the number of species(NbSp) and FEs (NbFEs), relative percentages (NbSpP,NbFEsP ), and the 4D convex hull volume (Frich) among the 3 temporal points for each habitat
Fric <- do.call(rbind, Fric)

colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Frich (Vol4D)")


######### plot convex hull

cols <- c("#CD2626", "#F9D71C", "#3A5FCD","#CD2626", "#F9D71C", "#3A5FCD","#CD2626", "#F9D71C", "#3A5FCD","#CD2626","#F9D71C","#3A5FCD","#CD2626", "#F9D71C", "#3A5FCD")
colstr <- c("#CD262670", "#FFFF0070", "#3A5FCD70", "#CD262670", "#FFFF0070", "#3A5FCD70", "#CD262670", "#FFFF0070", "#3A5FCD70",
            "#CD262670","#FFFF0070","#3A5FCD70","#CD262670", "#FFFF0070", "#3A5FCD70")

names(cols) <- c("Pzzu_cor_2003", "Pzzu_cor_2011" , "Pzzu_cor_2018", "Pzzu_par_2006",
                 "Pzzu_par_2011", "Pzzu_par_2018", "Gabin_par_1999", "Gabin_par_2007", 
                 "Gabin_par_2009", "Pzzinu_par_2006","Pzzinu_par_2011","Pzzinu_par_2016", 
                 "Passe_cor_2006","Passe_cor_2011","Passe_cor_2018")
names(colstr) <- c("Pzzu_cor_2003", "Pzzu_cor_2011" , "Pzzu_cor_2018", "Pzzu_par_2006",
                   "Pzzu_par_2011", "Pzzu_par_2018", "Gabin_par_1999", "Gabin_par_2007", 
                   "Gabin_par_2009", "Pzzinu_par_2006","Pzzinu_par_2011","Pzzinu_par_2016", 
                   "Passe_cor_2006","Passe_cor_2011","Passe_cor_2018")


###PLOTING FUNCTIONAL SPACES. Fig. 2a-o. 
dev.new()
n_axes = 4
 
tiff(filename="Figure2a_o.tif", height=17, width=11, units="cm", compression = c("lzw"), res=500, pointsize=8)

par(mfrow = c(5,3))

Fric <- Fric[,c(2,4,5)]



for (i in condition) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[i,] > 0)]
  
  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
  
  #computing the convex hull for each condition
  m <- fd.coord[rownames(fd.coord) %in% fes_cond, ]
  
  tr <-tri.mesh(m[,1],m[,2])
  ch <- convex.hull(tr)
  
  ##computing the global convex hull to be represented as background
  m2 <- fd.coord[rownames(fd.coord),]
  tr2 <-tri.mesh(m2[,1],m2[,2])
  ch2 <- convex.hull(tr2)   
  
  plot(fd.coord[,1], fd.coord[,2], xlab = "PCoA 1", ylab = "PCoA 2", type="n",main=paste0(i))
  
  polygon(ch2, col="#CCCCCC30", border=FALSE)
  polygon(ch, col=colstr[i], border=cols[i])
  points(m[,1:2], pch = 16, col=cols[i]) 
  
  
}#eo for convex

dev.off()

# The Site and Time specific values of Frich added to the plots (2a-o) in the final version of the figure can be found 
# in the Fric object (created above) as (Frich (Vol4D)). 



############################### 3) NULL MODELS OF FUNCTIONAL RICHNESS (Fig. 2p-t)

### We create separate models for each site to test the effect of the Time in each of them.

#################################### Passe_cor 2006-2018 (Fig. 2p)
condition <- c("Passe_cor_2006", "Passe_cor_2011" , "Passe_cor_2018") 
ab_Passe_cor = ab[289:360,]
# an other way is to filter if it contain "Passe_cor" in the rownames

# Data arrangements to c

ab.conditions <- lapply(condition, function(x) {
  
  quad <- rownames(sites[sites$Year == x,])
  
  colSums(ab_Passe_cor[rownames(ab_Passe_cor) %in% quad,])
  
})#eo lapply

ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition

#### Calculate convex hull (101 species in total for calculating the relative richness)

Fric <- lapply(condition, function (x) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
  
  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
  
  ch <- convhulln(m, options = "FA")
  
  chg <- convhulln(fd.coord, options = "FA")
  
  c(length(species), length(species)/111*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
  
})#eo lapply

names(Fric) = condition
Fric <- do.call(rbind, Fric)

colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol4D")

Fric <- Fric[,c(2,4,5)]

# null model 

n_perm = 100

spe_fes_r = spe_fes

Fric_perm <- lapply(condition, function (x) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
  
  perm <- sapply((1:n_perm), function (z) {
    
    
    spe_fes_r$FE <- sample(spe_fes$FE)      
    
    fes_cond <- spe_fes_r[rownames(spe_fes_r) %in% species, ]
    
    m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
    
    ch <- convhulln(m, options = "FA")
    
    chg <- convhulln(fd.coord, options = "FA")
    
    c(length(species), length(species)/111*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
    
    
    
  })#eo sapply
  
  rownames(perm) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol")
  
  
  perm 
  
})#eo lapply

names(Fric_perm) = condition



Fric_perm_Q <- lapply(Fric_perm, function (x) {
  
  rowQuantiles(x, probs=c(0.05, 0.95))
  
})#eo lapply

Fric = as.data.frame(Fric)

Fric$lowerFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,1] })#eo sapply
Fric$upperFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,2] })#eo sapply
Fric$lowerVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,1] })#eo sapply
Fric$upperVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,2] })#eo sapply
Fric$cond <- condition
condition <- factor(condition, levels = c("Passe_cor_2006", "Passe_cor_2011" , "Passe_cor_2018"))

Fric$cond <- as.factor(condition)
levels(Fric$cond)
colnames(Fric) <- c("NbSp", "NbFE", "Vol4D", "lowerFE", "upperFE", "lowerVol", "upperVol", "cond")


#Plot the null model (Passe_cor). Figure 2p

Time1 <- rbind(Frich[1,], Frich[1,])
Time1$Vol4D <- c(Time1$lowerVol[1], Time1$upperVol[1]) 
Time2 <- rbind(Fric[2,], Fric[2,])
Time2$Vol4D <- c(Time2$lowerVol[1], Time2$upperVol[1]) 
Time3<- rbind(Fric[3,], Fric[3,])
Time3$Vol4D <- c(Time3$lowerVol[1], Time3$upperVol[1]) 

p1 <- ggplot(Fric,  aes(x = cond, y=Vol4D)) +
  scale_y_continuous(breaks=c(20,40,60,80,100), limits = c(0,100))+
  ylab("Relative richness (%)") +
  xlab("Time")+
  geom_point(aes(x = cond, y=Vol4D), data=Fric, pch=16, col=cols[c(1,2,3)], cex=4)+
  geom_line(data=Time1, aes(x = cond, y=Vol4D),lwd=2, col=cols["Passe_cor_2006"]) +
  geom_line(data=Time2, aes(x = cond, y=Vol4D),lwd=2, col=cols["Passe_cor_2011"]) +
  geom_line(data=Time3, aes(x = cond, y=Vol4D),lwd=2, col=cols["Passe_cor_2018"])+
  ggtitle("Passe_cor")+
  theme(plot.title = element_text(hjust = 0.5,size=17))+
  coord_flip()+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=18,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.title.x = element_text(size=14))+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.text.y=element_text(colour="black", size=12), strip.text.x=element_text(size=18,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  scale_x_discrete(labels=c("T3","T2","T1"),limits = rev(levels(Fric$cond)))+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border= element_rect(colour = "black", fill=NA, size=0.8),
        panel.background=element_blank())


tiff(filename="Figure_2p.tif", height=10, width=10, units="cm", compression = c("lzw"), res=300, pointsize=8)
p1
dev.off()


####################### Pzzu_cor_2003_2018. Fig. 2q

condition <- c("Pzzu_cor_2003", "Pzzu_cor_2011" , "Pzzu_cor_2018") 
ab_Pzzu_cor = ab[1:72,]

# Data arrangements 
ab.conditions <- lapply(condition, function(x) {
  
  quad <- rownames(sites[sites$Year == x,])
  
  colSums(ab_Pzzu_cor[rownames(ab_Pzzu_cor) %in% quad,])
  
})#eo lapply

ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition

Fric <- lapply(condition, function (x) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
  
  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
  
  ch <- convhulln(m, options = "FA")
  
  chg <- convhulln(fd.coord, options = "FA")
  
  c(length(species), length(species)/111*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
  
})#eo lapply

names(Fric) = condition
Fric <- do.call(rbind, Fric)
colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol4D")
Fric <- Fric[,c(2,4,5)]

############### null model  

n_perm = 100

spe_fes_r = spe_fes

Fric_perm <- lapply(condition, function (x) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
  
  perm <- sapply((1:n_perm), function (z) {
    
    
    spe_fes_r$FE <- sample(spe_fes$FE)      
    
    fes_cond <- spe_fes_r[rownames(spe_fes_r) %in% species, ]
    
    m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
    
    ch <- convhulln(m, options = "FA")
    
    chg <- convhulln(fd.coord, options = "FA")
    
    c(length(species), length(species)/111*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
    
    
    
  })#eo sapply
  
  rownames(perm) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol")
  
  
  perm 
  
})#eo lapply

names(Fric_perm) = condition

Fric_perm_Q <- lapply(Fric_perm, function (x) {
  
  rowQuantiles(x, probs=c(0.05, 0.95))
  
})#eo lapply

Fric = as.data.frame(Fric)

Fric$lowerFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,1] })#eo sapply
Fric$upperFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,2] })#eo sapply
Fric$lowerVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,1] })#eo sapply
Fric$upperVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,2] })#eo sapply
Fric$cond <- condition
condition <- factor(condition, levels = c("Pzzu_cor_2003", "Pzzu_cor_2011" , "Pzzu_cor_2018"))

Fric$cond <- as.factor(condition)
levels(Fric$cond)
colnames(Fric) <- c("NbSp", "NbFE", "Vol4D", "lowerFE", "upperFE", "lowerVol", "upperVol", "cond")


#Plot the null model 
#Figure 2q. (Pzzu_cor)

Time1 <- rbind(Fric[1,], Fric[1,])
Time1$Vol4D <- c(Time1$lowerVol[1], Time1$upperVol[1]) 
Time2 <- rbind(Fric[2,], Fric[2,])
Time2$Vol4D <- c(Time2$lowerVol[1], Time2$upperVol[1]) 
Time3<- rbind(Fric[3,], Fric[3,])
Time3$Vol4D <- c(Time3$lowerVol[1], Time3$upperVol[1]) 

p2 <- ggplot(Fric,  aes(x = cond, y=Vol4D)) +
  scale_y_continuous(breaks=c(20,40,60,80,100), limits = c(0,100))+
  ylab("Relative richness (%)") +
  xlab("Time")+
  geom_point(aes(x = cond, y=Vol4D), data=Fric, pch=16, col=cols[c(1,2,3)], cex=4)+
  geom_line(data=Time1, aes(x = cond, y=Vol4D),lwd=2, col=cols["Pzzu_cor_2003"]) +
  geom_line(data=Time2, aes(x = cond, y=Vol4D),lwd=2, col=cols["Pzzu_cor_2011"]) +
  geom_line(data=Time3, aes(x = cond, y=Vol4D),lwd=2, col=cols["Pzzu_cor_2018"])+
  ggtitle("Pzzu_cor")+
  theme(plot.title = element_text(hjust = 0.5,size=17))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=18,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.title.x = element_text(size=14))+
  coord_flip()+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.text.y=element_text(colour="black", size=12), strip.text.x=element_text(size=18,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  scale_x_discrete(labels=c("T3","T2","T1"),limits = rev(levels(Fric$cond)))+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border= element_rect(colour = "black", fill=NA, size=0.8),
        panel.background=element_blank())




tiff(filename="Figure_2q.tif", height=10, width=10, units="cm", compression = c("lzw"), res=300, pointsize=8)
p2
dev.off()


####################### Pzzinu_par 2006-2016. Fig. 2r 


condition <- c("Pzzinu_par_2006","Pzzinu_par_2011","Pzzinu_par_2016")

ab_Pzzinu_par = ab[217:288,]

# Data arrangements 

ab.conditions <- lapply(condition, function(x) {
  
  quad <- rownames(sites[sites$Year == x,])
  
  colSums(ab_Pzzinu_par[rownames(ab_Pzzinu_par) %in% quad,])
  
})#eo lapply

ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition

#### Calculate convex hull (101 species in total for calculating the relative richness)

Fric <- lapply(condition, function (x) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
  
  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
  
  ch <- convhulln(m, options = "FA")
  
  chg <- convhulln(fd.coord, options = "FA")
  
  c(length(species), length(species)/111*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
  
})#eo lapply

names(Fric) = condition
Fric <- do.call(rbind, Fric)

colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol4D")

Fric <- Fric[,c(2,4,5)]

# null model of Functional richness 

n_perm = 100

spe_fes_r = spe_fes

Fric_perm <- lapply(condition, function (x) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
  
  perm <- sapply((1:n_perm), function (z) {
    
    
    spe_fes_r$FE <- sample(spe_fes$FE)      
    
    fes_cond <- spe_fes_r[rownames(spe_fes_r) %in% species, ]
    
    m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
    
    ch <- convhulln(m, options = "FA")
    
    chg <- convhulln(fd.coord, options = "FA")
    
    c(length(species), length(species)/111*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
    
    
    
  })#eo sapply
  
  rownames(perm) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol")
  
  
  perm 
  
})#eo lapply

names(Fric_perm) = condition

Fric_perm_Q <- lapply(Fric_perm, function (x) {
  
  rowQuantiles(x, probs=c(0.05, 0.95))
  
})#eo lapply


Fric = as.data.frame(Fric)

Fric$lowerFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,1] })#eo sapply
Fric$upperFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,2] })#eo sapply
Fric$lowerVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,1] })#eo sapply
Fric$upperVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,2] })#eo sapply
Fric$cond <- condition
condition <- factor(condition, levels = c("Pzzinu_par_2006","Pzzinu_par_2011","Pzzinu_par_2016"))

Fric$cond <- as.factor(condition)
levels(Fric$cond)
colnames(Fric) <- c("NbSp", "NbFE", "Vol4D", "lowerFE", "upperFE", "lowerVol", "upperVol", "cond")


#Plot the null model. (Pzzinu_par). Figure_2r

Time1 <- rbind(Fric[1,], Fric[1,])
Time1$Vol4D <- c(Time1$lowerVol[1], Time1$upperVol[1]) 
Time2 <- rbind(Fric[2,], Fric[2,])
Time2$Vol4D <- c(Time2$lowerVol[1], Time2$upperVol[1]) 
Time3<- rbind(Fric[3,], Fric[3,])
Time3$Vol4D <- c(Time3$lowerVol[1], Time3$upperVol[1]) 

p3 <- ggplot(Fric,  aes(x = cond, y=Vol4D)) +
  scale_y_continuous(breaks=c(20,40,60,80,100), limits = c(0,100))+
  ylab("Relative richness (%)") +
  xlab("Time")+
  geom_point(aes(x = cond, y=Vol4D), data=Fric, pch=16, col=cols[c(1,2,3)], cex=4)+
  geom_line(data=Time1, aes(x = cond, y=Vol4D),lwd=2, col=cols["Pzzinu_par_2006"]) +
  geom_line(data=Time2, aes(x = cond, y=Vol4D),lwd=2, col=cols["Pzzinu_par_2011"]) +
  geom_line(data=Time3, aes(x = cond, y=Vol4D),lwd=2, col=cols["Pzzinu_par_2016"])+
  ggtitle("Pzzinu_par")+
  theme(plot.title = element_text(hjust = 0.5,size=17))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=18,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.title.x = element_text(size=14))+
  coord_flip()+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.text.y=element_text(colour="black", size=12), strip.text.x=element_text(size=18,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  scale_x_discrete(labels=c("T3","T2","T1"),limits = rev(levels(Fric$cond)))+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border= element_rect(colour = "black", fill=NA, size=0.8),
        panel.background=element_blank())


tiff(filename="Figure_2r.tif", height=10, width=10, units="cm", compression = c("lzw"), res=300, pointsize=8)
p3
dev.off()


######################## Pzzu_par 2006-2018. Fig 2s 

condition <- c("Pzzu_par_2006", "Pzzu_par_2011" , "Pzzu_par_2018")
ab_Pzzu_par = ab[73:144,]

# Data arrangements 

ab.conditions <- lapply(condition, function(x) {
  
  quad <- rownames(sites[sites$Year == x,])
  
  colSums(ab_Pzzu_par[rownames(ab_Pzzu_par) %in% quad,])
  
})#eo lapply

ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition

#### Calculate convex hull (101 species in total for calculating the relative richness)

Fric <- lapply(condition, function (x) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
  
  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
  
  ch <- convhulln(m, options = "FA")
  
  chg <- convhulln(fd.coord, options = "FA")
  
  c(length(species), length(species)/111*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
  
})#eo lapply

names(Fric) = condition

# Fric contains the number of species(NbSp) and FEs (NbFEs), relative percentages (NbSpP,NbFEsP ) , and the volume among the 3 temporal points for each habitat
Fric <- do.call(rbind, Fric)

colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol4D")

Fric <- Fric[,c(2,4,5)]

############### null model 

n_perm = 100

spe_fes_r = spe_fes

Fric_perm <- lapply(condition, function (x) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
  
  perm <- sapply((1:n_perm), function (z) {
    
    
    spe_fes_r$FE <- sample(spe_fes$FE)      
    
    fes_cond <- spe_fes_r[rownames(spe_fes_r) %in% species, ]
    
    m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
    
    ch <- convhulln(m, options = "FA")
    
    chg <- convhulln(fd.coord, options = "FA")
    
    c(length(species), length(species)/111*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
    
    
    
  })#eo sapply
  
  rownames(perm) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol")
  
  
  perm 
  
})#eo lapply

names(Fric_perm) = condition



Fric_perm_Q <- lapply(Fric_perm, function (x) {
  
  rowQuantiles(x, probs=c(0.05, 0.95))
  
})#eo lapply

Fric = as.data.frame(Fric)

Fric$lowerFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,1] })#eo sapply
Fric$upperFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,2] })#eo sapply
Fric$lowerVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,1] })#eo sapply
Fric$upperVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,2] })#eo sapply
Fric$cond <- condition
condition <- factor(condition, levels = c("Pzzu_par_2006", "Pzzu_par_2011", "Pzzu_par_2018"))

Fric$cond <- as.factor(condition)
levels(Fric$cond)
colnames(Fric) <- c("NbSp", "NbFE", "Vol4D", "lowerFE", "upperFE", "lowerVol", "upperVol", "cond")


#Plot the null model (Pzzu_par). Figure 2s

Time1 <- rbind(Fric[1,], Fric[1,])
Time1$Vol4D <- c(Time1$lowerVol[1], Time1$upperVol[1]) 
Time2 <- rbind(Fric[2,], Fric[2,])
Time2$Vol4D <- c(Time2$lowerVol[1], Time2$upperVol[1]) 
Time3<- rbind(Fric[3,], Fric[3,])
Time3$Vol4D <- c(Time3$lowerVol[1], Time3$upperVol[1]) 

p4 <- ggplot(Fric,  aes(x = cond, y=Vol4D)) +
  scale_y_continuous(breaks=c(20,40,60,80,100), limits = c(0,100))+
  ylab("Relative richness (%)") +
  xlab("Time")+
  geom_point(aes(x = cond, y=Vol4D), data=Fric, pch=16, col=cols[c(1,2,3)], cex=4)+
  geom_line(data=Time1, aes(x = cond, y=Vol4D),lwd=2, col=cols["Pzzu_par_2006"]) +
  geom_line(data=Time2, aes(x = cond, y=Vol4D),lwd=2, col=cols["Pzzu_par_2011"]) +
  geom_line(data=Time3, aes(x = cond, y=Vol4D),lwd=2, col=cols["Pzzu_par_2018"])+
  ggtitle("Pzzu_Par")+
  theme(plot.title = element_text(hjust = 0.5,size=17))+
  coord_flip()+ 
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=18,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.title.x = element_text(size=14))+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.text.y=element_text(colour="black", size=12), strip.text.x=element_text(size=18,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  scale_x_discrete(labels=c("T3","T2","T1"),limits = rev(levels(Fric$cond)))+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border= element_rect(colour = "black", fill=NA, size=0.8),
        panel.background=element_blank())


tiff(filename="Figure_2s.tif", height=10, width=10, units="cm", compression = c("lzw"), res=300, pointsize=8)
p4
dev.off()


################################## Gabin_par 1999-2009. Fig 2t

condition <- c("Gabin_par_1999", "Gabin_par_2007" , "Gabin_par_2009")
ab_Gabin_par = ab[145:216,]

# Data arrangements 

ab.conditions <- lapply(condition, function(x) {
  
  quad <- rownames(sites[sites$Year == x,])
  
  colSums(ab_Gabin_par[rownames(ab_Gabin_par) %in% quad,])
  
})#eo lapply

ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition

#### Calculate convex hull (101 species in total for calculating the relative richness)

Fric <- lapply(condition, function (x) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
  
  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
  
  ch <- convhulln(m, options = "FA")
  
  chg <- convhulln(fd.coord, options = "FA")
  
  c(length(species), length(species)/111*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
  
})#eo lapply

names(Fric) = condition
Fric <- do.call(rbind, Fric)

colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol4D")

Fric <- Fric[,c(2,4,5)]

############### null model 

n_perm = 100

spe_fes_r = spe_fes

Fric_perm <- lapply(condition, function (x) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
  
  perm <- sapply((1:n_perm), function (z) {
    
    
    spe_fes_r$FE <- sample(spe_fes$FE)      
    
    fes_cond <- spe_fes_r[rownames(spe_fes_r) %in% species, ]
    
    m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
    
    ch <- convhulln(m, options = "FA")
    
    chg <- convhulln(fd.coord, options = "FA")
    
    c(length(species), length(species)/111*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
    
    
    
  })#eo sapply
  
  rownames(perm) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol")
  
  
  perm 
  
})#eo lapply

names(Fric_perm) = condition



Fric_perm_Q <- lapply(Fric_perm, function (x) {
  
  rowQuantiles(x, probs=c(0.05, 0.95))
  
})#eo lapply

Fric = as.data.frame(Fric)

Fric$lowerFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,1] })#eo sapply
Fric$upperFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,2] })#eo sapply
Fric$lowerVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,1] })#eo sapply
Fric$upperVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,2] })#eo sapply
Fric$cond <- condition
condition <- factor(condition, levels = c("Gabin_par_1999", "Gabin_par_2007" , "Gabin_par_2009"))

Fric$cond <- as.factor(condition)
levels(Fric$cond)
colnames(Fric) <- c("NbSp", "NbFE", "Vol4D", "lowerFE", "upperFE", "lowerVol", "upperVol", "cond")


#Plot the null model (Gabin_par) Fig. 2t


Time1 <- rbind(Fric[1,], Fric[1,])
Time1$Vol4D <- c(Time1$lowerVol[1], Time1$upperVol[1]) 
Time2 <- rbind(Fric[2,], Fric[2,])
Time2$Vol4D <- c(Time2$lowerVol[1], Time2$upperVol[1]) 
Time3<- rbind(Fric[3,], Fric[3,])
Time3$Vol4D <- c(Time3$lowerVol[1], Time3$upperVol[1]) 

p5 <- ggplot(Fric,  aes(x = cond, y=Vol4D)) +
  scale_y_continuous(breaks=c(20,40,60,80,100), limits = c(0,100))+
  ylab("Relative richness (%)") +
  xlab("Time")+
  geom_point(aes(x = cond, y=Vol4D), data=Fric, pch=16, col=cols[c(1,2,3)], cex=4)+
  geom_line(data=Time1, aes(x = cond, y=Vol4D),lwd=2, col=cols["Gabin_par_1999"]) +
  geom_line(data=Time2, aes(x = cond, y=Vol4D),lwd=2, col=cols["Gabin_par_2007"]) +
  geom_line(data=Time3, aes(x = cond, y=Vol4D),lwd=2, col=cols["Gabin_par_2009"])+
  ggtitle("Gabin_par")+
  coord_flip()+
  theme(plot.title = element_text(hjust = 0.5,size=17))+
  theme(axis.text.x=element_text(colour="black", size=12), strip.text.x=element_text(size=18,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  theme(axis.title.x = element_text(size=14))+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.text.y=element_text(colour="black", size=12), strip.text.x=element_text(size=18,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",size=1.1))+
  scale_x_discrete(labels=c("T3","T2","T1"),limits = rev(levels(Fric$cond)))+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border= element_rect(colour = "black", fill=NA, size=0.8),
        panel.background=element_blank())


tiff(filename="Figure_2t.tif", height=10, width=10, units="cm", compression = c("lzw"), res=300, pointsize=8)
p5
dev.off()

##################### 4) TRAIT VECTORS DIRECTIONS AND POSITION OF CATEGORICAL TRAITS IN FUNCTIONAL SPACE 

##We use the cdmscale() and envfit() to re-ecreate our trait-space and find vector directions and categories positions in functional space

fit <- cmdscale(gower,eig=TRUE, k=4) # PCoA


###Morphology

efit <- envfit(fit, fes[1], na.rm=TRUE) 

###Coloniality 

efit1 <- envfit(fit, fes[2], na.rm=TRUE) 


###Longevity 

efit2 <- envfit(fit, fes[3], na.rm=TRUE) 

###Height (cm)

efit3 <- envfit(fit, fes[4], na.rm=TRUE) 

###Width (cm)

efit4 <- envfit(fit, fes[5], na.rm=TRUE) 

###No-epibiosis 

efit5 <- envfit(fit, fes[6], na.rm=TRUE) 

###Heterotrophy

efit6 <- envfit(fit, fes[7], na.rm=TRUE) 

###Photosynthetic pigment

efit7 <- envfit(fit, fes[8], na.rm=TRUE) 

###Feeding strategy

efit8 <- envfit(fit, fes[9], na.rm=TRUE) 
###Age at reproduction

efit9 <- envfit(fit, fes[10], na.rm=TRUE) 

###Growth rate

efit11 <- envfit(fit, fes[11], na.rm=TRUE) 

###Physical defense

efit12 <- envfit(fit, fes[12], na.rm=TRUE) 


######### Plotting vectors and categories in trait space. Used to create Fig. 3h-l

###We first creat the background polygon (functional space)
m2 <- fd.coord[rownames(fd.coord),]
tr2 <-tri.mesh(m2[,1],m2[,2])
ch2 <- convex.hull(tr2)

n_axes = 4

pdf("Figure_3h_l.pdf", height=4.3, width= 6, useDingbats=FALSE)

par(mfrow=c(2,3))

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Semi-quantitative traits", col.main="#CD2626")
polygon(ch2, col="#CCCCCC30", border=FALSE)
plot(efit1,col =  "red", cex = 0.3)
plot(efit2, col = "red", cex = 0.3)
plot(efit3, col = "red", cex = 0.3)
plot(efit4, col = "red", cex = 0.3)
plot(efit5, col = "red", cex = 0.3)
plot(efit6, col = "red", cex = 0.3)
plot(efit9, col = "red", cex = 0.3)
plot(efit11, col = "red", cex = 0.3)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Morphology", col.main="#CD2626")
polygon(ch2, col="#CCCCCC30", border=FALSE)
plot(efit, col = "red", cex = .4)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Photosynthetic pigment", col.main="#CD2626")
polygon(ch2, col="#CCCCCC30", border=FALSE)
plot(efit7, col = "red", cex = .4)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Feeding strategy", col.main="#CD2626")
polygon(ch2, col="#CCCCCC30", border=FALSE)
plot(efit8, col = "red", cex = 0.4)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Physical defences", col.main="#CD2626")
polygon(ch2, col="#CCCCCC30", border=FALSE)
plot(efit12, col = "red", cex = .4)

dev.off()

