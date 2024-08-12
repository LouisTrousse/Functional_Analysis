## Code to compute: 
#1) Broad functional classification into functional clusters (with the overall pool).
#2) Broad functional classification into functional clusters (per site).
#3) Number of species per functional claster in each site and Time point.
#4) Temporal abundance changes in functional clusters per site.

#### 1. Broad Classification in Functional Clusters ####
# Convert the distanc list into a matrix
gower_mat = as.matrix(gower)

# First we retrieve the number of different FE, number of columns in the gower matrix
n_FE = dim(gower_mat)[1]

# How many clusters? We will use the "Silhouette Width" parameter to determine the optimal number of clusters
# We will iterate through 2-n_FE-1 clusters (because we need at least 2 clusters and n_FE clusters would not be a clustering) and see which shows the maximum separation between centres
# The "Silhoutte Width" parameter is used by default for PAM clusters

# We iterate through the number of clusters
sil_width <- lapply(2:(n_FE-1), function(i) {
  pam_fit = pam(gower, diss = TRUE, k=i)
  return(pam_fit$silinfo$avg.width)
})


# set the limit of the number of clusters to 20 for a broad classification
limit = 20
# We use a conditionnal loop to enshure that the limit set by the user is not exceeded the number of clusters possible 
if (n_FE > limit) {
  #unlist the sil_width list
  sil_width <- unlist(sil_width)
  # We associate the number of clusters with the silhouette width value in the same vector
  sil_width_v <- cbind(2:(n_FE-1), sil_width)
  # We discard values with a number of clusters higher than 20 because we are interested in a broad classification
  sil_width_filtered <- sil_width_v[sil_width_v[,1] < 20,]
  # We sort the values
  sil_width_ordered <- sil_width_filtered[order(sil_width_filtered[,2], decreasing = TRUE),]
  # We select the maximum value of silhouette width
  max_sil_width <- sil_width_ordered[1,]
  n_cluster <- max_sil_width[1]
  cat("The optimal number of clusters is",n_cluster,"\n with a silhouette width of",max_sil_width[2])
} else {
  cat("A limit of",limit,"culsters is higher than the number of possible clusters.\nPlease set a limit beetwen 2 and ",n_FE-1)
}

# Create the clusters using PAM for 8 clusters
pam_fit = pam(gower, diss=TRUE, k= n_cluster)

# Plotting
# We add the "cluster N" column to our data with coordinates. This will allow us to know to which cluster has been assigned to every FE. 
# We will use these data later to calculate the N of sp per cluster (functional redundancy) and if there are temporal changes in the relative abundance of functional clusters in the different assemblages. 

# We need fd.coord as before
fd.coord <- qfs$details_funct_space$mat_coord[,1:4]

pam_fit$clustering = as.data.frame(pam_fit$clustering)
binded <- cbind(pam_fit$clustering, fd.coord)
colnames(binded)[which(names(binded) == "pam_fit$clustering")] <- "Cluster"

#####plotting with colors 
#first we create the convex hulls of the clusters to plot them in the figure
# Retrieve the global convex hull (functionnal space)
m2 <- fd.coord[rownames(fd.coord),]
tr2 <-tri.mesh(m2[,1],m2[,2])
ch2 <- convex.hull(tr2) 
cluster <- "global"
ch_global <- list(global = c(ch2, Cluster = cluster))


# For each cluster we create the convex hull
Clusters = unique(binded$Cluster)

cluster_hull <- lapply(Clusters, function(cluster) {
  # retrieve coordinates of the cluster
  m <- binded %>% 
    filter(Cluster == cluster)
  # if there is only two or less point in the cluster, we can't create a convex hull
  # So we return coordinates of the points in the same form than the convex hull function
  if (nrow(m) <= 2) {
    ch = list(x = m[,"PC1"], y = m[,"PC2"], i = "NA", Cluster = cluster)
  } else {
    # create the triangulation
    tr <-tri.mesh(m[,"PC1"],m[,"PC2"])
    # create the convex hull
    ch <- convex.hull(tr)
    # add an entry on the c list to know to which cluster the convex hull belongs
    ch <- c(ch, Cluster = cluster)
  } #eo if 
}) #eo lapply

# We take the global convex hull and the convex hulls of the clusters and merge them in a data frame
# add global convex hull to the list with the name of the cluster as list entry (global) for the global convex hull
all_hulls <- c(cluster_hull, global = ch_global)

# convert to data frame each convex hull
all_hulls_df <- lapply(all_hulls, function(hull) {
  
  #retrieve hull index to know to which cluster it belong
  data.frame(x = hull$x, # x coordinates
             y = hull$y, # y coordinates
             i = hull$i,
             Cluster = as.character(hull$Cluster)) # index of this point
})

# merge all the convex hulls in a single data frame
all_hulls_df <- do.call(rbind, all_hulls_df)

# We specify the color of the clusters
Colors = c("1"= "black", "2" = "red", "3" = "green", "4" = "yellow", "5" = "#66FFFF", "6" = "cadetblue", "7" = "orange", "8" = "purple")
global_hull_color = "#CCCCCC30"

# If the number of clusters is egual to the number of cluster colors we use the colors specified above, otherwise we use a default palette (rainbow)
if (n_cluster == length(Colors)) {
  binded$color <- Colors[match(binded$Cluster, names(Colors))]
} else {
  # print a warning message
  cat("The number of clusters is higher than the number of colors specified. A default color palette has been used.\n")
  # We create a color palette based on the default rainbow and assign it to the clusters
  colors<-rainbow(n_cluster)
  binded$color <- colors[match(binded$Cluster, as.character(1:n_cluster))]
}

# we do the same for the global convex hull and add a color for the global convex hull
if (n_cluster == length(Colors)) {
  #add the color of each cluster
  all_hulls_df$color <- Colors[match(all_hulls_df$Cluster, names(Colors))]
  # add a color for the global convex hull
  all_hulls_df$color[all_hulls_df$Cluster == "global"] <- global_hull_color
} else {
  colors<-rainbow(n_cluster) 
  all_hulls_df$color <- colors[match(all_hulls_df$Cluster, as.character(1:n_cluster))]
  all_hulls_df$color <- colors[match(all_hulls_df$Cluster, as.character(1:n_cluster))]
}

# We plot using the ggplot2 package
ggplot(binded, aes(x = PC1, y = PC2)) +
  geom_polygon(data = all_hulls_df, aes(x = x, y = y, group = Cluster, fill = color), alpha = 0.7) +
  geom_point(aes(color = color, fill = color)) +
  labs(title = paste0(Title," Functional Clusters"), x = "PCoA1", y = "PCoA2")+
  scale_color_identity() + # color as specified in the dataset 
  scale_fill_identity()+ # color as specified in the dataset
  theme_classic()+ # classic theme
  theme(panel.border = element_rect(colour = "black", fill=NA))+ # border of the plot
  theme(axis.text = element_text(face = "bold", size = 12))+ # bold axis text
  theme(axis.title = element_text(face = "bold"))+ # bold axis title
  theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid"))+
  theme(legend.position = "none")+ # remove legend
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))# centered bold

# We can resume it in a function. 
create_convex_hull_plot <- function(fd.coord, binded, n_cluster, Colors = NULL, Title = "Global") {
  
  # Load necessary libraries
  library(ggplot2)
  library(tripack) # for tri.mesh and convex.hull
  library(dplyr)
  
  # Retrieve the global convex hull (functional space)
  m2 <- fd.coord[rownames(fd.coord),]
  tr2 <- tri.mesh(m2[,1], m2[,2])
  ch2 <- convex.hull(tr2)
  cluster <- "global"
  ch_global <- list(global = c(ch2, Cluster = cluster))
  
  # For each cluster, create the convex hull
  Clusters <- unique(binded$Cluster)
  
  cluster_hull <- lapply(Clusters, function(cluster) {
    # Retrieve coordinates of the cluster
    m <- binded %>% filter(Cluster == cluster)
    
    # If there are two or fewer points in the cluster, we can't create a convex hull
    if (nrow(m) <= 2) {
      ch <- list(x = m[,"PC1"], y = m[,"PC2"], i = "NA", Cluster = cluster)
    } else {
      # Create the triangulation
      tr <- tri.mesh(m[,"PC1"], m[,"PC2"])
      # Create the convex hull
      ch <- convex.hull(tr)
      # Add an entry on the list to know to which cluster the convex hull belongs
      ch <- c(ch, Cluster = cluster)
    }
    return(ch)
  })
  
  # Add global convex hull to the list
  all_hulls <- c(cluster_hull, global = ch_global)
  
  # Convert to data frame each convex hull
  all_hulls_df <- lapply(all_hulls, function(hull) {
    data.frame(
      x = hull$x, # x coordinates
      y = hull$y, # y coordinates
      i = hull$i, # index of this point
      Cluster = as.character(hull$Cluster)
    )
  })
  
  # Merge all the convex hulls in a single data frame
  all_hulls_df <- do.call(rbind, all_hulls_df)
  
  # Default color for global convex hull
  global_hull_color <- "#CCCCCC30"
  
  # Assign colors to clusters
  if (is.null(Colors)) {
    # Default colors based on the number of clusters using rainbow palette
    Colors <- setNames(rainbow(n_cluster), as.character(1:n_cluster))
  }
  
  if (n_cluster <= length(Colors)) {
    binded$color <- Colors[match(binded$Cluster, names(Colors))]
  } else {
    # Warning message
    cat("The number of clusters is higher than the number of colors specified. A default color palette has been used.\n")
    # Create a color palette based on the default rainbow and assign it to the clusters
    colors <- rainbow(n_cluster)
    binded$color <- colors[match(binded$Cluster, as.character(1:n_cluster))]
  }
  
  # Assign colors to the convex hulls and add color for the global convex hull
  if (n_cluster <= length(Colors)) {
    all_hulls_df$color <- Colors[match(all_hulls_df$Cluster, names(Colors))]
    all_hulls_df$color[all_hulls_df$Cluster == "global"] <- global_hull_color
  } else {
    colors <- rainbow(n_cluster)
    all_hulls_df$color <- colors[match(all_hulls_df$Cluster, as.character(1:n_cluster))]
    all_hulls_df$color[all_hulls_df$Cluster == "global"] <- global_hull_color
  }
  
  # Plot using ggplot2
  ggplot(binded, aes(x = PC1, y = PC2)) +
    geom_polygon(data = all_hulls_df, aes(x = x, y = y, group = Cluster, fill = color), alpha = 0.7) +
    geom_point(aes(color = color, fill = color)) +
    labs(title = paste0(Title, " Functional Clusters"), x = "PCoA1", y = "PCoA2")+
    scale_color_identity() + # color as specified in the dataset
    scale_fill_identity() + # color as specified in the dataset
    theme_classic() + # classic theme
    theme(panel.border = element_rect(colour = "black", fill=NA)) + # border of the plot
    theme(axis.text = element_text(face = "bold", size = 12)) + # bold axis text
    theme(axis.title = element_text(face = "bold")) + # bold axis title
    theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid")) +
    theme(legend.position = "none") + # remove legend
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) # centered bold title
}

# Example usage (assuming fd.coord, binded, and n_cluster are defined)
colors <- c("1"= "black", "2" = "red", "3" = "green", "4" = "yellow", "5" = "#66FFFF", "6" = "cadetblue", "7" = "orange", "8" = "purple")

create_convex_hull_plot(fd.coord, binded, n_cluster, Colors = colors, Title = "Global")


# We can do the same for each Site 
# First we retrieve the number of different Site 
Sites = unique(data$Site)

# Then we create a binded data frame for all site with information of abundance of FE that we can filter to retrieve only FE present in a given site
# We retrieve abundance of FE for each site
fe.conditions = t(ab.fe.conditions)
# we merge it to the binded data frame
fe_binded<-cbind(binded,fe.conditions)
# We transform it into a tidy format
fe_binded_tidy <- fe_binded %>%
  rownames_to_column(var = "FE") %>%
  gather(key = "Year", value = "Abundance", -Cluster, -PC1, -PC2, -PC3, -PC4, -FE)

# We iterate for each Site 
Site_Functionnal_Clusters <- lapply(Sites,function(site){
  
  # We retrieve the years in which the site is present
  position <-which(str_detect(unique(fe_binded_tidy$Year), site))
  Years= unique(fe_binded_tidy$Year)[position]
  
  #filter the data to retrieve only the FE present in the site
  Site_Fe_binded <- fe_binded_tidy %>%
    filter(Year %in% Years) %>% 
    filter(Abundance > 0) %>% 
    select(-Abundance)
  
  # We transform this data to have to right format for the convex hull plot
  Site_binded <- Site_Fe_binded %>%
    select(-Year) %>% # supress the year column 
    distinct(FE, .keep_all = TRUE) %>% # keep only one row per FE
    column_to_rownames(var = "FE") # set the FE as rownames
    
  # We apply the same method as before to create the clusters
  create_convex_hull_plot(fd.coord = fd.coord, Site_binded, n_cluster, Colors = colors, Title = site)
}) 

Site_Functionnal_Clusters

# Add a title to each plot
Site_Functionnal_Clusters <- lapply(1:length(Site_Functionnal_Clusters), function(i) {
  Site_Functionnal_Clusters[[i]] + ggtitle(paste0(Sites[i]))
})

# For each element of the list, we add the Site as entry
names(Site_Functionnal_Clusters) <- Sites

# We draw the plot in one using the gridExtra package
library(gridExtra)
grid.arrange(grobs = Site_Functionnal_Clusters, ncol = 5)

# We can save the plots in a jpg and a svg file 
ggsave("Site_Functionnal_Clusters.jpg", grid.arrange(grobs = Site_Functionnal_Clusters, ncol = 5), width = 20, height = 20, units = "cm")
ggsave("Site_Functionnal_Clusters.svg", grid.arrange(grobs = Site_Functionnal_Clusters, ncol = 5), width = 20, height = 20, units = "cm")


#### 2. Number of species per functional cluster in each site and Time point ####

### We load the file with the N of species per cluster previously calculated and obtained in STEP 1. 

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
