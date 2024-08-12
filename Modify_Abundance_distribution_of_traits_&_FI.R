#This script needs to be run after 1.Functional space & Frich 
#______________________________________________________________# 

##### Load required packages #####
# install.packages("FD") # FD is a package for functional diversity measures
library('FD')
# install.packages("tripack") # Triangulation of irregularly spaced data
library('tripack')
# install.packages("geometry") # Geometric operations on point sets
library('geometry')
# install.packages("matrixStats") # Functions that operate on rows and columns of matrices
library('matrixStats')
# install.packages("tidyverse") # Tidyverse packages : ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats
library('tidyverse')
# install.packages("gridExtra") # Arrange multiple grid-based grobs on a page
library('gridExtra')
# install.packages("ggrepel") # Repel overlapping text labels
library('ggrepel')
# install.packages("vegan") # Community Ecology Package
library('vegan')
# install.packages("cluster") # Cluster Analysis Extended
library('cluster') # For the Gower distance matrix and PAM

##### Load Data #####
# Load FEs data
fes <- read.csv2("./Raw_files/FE_ordered.csv", sep=";", dec=",", row.names=1)
#Load Species and FEs data
spe_fes <- read.csv2("./Raw_files/Species_FE.csv", sep=";", dec=",", row.names=1)
#Load Abundance data
ab <- read.csv2("./Raw_files/Abundances.csv", sep=";", dec=",", row.names=1)
#Load sites and quadrats metadata
sites <- read.csv2("Sites.csv", sep=";", dec=",", row.names=1)

#### 1. Abundance distributions of traits in the functional space ####
# Retrieve all conditions
conditions <- unique(sites$Year) # change by your own condition if other in your metadataset
# or specify conditions manually but be careful to use the same name as in the sites metadataset
# conditions <- c("condition_1","condition_2","condition_3",...,"condition_n")

# Data manipulation and arrangements
# Compute the abundance of the species in the different conditions (mean cover)
ab.conditions <- lapply(conditions, function(x) {
  # get the sites coresponding to the condition
  quad <- rownames(sites[sites$Year == x,])
  # get the abundance of the species in the condition as the mean of the abundance of the species in the sites
  colMeans(ab[rownames(ab) %in% quad,])
})#eo lapply

# Transform the list into a matrix
ab.conditions <- do.call(rbind, ab.conditions)

# Add the condition name as row name
rownames(ab.conditions) = conditions

# Transform the matrix into a data frame with species as row names and conditions as columns
Species_Weights <- t(ab.conditions) %>%  # Transpose the matrix (invert rows and columns)
  as.data.frame() %>% # Transform the matrix into a data frame
  rownames_to_column(var = "Species") # Add the species as a collumn

# Compute the standard deviation of the abundance of the species in the different conditions (individual variability of each condition)
ab.conditions_Sd <- lapply(conditions, function(x) {
  # get the sites coresponding to the condition
  quad <- rownames(sites[sites$Year == x,])
  # get the Sd of these abundances
  colSds(as.matrix(ab[rownames(ab) %in% quad,])) #
})#eo lapply

# Transform the list into a matrix
ab.conditions_Sd <- do.call(rbind, ab.conditions_Sd)

# Add the condition name as row name
rownames(ab.conditions_Sd) = conditions

# get the levels of the FEs : the different FEs
fes <- levels(as.factor(spe_fes$FE))

# compute the abundance of each FEs in the different conditions
ab.fe.conditions <- lapply(conditions, function (z) {
  # In each condition, compute the abundance of each FEs
  abund.fes <-  sapply(fes, function (x) {
    # get row names (Species) for the specific Fe
    species <- rownames(spe_fes)[which(spe_fes$FE == x)]
    # Calculate the abundance of the FEs in fonction of the abundance of the corresponding species
    sum(ab.conditions[z,species])
  })#eo sapply
  # return the abundance of the FEs in the condition
  abund.fes
})#eo lapply

# Add the condition name as row name
names(ab.fe.conditions) = conditions

# Transform the list into a matrix
ab.fe.conditions <- do.call(rbind, ab.fe.conditions)

#### 2. Functional identity (FI) trends ####
##### 2.1. Weighted Functional identity (FI) of the assemblages #####
# Retrieve the weighted centroids of the assemblages in trait space (Functional identity)

# load coordinates previously generated when creating the functional space
# read fd coord as a data frame with first collum = rownames and precise that data are numeric
fd.coord <- read.csv2("FE_4D_coord.csv", sep=",", dec=".", row.names=1, header=TRUE)
fd.coord.FE <- fd.coord %>% 
  rownames_to_column(var = "FE") # add the FE as a column

# retrieve the FE of the species in a data frame
spe_fes_df <- rownames_to_column(spe_fes, var = "Species")

# merge the coordinates and the species FEs to get the coordinates of the species in the functional space
FE_4D_SP_coord <- merge(fd.coord.FE, spe_fes_df, by = "FE") %>% 
  select(-FE) %>% # remove the first column
  column_to_rownames(var = "Species") # set the row names to the species

# Save this file to use it for further analyses
write.csv(FE_4D_SP_coord, file="FE_4D_SP_coord.csv")

# Obtaining weighted centroids for each assemblage and Time point 
# The weighted centroid is the sum of the product of the coordinates of the species in the functional space by the weight of the species. This species weight is the abundance of the species in the assemblage. 
# For each Condition, we calculate the weighted centroid of the assemblage in the functional space
# To this we multiply the coordinates of the species in the functional space by the abundance of the species in the assemblage and sum the results. 
# Then we average the results by the total abundance of the assemblage to get the weighted centroid of the assemblage in the functional space.

# First we merge The weight dataset to the species coordinates dataset using bindrows. 
# By doing this we will not have any doubt about the order of the species in the two datasets.
Species_Weights_coord <- left_join(Species_Weights, FE_4D_SP_coord %>% 
                                     rownames_to_column(var= "Species"), by = "Species")

# Then for each condition we calculate weighted species coordinates in the functional space
Weighted_centroid <- lapply(conditions, function (x) {
  # get the weighted coordinates of the species in the functional space for the condition z
  Condition_Weighted_coord <- Species_Weights_coord %>% 
    select(c(x, PC1, PC2)) %>% # select the columns of interest
    mutate(weighted_PC1 = PC1 * Species_Weights_coord[,x], # multiply the coordinates of the species in the functional space by the abundance of the species in the assemblage
           weighted_PC2 = PC2 * Species_Weights_coord[,x]) # multiply the coordinates of the species in the functional space by the abundance of the species in the assemblage
  
  # Sum the results to get the weighted centroid of the assemblage in the functional space and divide by the total abundance of the assemblage
  PC1 <- sum(Condition_Weighted_coord$weighted_PC1) / sum(Species_Weights_coord[,x])
  PC2 <- sum(Condition_Weighted_coord$weighted_PC2) / sum(Species_Weights_coord[,x])
  # return the weighted centroid of the assemblage in the functional space
  c(PC1, PC2)
}) #eo lapply
    
# Save the results in a data frame
Weighted_centroid <- do.call(rbind, Weighted_centroid) # Transform the list into a matrix
colnames(Weighted_centroid) <- c("PC1", "PC2") # Add column names
rownames(Weighted_centroid) <- conditions # Add row names

Weighted_centroid

##### 2.2. Plotting of the FI trends through time ####
# Plot the weighted centroids of the assemblages in the functional space

# First, we retrieve the global potential functional space (background polygon)
m2 <- fd.coord[rownames(fd.coord),] # Get the coordinates of the functional space
tr2 <-tri.mesh(m2[,1],m2[,2]) # Triangulate the functional space
ch2 <- convex.hull(tr2) # Compute the convex hull of the functional space

# Create a data frame with the convex hull coordinates
hull_df_global <- data.frame(x = ch2$x, # x coordinates
                             y = ch2$y, # y coordinates
                             i = ch2$i) # index of this point

# Then, for each condition, we plot with ggplot the weighted centroid of the assemblage in the functional space
FI_plots <- lapply(conditions, function (x) {
  # filter to keep only species that have an abundance > 1 in the condition
  # Set the threshold to the percentage of abundance you want to consider
  threshold <- 0
  
  # Get the species with an abundance > threshold
  ab.fe.conditions.temp <- ab.fe.conditions[x,] %>% 
    as.data.frame() %>% 
    filter(ab.fe.conditions[x,]>threshold) %>% 
    mutate(size = as.numeric(.)) # transform the abundance into numeric
  
  fd.coord.temp <- fd.coord[rownames(fd.coord) %in% rownames(ab.fe.conditions.temp),] # Get the coordinates of the species with an abundance > threshold
  
  # Plot the weighted centroid of the assemblage in the functional space
  ggplot(data = fd.coord, aes(x = PC1, y = PC2)) +
    coord_fixed() + # Set the aspect ratio
    geom_polygon(data = hull_df_global, aes(x = x, y = y), fill = "#CCCCCC30", color = NA) + # plot the background polygon
    geom_point(data = fd.coord.temp, aes(x = PC1, y = PC2, size = ab.fe.conditions.temp$size), shape = 21, fill = as.data.frame(cols)[x,], color = as.data.frame(colstr)[x,], alpha = 0.5)+ # plot the FE in the functional space
    scale_size_continuous(range = c(1, max(ab.fe.conditions.temp$size)*0.5)) + # set the range of the size of the points as a function of the abundance of the FEs (size of the points is proportional to the abundance of the FEs) 0.5 is a factor to adjust the size of the points
    # plot the centroid of the assemblage in the functional space
    annotate("text", x = Weighted_centroid[x,"PC1"], y = Weighted_centroid[x,"PC2"], label = "+", color = as.data.frame(cols)[x,], size = 14)+ # plot the centroid of the assemblage in the functional space
    labs(title = x, x = "PCoA 1", y = "PCoA 2") + #axis label 
    theme_classic()+ # classic theme
    theme(panel.border = element_rect(colour = "black", fill=NA))+ # border of the plot
    theme(axis.text = element_text(face = "bold", size = 12))+ # bold axis text
    theme(axis.title = element_text(face = "bold"))+ # bold axis title
    theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid"))+
    theme(legend.position = "none")+ # remove legend
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))# centered bold
}) #eo lapply

# FI_plots # please uncomment this line to see the plots

names(FI_plots) = conditions #Add the condition name as row name

# We use the previous code to plot the FI trends for all conditions
FI_plot <- display_multiple_plots(FI_plots, layout = "column") # Display the plots in one

# We can save the plot in a svg and png file
ggsave("FI_plot.svg", plot = FI_plot, device = "svg", width = 40, height = 20, units = "cm")
ggsave("FI_plot.png", plot = FI_plot, device = "png", width = 40, height = 20, units = "cm")


# We can want to have all plots of one site but the differents years in the same plot 
# First we retrieve the sites
unique_sites <- unique(sites$Site)

# For each site, we retrieve the years and plot the FI trends for each year
Plot_Sites_all_years <- lapply (unique_sites, function(x){
  # retrieve the years for the site
  years <- unique(sites[sites$Site == x,]$Year)
  
  # retieve abundance of the FE in the different conditions for the site and filter the FEs with an abundance > threshold
  threshold <- 0
  ab.conditions_site <- lapply(years, function (z) {
    ab.conditions_year <- ab.fe.conditions[z,] %>% # get the abundance of the FEs in the condition z
      as.data.frame() # transform the matrix into a data frame
  }) # eo lapply
  
  # merge the data frames into one and rename each column with the corresponding year
  ab.conditions_site <- do.call(cbind, ab.conditions_site)
  colnames(ab.conditions_site) <- years
  
  # Tranform the data frame to have all values in one collum named abundance and one more variable call year
  data_temp <- ab.conditions_site %>%
    rownames_to_column(var = "FE") %>% # add the FE as a column
    gather(key = "Year", value = "abundance", -FE) %>%  # transform the data frame to have all values in one collum named abundance and one more variable call year
    filter(abundance > threshold) %>% # filter the FEs with an abundance > threshold
    # merge the data frame with the coordinates of the FEs in the functional space
    left_join(as.data.frame(fd.coord) %>% rownames_to_column(var = "FE"), by = "FE")
  
  # We can now plot the FI trends for the site with all years in the same plot separated by a different color as previously using colors from the cols and colstr data frames
  # dependind on the length of the number of years, we add more geom-point layer to the plot to have all years in the same plot
  # we retrieve the number of years
  n_years <- length(years)
  
  # we define non conditional treatment
  plot_base <- ggplot(data = fd.coord, aes(x = PC1, y = PC2)) +
    coord_fixed() + # Set the aspect ratio
    geom_polygon(data = hull_df_global, aes(x = x, y = y), fill = "#CCCCCC30", color = NA)+ # plot the background polygon
    labs(title = x, x = "PCoA 1", y = "PCoA 2") + #axis label 
    theme_classic()+ # classic theme
    theme(panel.border = element_rect(colour = "black", fill=NA))+ # border of the plot
    theme(axis.text = element_text(face = "bold", size = 12))+ # bold axis text
    theme(axis.title = element_text(face = "bold"))+ # bold axis title
    theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid"))+
    theme(legend.position = "none")+ # remove legend
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))# centered bold
  
  #then we use conditional treatment : 
  if (n_years == 1){
    plot_base +
      geom_point(data = subset(data_temp, Year == years[1]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      scale_color_manual(values = cols[years]) + # set the color of the points as a function of the year
      scale_fill_manual(values = colstr[years]) + # set the fill color of the points as a function of the year
      scale_size_continuous(range = c(1, max(data_temp$abundance)*0.5)) # set the range of the size of the points as a function of the abundance of the FEs (size of the points is proportional to the abundance of the FEs) 0.5 is a factor to adjust the size of the points
  } else if (n_years == 2){
    plot_base + 
      geom_point(data = subset(data_temp, Year == years[1]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[2]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      scale_color_manual(values = cols[years]) +
      scale_fill_manual(values = colstr[years]) +
      scale_size_continuous(range = c(1, max(data_temp$abundance)*0.5))+
      annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols)[years,], size = 14)# plot the centroid of the assemblage in the functional space
  } else if (n_years == 3) {
    plot_base + 
      geom_point(data = subset(data_temp, Year == years[1]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[2]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[3]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      scale_color_manual(values = cols[years]) +
      scale_fill_manual(values = colstr[years]) +
      scale_size_continuous(range = c(1, max(data_temp$abundance)*0.5))+
      annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols)[years,], size = 14)
  } else if (n_years == 4) {
    plot_base + 
      geom_point(data = subset(data_temp, Year == years[1]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[2]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[3]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[4]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      scale_color_manual(values = cols[years]) +
      scale_fill_manual(values = colstr[years]) +
      scale_size_continuous(range = c(1, max(data_temp$abundance)*0.5))+
      annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols)[years,], size = 14)
  } else if (n_years == 5) {
    plot_base + 
      geom_point(data = subset(data_temp, Year == years[1]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[2]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[3]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[4]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[5]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      scale_color_manual(values = cols[years]) +
      scale_fill_manual(values = colstr[years]) +
      scale_size_continuous(range = c(1, max(data_temp$abundance)*0.5))+
      annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols)[years,], size = 14)
  } else if (n_years == 6) {
    plot_base + 
      geom_point(data = subset(data_temp, Year == years[1]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[2]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[3]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[4]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[5]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[6]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      scale_color_manual(values = cols[years]) +
      scale_fill_manual(values = colstr[years]) +
      scale_size_continuous(range = c(1, max(data_temp$abundance)*0.5))+
      annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols)[years,], size = 14)
  } else if (n_years == 7) {
    plot_base + 
      geom_point(data = subset(data_temp, Year == years[1]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[2]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[3]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[4]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[5]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[6]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[7]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      scale_color_manual(values = cols[years]) +
      scale_fill_manual(values = colstr[years]) +
      scale_size_continuous(range = c(1, max(data_temp$abundance)*0.5))+
      annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols)[years,], size = 14)
  } else if (n_years == 8) {
    plot_base + 
      geom_point(data = subset(data_temp, Year == years[1]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[2]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[3]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[4]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[5]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[6]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[7]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[8]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      scale_color_manual(values = cols[years]) +
      scale_fill_manual(values = colstr[years]) +
      scale_size_continuous(range = c(1, max(data_temp$abundance)*0.5))+
      annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols)[years,], size = 14)
  } else if (n_years == 9) {
    plot_base + 
      geom_point(data = subset(data_temp, Year == years[1]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[2]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[3]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[4]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[5]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[6]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[7]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[8]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[9]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      scale_color_manual(values = cols[years]) +
      scale_fill_manual(values = colstr[years]) +
      scale_size_continuous(range = c(1, max(data_temp$abundance)*0.5))+
      annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols)[years,], size = 14)
  } else if (n_years == 10) {
    plot_base + 
      geom_point(data = subset(data_temp, Year == years[1]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[2]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[3]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[4]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[5]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[6]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[7]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[8]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[9]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      geom_point(data = subset(data_temp, Year == years[10]),
                 aes(x = PC1, y = PC2, size = abundance, color = Year, fill = Year),
                 shape = 21,
                 alpha = 0.5)+ # plot the FE in the functional space
      scale_color_manual(values = cols[years]) +
      scale_fill_manual(values = colstr[years]) +
      scale_size_continuous(range = c(1, max(data_temp$abundance)*0.5))+
      annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols)[years,], size = 14)
  } else {
    print("Too many years to plot, the maxium is 10 years, if you want to plot more years, please refer to the code to add more year layers")
  } #eo if else
}) #eo lapply

# Plot_Sites_all_years # Display the plots , uncomment this line to see the plots

# add the site name to the list
names(Plot_Sites_all_years) <- unique_sites

Plot_FI_all_years <- do.call(grid.arrange, c(Plot_Sites_all_years, nrow = 1))

# We can save the plot in a svg and png file
ggsave("FI_plot_all_years.svg",Plot_FI_all_years, width = 40, height = 20, units = "cm")
ggsave("FI_plot_all_years.png",Plot_FI_all_years, width = 40, height = 20, units = "cm")

#### 3. Permutational Analysis of Variance (PERMANOVA) of Functional Identity ####
# We use the weighted centroids (based on weighted coordinates) of the community to conduct the PERMANOVA analysis
# This PERMANOVA analysis will test the effect of the time on the weighted centroids of the community 
# As weighted centroids are based on the weighted coordinates of the species, we will use the Canberra distance to measure the distance between the weighted centroids
# We use only first two dimensions of the functional space to calculate the weighted centroids, because the first two dimensions explain the most of the variance in the functional space and it will be easier to interpret the results on the 2D plot

# We calculate the weighted centroids of the community for each quadrat
ab_P <- t(ab) # transpose the abundance matrix
weighted_axis1 = FE_4D_SP_coord$PC1 * ab_P
Axis1sum <- colSums(weighted_axis1)
weighted_axis1_coordinates <- Axis1sum/100

weighted_axis2 = FE_4D_SP_coord$PC2 * ab_P
Axis2sum <- colSums(weighted_axis2)
weighted_axis2_coordinates <- Axis2sum/100

weighted_axis3 = FE_4D_SP_coord$PC3 * ab_P
Axis3sum <- colSums(weighted_axis3)
weighted_axis3_coordinates <- Axis3sum/100

weighted_axis4 = FE_4D_SP_coord$PC4 * ab_P
Axis4sum <- colSums(weighted_axis4)
weighted_axis4_coordinates <- Axis4sum/100


### We creat a data frame that contains the obtained weighted "x" "y" coordinates. for each quadrat. 
weighted_centroids <- cbind(weighted_axis1_coordinates, weighted_axis2_coordinates)
weighted_centroids = as.data.frame(weighted_centroids)

#We complete the data frame with more useful info for the analysis.
#Site
data <-  weighted_centroids %>%
  rownames_to_column(var = "rownames") %>%
  separate(rownames, c("Year", "site", "Type", "Quadrat"), sep = "_") %>%
  # merge collums Year and Site
  mutate(Site = paste(site,Type, sep = "_")) %>%
  mutate(Years = paste(Year,site,Type, sep = "_")) %>% 
  select(-Year, -site, -Type, -Quadrat)

### We use this new data frame to conduct the PERMANOVA ANALYSIS 
# we apply for each site the adonis2 function (PERMANOVA) to test changes on the weighted centroids of the community over time
# we use the Canberra distance to measure the distance between the weighted centroids as the canberra distance is more appropriate for compositional data than the Euclidean distance see legendre et al. 2012
# we use the Time as a factor to test the effect of the time on the weighted centroids of the community

n_perm <- 9999 # number of permutations

PERMANOVA_test <- lapply(unique(data$Site), function(x){
  data_site <- data[data$Site == x,]
  per_site <- adonis2(cbind(data_site$weighted_axis1_coordinates,
                            data_site$weighted_axis2_coordinates, 
                            data_site$weighted_axis3_coordinates, 
                            data_site$weighted_axis4_coordinates) ~ data_site$Years, method="canberra", permutations = n_perm)
  per_site
})

# we can rename the results with the site names
names(PERMANOVA_test) <- unique(data$Site)

# We can merge the results in a data frame
PERMANOVA_test_df <- do.call(rbind, PERMANOVA_test)
# We can add the site names to replace the names data_site$Years

PERMANOVA_test_df



