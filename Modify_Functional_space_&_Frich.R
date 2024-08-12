#### Workspace ####
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


##### Load Data #####
# Load FEs data
fes <- read.csv2("./Raw_files/FE_ordered.csv", sep=";", dec=",", row.names=1)
#Load Species and FEs data
spe_fes <- read.csv2("./Raw_files/Species_FE.csv", sep=";", dec=",", row.names=1)
#Load Abundance data
ab <- read.csv2("./Raw_files/Abundances.csv", sep=";", dec=",", row.names=1)
#Load sites and quadrats metadata
sites <- read.csv2("Sites.csv", sep=";", dec=",", row.names=1)

#### 1. Create functional Space ####
#load "Quality functionnal space" function
#Please read functional_trait_analysis.qmd to further informations. 
source("./Raw_files/quality_funct_space.R") # change by current function directory 

#Calculate Quality space, please change information to adapt to your own dataset
qfs <- quality_funct_space(mat_funct = fes,
                           traits_weights=NULL, 
                           nbdim=11,
                           metric="Gower",
                           dendro=FALSE, 
                           plot="quality_funct_space")

# Verify space quality
# quality of spaces (low meanSD = high quality)
round( qfs$meanSD , 4)

# keeping Functional entities coordinates on the 4 first dimensions, meanSD<0.004
fd.coord <- qfs$details_funct_space$mat_coord[,1:4]

# write result for further analyses
write.csv(fd.coord, file="FE_4D_coord.csv") #to use it for further analyses

#see variance explained by the PCoA axes
gower<-qfs$details_funct_space$mat_dissim # store matrix 
fit <- cmdscale(gower,eig=TRUE, k=4) # PCoA
#variance explained by the axes
cumsum(fit$eig[fit$eig>=0]) / sum(fit$eig[fit$eig>0])


#### 2. Explore trends in functional richness ####
# Retrieve all conditions
conditions <- unique(sites$Year) # change by your own condition if other in your metadataset
# or specify conditions manually but be careful to use the same name as in the sites metadataset
# conditions <- c("condition_1","condition_2","condition_3",...,"condition_n")

# rearange data to retrieve abundances for each specific conditions
ab.conditions <- lapply(conditions, function(x) {
  # filter the quadrats for the specific condition
  quad <- rownames(sites[sites$Year == x,])
  # filter the abundance matrix for the specific condition
  colSums(ab[rownames(ab) %in% quad,])
})#eo lapply

# merge the list into a matrix
ab.conditions <- do.call(rbind, ab.conditions)
# set the rownames to the conditions
rownames(ab.conditions) = conditions

# Calculation of relative richness as global functional space are calculate under the condition that all species are present in a theoricall space 
# To calculate relative richness, we need to calculate the total number of species in the dataset. 
NbSp_tot = nrow(spe_fes) #number of species specied in the species_FE dataset

#Calculate convex hull for the global functional space (111 taxonomic units in total for calculating the relative richness)

#### Calculate convex hull (111 taxonomic units in total for calculating the relative richness)
chg <- convhulln(fd.coord, options = "FA") # Calculate convex hull for the global functional space

Frich <- lapply(conditions, function (x) {
  #For each condition ... 
  # Identify species where abundance > 0
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  # Subset spe_fes to only include FE of present species (rows with names in species)
  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
  # Subset fd.coord to only include rows with names in fes_cond (only coordinates of present FE)
  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
  # Compute the convex hull of m and fd.coord using the convhulln function
  ch <- convhulln(m, options = "FA") # For actual condition
  # Return a vector with several calculated values
  c(length(species), length(species)/NbSp_tot*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol)
}) # End of lapply

# Name the list with each conditions
names(Frich) = conditions

# Fric contains the number of species(NbSp) and FEs (NbFEs), relative percentages (NbSpP,NbFEsP ), and the 4D convex hull volume (Frich (Vol4D)) among the 3 temporal points for each habitat
# Merge the list into a matrix
Frich <- do.call(rbind, Frich)

# set the column names
colnames(Frich) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Frich (Vol4D)")

# Plot convex hull volume for each condition with ggplot2
# Choose manually the colors for each condition, the number of conditions should be the same as the number of colors
# Choose colors manualy to fill convex hull polygons 
cols <- c("Pzzu_cor_2003" = "#CD2626",
          "Pzzu_cor_2011" = "#F9D71C", 
          "Pzzu_cor_2018" = "#3A5FCD",
          "Pzzu_par_2006" = "#CD2626",
          "Pzzu_par_2011" = "#F9D71C",
          "Pzzu_par_2018" = "#3A5FCD",
          "Gabin_par_1999" = "#CD2626",
          "Gabin_par_2007" = "#F9D71C", 
          "Gabin_par_2009" ="#3A5FCD",
          "Pzzinu_par_2006" = "#CD2626",
          "Pzzinu_par_2011" = "#F9D71C",
          "Pzzinu_par_2016" ="#3A5FCD", 
          "Passe_cor_2006" = "#CD2626",
          "Passe_cor_2011" = "#F9D71C",
          "Passe_cor_2018"= "#3A5FCD")
colstr <- c("Pzzu_cor_2003" = "#CD2626",
            "Pzzu_cor_2011" = "#F9D71C", 
            "Pzzu_cor_2018" = "#3A5FCD",
            "Pzzu_par_2006" = "#CD2626",
            "Pzzu_par_2011" = "#F9D71C",
            "Pzzu_par_2018" = "#3A5FCD",
            "Gabin_par_1999" = "#CD2626",
            "Gabin_par_2007" = "#F9D71C", 
            "Gabin_par_2009" ="#3A5FCD",
            "Pzzinu_par_2006" = "#CD2626",
            "Pzzinu_par_2011" = "#F9D71C",
            "Pzzinu_par_2016" ="#3A5FCD", 
            "Passe_cor_2006" = "#CD2626",
            "Passe_cor_2011" = "#F9D71C",
            "Passe_cor_2018"= "#3A5FCD")

#computing the global convex hull to be represented as background of each plot.
#The global trait space for all species (100% of space occupied / Frich = 1)
m2 <- fd.coord[rownames(fd.coord),]#select all coordinates of the functional entities (FE)
tr2 <-tri.mesh(m2[,1],m2[,2]) #Create Delaunay triangulation 
ch2 <- convex.hull(tr2) #Compute convex hull

# Create a data frame with the convex hull coordinates
hull_df_global <- data.frame(x = ch2$x, 
                             y = ch2$y, 
                             i = ch2$i)
# add a group column to the data frame
hull_df_global$group <- "Global"
conditions
# Define a function to calculate the convex hull for each condition and create a plot
Plots <- lapply(conditions,function(x){
  # Identify species where abundance > 0
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  # Subset spe_fes to only include FE of present species (rows with names in species)
  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
  # Subset fd.coord to only include rows with names in fes_cond (only coordinates of present FE)
  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]#subset coordinates of present FEs
  tr <-tri.mesh(m[,1],m[,2]) #Create Delaunay triangulation of a set of point m[,1] is the first colum containing PC1 coordinates
  ch <- convex.hull(tr) #Computes convex hull
  # Create a data frame with the convex hull coordinates
  hull_df <- data.frame(x = ch$x, 
                        y = ch$y, 
                        i = ch$i)
  # add a group column to the data frame
  hull_df$group <- x
  
  # merge dataset into a single dataframe
  hull_df_combined <- bind_rows(hull_df_global, hull_df)
  
  # Create the plot
  Plot <- ggplot(hull_df_combined, aes(x = x, y= y)) +
    geom_polygon(data = hull_df_global, aes(group = group, fill = group, color = group), linewidth = 1, alpha = 0.5) + # global convex hull
    geom_polygon(data = hull_df, aes(group = group, fill = group, color = group),linewidth = 1, alpha = 0.5) + # specific convex hull
    labs(title = x, x = "PCoA 1", y = "PCoA 2") + #axis label 
    theme_classic()+ # classic theme
    theme(panel.border = element_rect(colour = "black", fill=NA))+ # border of the plot
    theme(axis.text = element_text(face = "bold", size = 12))+ # bold axis text
    theme(axis.title = element_text(face = "bold"))+ # bold axis title
    theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid"))+
    scale_fill_manual(values = cols)+ # set fill colors
    scale_color_manual(values = colstr)+ # set border colors
    theme(legend.position = "none")+ # remove legend
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+# centered bold title
    annotate("text", x = 0.3, y = 0.3, label = paste("Frich = ", round(Frich[x,5],2)), size = 5, color = "black", fontface = "bold")+
    geom_point(data= m[,1:2], aes( x= PC1 , y = PC2), colour = colstr[x])
  return(Plot)
})

# Plots #uncomment to draw individual plots

# Add the name of the condition to the list
names(Plots) <- conditions

# display all plots corresponding to a same grouping factor (Site) in a single line or collum of the plot

display_multiple_plots <- function(plots,layout = "column"){
  plots = FI_plots
  # Retrieve Sites in the name of each condition (condition less the year)
  Sites <- unique(sapply(names(plots), function(x) {substr(x, 1, nchar(x)-5)}))
  
  # Associate each condition to a site from the names of the plot 
  sites <- sapply(names(plots), function(x) {substr(x, 1, nchar(x)-5)})
  
  # create a list of plots for each site
  plots_site <- lapply(Sites, function(x) {
    # filter the plots for the specific site
    plots[which(sites == x)]
  })#eo lapply
  
  # retrieve the number of Sites 
  nSites = length(Sites)
  
  # retrieve wich is the maximum number of plots for a site
  max_plots = max(sapply(plots_site, length))
  
  # In function of the layout, we will draw all plots with one line per site if layout = "row" or one column per site if layout = "column". 
  # To do so we adjust the position of each plot before running the grid.arrange function depending on the layout 
  
    # The position of each plot is adjusted in the Plots list to have the same number of plots for each site by adding empty plots if necessary
  Plots_extended <- lapply(plots_site, function(x) {
    
    # add empty plots to the list of plots for each site to have the same number of plots for each site
    if (length(x) < max_plots) {
      # add empty plots
      empty_plots <- ggplot(data = data.frame(x = 1, y = 1), aes(x = x, y = y)) +
        geom_point(aes(x = x, y = y), pch = 16, col = "white", cex = 4) + # plot empty plot
        theme_void() # remove all theme
      # transform this empty plot into a list
      empty_plots <- list(empty_plots)
      x <- c(x, empty_plots)
    }
    x # return the list of plots for each site
  }) #eo lapply
  
  # Merge and arrange the list to have all plots in the same list with all first poistion of the list corresponding to the first plot of each site, all second position of the list corresponding to the second plot of each site, etc.
  if (tolower(layout) == "column") {
    plots_drawing <- do.call(c, lapply(1:max_plots, function(x) {
      lapply(Plots_extended, function(y) {
        y[[x]]
      })
    }))#eo lapply
  } else if (tolower(layout) == "row") {
    plots_drawing <- do.call(c, Plots_extended)
  } else {
    print("error, layout should be 'row' or 'column'")
  }
  
  # Draw all plots of each site in a single plot containing all sites in a single line or column for each list of plots
  # draw all plot in one (don't forget to change gridextra)
  if (tolower(layout) == "column") {
    draw_all_plot <- do.call(grid.arrange, c(plots_drawing, ncol = nSites)) # addapt to your dataset
  } else if (tolower(layout) == "row") {
    draw_all_plot <- do.call(grid.arrange, c(plots_drawing, nrow = nSites)) # addapt to your dataset
  } else {
    print("error, layout should be 'row' or 'column'")
  }
}

Frich_Plot <- display_multiple_plots(Plots, layout = "Column")

# save the plot in a SVG file and a jpeg file specify dimension for your own dataset
ggsave("Frich_Plots.svg", Frich_Plot, width = 40, height = 20, units = "cm")
ggsave("Frich_Plots.jpeg", Frich_Plot, width = 40, height = 20, units = "cm")

#### 3. Functional Richness : Null hypothesis ####
# Retrieve all conditions 
Conditions <- unique(sites$Year)
# Or specify conditions manually
# Conditions <- c("condition_1","condition_2","condition_3",...,"condition_n")

# Retrieve list of sites on which the analysis will be performed
Sites = unique(sites$Site)
# Or specify the site manually
#Sites = c("site_1","site_2","site_3",...,"site_n")

Sites_Frich <- lapply(Sites, function(x){
  # filter conditions for the site of interest
  Specific_Conditions <- Conditions[grep(x, Conditions)]
  # filter the abundance matrix for the site of interest
  ab_specific_conditions = ab[grep(x, rownames(ab)),]
  
  # Data arrangements to retrieve abundances for each specific conditions
  ab.specific_condition <- lapply(Specific_Conditions, function(x) {
      # filter the quadrats for the specific condition
    quad <- rownames(sites[sites$Year == x,])
    # filter the abundance matrix for the specific condition and calculate the sum of abundance
    colSums(ab_specific_conditions[rownames(ab_specific_conditions) %in% quad,])
  })#eo lapply
  
  #replace the rownames of the abundance matrix with the conditions
  ab.specific_condition <- do.call(rbind, ab.specific_condition) #merge the list into a matrix
  rownames(ab.specific_condition) = Specific_Conditions #set the rownames to the conditions
  
  #Calculate convex hull for the global functional space
  Global_convexHull<- convhulln(fd.coord, options = "FA")
  
  #Calculate convex hull 
  #retrieve total number of species in the dataset
  NbSp_tot = nrow(spe_fes) #number of species specied in the species_FE dataset
  
  #Calculate convex hull for the global functional space (111 taxonomic units in total for calculating the relative richness)
  Frich <- lapply(Specific_Conditions, function (x) {
    # filter the species for the specific condition (species present in these samples)
    species <- colnames(ab.specific_condition)[which(ab.specific_condition[x,] > 0)]
    # filter the functional traits for this species
    fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
    # filter the functional space for these functional traits
    m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
    # calculate the convex hull for the functional space
    ch <- convhulln(m, options = "FA")
    # return the number of species, the relative number of species, the number of FEs, the relative number of FEs, and the volume of the convex hull
    c(length(species), length(species)/NbSp_tot*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/Global_convexHull$vol*100)
    
  })#eo lapply Frich
  
  #name the list with the conditions
  names(Frich) = Specific_Conditions
  # merge the list into a matrix
  Frich <- do.call(rbind, Frich)
  # name the columns
  colnames(Frich) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol4D")
  
  # Definition of the null model 
  # null model is based on the random selection of species in the dataset
  # Specify number of permutations
  n_perm = 100
  
  # create a copy of FEs to be used in the permutations
  spe_fes_r = spe_fes
  
  # Calculate permuted Frich 
  Frich_permuted <- lapply(Specific_Conditions, function (x) {
    # filter the species for the specific condition (species present in these samples)
    species <- colnames(ab.specific_condition)[which(ab.specific_condition[x,] > 0)]
    #create n_perm aleatory permutation of the functional traits for specified condition
    perm <- sapply((1:n_perm), function (z) {
      # shuffle the functional traits
      spe_fes_r$FE <- sample(spe_fes$FE)      
      # Filter the functional traits for the species present in specific condition
      fes_cond <- spe_fes_r[rownames(spe_fes_r) %in% species, ]
      # Filter the functional space for these functional traits
      m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
      # Calculate the convex hull for this random functional space
      ch <- convhulln(m, options = "FA")
      # Return the number of species, the relative number of species, the number of FEs, the relative number of FEs, and the volume of the convex hull
      c(length(species), length(species)/NbSp_tot*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
    })#eo sapply
    # Return the permuted Frich
    rownames(perm) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol")
    perm
  })#eo lapply Frich_permuted
  
  # Name the list with the conditions considered
  names(Frich_permuted) = Specific_Conditions
  
  Fric_perm_Quantiles <- lapply(Frich_permuted, function (x) {
    # Calculate the quantiles 0.05 and 0.95 of the permuted Frich
    rowQuantiles(x, probs=c(0.05, 0.95))
  })#eo lapply Fric_perm_Quantiles
  
  # retrieve actual Frich values of each condition
  Frich = as.data.frame(Frich)
  # Add the lower (0.05) and upper quantiles 0.95 to the Frich data frame
  Frich$lowerFE <- sapply(Specific_Conditions, function (x) { Fric_perm_Quantiles[[x]][3,1] })#eo sapply
  Frich$upperFE <- sapply(Specific_Conditions, function (x) { Fric_perm_Quantiles[[x]][3,2] })#eo sapply
  Frich$lowerVol <- sapply(Specific_Conditions, function (x) { Fric_perm_Quantiles[[x]][5,1] })#eo sapply
  Frich$upperVol <- sapply(Specific_Conditions, function (x) { Fric_perm_Quantiles[[x]][5,2] })#eo sapply
  
  # Add the condition as a factor to the Frich data frame
  Frich$cond <- as.factor(Specific_Conditions)
  
  # set the column names
  colnames(Frich) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol4D", "lowerFE", "upperFE", "lowerVol", "upperVol", "cond")
  
  # Add the condition as a factor to the Frich data frame
  Frich$year <- as.factor(as.character(str_sub(Specific_Conditions, start = -4)))
  
  # Plot the Frich values with the 95% confidence interval
  Plot <-ggplot(data = Frich, x = year , y = Vol4D)+
    geom_point(aes(x = year, y = Vol4D), pch = 16, col = cols[c(1:length(Frich$year))], cex = 4)+ # plot actual Frich values
    geom_linerange(aes(x = year, ymin = lowerVol, ymax = upperVol), col= cols[c(1:length(Frich$year))],linewidth = 2) + # plot 95% confidence interval
    scale_y_continuous(breaks=c(0,20,40,60,80,100), limits = c(0,100))+ # set y axis limits
    ylab("Relative richness (%)") + # y axis label
    xlab("Time")+ # x axis label
    ggtitle(x)+ # title
    coord_flip()+ # flip x and y axis
    scale_x_discrete(limits = as.character(rev(levels(Frich$year))))+
    theme_classic()+ # classic theme
    theme(panel.border = element_rect(colour = "black", fill=NA))+ # border of the plot
    theme(axis.text = element_text(face = "bold", size = 12))+ # bold axis text
    theme(axis.title = element_text(face = "bold", size = 12))+ # bold axis title
    theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid"))+
    theme(legend.position = "none")+ # remove legend
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+# centered bold title
    theme(axis.title.x = element_text(size=14))+ # set x axis title size
    theme(axis.title.y = element_text(size=14))+ # set y axis title size
    theme(axis.line.y = element_line(), # y axis line
          axis.line.x=element_line(), # x axis line
          panel.grid.major=element_blank(), # remove major grid
          panel.border= element_rect(colour = "black", fill=NA),
          panel.background=element_blank())
  
  return(Plot)
} )#eo lapply Sites_Frich

Sites_Frich

ncol = length(Sites_Frich)

# display the plots
Plots <- do.call(grid.arrange, c(Sites_Frich, ncol = ncol))

# Save the plot in a SVG file and a jpeg file specify dimension for your own dataset
ggsave("Frich_Plots_Sites.svg", Plots, width = 40, height = 5, units = "cm")
ggsave("Frich_Plots_Sites.jpeg", Plots, width = 40, height = 5, units = "cm")

#### 4. Trait vectors direction and traits in functional Space ####
#We use the cdmscale() and envfit() to re-ecreate our trait-space and find vector directions and categories positions in functional space

# We first create the background polygon (global functional space)
m2 <- fd.coord[rownames(fd.coord),] # select all coordinates of the functional entities (FE)
tr2 <-tri.mesh(m2[,1],m2[,2]) # Create Delaunay triangulation
ch2 <- convex.hull(tr2) # Compute convex hull

# Create a data frame with the convex hull coordinates
hull_df_global <- data.frame(x = ch2$x, # x coordinates
                             y = ch2$y, # y coordinates
                             i = ch2$i) # index of this point

# Based on the gower distance matrix, we calculate the PCoA axes : the ordination object with the eigenvectors and eigenvalues and the coordinates of the functional entities in the PCoA space.
fit <- cmdscale(gower,eig=TRUE, k=4) # PCoA

# For our variables, we calculate the vector direction in the functional space or the position of each categories in the functional space
efit <- envfit(fit, fes, na.rm=TRUE) # envfit

##### 4.1. PCoA with Traits Vectors ####
# We can plot the functional space with the vector direction for each trait (numerical)

# To enshure that the vector direction and the categories positions are in the same space, we use the same PCoA axes for all the variables. 
# To enshure that vectors have right length, we use the arrow.mul argument in the scores() function. But, we need to adjust the length of the vectors to fit the plot. 
# We can calculate the maximum and minimum values of the PCoA axes and adjust the length of the vectors to fit the plot consequently.
# Calculate max and min values of the PCoA axes
Max_x <- max(fd.coord[,1])
Min_x <- min(fd.coord[,1])
Max_y <- max(fd.coord[,2])
Min_y <- min(fd.coord[,2])

# We select the maximum distance from the origin to the maximum and minimum values of the PCoA axes
Max_dist <- max(c(abs(Max_x), abs(Min_x), abs(Max_y), abs(Min_y)))

# We can calculate the coefficient to adjust the length of the vectors based on the maximum distance from the origin to the maximum and minimum values of the PCoA axes. By doing this we can adjust the length of the vectors to fit the plot, as the longest vector will be equal to the maximum distance from the origin to the maximum and minimum values of the PCoA axes. To be sure, we use the round() function to round the coefficient to the nearest superior integer.
# As in R there is no function to round to the nearest superior float, we use the ceiling() function to round to the nearest superior integer. but made it a custom function to work with float
# Custom function to round to the nearest superior float
round_up <- function(x, digits = 0) {
  ### Round up to the nearest multiple of n
  n <- 10^digits # 10^0 = 1, 10^1 = 10, 10^2 = 100, etc.
  return (ceiling(x * n) / n) # Round up to the nearest upper integer and divide by n. ex 3.14 * 10 = 31.4, round up to 32 by ceiling, and 32/10 = 3.2
}

# Calculate the coefficient to adjust the length of the vectors
Adjust_coef <- round_up(Max_dist, 1)

# Extract the vector coordinates for continuous traits
vector_data <- scores(efit, "vector", arrow.mul = Adjust_coef) %>%  # scores ajusted to the coefficient
  as.data.frame() %>% # transform to data frame
  mutate(trait = rownames(.)) # add the trait name

# Create a plot with the functional space and the vector direction for each numerical trait
pcoa_vector_plot <- ggplot(fd.coord, aes(x = PC1, y = PC2)) + # Create the plot based on the PCoA coordinates
  coord_fixed() + # Set the aspect ratio
  geom_polygon(data = hull_df_global, aes(x = x, y = y), fill = "#CCCCCC30", color = NA) + # Add the global convex hull
  geom_segment(data = vector_data, aes(x = 0, y = 0, xend = Dim1 , yend = Dim2), arrow = arrow(length = unit(0.2, "cm")), color = "black") + # Add the vector direction or coordinates
  geom_label_repel(data= vector_data, aes(label = trait, x = Dim1, y = Dim2), box.padding = 0.5) + # Add the trait labels
  labs(title = "PCoA with Traits Vectors", x = "PCoA1", y = "PCoA2") + # Add the title and axis labels
  theme_classic()+# Apply the minimal theme
  theme(panel.border = element_rect(colour = "black", fill=NA)) # border of the plot

pcoa_vector_plot

# Save the plot in a SVG file and a jpeg file
ggsave("FES_Traits_Vectors.svg", pcoa_vector_plot, width = 15, height = 15, units = "cm")
ggsave("FES_Traits_Vectors.jpeg", pcoa_vector_plot, width = 15, height = 15, units = "cm")

##### 4.2. PCoA with Traits Categories ####
# We can plot the functional space with the categories positions for each trait (categorical)

# Extract the factor coordinates for categorical traits
factor_data <- scores(efit, "factor") %>%  # scores
  as.data.frame() %>% # transform to data frame
  mutate(trait =  substr(rownames(.), 1, nchar(rownames(.))-1)) # add the trait name but only the trait not the level (last letter)

# Create different plot with the functional space and the categories positions for each categorical trait
# To do this we can retrieve all the categorical traits and create a plot for each of them
Factorial_traits <- unique(factor_data$trait) # retrieve all the categorical traits

# Create a plot for each categorical trait
factorial_plots <- lapply(Factorial_traits, function(x) {
  data <- factor_data[factor_data$trait == x,] # Extract the factor coordinates for the specific trait
  ggplot(fd.coord, aes(x = PC1, y = PC2)) + # Create the plot based on the PCoA coordinates
    geom_polygon(data = hull_df_global, aes(x = x, y = y), fill = "#CCCCCC30", color = NA) + # Add the global convex hull
    geom_point(data = data, aes(x = Dim1, y = Dim2), color = "black") + # Add the factor coordinates
    geom_label_repel(data= data, aes(label = rownames(data), x = Dim1, y = Dim2), box.padding = 0.5) + # Add the trait labels
    labs(title = paste("Trait:", x), x = "PCoA1", y = "PCoA2") + # Add the title and axis labels
    theme_classic()+# Apply the minimal theme
    theme(panel.border = element_rect(colour = "black", fill=NA) ) # border of the plot
})

factorial_plots

# Display the plots all together in a grid with vector first and then categories
# Calculate the number of plots to draw 
n_plots <- length(Factorial_traits) + 1 # number of categorical traits + 1 for the vector plot

# SELECT MANUALY NUMBER OF LINES
plot_row <- 1 # set the number of lines and columns for the grid

# Draw the plots
do.call(grid.arrange, c(factorial_plots))

# Save the plot in a SVG file and a jpeg file
ggsave("FES_Traits_Categories.svg", do.call(grid.arrange, c(factorial_plots, nrow = plot_row)), width = 40, height = 10, units = "cm")
ggsave("FES_Traits_Categories.jpeg", do.call(grid.arrange, c(factorial_plots, nrow = plot_row)), width = 40, height = 10, units = "cm")


