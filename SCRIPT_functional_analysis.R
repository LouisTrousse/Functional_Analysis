# Author : Troussé Louis
# Created : 2024-08-15
# Affiliation : Septentrion Environnement

#### Foreword ####
# This work is based on the work of the following authors (please cite the following articles if you use this script) :
# Gómez-Gras, D., Linares, C., Dornelas, M., Madin, J., Brambilla, V., Ledoux, J.-B., López-Sendino, P., Bensoussan, N., & Garrabou, J. (2021). Climate change transforms the functional identity of Mediterranean coralligenous assemblages. Ecology Letters. https://doi.org/10.1111/ele.13718
# Teixidó, N., Gambi, M. C., Parravacini, V., Kroeker, K., Micheli, F., Villéger, S., & Ballesteros, E. (2018). Functional biodiversity loss along natural CO2 gradients. Nature Communications, 9(1), 5149. https://doi.org/10.1038/s41467-018-07592-1

# The purpose of this script is to perform a functional analysis on a dataset containing informations about the abundances of different species classified in different groups based on their functional traits.
# To work properly, the script need to be in the same folder as the data file and this script and the one that contain the functions need to be in the same folder. 

# This work is a reproduction and an improvement of the script profided by Gomez-Gras et al. (2021) and Teixidó et al. (2018) to make it easier to use and to make it more flexible. All treatments are similar to thoses performed by the authors, but are here done with functions to make them easier to reproduce with different datasets. Details about the functions are here hidden to make the script easier to read. All details can be found in the script of the function itself.

#### Data Importation ####
# Please set the working directory to the folder containing your data files. Ensure that the data files are in the same folder as this script and the scripts containing the different functions.
# setwd("/YourFolder") # Please enter the path to your folder (if you are familiar with Rproject, no need to set the working directory, the script will find the data files in the same folder as the script). 

##### Load dataset ##### 
# Please enter your data file names
Functional_Entities <- read.csv2("./Functional_Entities.csv", sep=",", dec=".")# Load information about species and their functional traits
Abundance_data <- read.csv2("./Abundance_data.csv", sep=",", dec=".", row.names=1)# Load information about the abundances of the species in different samples (e.g. quadrats)
Sample_metadata <- read.csv2("./Sample_metadata.csv", sep=",", dec=".") # Load information about the samples (e.g. location, site, year, treatment, etc.)

##### Load Functions ####
# CAREFUL : The functions scripts need to be in the same folder as this script
source("functionalanalysis.R")

#### Data Preprocessing ####
##### Data Filtering #####
# These steps are not mandatory but by filtering the data you can choose the site you want to analyse, the species you want to keep, etc.

# Filter the data based on the metadata
# Filter Sample Metadata to keep only Pallazu and Pallazinu 
Sample_metadata <- Sample_metadata[Sample_metadata$Site %in% c("Pzzu_cor","Pzzu_par","Pzzinu_par"),]
# Filter Abundances data set to keep Pallazinu and Pallazu samples
Abundance_data <- Abundance_data[rownames(Abundance_data) %in% Sample_metadata$Quadrat,] %>%
  select(which(colSums(Abundance_data) > 0)) #Filter to keep only present taxa
# Keep present Functional Entities
Functional_Entities <- Functional_Entities[Functional_Entities$Species %in% colnames(Abundance_data),]

##### Functional Space ##### 
# CAREFUL : This function need your intervention to choose the number of dimensions to keep. Please enter the number of dimensions you want to keep in the console when asked. 
# The mean squared_deviation index is display and will help you to choose the number of dimensions to keep.
Functional_space <- create_functional_space(mat_funct = Functional_Entities,
                                            nbdim = 11,
                                            metric = "Gower",
                                            dendro = FALSE,
                                            plot = "quality_funct_space",
                                            output_file = "FE_coord.csv")

# Pass the output of the function to the workspace
fd.coord <- Functional_space$fd.coord
explained_variance <- Functional_space$explained_variance
n_dims <- Functional_space$n_dims
gower <- Functional_space$gower
fit <- Functional_space$fit

#### Data Analysis ####
##### 1. Define your argments ##### 
# In order to perform the analysis, you need to define the arguments of the function based on your datasets. 

###### Conditions ######
# Condition to compare among samples sites.
condition_column <- "Year" # Please enter the name of the column in the Sample_metadata file that contain the condition you want to compare (e.g. Year, Site, Treatment, etc.)

###### Visual parameters ######
# Layout of the plots : Compare different conditions in columns or rows
layout <- "column" # Choose your layout : "column" or "row" as needed (visualisation of different conditions in columns or rows)

# Colors for the different conditions
colors <- c("Pzzu_cor_2003" = "cadetblue",
            "Pzzu_cor_2011" = "#F9D71C",
            "Pzzu_cor_2018" = "coral",
            "Pzzu_par_2006" = "cadetblue",
            "Pzzu_par_2011" = "#F9D71C",
            "Pzzu_par_2018" = "coral",
            "Gabin_par_1999" = "cadetblue",
            "Gabin_par_2007" = "#F9D71C",
            "Gabin_par_2009" ="coral",
            "Pzzinu_par_2006" = "cadetblue",
            "Pzzinu_par_2011" = "#F9D71C",
            "Pzzinu_par_2016" ="coral",
            "Passe_cor_2006" = "cadetblue",
            "Passe_cor_2011" = "#F9D71C",
            "Passe_cor_2018"= "coral")
edges_colors <- c("Pzzu_cor_2003" = "cadetblue",
                  "Pzzu_cor_2011" = "#F9D71C",
                  "Pzzu_cor_2018" = "coral",
                  "Pzzu_par_2006" = "cadetblue",
                  "Pzzu_par_2011" = "#F9D71C",
                  "Pzzu_par_2018" = "coral",
                  "Gabin_par_1999" = "cadetblue",
                  "Gabin_par_2007" = "#F9D71C",
                  "Gabin_par_2009" ="coral",
                  "Pzzinu_par_2006" = "cadetblue",
                  "Pzzinu_par_2011" = "#F9D71C",
                  "Pzzinu_par_2016" ="coral",
                  "Passe_cor_2006" = "cadetblue",
                  "Passe_cor_2011" = "#F9D71C",
                  "Passe_cor_2018"= "coral")


# Number of clusters to compute during broad clustering step
# This number is the maximum number of functional clusters that you want to compute. Function will compute the optimal number of clusters based on the data.
Cluster_limit <- 20

# Colors for the different clusters. If you have more than 20 clusters, please add more colors. If less, only the first clusters colors will be used.
cluster_colors <- c("1"= "black", 
                    "2" = "red",
                    "3" = "green", 
                    "4" = "yellow", 
                    "5" = "#66FFFF",
                    "6" = "cadetblue",
                    "7" = "orange", 
                    "8" = "purple", 
                    "9" = "brown",
                    "10" = "blue",
                    "11" = "pink",
                    "12" = "grey",
                    "13" = "darkgreen",
                    "14" = "darkred",
                    "15" = "darkblue",
                    "16" = "darkorange",
                    "17" = "orange4",
                    "18" = "grey20",
                    "19" = "magenta", 
                    "20" = "coral")


##### 2. Run functional analysis #####
# In this part, we will run the different functions to analyse the data and display the results. To conduct the analysis, please run the following code. 
# If needed, some parameters can be changed in this part of the script, but the previous part is the one where you should enter the main parameters. Some export parameters are set to save the results in .svg and .jpeg format, please change the format and the dimensions of the plots as needed.

###### Functional Richness in Functional Space ###### 
# Launch function
Frich_in_Functional_space <- frich_in_functional_space(Site_Data = Sample_metadata,
                                                         Functional_diversity_coord = fd.coord,
                                                         Species_Functional_Entities = Functional_Entities,
                                                         Abundance_matrix = Abundance_data,
                                                         condition_column = condition_column,
                                                         colors = colors,
                                                         edges_colors = edges_colors, 
                                                         layout = layout, 
                                                         compute_null_hypothesis = TRUE,
                                                         n_perm = 100)

# Pass results to workspace
Frich_in_space <- Frich_in_Functional_space[[1]]
Null_Frich <- Frich_in_Functional_space[[2]]

# Save results 
# In .SVG and .jpeg format (choose the format you want)
# Dimensions are set for the example, please specify dimension for your own dataset
# For Functional richness in space
ggsave("Frich_Plots.svg", Frich_in_space, width = 40, height = 20, units = "cm")
ggsave("Frich_Plots.jpeg", Frich_in_space, width = 40, height = 20, units = "cm")
# For Functional richness under null hypothesis
ggsave("Frich_Plots_Sites.svg", Null_Frich, width = 40, height = 5, units = "cm")
ggsave("Frich_Plots_Sites.jpeg", Null_Frich, width = 40, height = 5, units = "cm")

###### Traits in functional Space (vector and organisation) ######
# Launch function
space_traits <- space_traits(Functional_diversity_coord = fd.coord,
                             Species_functional_traits = Functional_Entities,
                             gower = gower,
                             fit = fit)

# Pass results to workspace
pcoa_vector_plot <- space_traits$pcoa_vector_plot
Factorial_plots <- space_traits$Factorial_plots

# Save results
# In .SVG and .jpeg format (choose the format you want)
# Dimensions are set for the example, please specify dimension for your own dataset

# For the vector direction and traits in functional space
ggsave("FES_Traits_Vectors.svg", pcoa_vector_plot, width = 20, height = 20, units = "cm")
ggsave("FES_Traits_Vectors.jpeg", pcoa_vector_plot, width = 20, height = 20, units = "cm")
# For the categories positions in functional space
ggsave("FES_Traits_Categories.svg", Factorial_plots, width = 40, height = 20, units = "cm")
ggsave("FES_Traits_Categories.jpeg", Factorial_plots , width = 40, height = 20, units = "cm")


###### Trait distribution in Functional space ######
# Launch function
traits_distribution <- traits_distribution(Site_Data = Sample_metadata,
                                   condition_column = condition_column,
                                   Abundance_matrix = Abundance_data,
                                   colors = colors,
                                   edges_colors = colors,
                                   Species_Functional_Entities = Functional_Entities,
                                   Functional_diversity_coord = fd.coord,
                                   threshold = 0, 
                                   layout = layout,
                                   all_in_one = TRUE, 
                                   PERMANOVA = TRUE, 
                                   n_perm = 9999,
                                   n_dim = 4, 
                                   trait_relative_abundances = TRUE)

# Pass results to workspace
FI_plot <- traits_distribution$FI_plot
Plot_FI_all_years <- traits_distribution$Plot_FI_all_years
PERMANOVA_test_df <- traits_distribution$PERMANOVA_test_df
Trait_abundances <- traits_distribution$Trait_abundances

# Save results 
# In .SVG and .jpeg format (choose the format you want)
# Dimensions are set for the example, please specify dimension for your own dataset
ggsave("FI_plot.svg", plot = FI_plot, device = "svg", width = 40, height = 20, units = "cm")
ggsave("FI_plot.jpeg", plot = FI_plot, device = "jpeg", width = 40, height = 20, units = "cm")

ggsave("Plot_FI_all_years.svg", plot = Plot_FI_all_years, device = "svg", width = 40, height = 10, units = "cm")
ggsave("Plot_FI_all_years.jpeg", plot = Plot_FI_all_years, device = "jpeg", width = 40, height = 10, units = "cm")

#install.packages('writexl') # Uncomment if you don't have the package installed
require('writexl')
write_xlsx(PERMANOVA_test_df,"PERMANOVA_results.xlsx")

# The last element is a list of plots, that we can save using 

lapply(1:length(Trait_abundances), function(x){
  title <- names(Trait_abundances[x])
  plot <- Trait_abundances[[x]]
  ggsave(paste0(title,"_Trait_abundances.svg"), plot, device = "svg", width = 40, height = 30, units = "cm")
  ggsave(paste0(title,"_Trait_abundances.jpeg"), plot, device = "jpeg", width = 40, height = 30, units = "cm")
})

###### Functional traits clustering ######
# Launch function
Functional_clusters <- functional_clustering(gower = gower,
                                               Cluster_limit = Cluster_limit,
                                               Functional_diversity_coord = fd.coord,
                                               Species_Functional_Entities = Functional_Entities,
                                               Abundance_matrix = Abundance_data,
                                               Site_Data = Sample_metadata,
                                               Colors = cluster_colors, 
                                               Title ="Global")

# Pass results to workspace
Cluster_plot <- Functional_clusters$Cluster_plot
Sites_functional_cluster <- Functional_clusters$Sites_functional_cluster
Sites_functional_Abundance <- Functional_clusters$Sites_functional_Abundance
plot(Sites_functional_Abundance)
Sites_Functional_Clusters_abundances <- Functional_clusters$Sites_Functional_Clusters_abundances
Tab_Site_cluster_abundances <- Functional_clusters$Tab_Site_cluster_abundances

# Save results
# In .SVG and .jpeg format (choose the format you want)
# Dimensions are set for the example, please specify dimension for your own dataset

# Global clustering 
ggsave("Cluster_plot.svg", plot = Cluster_plot, device = "svg", width = 20, height = 20, units = "cm")
ggsave("Cluster_plot.jpeg", plot = Cluster_plot, device = "jpeg", width = 20, height = 20, units = "cm")

# Clustering by Site
ggsave("Sites_functional_cluster.svg", plot = Sites_functional_cluster, device = "svg", width = 40, height = 10, units = "cm")
ggsave("Sites_functional_cluster.jpeg", plot = Sites_functional_cluster, device = "jpeg", width = 40, height = 10, units = "cm")

# Number of species per cluster in each site 
ggsave("Sites_functional_Abundance.svg", plot = Sites_functional_Abundance, device = "svg", width = 40, height = 10, units = "cm")
ggsave("Sites_functional_Abundance.jpeg", plot = Sites_functional_Abundance, device = "jpeg", width = 40, height = 10, units = "cm")

# Abundances of each cluster
ggsave("Sites_Functional_Clusters_abundances.svg", plot = Sites_Functional_Clusters_abundances, device = "svg", width = 40, height = 10, units = "cm")
ggsave("Sites_Functional_Clusters_abundances.jpeg", plot = Sites_Functional_Clusters_abundances, device = "jpeg", width = 40, height = 10, units = "cm")

# Table of abundances of each cluster
write_xlsx(Tab_Site_cluster_abundances,"Tab_Site_cluster_abundances.xlsx")