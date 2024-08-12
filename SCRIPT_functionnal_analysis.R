# Author : Troussé Louis
# Created : 2024-07-30
# Affiliation : Septentrion Environnement

# This work is based on the following articles : 


# The following script is meant to be use to perform the different functionnal analysis on data about Abundances of differents species classifed in different groups based on their functional traits.
# To work properly, the script need to be in the same folder as the data file and this script and the one that contain the functions need to be in the same folder. 

#### Data Importation ####

#### Data PreProcessing ####
##### Functionnal Space ##### 
Functional_space <- create_functional_space(mat_funct = FES,
                                            nbdim = 11,
                                            metric = "Gower",
                                            dendro = FALSE,
                                            plot = "quality_funct_space",
                                            output_file = "FE_4D_coord.csv")

# here we choose to keep 4 dimensions : PLEASE ENTER 4


# We can now pass the output of the function to the workspace
fd.coord <- Functional_space$fd.coord
variance_explained <- Functional_space$variance_explained
n_dims <- Functional_space$n_dims
gower <- Functional_space$gower
fit <- Functional_space$fit