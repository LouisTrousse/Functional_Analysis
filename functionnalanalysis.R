# This file contain functions to help to conduct a Functionnal Analysis
# Theses functions can be source in your dataset to be used in your analysis

#### GENERAL FUNCTIONS ####
# These functions are not meant to be used directly but are used by the other functions

##### Assign Colors to Clusters #####
assign_colors <- function(dataset, cluster_column, colors) {
  ### This function assigns colors to clusters based on the provided colors ###
  #Load necessary packages 
  require(tidyverse)
  
  # Check if the dataset is correctly provided
  if (!"data.frame" %in% class(dataset)) {
    stop("The provided dataset is not a data frame.")
  }
  # Check if the cluster columns is correctly provided
  if (!cluster_column %in% names(dataset)) {
    stop(paste("The specified cluster column", cluster_column, "does not exist in the dataset."))
  }
  
  # Determine the number of clusters
  n_cluster <- length(unique(dataset[[cluster_column]]))
  
  # If the number of clusters is equal to the number of cluster colors, assign the colors to the clusters
  if (n_cluster == length(colors)) {
    dataset$color <- colors[match(dataset[[cluster_column]], names(colors))]
  } else {
    # Print a warning message
    warning("The number of clusters is higher than the number of colors specified. A default color palette has been used.")
    # Create a color palette based on the default rainbow and assign it to the clusters
    default_colors <- rainbow(n_cluster)
    dataset$color <- default_colors[match(dataset[[cluster_column]], as.character(1:n_cluster))]
  }
  # Determine the number of clusters
  # Determine the unique clusters
  unique_clusters <- unique(dataset[[cluster_column]])
  n_cluster <- length(unique_clusters)
  
  # If the number of clusters is equal to the number of cluster colors, assign the colors to the clusters
  if (n_cluster == length(colors)) {
    dataset$color <- colors[match(dataset[[cluster_column]], names(colors))]
  } else {
    # Print a warning message
    warning("The number of clusters is different of the number of colors specified. A default color palette has been used.")
    # Create a color palette based on the default rainbow and assign it to the clusters
    default_colors <- SetNames(rainbow(n_cluster), unique_clusters)
    # Assign the colors to the clusters based on the cluster column and the poistion of each unique cluster in the list of clusters, to enshure to associate a color if the cluster is not a number 
    dataset$color <- default_colors[dataset[[cluster_column]]]
  }
  return(dataset) # Return the dataset with the assigned colors
}

##### Extract Year from character string #####
extract_year <- function(x) {
  ### Function to extract the year from a string that contains a 4-digit year in a site name that contains the year and other characters as in the example "Site1_2010" ###
  
  # load necessary packages 
  require(tidyverse)
  
  # Regular expression to match a 4-digit year in a string
  match <- regmatches(x, gregexpr("\\d{4}", x))
  
  if (length(match[[1]]) > 0) {
    return(match[[1]][1])  # Return the first 4-digit number found
  } else {
    return(NA)  # Return NA if no 4-digit number is found
  }
}

##### Display Multiple plots #####
display_multiple_plots <- function(plots,layout = "column"){
  ### Function to display multiple plots in a single plot with one line or column per group of plots merged by their common factor (e.g. Site) ###
  
  # Charge required libraries
  require(gridExtra)
  require(tidyverse)
  
  # Retrieve Sites in the name of each condition (condition less the year)
  Sites <- unique(sapply(names(plots), function(x) {
    year <- extract_year(x)
    
    if (!is.na(year)) {
      # Remove the year from the string
      site <- sub(year, "", x)
      # Replace whitespace with underscore
      site <- gsub("\\s", "_", site)
      # Suppress the last underscore if it exists
      site <- gsub("_$", "", site)
      return(site)
    } else {
      # If no year is found, replace whitespace with underscore
      site<-(gsub("\\s", "_", x))
      # Suppress the last underscore if it exists
      return(site)
    }
  }))
  
  # Associate each condition to a site from the names of the plot 
  sites <- sapply(names(plots), function(x) {
    year <- extract_year(x)
    
    if (!is.na(year)) {
      # Remove the year from the string
      site <- sub(year, "", x)
      # Replace whitespace with underscore
      site <- gsub("\\s", "_", site)
      # Suppress the last underscore if it exists
      site <- gsub("_$", "", site)
      return(site)
    } else {
      # If no year is found, replace whitespace with underscore
      site<-(gsub("\\s", "_", x))
      # Suppress the last underscore if it exists
      site <- gsub("_$", "", site)
      return(site)
    }
  })
  
  # create a list of plots for each Site
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
      # add any empty plots that needed 
      #retrieve the number of empty plots to add 
      n_empty_plots <- max_plots - length(x)
      
      # create the number of empty plots needed 
      empty_plots <- lapply(1:n_empty_plots, function(x) {
        ggplot(data = data.frame(x = 1, y = 1), aes(x = x, y = y)) +
          geom_point(aes(x = x, y = y), pch = 16, col = "white", cex = 4) + # plot empty plot
          theme_void() # remove all theme
      })#eo lapply
      
      # add this empty plots into the list
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
  # draw all plot in one (don't forget to charge gridextra)
  if (tolower(layout) == "column") {
    draw_all_plot <- do.call(grid.arrange, c(plots_drawing, ncol = nSites)) # addapt to your dataset
  } else if (tolower(layout) == "row") {
    draw_all_plot <- do.call(grid.arrange, c(plots_drawing, nrow = nSites)) # addapt to your dataset
  } else {
    print("error, layout should be 'row' or 'column'")
  }
  
  return(draw_all_plot)
}

##### Round Up any number #####
round_up <- function(x, digits = 0) {
  ### Round up to the nearest multiple of n ###
  n <- 10^digits # 10^0 = 1, 10^1 = 10, 10^2 = 100, etc.
  return (ceiling(x * n) / n) # Round up to the nearest upper integer and divide by n. ex 3.14 * 10 = 31.4, round up to 32 by ceiling, and 32/10 = 3.2
}

##### Create Convex Hull for multiple conditions #####
create_convex_hull_plot <- function(Fonctional_diversity_coord, binded, n_cluster, Colors = NULL, Title = "Global") {
  ### Function to create a convex hull plot for functional diversity ###
  
  # Charge required packages 
  require(tripack)
  require(tidyverse)
  
  #change the name of the functionnal diversity coordinates
  fd.coord = Fonctional_diversity_coord
  
  # Load necessary libraries
  require(ggplot2)
  require(tripack) # for tri.mesh and convex.hull
  require(dplyr)
  
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
  
  # add colors to the binded dataset to be abble to retrive them 
  binded <- assign_colors(dataset = binded, 
                          cluster_column = "Cluster", 
                          colors = Colors)
  
  # Default color for global convex hull
  global_hull_color <- "#CCCCCC30"
  
  
  # add the color to the named vector Colors 
  Colors <- c(Colors, "global" = global_hull_color)
  
  # Assign colors to clusters
  all_hulls_df <- assign_colors(dataset = all_hulls_df,
                                cluster_column = "Cluster" , 
                                colors = Colors)
  
  # Plot using ggplot2
  ggplot(binded, aes(x = PC1, y = PC2)) +
    geom_polygon(data = all_hulls_df, aes(x = x, y = y, group = Cluster, fill = color), alpha = 0.7) +
    geom_point(data = all_hulls_df, aes(x = x, y = y, color = color, fill = color)) +
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

#### FUNCTIONAL ANALYSIS FUNCTIONS ####
# These functions are meant to be used directly in your analysis and allow you to conduct analysis as in the following papers: 

##### Function to create a functional space #####
create_functional_space <- function(mat_funct, traits_weights = NULL, nbdim = 11, metric = "Gower", dendro = FALSE, plot = "quality_funct_space", output_file = "FE_4D_coord.csv") {
  ### Function to create a functional space based on functional entities and their traits ###
  
  # Load necessary packages
  require(tidyverse)
  
  # Load "Quality functionnal space" function
  source("./Raw_files/quality_funct_space.R")  # Change to the current function directory
  
  #retrieve mat_funct as required by the function quality_funct_space
  mat_funct <- mat_funct %>% select(-Species)%>%
    distinct() %>%
    column_to_rownames(var="FE") # If any duplicate FE, an error message will be printed
  
  # Calculate Quality Space
  qfs <- quality_funct_space(mat_funct = mat_funct,
                             traits_weights = traits_weights,
                             nbdim = nbdim,
                             metric = metric,
                             dendro = dendro,
                             plot = plot)
  
  # Verify space quality
  mean_sd <- round(qfs$meanSD, 4)
  cat("Mean SD of the quality space:", mean_sd, "\n")
  
  # Ask user for the number of dimensions to keep
  n_dims <- as.integer(readline(prompt = "Enter the number of dimensions to keep: "))
  
  # Ensure the user input is valid
  if (is.na(n_dims) || n_dims <= 0 || n_dims > nbdim) {
    stop("Invalid number of dimensions entered. Please enter a positive integer up to the number of calculated dimensions.")
  }
  
  # Keep Functional Entities coordinates on the specified dimensions
  fd.coord <- qfs$details_funct_space$mat_coord[, 1:n_dims]
  
  # Write result for further analyses
  write.csv(fd.coord, file = output_file)
  cat("Functional entities coordinates written to", output_file, "\n")
  
  # See variance explained by the PCoA axes
  gower <- qfs$details_funct_space$mat_dissim  # Store matrix
  fit <- cmdscale(gower, eig = TRUE, k = n_dims)  # PCoA
  
  # Variance explained by the axes
  variance_explained <- cumsum(fit$eig[fit$eig >= 0]) / sum(fit$eig[fit$eig > 0])
  cat("Variance explained by the axes:\n")
  print(variance_explained)
  
  invisible(list(fd.coord = fd.coord, variance_explained = variance_explained, n_dims = n_dims, gower = gower, fit = fit))
}

##### Functional Richness in Functional Space #####
frich_in_functionnal_space <- function(Site_Data, Abundance_matrix, Fonctional_diversity_coord, Species_Functional_Entities, condition_column, colors, edges_colors = colors, layout = "column") {
  ### Function to plot the functional richness of each condition of a dataset ###
  
  # load necessary packages
  require(tripack)
  require(tidyverse)
  
  #retrieve all the conditions in the dataset
  conditions <- unique(Site_Data[[condition_column]])
  
  # Assign colors to each condition for later plotting
  Site_Data <- assign_colors(dataset = Site_Data, cluster_column = condition_column, colors = colors)
  Site_Data2 <- assign_colors(dataset = Site_Data, cluster_column = condition_column, colors = edges_colors)
  
  # Create a named vector of colors for conditions
  condition_colors <- setNames(Site_Data$color, Site_Data[[condition_column]])
  edges_colors <- setNames(Site_Data2$color, Site_Data2[[condition_column]])
  
  # Rearrange data to retrieve abundances for each specific condition
  ab = Abundance_matrix
  ab.conditions <- lapply(conditions, function(x) {
    # Filter the quadrats for the specific condition
    if ("Quadrat" %in% colnames(Site_Data)){
      quad <- Site_Data$Quadrat[Site_Data[[condition_column]] == x]
    } else {
      quad <- rownames(Site_Data[Site_Data[[condition_column]] == x,])
    }
    
    # Filter the abundance matrix for the specific condition and calculate the sum of abundances because we are interrested in the presence of species for the condition
    colSums(ab[rownames(ab) %in% quad, ])
  })
  
  # Merge the list into a matrix
  ab.conditions <- do.call(rbind, ab.conditions)
  # Set the rownames to the conditions
  rownames(ab.conditions) <- conditions
  
  # Calculation of relative richness
  NbSp_tot <- length(unique(Species_Functional_Entities$Species)) # Number of species in the species_FE dataset
  
  # Calculate convex hull for the global functional space
  chg <- convhulln(Fonctional_diversity_coord, options = "FA") # Calculate convex hull for the global functional space
  
  Frich <- lapply(conditions, function(x) {
    # Identify species where abundance > 0
    species <- colnames(ab.conditions)[which(ab.conditions[x, ] > 0)]
    # Subset Species_Functional_Entities to only include FE of present species
    fes_cond <- Species_Functional_Entities %>% 
      subset(Species %in% species)
    
    # Subset Fonctional_diversity_coord to only include rows with names in fes_cond
    m <- Fonctional_diversity_coord[rownames(Fonctional_diversity_coord) %in% fes_cond$FE, ]
    
    # Compute the convex hull of m
    ch <- convhulln(m, options = "FA")
    # Return a vector with several calculated values
    c(length(species), length(species) / NbSp_tot * 100, dim(m)[1], dim(m)[1] / dim(Fonctional_diversity_coord)[1] * 100, ch$vol / chg$vol)
  })
  
  # Name the list with each condition
  names(Frich) <- conditions
  
  # Merge the list into a matrix
  Frich <- do.call(rbind, Frich)
  
  # Set the column names
  colnames(Frich) <- c("NbSp", "NbSpP", "NbFEs", "NbFEsP", "Frich (Vol4D)")
  
  # Create a data frame with the convex hull coordinates for the global space
  m2 <- Fonctional_diversity_coord
  tr2 <- tri.mesh(m2[, 1], m2[, 2])
  ch2 <- convex.hull(tr2)
  hull_df_global <- data.frame(x = ch2$x, y = ch2$y, i = ch2$i, group = "Global")
  
  # Apply a function to calculate the convex hull for each condition and create the corresponding plot
  Plots <- lapply(conditions, function(x) {
    # Identify species where abundance > 0
    species <- colnames(ab.conditions)[which(ab.conditions[x, ] > 0)]
    # Subset Species_Functional_Entities to only include FE of present species
    fes_cond <- Species_Functional_Entities %>% 
      subset(Species %in% species, )
    # Subset Fonctional_diversity_coord to only include rows with names in fes_cond
    m <- Fonctional_diversity_coord[rownames(Fonctional_diversity_coord) %in% fes_cond$FE, ]
    # Compute the convex hull of m
    tr <- tri.mesh(m[, 1], m[, 2])
    ch <- convex.hull(tr)
    # Create a data frame with the convex hull coordinates
    hull_df <- data.frame(x = ch$x, y = ch$y, i = ch$i, group = x)
    
    # Merge datasets into a single dataframe
    hull_df_combined <- bind_rows(hull_df_global, hull_df)
    
    # Create the plot
    Plot <- ggplot(hull_df_combined, aes(x = x, y = y)) +
      geom_polygon(data = hull_df_global, aes(group = group, fill = group, color = group), linewidth = 1, alpha = 0.5) + 
      geom_polygon(data = hull_df, aes(group = group, fill = group, color = group), linewidth = 1, alpha = 0.5) +
      labs(title = x, x = "PCoA 1", y = "PCoA 2") +
      theme_classic() +
      theme(panel.border = element_rect(colour = "black", fill = NA)) +
      theme(axis.text = element_text(face = "bold", size = 12)) +
      theme(axis.title = element_text(face = "bold")) +
      theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid")) +
      scale_fill_manual(values = condition_colors) +
      scale_color_manual(values = edges_colors) +
      theme(legend.position = "none") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
      annotate("text", x = 0.3, y = 0.3, label = paste("Frich = ", round(Frich[x, 5], 2)), size = 5, color = "black", fontface = "bold") +
      geom_point(data = data.frame(Fonctional_diversity_coord[rownames(Fonctional_diversity_coord) %in% species, ]), aes(x = Fonctional_diversity_coord[rownames(Fonctional_diversity_coord) %in% species, 1], y = Fonctional_diversity_coord[rownames(Fonctional_diversity_coord) %in% species, 2]), color = colors[x])
    
    return(Plot)
  })
  
  # Add the name of the condition to the list of plots
  names(Plots) <- conditions
  
  # Display the plots in a single plot using the display_multiple_plots function
  Frich_plots <- display_multiple_plots(plots = Plots, layout = layout)
  
  return(list(Frich_plots = Frich_plots, hull_df_global = hull_df_global, conditions = conditions, ab.conditions = ab.conditions))
}

##### Functional Richness against Null Hypothesis #####
plot_frich <- function(Site_Data, Abundance_matrix, Fonctional_diversity_coord, Species_Functional_Entities, condition_column, colors, edges_colors = colors, layout = "column", compute_null_hypothesis = FALSE, n_perm = 100) {
  ### Function to plot the functional richness of each condition of a dataset ###
  
  # load necessary packages
  require(tripack)
  require(tidyverse)
  require(gridExtra)
  
  #retrieve all the conditions in the dataset
  conditions <- unique(Site_Data[[condition_column]])
  
  # Assign colors to each condition for later plotting
  Site_Data <- assign_colors(dataset = Site_Data, cluster_column = condition_column, colors = colors)
  Site_Data2 <- assign_colors(dataset = Site_Data, cluster_column = condition_column, colors = edges_colors)
  
  # Create a named vector of colors for conditions
  condition_colors <- setNames(Site_Data$color, Site_Data[[condition_column]])
  edges_colors <- setNames(Site_Data2$color, Site_Data2[[condition_column]])
  
  # Rearrange data to retrieve abundances for each specific condition
  ab = Abundance_matrix
  ab.conditions <- lapply(conditions, function(x) {
    # Filter the quadrats for the specific condition
    if ("Quadrat" %in% colnames(Site_Data)){
      quad <- Site_Data$Quadrat[Site_Data[[condition_column]] == x]
    } else {
      quad <- rownames(Site_Data[Site_Data[[condition_column]] == x,])
    }
    # Filter the abundance matrix for the specific condition
    colSums(ab[rownames(ab) %in% quad, ])
  })
  
  # Merge the list into a matrix
  ab.conditions <- do.call(rbind, ab.conditions)
  # Set the rownames to the conditions
  rownames(ab.conditions) <- conditions
  
  # Calculation of relative richness
  NbSp_tot <- length(Species_Functional_Entities$Species) # Number of species in the species_FE dataset
  
  # Calculate convex hull for the global functional space
  chg <- convhulln(Fonctional_diversity_coord, options = "FA") # Calculate convex hull for the global functional space
  
  Frich <- lapply(conditions, function(x) {
    # Identify species where abundance > 0
    species <- colnames(ab.conditions)[which(ab.conditions[x, ] > 0)]
    # Subset Species_Functional_Entities to only include FE of present species
    fes_cond <- Species_Functional_Entities %>% 
      subset(Species %in% species, )
    # Subset Fonctional_diversity_coord to only include rows with names in fes_cond
    m <- Fonctional_diversity_coord[rownames(Fonctional_diversity_coord) %in% fes_cond$FE, ]
    # Compute the convex hull of m
    ch <- convhulln(m, options = "FA")
    # Return a vector with several calculated values
    c(length(species), length(species) / NbSp_tot * 100, dim(m)[1], dim(m)[1] / dim(Fonctional_diversity_coord)[1] * 100, ch$vol / chg$vol)
  })
  
  # Name the list with each condition
  names(Frich) <- conditions
  
  # Merge the list into a matrix
  Frich <- do.call(rbind, Frich)
  
  # Set the column names
  colnames(Frich) <- c("NbSp", "NbSpP", "NbFEs", "NbFEsP", "Frich (Vol4D)")
  
  # Create a data frame with the convex hull coordinates for the global space
  m2 <- Fonctional_diversity_coord
  tr2 <- tri.mesh(m2[, 1], m2[, 2])
  ch2 <- convex.hull(tr2)
  hull_df_global <- data.frame(x = ch2$x, y = ch2$y, i = ch2$i, group = "Global")
  
  # Apply a function to calculate the convex hull for each condition and create the corresponding plot
  Plots <- lapply(conditions, function(x) {
    # Identify species where abundance > 0
    species <- colnames(ab.conditions)[which(ab.conditions[x, ] > 0)]
    # Subset Species_Functional_Entities to only include FE of present species
    fes_cond <- Species_Functional_Entities %>% 
      subset(Species %in% species, )
    # Subset Fonctional_diversity_coord to only include rows with names in fes_cond
    m <- Fonctional_diversity_coord[rownames(Fonctional_diversity_coord) %in% fes_cond$FE, ]
    # Compute the convex hull of m
    tr <- tri.mesh(m[, 1], m[, 2])
    ch <- convex.hull(tr)
    # Create a data frame with the convex hull coordinates
    hull_df <- data.frame(x = ch$x, y = ch$y, i = ch$i, group = x)
    
    # Merge datasets into a single dataframe
    hull_df_combined <- bind_rows(hull_df_global, hull_df)
    
    # Create the plot
    Plot <- ggplot(hull_df_combined, aes(x = x, y = y)) +
      geom_polygon(data = hull_df_global, aes(group = group, fill = group, color = group), linewidth = 1, alpha = 0.5) + 
      geom_polygon(data = hull_df, aes(group = group, fill = group, color = group), linewidth = 1, alpha = 0.5) +
      labs(title = x, x = "PCoA 1", y = "PCoA 2") +
      theme_classic() +
      theme(panel.border = element_rect(colour = "black", fill = NA)) +
      theme(axis.text = element_text(face = "bold", size = 12)) +
      theme(axis.title = element_text(face = "bold")) +
      theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid")) +
      scale_fill_manual(values = condition_colors) +
      scale_color_manual(values = edges_colors) +
      theme(legend.position = "none") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
      annotate("text", x = 0.3, y = 0.3, label = paste("Frich = ", round(Frich[x, 5], 2)), size = 5, color = "black", fontface = "bold") +
      geom_point(data = data.frame(Fonctional_diversity_coord[rownames(Fonctional_diversity_coord) %in% species, ]), aes(x = Fonctional_diversity_coord[rownames(Fonctional_diversity_coord) %in% species, 1], y = Fonctional_diversity_coord[rownames(Fonctional_diversity_coord) %in% species, 2]), color = colors[x])
    
    return(Plot)
  })
  
  # Add the name of the condition to the list of plots
  names(Plots) <- conditions
  
  # Display the plots in a single plot using the display_multiple_plots function
  Frich_plots <- display_multiple_plots(plots = Plots, layout = layout)
  
  # If the null hypothesis is set to TRUE, run the null hypothesis test
  if (compute_null_hypothesis == TRUE){
    
    # Retrieve list of sites on which the analysis will be performed
    Sites = unique(Site_Data$Site)
    
    Sites_Frich <- lapply(Sites, function(x){
      # filter conditions for the site of interest
      Specific_Conditions <- conditions[grep(x, conditions)]
      # filter the abundance matrix for the site of interest
      ab_specific_conditions = ab[grep(x, rownames(ab)),]
      
      # Data arrangements to retrieve abundances for each specific conditions
      ab.specific_condition <- lapply(Specific_Conditions, function(x) {
        # filter the quadrats for the specific condition
        if ("Quadrat" %in% colnames(Site_Data)){
          quad <- Site_Data$Quadrat[Site_Data[[condition_column]] == x]
        } else {
          quad <- rownames(Site_Data[Site_Data[[condition_column]] == x,])
        } #eo if statement
        
        # filter the abundance matrix for the specific condition and calculate the sum of abundance
        colSums(ab_specific_conditions[rownames(ab_specific_conditions) %in% quad,])
      })#eo lapply
      
      #replace the rownames of the abundance matrix with the conditions
      ab.specific_condition <- do.call(rbind, ab.specific_condition) #merge the list into a matrix
      rownames(ab.specific_condition) = Specific_Conditions #set the rownames to the conditions
      
      #Calculate convex hull for each specific condtion (111 taxonomic units in total for calculating the relative richness)
      Frich <- lapply(Specific_Conditions, function (x) {
        # filter the species for the specific condition (species present in these samples)
        species <- colnames(ab.specific_condition)[which(ab.specific_condition[x,] > 0)]
        # filter the functional traits for this species
        fes_cond <- Species_Functional_Entities %>% 
          subset(Species %in% species, )
        # Subset Fonctional_diversity_coord to only include rows with names in fes_cond
        m <- Fonctional_diversity_coord[rownames(Fonctional_diversity_coord) %in% fes_cond$FE, ]
        # calculate the convex hull for the functional space
        ch <- convhulln(m, options = "FA")
        # return the number of species, the relative number of species, the number of FEs, the relative number of FEs, and the volume of the convex hull
        c(length(species), length(species)/NbSp_tot*100, dim(m)[1], dim(m)[1]/dim(Fonctional_diversity_coord)[1]*100, ch$vol/chg$vol*100)
      })#eo lapply Frich
      
      #name the list with the conditions
      names(Frich) = Specific_Conditions
      # merge the list into a matrix
      Frich <- do.call(rbind, Frich)
      # name the columns
      colnames(Frich) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol4D")
      
      # Definition of the null model 
      # null model is based on the random selection of species in the dataset
      
      # create a copy of FEs to be used in the permutations
      spe_fes_random = Species_Functional_Entities
      
      # Calculate permuted Frich 
      Frich_permuted <- lapply(Specific_Conditions, function (x) {
        # filter the species for the specific condition (species present in these samples)
        species <- colnames(ab.specific_condition)[which(ab.specific_condition[x,] > 0)]
        #create n_perm aleatory permutation of the functional traits for specified condition
        perm <- sapply((1:n_perm), function (i) {
          # shuffle the functional traits to associate them to random species (create a random functional space)
          spe_fes_random$FE <- sample(Species_Functional_Entities$FE)    
          # Filter the functional traits for the species present in specific condition
          random_fes_cond <- spe_fes_random %>% 
            subset(Species %in% species, )
          # Filter the functional space for these functional traits
          random_m <- fd.coord[rownames(fd.coord) %in% random_fes_cond,]
          # Calculate the convex hull for this random functional space
          random_ch <- convhulln(random_m, options = "FA")
          # Return the number of species, the relative number of species, the number of FEs, the relative number of FEs, and the volume of the convex hull
          c(length(species), length(species)/NbSp_tot*100, dim(random_m)[1], dim(random_m)[1]/dim(fd.coord)[1]*100, random_ch$vol/chg$vol*100)
        })#eo sapply
        # rename the rownames of the permuted Frich
        rownames(perm) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol")
        # Return perm
        return(perm)
      })#eo lapply Frich_permuted
      
      # Name the list with the conditions considered
      names(Frich_permuted) = Specific_Conditions
      
      Frich_perm_Quantiles <- lapply(Frich_permuted, function (x) {
        # Calculate the quantiles 0.05 and 0.95 of the permuted Frich
        rowQuantiles(x, probs=c(0.05, 0.95))
      })#eo lapply Frich_perm_Quantiles
      
      # retrieve actual Frich values of each condition
      Frich = as.data.frame(Frich)
      # Add the lower (0.05) and upper quantiles 0.95 to the Frich data frame
      Frich$lowerFE <- sapply(Specific_Conditions, function (x) { Frich_perm_Quantiles[[x]][3,1] })#eo sapply
      Frich$upperFE <- sapply(Specific_Conditions, function (x) { Frich_perm_Quantiles[[x]][3,2] })#eo sapply
      Frich$lowerVol <- sapply(Specific_Conditions, function (x) { Frich_perm_Quantiles[[x]][5,1] })#eo sapply
      Frich$upperVol <- sapply(Specific_Conditions, function (x) { Frich_perm_Quantiles[[x]][5,2] })#eo sapply
      
      # Add the year as a factor to the Frich data frame 
      Frich$year <- as.factor(sapply(Specific_Conditions, function(x) {extract_year(x)}))
      
      # set the column names
      colnames(Frich) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol4D", "lowerFE", "upperFE", "lowerVol", "upperVol", "year")
      
      # in the named vector, extract from each name the year 
      year_colors <- setNames(condition_colors,lapply(names(condition_colors), extract_year))
      # Plot the Frich values with the 95% confidence interval
      Plot_null <-ggplot(data = Frich, x = year , y = Vol4D)+
        geom_point(aes(x = year, y = Vol4D, color = year), pch = 16, cex = 4)+ # plot actual Frich values
        geom_linerange(aes(x = year, ymin = lowerVol, ymax = upperVol, color = year), linewidth = 2) + # plot 95% confidence interval
        scale_y_continuous(breaks=c(0,20,40,60,80,100), limits = c(0,100))+ # set y axis limits
        ylab("Relative richness (%)") + # y axis label
        xlab("Time")+ # x axis label
        ggtitle(x)+ # title
        coord_flip()+ # flip x and y axis
        scale_x_discrete(limits = as.character(rev(levels(Frich$year))))+
        theme_classic()+ # classic theme
        scale_color_manual(values = year_colors) +
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
      
      return(Plot_null)
    })#eo lapply Sites_Frich
    
    #Retrieve the number of Sites
    ncol = length(Sites)
    # display the plots
    Plots_Frich_null <- do.call(grid.arrange, c(Sites_Frich, ncol = ncol))
  } else {
    Plots_Frich_null = NA
  }
  return(list(Frich_plots, Plots_Frich_null))
} #eo Plot_functional_richness

##### Trait vectors direction and traits in functional Space #####
space_traits<- function(Fonctional_diversity_coord, Species_functionnal_traits, gower, fit) {
  ### Function to plot the traits (categorical and numerical in the functional space ###
  
  # load necessary packages
  require(tidyverse)
  require(ggrepel)
  require(gridExtra)
  require(vegan)
  
  #We use the cdmscale() and envfit() to re-ecreate our trait-space and find vector directions and categories positions in functional space
  
  # We first create the background polygon (global functional space)
  M_global <- fd.coord[rownames(fd.coord),] # select all coordinates of the functional entities (FE)
  TR_global <-tri.mesh(M_global[,1],M_global[,2]) # Create Delaunay triangulation
  CH_global <- convex.hull(TR_global) # Compute convex hull
  
  # Create a data frame with the convex hull coordinates
  hull_df_global <- data.frame(x = CH_global$x, # x coordinates
                               y = CH_global$y, # y coordinates
                               i = CH_global$i) # index of this point
  
  # For our variables, we calculate the vector direction in the functional space or the position of each categories in the functional space
  
  # We retrieve fes without species from the Species functionnal trait
  fes = Species_functionnal_traits %>% 
    select(-Species)%>% 
    distinct() %>% # remove duplicates
    column_to_rownames(var="FE") # If any duplicate FE, an error message will be printed
  
  efit <- envfit(fit, fes, na.rm=TRUE) # envfit fit the environmental variables to the functional space
  
  # TRAIT VECTOR PLOT
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
  
  # We can calculate the coefficient to adjust the length of the vectors based on the maximum distance from the origin to the maximum and minimum values of the PCoA axes. By doing this we can adjust the length of the vectors to fit the plot, as the longest vector will be equal to the maximum distance from the origin to the maximum and minimum values of the PCoA axes. 
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
    theme(panel.border = element_rect(colour = "black", fill=NA))+ # border of the plot
    theme(axis.text = element_text(face = "bold", size = 12))+ # bold axis text
    theme(axis.title = element_text(face = "bold", size = 12))+ # bold axis title
    theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid"))+
    theme(legend.position = "none")+ # remove legend
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +# centered bold title
    theme(axis.title.x = element_text(size=14))+ # set x axis title size
    theme(axis.title.y = element_text(size=14)) # set y axis title size
  
  # CATEGORY PLOT
  # We can also visualize the position of the categories in the functional space for Factorial traits
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
      theme(panel.border = element_rect(colour = "black", fill=NA)) +# border of the plot
      theme(axis.text = element_text(face = "bold", size = 12))+ # bold axis text
      theme(axis.title = element_text(face = "bold", size = 12))+ # bold axis title
      theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid"))+
      theme(legend.position = "none")+ # remove legend
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +# centered bold title
      theme(axis.title.x = element_text(size=14))+ # set x axis title size
      theme(axis.title.y = element_text(size=14)) # set y axis title size
  })
  
  # Draw the plots in one plot
  Factorial_plots <- do.call(grid.arrange, c(factorial_plots))
  
  return(list(pcoa_vector_plot = pcoa_vector_plot, Factorial_plots = Factorial_plots))
}

##### Trait distribution in functional space #####
abund_traits_distribution <- function (Site_Data, condition_column, Abundance_matrix, colors = "black", edges_colors = colors, Species_Functional_Entities, Fonctional_diversity_coord, threshold = 0, layout = "column", all_in_one = FALSE, PERMANOVA = TRUE, n_dim = 2, n_perm = 9999) {
  ### Function to plot the abundance and distribution of traits in the functional space ###
  
  # load necessary packages
  require(tripack)
  require(tidyverse)
  
  #retrieve all the conditions in the dataset
  conditions <- unique(Site_Data[[condition_column]])
  
  cols <- colors 
  colstr <- edges_colors
  # Assign colors to each condition for later plotting
  Site_Data <- assign_colors(dataset = Site_Data, cluster_column = condition_column, colors = cols)
  Site_Data2 <- assign_colors(dataset = Site_Data, cluster_column = condition_column, colors = colstr)
  
  # Create a named vector of colors for conditions
  cols <- setNames(Site_Data$color, Site_Data[[condition_column]])
  colstr <- setNames(Site_Data2$color, Site_Data2[[condition_column]])
  
  # Data manipulation and arrangements
  # Compute the abundance of the species in the different conditions (mean cover)
  ab = Abundance_matrix
  ab.conditions <- lapply(conditions, function(x) {
    # get the quadrats coresponding to the condition
    # if Site_Data have a column named Quadrat, use it to filter the quadrats
    if ("Quadrat" %in% colnames(Site_Data)){
      quad <- Site_Data$Quadrat[Site_Data[[condition_column]] == x]
    } else {
      quad <- rownames(Site_Data[Site_Data[[condition_column]] == x,])
    }
    
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
  
  # get the levels of the FEs : the different FEs
  fes <- levels(as.factor(Species_Functional_Entities$FE))
  
  # compute the abundance of each FEs in the different conditions
  ab.fe.conditions <- lapply(conditions, function (z) {
    # In each condition, compute the abundance of each FEs
    abund.fes <-  sapply(fes, function (x) {
      # get row names (Species) for the specific Fe
      species <- Species_Functional_Entities$Species[which(Species_Functional_Entities$FE == x)]
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
  
  ## Functional identity (fI) trends
  fd.coord.FE <- Fonctional_diversity_coord %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "FE") # add the FE as a column
  
  # merge the coordinates and the species FEs to get the coordinates of the species in the functional space
  FE_SP_coord <- merge(fd.coord.FE, Species_Functional_Entities %>%
                        select(Species, FE)# select only wished columns
                       , by = "FE") %>% 
    select(-FE) # remove FE column
  
  # Obtaining weighted centroids for each assemblage and Time point
  # The weighted centroid is the sum of the product of the coordinates of the species in the functional space by the weight of the species. This species weight is the abundance of the species in the assemblage. 
  # For each Condition, we calculate the weighted centroid of the assemblage in the functional space
  # To this we multiply the coordinates of the species in the functional space by the abundance of the species in the assemblage and sum the results. 
  # Then we average the results by the total abundance of the assemblage to get the weighted centroid of the assemblage in the functional space.
  
  # First we merge The weight dataset to the species coordinates dataset using bindrows. 
  # By doing this we will not have any doubt about the order of the species in the two datasets.
  Species_Weights_coord <- left_join(Species_Weights, FE_SP_coord, by = "Species")
  
  # Then for each condition we calculate weighted species coordinates in the functional space
  Weighted_centroid <- lapply(conditions, function (x) {
    # get the weighted coordinates of the species in the functional space for the condition z
    Condition_Weighted_coord <- Species_Weights_coord %>% 
      select(x, PC1, PC2) %>% # select the columns of interest
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
    # filter to keep only species that have an abundance > 0 in the condition
    # Set the threshold to the percentage of abundance you want to consider
    threshold <- threshold
    
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
      geom_point(data = fd.coord.temp, aes(x = PC1, y = PC2, size = ab.fe.conditions.temp$size), shape = 21, fill = as.data.frame(cols, row.names= names(cols))[x,], color = as.data.frame(edges_colors, row.names= names(edges_colors))[x,], alpha = 0.5)+ # plot the FE in the functional space
      scale_size_continuous(range = c(1, max(ab.fe.conditions.temp$size)*0.5)) + # set the range of the size of the points as a function of the abundance of the FEs (size of the points is proportional to the abundance of the FEs) 0.5 is a factor to adjust the size of the points
      # plot the centroid of the assemblage in the functional space
      annotate("text", x = Weighted_centroid[x,"PC1"], y = Weighted_centroid[x,"PC2"], label = "+", color = as.data.frame(cols, row.names = names(cols))[x,], size = 14)+ # plot the centroid of the assemblage in the functional space
      labs(title = x, x = "PCoA 1", y = "PCoA 2") + #axis label 
      theme_classic()+ # classic theme
      theme(panel.border = element_rect(colour = "black", fill=NA))+ # border of the plot
      theme(axis.text = element_text(face = "bold", size = 12))+ # bold axis text
      theme(axis.title = element_text(face = "bold"))+ # bold axis title
      theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid"))+
      theme(legend.position = "none")+ # remove legend
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))# centered bold
  }) #eo lapply
  
  names(FI_plots) = conditions #Add the condition name as row name
  
  # We use the previous code to plot the FI trends for all conditions
  FI_plot <- display_multiple_plots(FI_plots, layout = layout) # Display the plots in one
  
  # We can want to have all plots of one site but the differents years in the same plot 
  if (all_in_one == TRUE) { 
    # First we retrieve the sites
    unique_sites <- unique(Site_Data$Site)
    
    # For each site, we retrieve the conditions and plot the FI trends for each condition
    Plot_Sites_all_years <- lapply (unique_sites, function(x){
      # retrieve the years for the site
      years <- unique(Site_Data[Site_Data$Site == x,][[Condition_column]])
      years <- unique(Site_Data[Site_Data$Site == x,]$Year) ##OLD 
      
      # retieve abundance of the FE in the different conditions for the site and filter the FEs with an abundance > threshold
      threshold <- threshold
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
      
      # we define non conditional treatment part
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
          scale_color_manual(values = cols[Year]) + # set the color of the points as a function of the year
          scale_fill_manual(values = colstr[Year]) + # set the fill color of the points as a function of the year
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
          annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols, row.names=names(cols))[years,], size = 14)# plot the centroid of the assemblage in the functional space
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
          annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols, row.names=names(cols))[years,], size = 14)
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
          annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols, row.names=names(cols))[years,], size = 14)
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
          annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols, row.names=names(cols))[years,], size = 14)
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
          annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols, row.names=names(cols))[years,], size = 14)
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
          annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols, row.names=names(cols))[years,], size = 14)
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
          annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols, row.names=names(cols))[years,], size = 14)
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
          annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols, row.names=names(cols))[years,], size = 14)
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
          annotate("text", x = Weighted_centroid[years,"PC1"], y = Weighted_centroid[years,"PC2"], label = "+", color = as.data.frame(cols, row.names=names(cols))[years,], size = 14)
      } else {
        print("Too many years to plot, the maxium is 10 years, if you want to plot more years, please refer to the code to add more year layers")
      } #eo if else
    }) #eo lapply
    
    # add the site name to the list
    names(Plot_Sites_all_years) <- unique_sites
    
    Plot_FI_all_years <- do.call(grid.arrange, c(Plot_Sites_all_years, nrow = 1))
  } else {
    Plot_FI_all_years = NA
  }
  
  # Lauch a PERMANOVA to test the effect of the condition on the functional space
  if (PERMANOVA == TRUE) { 
    # We use the weighted centroids (based on weighted coordinates) of the community to conduct the PERMANOVA analysis
    # This PERMANOVA analysis will test the effect of the condition on the weighted centroids of the community 
    # As weighted centroids are based on the weighted coordinates of the species, we will use the Canberra distance to measure the distance between the weighted centroids
    # We calculate the weighted centroids of the community for each quadrat
    ab_P <- t(ab) # transpose the abundance matrix
    
    #apply for the number of dimensions specified for this PERMANOVA
    # Calculate the coordinates of each quadrats allong axis 
    Weighted_axis <- lapply(1:n_dims, function(x){
      weighted_axis = FE_SP_coord[,x] * ab_P
      Axis_sum <- colSums(weighted_axis)
      weighted_axis_coordinates <- Axis_sum/100
    })
    
    # Binding all axis and calculating the sum of the weighted coordinates
    #firstwe create a dataframe containing the first axis
    weighted_centroids <- as.data.frame(Weighted_axis[[1]])
    for (i in 2:n_dims){
      i = i
      weighted_centroids <- cbind(weighted_centroids, Weighted_axis[[i]]) %>% 
        #rename the variable with the axis number
        setNames(paste0("PC",1:i))
    }
    
    #We complete the data frame with more useful info for the analysis.
    #Site
    data <-  weighted_centroids %>%
      rownames_to_column(var = "Rownames") %>%
      separate(Rownames, c(condition_column, "site", "Type", "Quadrat"), sep = "_") %>%
      # merge collums Year and Site
      mutate(Site = paste(site,Type, sep = "_")) %>%
      select(-site, -Type, -Quadrat)
    
    ### We use this new data frame to conduct the PERMANOVA ANALYSIS 
    # we apply for each site the adonis2 function (PERMANOVA) to test changes on the weighted centroids of the community over the condition
    # we use the Canberra distance to measure the distance between the weighted centroids as the canberra distance is more appropriate for compositional data than the Euclidean distance see legendre et al. 2012
    # we use the condition as a factor to test the effect of the condition on the weighted centroids of the community
    
    PERMANOVA_test <- lapply(unique(data$Site), function(x){
      data_site <- data[data$Site == x,]
      paste(condition_column)
      data_site_num <- data_site %>% 
        select(-Site, - paste(condition_column))
      per_site <- adonis2(data_site_num ~ data_site[,condition_column], method="canberra", permutations = n_perm)
      per_site
    })
    
    # we can rename the results with the site names
    names(PERMANOVA_test) <- unique(data$Site)
    
    # We can merge the results in a data frame and adding an empty column to add the site names
    PERMANOVA_test_df <- do.call(rbind, PERMANOVA_test) %>%
      as.data.frame() %>% 
      rownames_to_column(var = "Rownames") %>% 
      as_tibble() %>% 
      separate_wider_delim(cols = Rownames,names= c("Site","Attr"), delim = ".") %>% 
      mutate(Attr = ifelse(Attr == "data_site[, condition_column]", paste(condition_column), Attr))# Modify Data_site[, condition_column] into paste(condition_column)
  } else {
    PERMANOVA_test_df = NA # Define an empty result
  }
  
  # _______________________________________
  
  # We Can plot changes in relative abundances of functional trait categories face to the chosen condition 
  # This part is based on the work of Nuria Teixidó et al., 2018 : Functional biodiversity loss along natural CO2 (DOI: 10.1038/s41467-018-07592-1)
  
  #First we link the traits to the abundance matrix
  # We merge species functionnal entities, their character detail and their abundances in a long format data frame
  if 
  
  
  
  
  
  
  Species_traits_abund <- merge(Species_Functional_Entities %>% 
                                  pivot_longer(cols = -c(FE, Species), names_to = "Trait", values_to = "Value",values_transform = as.character),
                                Species_Weights %>% 
                                  pivot_longer(cols = - Species , names_to = "Condition", values_to = "Abundance", values_transform = as.double), 
                                by = "Species") %>% 
    filter(Abundance > 0) %>% # we filter the species with an abundance higher than 0
    left_join(Site_Data %>% 
                rownames_to_column(var= "Quadrat")%>% 
                select(-Quadrat) %>% 
                distinct(), join_by("Condition"== "Year"))# we merge with the Sites data to have the condition of each site
  
  # We calculate for each site the relative abundance of each functional trait category
  Site_trait_cat_abundances <- Species_traits_abund %>% 
    group_by(Site, Condition, Trait) %>% 
    mutate(Abundance_tot = sum(Abundance)) %>%
    ungroup() %>% 
    mutate(Relativ_Abundance = sum(Abundance)/Abundance_tot)
    summarise(Relativ_Abundance = sum(Abundance)/Abundance_tot)
  
  view(Site_trait_cat_abundances)
  
  
  
  # We retrieve the unique sites
  Sites <- unique(Site_Data$Site)

  # We calculates the abundances
  Site_trait_cat_abundances <- lapply(Sites, function (site){ 
    
    # We filter the data for the site (condition that contain the name of the site)
    data_site <- Species_traits_abund %>% 
      filter(Site == site)
    
    # Then, we iterate through the traits to calculate the relative abundance of each trait
    trait_cat_abundances <- lapply(unique(data_site$Trait), function(trait){
      
      # We filter the data for the trait
      data_trait <- data_site %>% 
        filter(Trait == trait)
      
      # We calculate the relative abundance of the trait category
      cat_abundances <- data_trait %>% 
        summarise(Relative_abundance = sum(Abundance)/sum(data_site$Abundance))
      
      # We add the trait name to the data frame
      trait_abundance$Trait <- trait
      
      return(trait_abundance)
    })
      
    })
  
  return(list(FI_plot = FI_plot ,Plot_FI_all_years = Plot_FI_all_years , PERMANOVA_test_df = PERMANOVA_test_df))
  
} # eo function

##### Functionnal traits clustering #####
functionnal_clustering <- function(gower, Cluster_limit = 20, Fonctional_diversity_coord,Species_Functional_Entities, Site_Data, Colors = NULL, Title ="Broad clustering"){
  ### This function compute: 
  # Broad functional classification into functional clusters (with the overall pool).
  # Broad functional classification into functional clusters (per site).
  # Number of species per functional cluster in each site and Time point.
  # Temporal abundance changes in functional clusters per site. ###
  
   # load the required libraries
  library(cluster)
  library(tidyverse)
  library(gridExtra)
  library(geometry)
  
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
  
  # We use a conditionnal loop to enshure that the limit set by the user is not exceeded the number of clusters possible 
  if (n_FE > Cluster_limit) {
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
    cat("The optimal number of clusters is",n_cluster,"\n with a silhouette width of",max_sil_width[2], "\n")
  } else {
    cat("A limit of",Cluster_limit,"culsters is higher than the number of possible clusters.\nPlease set a limit beetwen 2 and ",n_FE-1, "\n")
    stop()
  }
  
  # Create the clusters using PAM for n clusters
  pam_fit = pam(gower, diss=TRUE, k= n_cluster)
  
  # We add the "cluster N" column to our data with coordinates. This will allow us to know to which cluster has been assigned to every FE. 
  # We will use these data later to calculate the N of sp per cluster (functional redundancy) and if there are temporal changes in the relative abundance of functional clusters in the different assemblages
  pam_fit$clustering = as.data.frame(pam_fit$clustering)
  binded <- cbind(pam_fit$clustering, fd.coord)
  colnames(binded)[which(names(binded) == "pam_fit$clustering")] <- "Cluster"
  
  # Plotting
  # For each cluster we create the corresponding convex hull
  Cluster_plot <- create_convex_hull_plot(Fonctional_diversity_coord = Fonctional_diversity_coord, n_cluster = n_cluster, binded = binded, Colors = Colors, Title = Title)
  
  # We retrieve each color of the clusters 
  colors = setNames(unique(Cluster_plot$data$color), unique(Cluster_plot$data$Cluster))
  
  #___________________________________
  # We can do the same for each site
  # First we retrieve the number of different Site 
  Sites = unique(Site_Data$Site)
  
  # Then we create a binded data frame for all site with information of abundance of FE that we can filter to retrieve only FE present in a given site
  # We retrieve abundance of FE for each site
  # get the levels of the FEs : the different FEs
  fes <- levels(as.factor(Species_Functional_Entities$FE))
  
  # compute the abundance of each FEs in the different conditions
  ab.fe.conditions <- lapply(conditions, function (z) {
    # In each condition, compute the abundance of each FEs
    abund.fes <-  sapply(fes, function (x) {
      # get row names (Species) for the specific Fe
      species <- rownames(Species_Functional_Entities)[which(Species_Functional_Entities$FE == x)]
      # Calculate the abundance of the FEs in fonction of the abundance of the corresponding species
      sum(ab.conditions[z,species])
    })#eo sapply
    # return the abundance of the FEs in the condition
  })#eo lapply
  
  # Add the condition name as row name
  names(ab.fe.conditions) = conditions
  
  # Transform the lists into matrix
  ab.fe.conditions <- do.call(rbind, ab.fe.conditions)
  
  #Transpose it
  fe.conditions = t(ab.fe.conditions)
  
  # we merge it to the binded data frame
  fe_binded<-cbind(binded,fe.conditions)
  
  # We identify the number of axis PC
  n_axis <- dim(fd.coord)[2]
  #and create a lsit of PC in function of this 
  PCs <- paste0("PC",1:n_axis)
  
  # We transform it into a tidy format
  fe_binded_tidy <- fe_binded %>%
    rownames_to_column(var = "FE") %>%
    gather(key = "Year", value = "Abundance", -Cluster, -all_of(PCs), -FE) # Regroup all the years in one column
  
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
      distinct(FE, .keep_all = TRUE)# keep only one row per FE
    
    # We filter the list of colors with only colors in the Site_binded$cluster
    colors <- colors[as.character(unique(Site_binded$Cluster))]
    
    # We apply the same method as before to create the clusters
    create_convex_hull_plot(Fonctional_diversity_coord = fd.coord, binded = Site_binded, n_cluster = n_cluster, Colors = colors, Title = site)
  }) 
  
  # Add a title to each plot
  Site_Functionnal_Clusters <- lapply(1:length(Site_Functionnal_Clusters), function(i) {
    Site_Functionnal_Clusters[[i]] + ggtitle(paste0(Sites[i]))
  })
  
  # For each element of the list, we add the Site as entry
  names(Site_Functionnal_Clusters) <- Sites
  
  # We draw the plot in one using the gridExtra package
  n_col = length(Site_Functionnal_Clusters) # number of column = number of year
  
  require(gridExtra)
  Sites_functionnal_cluster <- grid.arrange(grobs = Site_Functionnal_Clusters, ncol = n_col)
  
  #_____________________________
  #Retrieve the number of species per cluster in each site and Time point
  
  # We retrieve abundances of each species for each site
  # get the levels of the species : the different species
  sp <- levels(as.factor(rownames(Species_Functional_Entities)))
  
  # compute the abundance of each species in the different conditions
  ab.sp.conditions <- lapply(conditions, function (z) {
    # In each condition, compute the abundance of each species 
    abund.sp <-  sapply(sp, function (x) {
      # Calculate the abundance of the species in the condition
      sum(ab.conditions[z,x])
    })#eo sapply
    # return the abundance of the FEs in the condition
  })#eo lapply
  
  # Add the condition name as row name
  names(ab.sp.conditions) = conditions
  
  # Transform the lists into matrix
  ab.sp.conditions <- do.call(rbind, ab.sp.conditions)
  
  #Transpose it
  sp.conditions = t(ab.sp.conditions)
  
  # we merge it to the FEs informations
  sp.fe.conditions <- cbind(Species_Functional_Entities,sp.conditions) # merge the abundance of the species with the FEs information by theur names (rownames)
  
  # we pass the data into a tidy format and add cluster information
  sp.fe.conditions.tidy <- sp.fe.conditions %>%
    rownames_to_column(var = "Species") %>%
    gather(key = "Year", value = "Abundance", -Species, -FE) %>%  # Regroup all the years in one column
    merge(fe_binded_tidy %>% 
            select(FE, Cluster, Year), by = c("FE", "Year")) # We add the cluster information to the data
  
  # We iterate for each site
  Site_Cluster_abundances <- lapply(Sites,function(site){
    
    # We retrieve the years in which the site is present
    position <-which(str_detect(unique(sp.fe.conditions.tidy$Year), site))
    Years= unique(sp.fe.conditions.tidy$Year)[position]
    
    #filter the data to retrieve only the species present in the site and count the number of species per cluster and year
    Site_Fe_sp_binded <- sp.fe.conditions.tidy %>%
      filter(Year %in% Years) %>% 
      filter(Abundance > 0) %>% 
      group_by(Cluster, Year) %>%
      mutate(Year = extract_year(Year)) %>%
      summarise(N_sp = n_distinct(Species)) %>% 
      merge(as.data.frame(colors) %>% 
              rownames_to_column(var = "Cluster"), by= "Cluster") %>%  # add color to each cluster 
      mutate(Site= site)
  }) #eo lapply
  
  # We merge all in one data frame and add a column site 
  Sites_Fe_sp_binded <- do.call(rbind, Site_Cluster_abundances)
  
  # We retrieve the maximum number of species per cluster
  max_sp <- round_up(max(Sites_Fe_sp_binded$N_sp), digits = -1) # roud max value to the nearest 10 for the y axis
  
  # We iterate for each site to draw plot
  Site_Cluster_abundances_plot <- lapply(Sites,function(site){
    
    # We draw a plot that represente number of species as point in function of the year,plot  attribute color of the cluster to the number of species (color) 
    Cluster_species_plot <- ggplot(Sites_Fe_sp_binded %>% 
                                     filter(Site== site), aes(x = Year, y = N_sp)) +
      scale_y_continuous(limits = c(0, max_sp), breaks = seq(0, max_sp, by = 2)) +
      geom_point(aes(color = colors, fill = colors), size = 4, position=position_jitter(w=0.2, h=0, seed=NULL)) +
      labs(title = paste0(site), x = " ", y = "Number of species per cluster")+
      scale_color_identity() + # color as specified in the dataset
      scale_fill_identity() + # color as specified in the dataset
      theme_classic() +
      theme(panel.border = element_rect(colour = "black", fill = NA)) +
      theme(axis.text = element_text(face = "bold", size = 12)) +
      theme(axis.title = element_text(face = "bold")) +
      theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid")) +
      theme(legend.position = "none") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  }) #eo lapply
  
  # For each element of the list, we add the Site as entry
  names(Site_Cluster_abundances_plot) <- Sites
  
  # We draw the plots in one using the gridExtra package
  n_col = length(Site_Cluster_abundances_plot) # number of column = number of year
  
  Sites_functionnal_Abundance <- grid.arrange(grobs = Site_Cluster_abundances_plot, ncol = n_col)
  
  #____________________________
  # We can now calculate the abundance of each cluster 
  # We retrieve the previously calculated abundances per FE linked to cluster
  # We iterate for each Site 
  
  Site_Functionnal_Clusters_abundances <- lapply(Sites,function(site){
    # We retrieve the years in which the site is present
    position <-which(str_detect(unique(fe_binded_tidy$Year), site))
    Years= unique(fe_binded_tidy$Year)[position]
    
    #filter the data to retrieve only the FE present in the site and calculate the abundance of each cluster
    Site_Fe_binded <- fe_binded_tidy %>%
      filter(Year %in% Years) %>% 
      filter(Abundance > 0) %>% 
      group_by(Cluster, Year) %>% 
      summarise(Abundance = sum(Abundance))%>% 
      # we calculate the relative abundance of each cluster
      group_by(Year) %>% 
      mutate(Abundance_rel= Abundance/sum(Abundance)*100) %>% 
      mutate(Year = extract_year(Year)) # extract only year value
    
    # We filter the list of colors with only colors in the Site_binded$cluster
    colors <- colors[as.character(unique(Site_Fe_binded$Cluster))]
    
    # We apply the color to each cluster 
    Site_Fe_binded <- assign_colors(Site_Fe_binded, cluster_column = "Cluster", colors = colors)
    
    # We can draw plots that represent the abundance of each cluster in function of the year and represent clusters from the top : from 1 to 8
    # First we need to precise the order that we want to display colors 
    Site_Fe_binded$color <- factor(Site_Fe_binded$color, levels = colors) # enshure that colors are display in the order precise by the named vector colors
    
    # Draw plot 
    Cluster_abundance_plot <- ggplot(Site_Fe_binded, aes(x = Year, y = Abundance_rel)) +
      geom_bar(stat = "identity", aes(fill = color), position = "Stack") +
      labs(title = paste0(site), x = " ", y = "Relative abundance (%)")+
      scale_fill_identity() + # color as specified in the dataset
      theme_classic() +
      theme(panel.border = element_rect(colour = "black", fill = NA)) +
      theme(axis.text = element_text(face = "bold", size = 12)) +
      theme(axis.title = element_text(face = "bold")) +
      theme(axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid")) +
      theme(legend.position = "none") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  }) #eo lapply
  
  # For each element of the list, we add the Site as entry
  names(Site_Functionnal_Clusters_abundances) <- Sites
  
  # retrieve table of abundances in the data part of each plot 
  Tab_Site_cluster_abundances <- lapply(Site_Functionnal_Clusters_abundances, function(x) {
    # add name of the list as Site and return it 
    Site_cluster_tab <- x$data
  })
  
  # For each element of the list, we add the Site as entry
  names(Tab_Site_cluster_abundances) <- Sites
  
  # We add the Site as a value in one clumn of each element of the list 
  Tab_Site_cluster_abundances <- lapply(1:length(Tab_Site_cluster_abundances), function(i) {
    Tab_Site_cluster_abundances[[i]]$Site <- Sites[i]
    Tab_Site_cluster_abundances[[i]]
  })
  
  # merge all in one data frame
  Tab_Site_cluster_abundances <- do.call(rbind, Tab_Site_cluster_abundances)
  
  # We select only variables in which we are interested in 
  Tab_Site_cluster_abundances <- Tab_Site_cluster_abundances %>% 
    select(Site, Year, Cluster, Abundance_rel) %>% 
    arrange(Site, Year, Cluster)
  
  # We draw the plots in one using the gridExtra package
  n_col = length(Site_Functionnal_Clusters_abundances) # number of column
  
  Sites_Functionnal_Clusters_abundances <- grid.arrange(grobs = Site_Functionnal_Clusters_abundances, ncol = n_col)
  
  return(list(Cluster_plot = Cluster_plot, Sites_functionnal_cluster= Sites_functionnal_cluster, Sites_functionnal_Abundance= Sites_functionnal_Abundance,Tab_Site_cluster_abundances = Tab_Site_cluster_abundances, Sites_Functionnal_Clusters_abundances= Sites_Functionnal_Clusters_abundances))
  
} # eo function


