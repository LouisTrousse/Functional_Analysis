#### Data preparation test ####

# First we merge Fes and Spe_fes by Fe to have one table containing all the information and from wich we can retrieve each dataset 
Species_functionnal_traits<-merge(spe_fes %>% 
                                    rownames_to_column(var="Species"), fes %>% 
                                    rownames_to_column(var= "FE"), by ="FE") %>% 
  relocate(Species, .before=FE)

# Save Data 
write.csv(Species_functionnal_traits, "Species_functionnal_traits.csv", row.names = FALSE)

# Based on a document of that type, we can retrieve the two datasets
# Species_functionnal_traits<-read.csv2("Species_functionnal_traits.csv") # uncomment to charge 
spe_fes <- Species_functionnal_traits %>% select(Species, FE) %>% column_to_rownames(var="Species")
fes <- Species_functionnal_traits %>% select(-Species)%>%
  distinct() %>%
  column_to_rownames(var="FE") # If any duplicate FE, an error message will be printed

# The following function will check that if Spe and Fes are given separatly they are merged, If only one file is given, it will check if the format is the right one and converse it 
Check_Species_functionnal_traits <- function (Species_functionnal_traits = NULL, Fonctionnal_traits = NULL){ 
  # If both are NULL print an error 
  if (is.null(Species_functionnal_traits) & is.null(Fonctionnal_traits)){
    stop("You must provide at least one of the two datasets")
  } else if (is.null(Species_functionnal_traits) & !is.null(Fonctionnal_traits)){
    # We check if it contains a column named Species and a collum named FE or it's rownames are FEs 
    if (is.character(rownames(Fonctionnal_traits)) == TRUE & "Species" %in%  colnames(Fonctionnal_traits))
      # If the rownames are characters, we assume that they are FE
      fes <- Fonctionnal_traits %>% 
        distinct() %>%
        column_to_rownames(var="FE") # If any duplicate FE, an error message will be printed
      # We then check if the dataset contains a collumn named Species
      if ("Species" %in% colnames(Fonctionnal_traits)){
        spe_fes <- Fonctionnal_traits %>% 
          select(Species, FE) %>% 
          column_to_rownames(var="Species")
      } else {
        stop("The dataset must contain a column named Species")
      }
    } else {
      # If the rownames are not FEs, we assume that the dataset is the one containing the Species
      spe_fes <- Fonctionnal_traits %>% 
        rownames_to_column(var="Species")
      # We then check if the dataset contains a collumn named FE
      if ("FE" %in% colnames(Fonctionnal_traits)){
        fes <- Fonctionnal_traits %>% 
          column_to_rownames(var="FE") # If any duplicate FE, an error message will be printed
      } else {
        stop("The dataset must contain a column named FE")
      }
    }
    
  }
  }




# If the abundance matrix contain more than one ID and one collum for each condition, we can retrieve the dataset ab as needed for our study. 
# We will use the following function to enshure that the dataset is in the right format.
Check_Abundance_matrix <- function(Abundance_matrix, ID, Species){ 
  # check if rownames of the given abundance_matrix is an integer or an identifier
  if (is.integer(rownames(Abundance_matrix)) == TRUE){
    Abundance_matrix <- Abundance_matrix %>% 
      column_to_rownames(var=ID) %>% 
      # We now enshure that collumns are only species based on their appartenance to the list 
      select(Species)
  } else { 
    # If rowanmes aren't integer, we assume that they are the ID of the sampling unit
    # Then we only need to check if the collumns are species
    Abundance_matrix <- Abundance_matrix %>% 
      select(Species)
  }
}


# computing relative abundances of trait categories among pH conditions

# First we need to retrieve FEs_Trait as a table containing on each row one unique FE and the value of each traits


# list to store results
ab.trait_condition<-list()


# ordering FEs in trait table as in abundance table
fes_traits<-fes_traits[ colnames(ab.fe.conditions), ]


# for each trait, computing total relative abundance for each modality among FEs among pH zones
for (t in colnames(fes_traits) )
{
  # levels of trait t
  levels_t<-as.character( sort( unique(fes_traits[,t]) ) )
  
  # empty vectors to store results
  trait_t<-c()
  condition_t<-c()
  rel_ab_t<-c()
  
  # computing relative abundance of each trait modality in each condition
  for (i in condition)
    for (j in levels_t)
    {
      trait_t<-c(trait_t, j)
      condition_t<-c(condition_t, i)
      rel_ab_t<-c(rel_ab_t, sum( ab.fe.conditions[ i , which(fes_traits[,t]==j) ] )  )
    }# end of i,j
  
  # setting correcdt order for levels of conditions
  condition_t <- factor(condition_t, levels = condition )
  
  # storing results in a dataframe
  ab.trait_condition[[t]]<-data.frame(  trait_val= trait_t,  condition= condition_t, rel_ab_FE=rel_ab_t )
  
}# end of t


# Figure 3. Change in relative abundance of functional trait categories along the pH gradient. 
#ggplot does not like list of dataframes; thus we code 15 individual plots and later we grouped together

# define color palette for 13, 7, 6, 5, 3, 2 trait categories

cc13<-c("#a6cee3","#1f78b4","#b2df8a", "#33a02c", "#fb9a99", "#e31a1c","#fdbf6f","#ff7f00", "#cab2d6",  "#6a3d9a", "#ffff99","#b15928", "#1c1b1b") 
cc7<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#ffb24f")
cc6<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")
cc5<-c("#e31a1c", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
cc3<-c("#ffb24f",  "#1f78b4","#33a02c")
cc2<-c( "#a6cee3","#1f78b4")


#plot1: Morphological.form, 13 trait values, cc13
plot1 <-ggplot(data=ab.trait_condition$Morphological.form, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  guides(fill=guide_legend(ncol=2))+
  labs( fill="Traits")+
  ggtitle("Morphological form")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc13)

plot1

#plot2: Solitary.Colonial, 3 trait values, cc3
plot2 <-ggplot(data=ab.trait_condition$Solitary.Colonial, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Solitary-Colonial")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc3)

plot2

#plot 3 Max.Longevity, 7 trait values, cc7 
plot3 <-ggplot(data=ab.trait_condition$Max.Longevity, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Max Longevity")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc7)

plot3

#plot 4 Height, 5 trait values, cc5

plot4 <-ggplot(data=ab.trait_condition$Height, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Height")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc5)

plot4


#plot 5 Width, 6 trait values, cc6

plot5 <-ggplot(data=ab.trait_condition$Width, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Width")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc6)

plot5

#plot 6 Epibiosis, 3 trait values, cc3
plot6 <-ggplot(data=ab.trait_condition$Epibiosis, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs( fill="Traits")+
  ggtitle("Epibiosis")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc3)

plot6

#plot7: Energetic.resource, 3 trait values, cc3
plot7 <-ggplot(data=ab.trait_condition$Energetic.resource, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(y="Relative Abundance of Functional Categories (%)", fill="Traits")+
  ggtitle("Energetic resource")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc3)

plot7


#plot 8 Major.photosynthetic.pigments, 7 trait values, cc7 
plot8 <-ggplot(data=ab.trait_condition$Major.photosynthetic.pigments, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Photosynthetic pigments")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc7)

plot8

#plot 9 Feeding, 5 trait values, cc5

plot9 <-ggplot(data=ab.trait_condition$Feeding, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Feeding ")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc5)

plot9

#plot 10 Age.reproductive.maturity, 6 trait values, cc6

plot10 <-ggplot(data=ab.trait_condition$Age.reproductive.maturity, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Age reproductive maturity")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc6)

plot10


#plot11: Asexual.Reproduction, 2 trait values, cc2

plot11<-ggplot(data=ab.trait_condition$Asexual.Reproduction, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(x="pH conditions",  fill="Traits")+
  ggtitle("Asexual Reproduction")+
  scale_y_continuous( expand = c(0,0))+
  scale_x_discrete(labels=c("Ambient"="Ambient", "Low"="Low", "Extreme Low"= "Extreme Low"))+
  scale_fill_manual(values = cc2)
plot11

#plot12: Growth.rates, 5 trait values, cc5

plot12<-ggplot(data=ab.trait_condition$Growth.rates, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs( fill="Traits")+
  ggtitle("Growth rates")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc5)
plot12

#plot13: Calcification, 5 trait values, cc5

plot13<-ggplot(data=ab.trait_condition$Calcification, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x=element_text(angle=60,hjust=1, size=12),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        axis.title.x = element_text(size=14),
        legend.title=element_blank())+
  labs(x="pH conditions",  fill="Traits")+
  ggtitle("Calcification")+
  scale_y_continuous( expand = c(0,0))+
  scale_x_discrete(labels=c("Ambient"="Ambient", "Low"="Low", "Extreme Low"= "Extreme Low"))+
  scale_fill_manual(values = cc5)
plot13


#plot14: Chemical.defenses, 2 trait values, cc2

plot14<-ggplot(data=ab.trait_condition$Chemical.defenses, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x=element_text(angle=60,hjust=1, size=12),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        axis.title.x = element_text(size=14),
        legend.title=element_blank())+
  labs(x="pH conditions",  fill="Traits")+
  ggtitle("Chemical defenses")+
  scale_y_continuous( expand = c(0,0))+
  scale_x_discrete(labels=c("Ambient"="Ambient ", "Low"="Low ", "Extreme Low"= "Extreme Low "))+
  scale_fill_manual(values = cc2)
plot14

#plot15: Mobility, 2 trait values, cc2

plot15<-ggplot(data=ab.trait_condition$Mobility , aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x=element_text(angle=60,hjust=1, size=12),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        axis.title.x = element_text(size=14),
        legend.title=element_blank())+
  labs(x="pH conditions",  fill="Trait Value")+
  ggtitle("Mobility")+
  scale_y_continuous( expand = c(0,0))+
  scale_x_discrete(labels=c("Ambient"="Ambient", "Low"="Low", "Extreme Low"= "Extreme Low"))+
  scale_fill_manual(values = cc2)
plot15
