#packages & publication theme: 
##############
library(factoextra)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sp)
library(randomForest)
library(caTools)
library(e1071)
library(GGally)
library(plyr)

theme_Publication <- function(base_size=13, base_family="Garamond") {
  library(grid)
  library(ggthemes)
  library(gridExtra)
  windowsFonts("Garamond" = windowsFont("Garamond"))
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}


###############
#Guilds <- read.csv("D:/00_Work/Panama_Microbiome_Project/Guilds.csv")
#Seedling_Metadata = read.csv("D:/00_Work/Panama_Microbiome_Project/Seedling_Metadata_05.27.2021.csv")
#SeedlingGuilds = merge(Guilds, Seedling_Metadata, by = "Sample_Name")

SeedlingGuilds = read.csv("D:/00_Work/Panama_Microbiome_Project/SeedlingGuilds.csv", header=TRUE)
SeedlingGuilds = SeedlingGuilds %>%
  filter(Status != "Natural_Germination")
############################
prcomp_Seedlings =  prcomp(SeedlingGuilds[,12:9186])
##############################
summary(prcomp_Seedlings)
################################
SCREE_PLOT_Seedlings  = fviz_eig(prcomp_Seedlings)
#############################

VARIABLES_Seedlings = fviz_pca_var(prcomp_Seedlings,
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
##########################
INDIVIDUALS_SeedlingSpecies  = fviz_pca_ind(prcomp_Seedlings, geom = "point", 
                    habillage = SeedlingGuilds$Seedling_Species, addEllipses = TRUE, 
                    ellipse.level = 0.8) +
                    ggtitle("Fungal_Community (Seedling Species)") +
                    theme_Publication() +
                    scale_colour_Publication()
######################
INDIVIDUALS_AdultSpecies = fviz_pca_ind(prcomp_Seedlings, geom = "point", 
                    habillage = SeedlingGuilds$Adult_Species, addEllipses = TRUE, 
                    ellipse.level = 0.8)+
                    ggtitle("Fungal_Community (Adult Host)") +
                    theme_Publication() +
                    scale_colour_Publication()
################
INDIVIDUALS_SeedlingCondition = fviz_pca_ind(prcomp_Seedlings, geom = "point", 
                    habillage = SeedlingGuilds$Status, addEllipses = TRUE, 
                    ellipse.level = 0.8,)+
                    ggtitle("Fungal_Community (Seedling Status)") +
                    theme_Publication() +
                    scale_colour_Publication()
####################

pca_SeedlingSpecies = data.frame(PC1 = prcomp_Seedlings$x[,1],
                        PC2 = prcomp_Seedlings$x[,2],
                        label = SeedlingsGuilds$Sample_Name,
                        classification = SeedlingsGuilds$Seedling_Species)

SeedlingSpeciesPCA = ggplot(pca_SeedlingSpecies, aes(x = PC1, y = PC2, label = label, col = classification)) +
  geom_point() +
  ggtitle("Fungal_Community (Seedling Species)")+
  ggrepel::geom_text_repel(cex=2.5)+
  theme_Publication() +
  scale_colour_Publication()

##############################  

pca_AdultHost= data.frame(PC1 = prcomp_Seedlings$x[,1],
                         PC2 = prcomp_Seedlings$x[,2],
                         label = SeedlingsGuilds$Sample_Name,
                         classification = SeedlingsGuilds$Adult_Species)

AdultHostPCa = ggplot(pca_AdultHost, aes(x = PC1, y = PC2, label = label, col = classification)) +
  geom_point() +
  ggtitle("Fungal_Community (Adult Host)")+
  ggrepel::geom_text_repel(cex=2.5)+
  theme_Publication() +
  scale_colour_Publication()

###############################

pca_SeedlingStatus= data.frame(PC1 = prcomp_Seedlings$x[,1],
                          PC2 = prcomp_Seedlings$x[,2],
                          label = SeedlingsGuilds$Sample_Name,
                          classification = SeedlingsGuilds$Status)

SeedlingStatusPCA = ggplot(pca_SeedlingStatus, aes(x = PC1, y = PC2, label = label, col = classification)) +
  geom_point() +
  ggtitle("Fungal_Community (Adult Host)")+
  ggrepel::geom_text_repel(cex=2.5)+
  theme_Publication() +
  scale_colour_Publication()

#################################


#MDS 


dist_Seedlings = dist(SeedlingGuilds[,12:9186], method = "canberra")
mds_Seedling = cmdscale(dist_Seedlings, eig = TRUE, k=2)

dist_Seedlings = bc_seedlings
#########
mds_SeedlingSpecies = data.frame(
  MDS1 = mds_Seedling$points[,1],
  MDS2 = mds_Seedling$points[,2],
  label= SeedlingGuilds$Sample_Name,
  SeedlingSpecies = SeedlingGuilds$Seedling_Species)

#SeedlingSpecieshulls = ddply(mds_SeedlingSpecies, "SeedlingSpecies", find_hull)

#######
mds_AdultHost = data.frame(
  MDS1 = mds_Seedling$points[,1],
  MDS2 = mds_Seedling$points[,2],
  label= SeedlingGuilds$Sample_Name,
  AdultHost = SeedlingGuilds$Adult_Species)
#########
mds_SeedlingStatus = data.frame(
  MDS1 = mds_Seedling$points[,1],
  MDS2 = mds_Seedling$points[,2],
  label= SeedlingGuilds$Sample_Name,
  Status = SeedlingGuilds$Status)
#########
MDS_SeedlingSpecies = ggplot(mds_SeedlingSpecies, aes(
x = MDS1, y = MDS2,
label = label, col = SeedlingSpecies)) +
geom_point() +
theme_Publication() +
scale_colour_Publication()+
ggrepel::geom_text_repel(cex = 2.5)
#######
MDS_AdultHost = ggplot(mds_AdultHost, aes(
  x = MDS1, y = MDS2,
  label = label, col = AdultHost)) +
  geom_point() +
  theme_Publication() +
  scale_colour_Publication()+
  ggrepel::geom_text_repel(cex = 2.5)
#########
MDS_SeedlingStatus = ggplot(mds_SeedlingStatus, aes(
  x = MDS1, y = MDS2,
  label = label, col = Status)) +
  geom_point() +
  theme_Publication() +
  scale_colour_Publication()+
  ggrepel::geom_text_repel(cex = 2.5)
#########




###NMDS, PCOA, BC dissimilarity
bc_seedlings = vegdist(SeedlingGuilds[,12:9186], method = "bray")

NMDS = metaMDS(SeedlingGuilds[,12:9186], k = 2, try = 20, trymax = 100)

stressplot(NMDS)

