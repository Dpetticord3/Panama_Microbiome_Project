#########################
##### PRELIMINARY #######
##### STRI/KU CSD #######
#########################


#########################
######## CLEAN ##########
###### THE METADATA #####
#########################

#load packages:
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(vegan)
library(tidyverse)
library(mvabund)

# import meta data
meta.seed<-read.table('data/Seedling_Metadata_05.27.2021.csv', header=TRUE,sep=',')
meta.soil.1<-read.table('data/Sequenced_soil_metadata_1.csv', header=TRUE,sep=',')
meta.soil.2<-read.table('data/Sequenced_soil_metadata_2.csv', header=TRUE,sep=',')
str(meta.soil.1)
str(meta.soil.2)
meta.soil<-meta.soil.1 %>% full_join(meta.soil.2, by ="Soil_ID") %>%
  select(-c("Tube.Sample.Name","Soil_ID","Tube_label","Extraction_date",
            "Sample_weight_mg","Notes","Adult_tag","Date_soil_collected",
            "Date_soil_processed","Soil_collector","Seedling_ID")) %>%
  rename(Adult_Species = Adult_species, Sample_Name= Sample.Name,
         Sample_Type=Sample_type, Seedling_Species = Exclosure)
meta.seed<-meta.seed %>%
  select(-c("Seedling_ID","AdultTag","Tube_Sample","Sample","Date_sample_collected","Status"))
# add Sample_Type:
meta.seed.root<-filter(meta.seed, grepl("_R",Sample_Name)) %>%
  mutate(Sample_Type="Root")
meta.seed.leaf<-filter(meta.seed, grepl("_L",Sample_Name)) %>%
  mutate(Sample_Type="Leaf")
meta.seed<-rbind(meta.seed.root,meta.seed.leaf)

#join meta.seed and meta.soil:

meta<-bind_rows(meta.seed,meta.soil)
meta<- column_to_rownames(meta, "Sample_Name")
str(meta)
meta$Sample_Type<- as.factor(meta$Sample_Type)
meta$Adult_Species2<-"HT"
meta$Adult_Species2[meta$Adult_Species=="SA"]<-"SA"
meta$Adult_Species2[meta$Adult_Species=="DG"]<-"DG"
meta$Adult_Species2[meta$Adult_Species=="HC"]<-"HC"
meta$Adult_Species<-as.factor(meta$Adult_Species2)
meta$Adult_Species2<-NULL
saveRDS(meta, file="data/meta.data.complete_5.28.21.RDS")