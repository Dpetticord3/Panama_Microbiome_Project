CamilleMetadata <- read_excel("D:/00_Work/Seedling Microbes/CamilleMetadata.xlsx")
CollectedSeedlingList <- read_excel("D:/00_Work/Seedling Microbes/CollectedSeedlingList.xlsx")

x= merge(CamilleMetadata, CollectedSeedlingList, by.x = "Seedling_ID", by.y = "Seedling_ID", all.x = TRUE, incomparables=TRUE)

new =separate(x, Seedling_ID, into =c("A", "B", "C", "D"), sep = "_", remove = FALSE)


write.csv(new, "Seedling_Metadata.csv")
