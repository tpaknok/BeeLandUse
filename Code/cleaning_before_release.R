Alpha <- read.csv("Data/Alpha.csv")
Beta <- read.csv("Data/Beta.csv")
Genus <- read.csv("Data/Genus.csv")

exclude <- c("Guti01","Land01","Liuy01","Liuy02","otie01")

library(tidyverse)

cleaned_alpha <- Alpha %>% filter(!studyID %in% exclude) %>% select("studyID",  
                                                                    "Lat",
                                                                    "Long",
                                                                    "Exclude.AM.Est.Shannon",
                                                                    "qPD.S",
                                                                    "Agr1",
                                                                    "Urb1",
                                                                    "abundanceMethod",
                                                                    "PC1",
                                                                    "PC2",
                                                                    "SC")

cleaned_beta <- Beta %>% filter(!study %in% exclude) %>% select("study",
                                                                "T_beta",
                                                                "mean.dist",
                                                                "metrics",
                                                                "LC",
                                                                "mean.spatial.dist",
                                                                "SC.diff.mean",
                                                                "abundanceMethod",
                                                                "mean_PC1",
                                                                "mean_PC2",
                                                                "site")

cleaned_genus<- Genus %>% filter(!studyID %in% exclude) %>% select("studyID", 
                                                                   "species",
                                                                   "Lat",
                                                                   "Long",
                                                                   "value",
                                                                   "Agriculture",
                                                                   "Urban",
                                                                   "abundanceMethod",
                                                                   "SamplingEffort",
                                                                   "PC1",
                                                                   "PC2",
                                                                   "LC.author.coarse")

colnames(cleaned_beta)[[3]] <- "P_beta"
colnames(cleaned_genus)[[5]] <- "abundance"
colnames(cleaned_genus)[[2]] <- "genus"

write.csv(cleaned_alpha,"Data/Cleaned_alpha.csv")
write.csv(cleaned_beta,"Data/Cleaned_beta.csv")
write.csv(cleaned_genus,"Data/Cleaned_genus.csv")
