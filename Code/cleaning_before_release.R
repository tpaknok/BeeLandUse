Alpha <- read.csv("Data/Alpha.csv")
Beta <- read.csv("Data/Beta.csv")
Genus <- read.csv("Data/Genus.csv")

exclude <- NULL

library(tidyverse)

cleaned_alpha <- Alpha %>% filter(!studyID %in% exclude) %>% select("studyID",  
                                                                    "Lat",
                                                                    "Long",
                                                                    "Exclude.AM.Est.Shannon",
                                                                    "qPD.S",
                                                                    "LC.author.coarse",
                                                                    "Agr1",
                                                                    "Urb1",
                                                                    "Nat1",
                                                                    "abundanceMethod",
                                                                    "PC1",
                                                                    "PC2",
                                                                    "SC")

cleaned_beta <- Beta %>% filter(!study %in% exclude) %>% select("study",
                                                                "T.gamma",
                                                                "T.gamma.SC",
                                                                "P.gamma",
                                                                "LC",
                                                                "extent",
                                                                "abundanceMethod",
                                                                "mean_PC1",
                                                                "mean_PC2",
                                                                "SamplingEffort.sd",
                                                                "SC.sd",
                                                                "mean_SC",
                                                                "site",
                                                                "mean_T_alpha",
                                                                "mean_P_alpha")

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
                                                                   "LC.author.coarse",
                                                                   "orig_value")

colnames(cleaned_genus)[[2]] <- "genus"

write.csv(cleaned_alpha,"Data/Cleaned_alpha.csv")
write.csv(cleaned_beta,"Data/Cleaned_beta.csv")
write.csv(cleaned_genus,"Data/Cleaned_genus.csv")

