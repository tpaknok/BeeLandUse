#analysis.df.PD <- read.csv("Data/Alpha.csv")
source("Code/get_R2.R")
analysis.df.PD <- read.csv("Data/Cleaned_alpha.csv")

analysis.df.PD$LC.author.coarse <- relevel(factor(analysis.df.PD$LC.author.coarse,ordered=F),ref="Natural")
analysis.df.PD1 <- analysis.df.PD
analysis.df.PD2 <- analysis.df.PD

analysis.df.PD1$LC.author.coarse <- relevel(factor(analysis.df.PD1$LC.author.coarse,ordered=F),ref="Agriculture")
analysis.df.PD2$LC.author.coarse <- relevel(factor(analysis.df.PD2$LC.author.coarse,ordered=F),ref="Urban")

### Text S3
library(geosphere)
analysis.df.PD$Long <- as.numeric(analysis.df.PD$Long)
dist <- as.matrix(distm(cbind(analysis.df.PD$Long,analysis.df.PD$Lat)))
study <- as.matrix(dist(as.numeric(as.factor(analysis.df.PD$studyID))))
SR <- as.matrix(dist(analysis.df.PD$Exclude.AM.Est.Shannon))

SAC_data <- data.frame(Dist=c(dist)[c(lower.tri(dist))],
                       Study=c(study)[c(lower.tri(study))],
                       SR=c(SR)[c(lower.tri(SR))])

m_SAC <- lm(SR~log(Dist+1),data=SAC_data)
summary(m_SAC)

SAC_result <- SAC_coef <- NULL
for (i in 1:length(unique(analysis.df.PD$studyID))) {
  subset_PD_df <- subset(analysis.df.PD,studyID == unique(analysis.df.PD$studyID)[[i]])
  
  if (nrow(subset_PD_df) >= 5){
    dist <- as.matrix(distm(cbind(subset_PD_df$Long,subset_PD_df$Lat)))
    study <- as.matrix(dist(as.numeric(as.factor(subset_PD_df$studyID))))
    SR <- as.matrix(dist(subset_PD_df$Exclude.AM.Est.Shannon))
    
    SAC_data <- data.frame(Dist=c(dist)[c(lower.tri(dist))],
                           Study=c(study)[c(lower.tri(study))],
                           SR=c(SR)[c(lower.tri(SR))])
    m_SAC <- lm(SR~log(Dist+1),data=SAC_data)
    summary(m_SAC)
    SAC_result[i] <-  summary(m_SAC)$r.squared
    SAC_coef[i] <- coefficients(m_SAC)[[2]]
  }
}

mean(SAC_result,na.rm=T)
max(SAC_result,na.rm=T)
mean(SAC_coef,na.rm=T)

### Obj 1a spaMM
library(spaMM)
library(DHARMa)
library(performance)
avail_thr <- parallel::detectCores(logical=FALSE) - 2L 

mTD2 <- fitme(log(Exclude.AM.Est.Shannon)~LC.author.coarse+abundanceMethod+PC1+PC2+(Agr1+Urb1||studyID)+Matern(1|Long + Lat %in% studyID),data=analysis.df.PD,
              family=gaussian,weights.form=~SC,method="REML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth")
              ,control.HLfit=list(NbThreads=max(avail_thr, 1L))) #note that this gives the same result as Agr1+Urb1 / LC.author.coarse in FE. Tested without SA terms in spaMM and glmmTMB
summary(mTD2,details=T)
aov_TD <- anova(mTD2,method="t.Chisq")

mTD2_no_SC <- fitme(log(Exclude.AM.Est.Shannon)~LC.author.coarse+abundanceMethod+PC1+PC2+(Agr1+Urb1||studyID),data=analysis.df.PD,
              family=gaussian,weights.form=~SC,method="REML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth")
              ,control.HLfit=list(NbThreads=max(avail_thr, 1L))) #note that this gives the same result as Agr1+Urb1 / LC.author.coarse in FE. Tested without SA terms in spaMM and glmmTMB

AIC(mTD2)[[2]]-AIC(mTD2_no_SC)[[2]]

mTD2_Agr <- fitme(log(Exclude.AM.Est.Shannon)~LC.author.coarse+abundanceMethod+PC1+PC2+(Nat1+Urb1||studyID)+Matern(1|Long + Lat %in% studyID),data=analysis.df.PD1,
              family=gaussian,weights.form=~SC,method="REML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth")
              ,control.HLfit=list(NbThreads=max(avail_thr, 1L))) #obtain p-value for contrast
summary(mTD2_Agr,details=T)
anova(mTD2_Agr,method="t.Chisq")
mTD2_Urb <- fitme(log(Exclude.AM.Est.Shannon)~LC.author.coarse+abundanceMethod+PC1+PC2+(Nat1+Agr1||studyID)+Matern(1|Long + Lat %in% studyID),data=analysis.df.PD2,
                  family=gaussian,weights.form=~SC,method="REML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth")
                  ,control.HLfit=list(NbThreads=max(avail_thr, 1L))) #obtain p-value for contrast
summary(mTD2_Urb,details=T) #results largely stable regardless of reference
anova(mTD2_Urb,method="t.Chisq")

mTD2_R2m <- get_R2(mTD2)[[1]]
mTD2_R2c <- get_R2(mTD2)[[2]]

plot(simulateResiduals(mTD2))
pair_comp_TD <- rbind(summary(mTD2,details=T,verbose=F)$beta_table[2:3,],
                   summary(mTD2_Agr,details=T,verbose=F)$beta_table[3,])

rownames(pair_comp_TD) <- c("Nat-Agr","Nat-Urb","Agr-Urb")
pair_comp_TD <- as.data.frame(pair_comp_TD)
pair_comp_TD$adjusted_p <- p.adjust(pair_comp_TD[,4],"fdr")

tableS1_1 <- data.frame(aov_TD,R2c=mTD2_R2c, R2m=mTD2_R2m)

###
mPD.S2 <- fitme(log(qPD.S)~LC.author.coarse+abundanceMethod+PC1+PC2+(Agr1+Urb1||studyID)+Matern(1|Long + Lat %in% studyID),data=analysis.df.PD,
                family=gaussian,weights.form=~SC,method="REML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth"))
summary(mPD.S2,details=T)

mPD.S2_no_SC <- fitme(log(qPD.S)~LC.author.coarse+abundanceMethod+PC1+PC2+(Agr1+Urb1||studyID),data=analysis.df.PD,
                family=gaussian,weights.form=~SC,method="REML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth"))
summary(mPD.S2,details=T)

AIC(mPD.S2)[[2]]- AIC(mPD.S2_no_SC)[[2]]
mPD.S2_Agr <- fitme(log(qPD.S)~LC.author.coarse+abundanceMethod+PC1+PC2+(Nat1+Urb1||studyID)+Matern(1|Long + Lat %in% studyID),data=analysis.df.PD1,
                  family=gaussian,weights.form=~SC,method="REML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth")
                  ,control.HLfit=list(NbThreads=max(avail_thr, 1L))) #note that this gives the same result as Agr1+Urb1 in FE. Tested without SA terms in spaMM and glmmTMB
summary(mTD2,details=T)

pair_comp_PD <- rbind(summary(mPD.S2,details=T,verbose=F)$beta_table[2:3,],
                   summary(mPD.S2_Agr,details=T,verbose=F)$beta_table[3,])

rownames(pair_comp_PD) <- c("Nat-Agr","Nat-Urb","Agr-Urb")
pair_comp_PD <- as.data.frame(pair_comp_PD)
pair_comp_PD$adjusted_p <- p.adjust(pair_comp_PD[,4],"fdr")

aov_PD <- anova(mPD.S2,method="t.Chisq")

plot(simulateResiduals(mPD.S2))
mPD.S2_R2m <- get_R2(mPD.S2)[[1]]
mPD.S2_R2c <- get_R2(mPD.S2)[[2]]

tableS1_2 <- data.frame(aov_PD,R2c=mPD.S2_R2c,R2m=mPD.S2_R2m)

###
analysis.df.PD$Covariate.Est.Shannon <- scale(log(analysis.df.PD$Exclude.AM.Est.Shannon))

mPD.S.controlled2.full <- fitme(log(qPD.S)~LC.author.coarse*Covariate.Est.Shannon+abundanceMethod+PC1+PC2+(Agr1*Covariate.Est.Shannon+Urb1*Covariate.Est.Shannon||studyID)+Matern(1|Long + Lat %in% studyID),data=analysis.df.PD,
                           family=gaussian,weights.form=~SC,method="REML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth"))

summary(mPD.S.controlled2.full,details=T)
aov_controlled <- anova(mPD.S.controlled2.full,method="t.Chisq")
plot(simulateResiduals(mPD.S.controlled2.full))

mPD.S.controlled2 <- fitme(log(qPD.S)~LC.author.coarse+Covariate.Est.Shannon+abundanceMethod+PC1+PC2+(Agr1+Covariate.Est.Shannon+Urb1||studyID)+Matern(1|Long + Lat %in% studyID),data=analysis.df.PD,
                           family=gaussian,weights.form=~SC,method="REML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth"))
summary(mPD.S.controlled2,details=T)
plot(simulateResiduals(mPD.S.controlled2))

aov_PD_controlled <- anova(mPD.S.controlled2,method="t.Chisq")
mPD.S.controlled2_R2m <- get_R2(mPD.S.controlled2)[[1]]
mPD.S.controlled2_R2c <-   get_R2(mPD.S.controlled2)[[2]]
tableS1_3 <- data.frame(aov_PD_controlled,R2c=mPD.S.controlled2_R2c,R2m=mPD.S.controlled2_R2m)
tableS1 <- rbind(tableS1_1,tableS1_2,tableS1_3)
tableS1 <- round(tableS1,3)
write.csv(tableS1,"Table/Alpha_table.csv")

write.csv(rbind(pair_comp_TD,pair_comp_PD),"Table/Alpha_table_pair_comp.csv")