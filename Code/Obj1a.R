#analysis.df.PD <- read.csv("Data/Alpha.csv") original file 
analysis.df.PD <- read.csv("Data/cleaned_Alpha.csv")

### Text S3
library(geosphere)
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
avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
mTD2 <- fitme(log(Exclude.AM.Est.Shannon)~Agr1+Urb1+abundanceMethod+PC1+PC2+(Agr1+Urb1||studyID)+Matern(1|Long + Lat %in% studyID),data=analysis.df.PD,
              family=gaussian,weights.form=~SC,method="ML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth")
              ,control.HLfit=list(NbThreads=max(avail_thr, 1L)))
summary(mTD2,details=T)
r2(mTD2)
plot(simulateResiduals(mTD2))
mTD2_R2c <- pseudoR2(mTD2)
mTD2_R2m <- pseudoR2(mTD2,nullform=.~1+(Agr1+Urb1||studyID)+Matern(1|Long + Lat %in% studyID))
tableS1_1 <- data.frame(summary(mTD2,details=T)$beta_table,R2c=mTD2_R2c, R2m=mTD2_R2m)
AIC(mTD2)

mPD.S2 <- fitme(log(qPD.S)~Agr1+Urb1+abundanceMethod+PC1+PC2+(Agr1+Urb1||studyID)+Matern(1|Long + Lat %in% studyID),data=analysis.df.PD,
                family=gaussian,weights.form=~SC,method="ML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth"))
summary(mPD.S2,details=T)
plot(simulateResiduals(mPD.S2))
mPD.S2_R2c <- pseudoR2(mPD.S2)
mPD.S2_R2m <- pseudoR2(mPD.S2,nullform=.~1+(Agr1+Urb1||studyID)+Matern(1|Long + Lat %in% studyID))
tableS1_2 <- data.frame(summary(mPD.S2,details=T)$beta_table,R2c=mPD.S2_R2c,R2m=mPD.S2_R2m)

AIC(mPD.S2)

analysis.df.PD$Covariate.Est.Shannon <- scale(log(analysis.df.PD$Exclude.AM.Est.Shannon))

mPD.S.controlled2 <- fitme(log(qPD.S)~Agr1*Covariate.Est.Shannon+Urb1*Covariate.Est.Shannon+abundanceMethod+PC1+PC2+(Agr1*Covariate.Est.Shannon+Urb1*Covariate.Est.Shannon||studyID)+Matern(1|Long + Lat %in% studyID),data=analysis.df.PD,
                           family=gaussian,weights.form=~SC,method="ML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth"))
summary(mPD.S.controlled2,details=T)
plot(simulateResiduals(mPD.S.controlled2))

mPD.S.controlled2_R2c <- pseudoR2(mPD.S.controlled2)
mPD.S.controlled2_R2m <- pseudoR2(mPD.S.controlled2,nullform=.~1+(Agr1+Urb1||studyID)+Matern(1|Long + Lat %in% studyID))
tableS1_3 <- data.frame(summary(mPD.S.controlled2,details=T)$beta_table,R2c=mPD.S.controlled2_R2c,R2m=mPD.S.controlled2_R2m)
tableS1 <- rbind(tableS1_1,tableS1_2,tableS1_3)
tableS1 <- round(tableS1,3)
write.csv(tableS1,"Table/tableS1.csv")

###
newdataTD = expand.grid(Agr1=c(0,1),Urb1=c(0,1),abundanceMethod="Abundance",PC1=mean(analysis.df.PD$PC1),PC2=mean(analysis.df.PD$PC2))
newdataTD = newdataTD[!(newdataTD$Agr1 == 1 & newdataTD$Urb1 == 1),]

newdataTD$predicted_y <- predict(mTD2,newdata=newdataTD,re.form=NA,type="response")
newdataTD <- cbind(newdataTD,get_intervals(object=mTD2,newdata=newdataTD,levels=0.95,re.form=NA,type="response"))
newdataTD$LC <- c("Natural","Agricultural","Urban")
newdataTD$LC <- factor(newdataTD$LC, levels = c("Natural","Agricultural","Urban"))
newdataTD[,c("predicted_y","predVar_0.025","predVar_0.975")] <- exp(newdataTD[,c("predicted_y","predVar_0.025","predVar_0.975")])
natural_ref <- c(newdataTD[newdataTD$LC == "Natural","predicted_y"])
newdataTD$percent_mean <- (c(newdataTD$predicted_y)-natural_ref)/natural_ref*100
newdataTD$percent_low <- (c(newdataTD$predVar_0.025)-natural_ref)/natural_ref*100
newdataTD$percent_high <- (c(newdataTD$predVar_0.975)-natural_ref)/natural_ref*100

newdataPD = expand.grid(Agr1=c(0,1),Urb1=c(0,1),abundanceMethod="Abundance",PC1=mean(analysis.df.PD$PC1),PC2=mean(analysis.df.PD$PC2))
newdataPD = newdataPD[!(newdataPD$Agr1 == 1 & newdataPD$Urb1 == 1),]
newdataPD$predicted_y <- predict(mPD.S2,newdata=newdataPD,re.form=NA,type="response")
newdataPD <- cbind(newdataPD,get_intervals(object=mPD.S2,newdata=newdataPD,levels=0.95,re.form=NA,type="response"))
newdataPD$LC <- c("Natural","Agricultural","Urban")
newdataPD$LC <- factor(newdataPD$LC, levels = c("Natural","Agricultural","Urban"))
newdataPD[,c("predicted_y","predVar_0.025","predVar_0.975")] <- exp(newdataPD[,c("predicted_y","predVar_0.025","predVar_0.975")])
natural_ref <- c(newdataPD[newdataPD$LC == "Natural","predicted_y"])
newdataPD$percent_mean <- (c(newdataPD$predicted_y)-natural_ref)/natural_ref*100
newdataPD$percent_low <- (c(newdataPD$predVar_0.025)-natural_ref)/natural_ref*100
newdataPD$percent_high <- (c(newdataPD$predVar_0.975)-natural_ref)/natural_ref*100

newdataPD_C = expand.grid(Agr1=c(0,1),Urb1=c(0,1),abundanceMethod="Abundance",PC1=mean(analysis.df.PD$PC1),PC2=mean(analysis.df.PD$PC2),Covariate.Est.Shannon=mean(analysis.df.PD$Covariate.Est.Shannon))
newdataPD_C = newdataPD_C[!(newdataPD_C$Agr1 == 1 & newdataPD_C$Urb1 == 1),]
newdataPD_C$predicted_y <- predict(mPD.S.controlled2,newdata=newdataPD_C,re.form=NA,type="response")
newdataPD_C <- cbind(newdataPD_C,get_intervals(object=mPD.S.controlled2,newdata=newdataPD_C,levels=0.95,re.form=NA,type="response"))
newdataPD_C$LC <- c("Natural","Agricultural","Urban")
newdataPD_C$LC <- factor(newdataPD_C$LC, levels = c("Natural","Agricultural","Urban"))
newdataPD_C[,c("predicted_y","predVar_0.025","predVar_0.975")] <- exp(newdataPD_C[,c("predicted_y","predVar_0.025","predVar_0.975")])
natural_ref <- c(newdataPD_C[newdataPD_C$LC == "Natural","predicted_y"])
newdataPD_C$percent_mean <- (c(newdataPD_C$predicted_y)-natural_ref)/natural_ref*100
newdataPD_C$percent_low <- (c(newdataPD_C$predVar_0.025)-natural_ref)/natural_ref*100
newdataPD_C$percent_high <- (c(newdataPD_C$predVar_0.975)-natural_ref)/natural_ref*100

pTD <- ggplot()+
  geom_hline(yintercept=0)+
  geom_point(data=newdataTD,aes(x=LC,y=percent_mean,colour=LC),position=position_dodge(width=0.5),size=2)+
  geom_errorbar(data=newdataTD,aes(x=LC,ymin=percent_low,ymax=percent_high,colour=LC),width=0.2,position=position_dodge(width=0.5),show.legend=T)+
  ylab(expression(paste(alpha,"-diversity difference (%)")))+
  xlab("")+
  ylim(-40,20)+
  annotate('text', label=' (b) Taxonomic', x=-Inf, y=Inf, hjust=0, vjust=1,size=4)+
  scale_colour_manual(values=c("green4","darkgoldenrod2","grey25"),name="Habitat")+
  theme_classic()+  theme(axis.text = element_text(size=10))

plot(pTD)

pPD <- ggplot()+
  geom_hline(yintercept=0)+
  geom_point(data=newdataPD,aes(x=LC,y=percent_mean,colour=LC),position=position_dodge(width=0.5),size=2)+
  geom_errorbar(data=newdataPD,aes(x=LC,ymin=percent_low,ymax=percent_high,colour=LC),width=0.2,position=position_dodge(width=0.5),show.legend=F)+
  ylab("")+
  xlab("")+
  ylim(-40,20)+
  annotate('text', label=' (c) Phylogenetic', x=-Inf, y=Inf, hjust=0, vjust=1,size=4)+
  scale_colour_manual(values=c("green4","darkgoldenrod2","grey25"))+
  theme_classic()+  theme(axis.text = element_text(size=10))

plot(pPD)

pPD_C <- ggplot()+
  geom_hline(yintercept=0)+
  geom_point(data=newdataPD_C,aes(x=LC,y=percent_mean,colour=LC,shape=LC),position=position_dodge(width=0.5),size=2)+
  geom_errorbar(data=newdataPD_C,aes(x=LC,ymin=percent_low,ymax=percent_high,colour=LC),width=0.2,position=position_dodge(width=0.5),show.legend=F)+
  ylab("")+
  xlab("")+
  ylim(-40,20)+
  annotate('text', label=' (d) Phylogenetic (Controlling \n      for taxonomic)', x=-Inf, y=Inf, hjust=0, vjust=1,size=4)+
  scale_colour_manual(values=c("green4","darkgoldenrod2","grey25"))+
  scale_shape_manual(values=c(19,4,4))+
  theme_classic()+
  theme( axis.text = element_text(size=10))

plot(pPD_C)

library(ggpubr)

p_alpha <- ggarrange(pTD,pPD,pPD_C,nrow=1,ncol=3,common.legend=T,legend="bottom")
plot(p_alpha)
ggsave("Figure/p_alpha.tiff",width=8,height=4,compression="lzw",dpi=600)

