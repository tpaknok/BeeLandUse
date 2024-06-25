library(interactions)
library(ggplot2)
library(lmerTest)
library(lme4)
library(ggeffects)
library(performance)
P_beta <- read.csv("Data/cleaned_Beta.csv")
all_results_reg <- all_graph_df <- all_table <- NULL
beta.metrics <- c("Btotal","Bturn","Bab")
result_reg <- graph_df<-NULL

for (i in 1:3) {
  for (j in 1:3) {
    message(beta.metrics[[j]]," ",i)
    subset.beta <- subset(P_beta,metrics == beta.metrics[[j]])
    subset.beta$LC <- factor(subset.beta$LC, ordered=F) 
    subset.beta$LC<-relevel(subset.beta$LC, ref = "Natural")
    
    subset.beta$site <- log(subset.beta$site)
    #subset.beta$scaled.SC.diff <- c(scale(subset.beta$SC.diff.mean)) #not using this. SC.diff.mean at zero has a true meaning
    subset.beta$scaled.spatial.dist <- c(scale(log(subset.beta$mean.spatial.dist)))
    subset.beta$Agriculture <- 0
    subset.beta$Agriculture[subset.beta$LC == "Agriculture"] <- 1
    subset.beta$Urban <- 0
    subset.beta$Urban[subset.beta$LC == "Urban"] <- 1
    
    results_reg  <- NULL
    graph_df <- NULL
    
    if (i == 1) {
      m <- lmer(T_beta~LC*scaled.spatial.dist+SC.diff.mean+mean_PC1+mean_PC2+abundanceMethod+(1|study),data=subset.beta)
    }
    if (i == 2) {
      m <- lmer(P_beta~LC*scaled.spatial.dist+SC.diff.mean+abundanceMethod+mean_PC1+mean_PC2+(1|study),data=subset.beta)
    }
    if (i==3){
      subset.beta$scaled.T.beta <- scale(subset.beta$T_beta)
      m <- lmer(P_beta~LC*scaled.spatial.dist+LC*scaled.T.beta+SC.diff.mean+abundanceMethod++mean_PC1+mean_PC2+(1|study),data=subset.beta)
      
      y_lab_corr <- ifelse(j==1,expression(Phylogenetic~beta[total]),ifelse(j==2,expression(Phylogenetic~beta[turn]),expression(Phylogenetic~beta[ab])))
      x_lab_corr <- ifelse(j==1,expression(Taxonomic~beta[total]),ifelse(j==2,expression(Taxonomic~beta[turn]),expression(Taxonomic~beta[ab])))
    }
    
    r2 <- r2(m)
    model_result <- summary(m)$coefficients
    no_SC_diff <- unique(subset.beta$scaled.SC.diff[subset.beta$SC.diff.mean == 0])
    predict.df <- as.data.frame(ggemmeans(m,terms=c("LC"),type="fe",ci.lvl=0.95,condition=c(abundanceMethod = "Abundance",SC.diff.mean = 0)))
    
    natural_ref <- predict.df$predicted[predict.df$x == "Natural"]
    predict.df$percent <- ((predict.df$predicted-natural_ref)/natural_ref)*100
    predict.df$low_CI <- ((predict.df$conf.low-natural_ref)/natural_ref)*100
    predict.df$high_CI <- ((predict.df$conf.high-natural_ref)/natural_ref)*100
    
    temp_results_reg <- data.frame(model=ifelse(i<3,"Observed","T_controlled"),LCAgriculture_p=model_result[2,5],LCUrban_p=model_result[3,5],LCAgr_coef=model_result[2,1],LCUrban_coef=model_result[3,1],
                                   DistxAgr_p = model_result[ifelse(i == 3,11,10),5],DistxUrb_p = model_result[ifelse(i == 3,12,11),5],
                                   DistxAgr_coef = model_result[ifelse(i == 3,11,10),1],DistxUrb_coef = model_result[ifelse(i == 3,12,11),1],
                                   Distance_p = model_result[4,5],Distance_coef = model_result[4,1])
    
    if (nrow(predict.df) == 3) {
      predict.df$p <- c(NA,temp_results_reg$LCAgriculture_p,temp_results_reg$LCUrban_p)
    } else {
      predict.df$p <- c(NA,temp_results_reg$LCAgriculture_p)
    }
    
    predict.df$p <- as.numeric(predict.df$p)
    predict.df$sig <- ifelse(predict.df$p < 0.05, "Sig","Insig")
    
    plot(simulateResiduals(m)) #can ignore the outer Newton message...from qgam https://github.com/florianhartig/DHARMa/issues/300
    
    temp_results_reg$metrics <- ifelse(j == 1, "Beta", ifelse(j==2,"Turnover","Abundance difference"))
    temp_results_reg$facet <- ifelse(i == 1, "Taxonomic", "Phylogenetic")
    predict.df$metrics <-  ifelse(j == 1, "Beta", ifelse(j==2,"Turnover","Abundance difference"))
    predict.df$facet <-  ifelse(i == 1, "Taxonomic", "Phylogenetic")
    predict.df$model <- unique(temp_results_reg$model)
    levels(predict.df$x) <- c("Natural","Agricultural","Urban")
    
    results_reg <- rbind(temp_results_reg,results_reg)
    graph_df <- rbind(graph_df,predict.df)
    
    ylim_min <- ifelse(j==1,-20,ifelse(j==2,-57,-35))
    ylim_max <- ifelse(j==1,20,ifelse(j==2,57,35))
    
    set <- ggemmeans(m,terms=c("LC","scaled.spatial.dist[-2.388183:2.918999796]"),type="fe",ci.lvl=0.95,condition=c(abundanceMethod = "Abundance",SC.diff.mean = 0))
    set$distance <- as.numeric(as.character(set$group))*sd(log(subset.beta$mean.spatial.dist))+mean(log(subset.beta$mean.spatial.dist))
    set$distance <- 10^set$distance
    set <- set[!(set$x == "Urban" & set$distance >180),]
    
    assign(paste0("p_",unique(graph_df$facet),unique(graph_df$metrics)) ,ggplot()+
             geom_line(data=graph_df,aes(x=distance,y=predicted,group=x,colour=x))+
             geom_ribbon(data=graph_df,aes(x=distance,ymin=conf.low,ymax=conf.high,fill=x),alpha=0.1)+
             ylab(ifelse(i==1,paste0(unique(graph_df$metrics)),paste0("")))+
             xlab(ifelse(i==2 & j == 3,paste0("Distance (km)"),paste0("")))+
             ggtitle(ifelse(j==1,paste0(unique(graph_df$facet)),paste0("")))+
             ylim(-0.1,1.05)+
             labs(colour = "Habitat",fill="Habitat")+
             scale_colour_manual(values=c("green4","gold1","grey25"))+
             scale_fill_manual(values=c("green4","gold1","grey25"))+
             theme_classic()+
             theme(plot.title=element_text(hjust=0.5)))
    
    if (length(unique(graph_df$sig)) == 4)
      shape <- c(4,1,19)
    if( length(unique(graph_df$sig)) == 3 & "Sig" %in% graph_df$sig & "Insig" %in% graph_df$sig)
      shape <- c(4,19)
    if(length(unique(graph_df$sig)) == 3 & "Sig" %in% graph_df$sig & "Marginal" %in% graph_df$sig)
      shape <- c(1,19)
    if(length(unique(graph_df$sig)) == 3 & "Insig" %in% graph_df$sig & "Marginal" %in% graph_df$sig)
      shape <- c(4,1)
    if(length(unique(graph_df$sig)) == 2 & "Insig" %in% graph_df$sig )
      shape <- 4
    
    graph_title <- ifelse(i==1, "Taxonomic \n", ifelse(i==2,"Phylogenetic \n","Phylogenetic \n (Controlling for taxonomic)"))
    if (j > 1) {
      graph_title <- ""
    }
    
    fig_lab <- letters[[j+3*(i-1)]]
    y_lab <- ifelse(j==1,expression(beta[total]~difference~"(%)"),ifelse(j==2,expression(beta[turn]~difference~"(%)"),expression(beta[ab]~difference~"(%)")))
    assign(paste0("p_",unique(graph_df$facet),unique(graph_df$metrics),unique(graph_df$model)) ,ggplot()+
             geom_hline(yintercept=0)+
             geom_point(data=graph_df,aes(x=x,y=percent,colour=x,shape=sig),position=position_dodge(width=0.5),size=2)+
             geom_errorbar(data=graph_df,aes(x=x,ymin=low_CI,ymax=high_CI,colour=x),width=0.2,position=position_dodge(width=0.5))+
             annotate('text', label=paste0(" (",fig_lab,")"), x=0, y=Inf, hjust=0, vjust=1,size=4)+
             guides(shape="none")+
             ylim(ylim_min,ylim_max)+
             ylab(ifelse(i==1,y_lab," "))+
             xlab("")+
             ggtitle(graph_title)+
             labs(colour = "Habitat",fill="Habitat")+
             scale_colour_manual(values=c("green4","darkgoldenrod2","grey25"))+
             scale_fill_manual(values=c("green4","darkgoldenrod2","grey25"))+
             scale_shape_manual(values=shape,na.value=19)+
             theme_classic()+
             theme(plot.title=element_text(size=10,hjust=0.5),
                   axis.text.x = element_text(angle = 10,size=7.5))
    )
    y_lab <- ifelse(j==1,expression(beta[total]),ifelse(j==2,expression(beta[turn]),expression(beta[ab])))
    
    assign(paste0("p_dist_",unique(graph_df$facet),unique(graph_df$metrics),unique(graph_df$model)),ggplot()+
             geom_line(data=set,aes(x=distance,y=predicted,group=x,colour=x))+
             geom_ribbon(data=set,aes(x=distance,ymin=conf.low,ymax=conf.high,fill=x),alpha=0.1)+
             ylab(ifelse(i==1,y_lab," "))+
             xlab(ifelse(i==1 & j==1,"Distance (km)"," "))+
             labs(colour = "Habitat",fill="Habitat")+
             annotate('text', label=paste0(" (",fig_lab,")"), x=0, y=Inf, hjust=0, vjust=1,size=4)+
             scale_colour_manual(values=c("green4","darkgoldenrod2","grey25"),labels=c("Natural","Agricultural","Urban"))+
             scale_fill_manual(values=c("green4","darkgoldenrod2","grey25"),labels=c("Natural","Agricultural","Urban"))+
             scale_x_continuous(trans="log")+
             theme_classic())
    
    
    all_results_reg <- rbind(all_results_reg,results_reg) 
    all_graph_df <- rbind(all_graph_df,graph_df)
    all_table <- rbind(all_table,round(data.frame(summary(m)$coefficients, 
                                                  R2c = r2$R2_conditional,
                                                  R2m = r2$R2_marginal),3))
  }
}

write.csv(all_table,"Table/Beta_table.csv")
library(ggpubr)

p <- ggarrange(p_TaxonomicBetaObserved,p_PhylogeneticBetaObserved,p_PhylogeneticBetaT_controlled,
               p_TaxonomicTurnoverObserved,p_PhylogeneticTurnoverObserved,p_PhylogeneticTurnoverT_controlled,
               `p_TaxonomicAbundance differenceObserved`,`p_PhylogeneticAbundance differenceObserved`,`p_PhylogeneticAbundance differenceT_controlled`,
               ncol=3,nrow=3,
               common.legend=T,legend="bottom")
plot(p)

ggsave("Figure/p_beta.tiff",dpi=600,compression="lzw",width=16.8,height=16.8,units="cm")


### Proportion
cor_df_Btotal <- data.frame(TBtotal=P_beta$T_beta[P_beta$metrics == "Btotal"],PBtotal=P_beta$P_beta[P_beta$metrics == "Btotal"])
cor_df_Bturn <- data.frame(TBturn=P_beta$T_beta[P_beta$metrics == "Bturn"],PBturn=P_beta$P_beta[P_beta$metrics == "Bturn"])
cor_df_Bab <- data.frame(TBab=P_beta$T_beta[P_beta$metrics == "Bab"],PBab=P_beta$P_beta[P_beta$metrics == "Bab"])

prop_df <- cbind(cor_df_Bab,cor_df_Bturn,cor_df_Btotal,LC=subset.beta$LC,study=subset.beta$study,scaled.spatial.dist=subset.beta$scaled.spatial.dist,SC.diff.mean=subset.beta$SC.diff.mean,abundanceMethod=subset.beta$abundanceMethod,mean_PC1=subset.beta$mean_PC1,mean_PC2=subset.beta$mean_PC2)
prop_df$Tprop <- prop_df$TBturn/prop_df$TBtotal
prop_df$Pprop <- prop_df$PBturn/prop_df$PBtotal

m_prop_T<- lmer(Tprop~LC*scaled.spatial.dist+SC.diff.mean+abundanceMethod+mean_PC1+mean_PC2+(1|study),data=prop_df)
summary(m_prop_T)

m_prop_P<- lmer(Pprop~LC*scaled.spatial.dist+SC.diff.mean+abundanceMethod+mean_PC1+mean_PC2+(1|study),data=prop_df)
summary(m_prop_P)

r2T <- r2(m_prop_T)
r2P <- r2(m_prop_P)

predict_df_prop_T <- ggpredict(m_prop_T,terms="LC",type="fe",ci.lvl=0.95,condition=c(abundanceMethod = "Abundance",SC.diff.mean = 0))
predict_df_prop_P <- ggpredict(m_prop_P,terms="LC",type="fe",ci.lvl=0.95,condition=c(abundanceMethod = "Abundance",SC.diff.mean = 0))

predict_df_prop <- rbind(predict_df_prop_T,predict_df_prop_P)
predict_df_prop$metrics <- c(rep("Taxonomic",3),rep("Phylogenetic",3))
predict_df_prop$predicted <- predict_df_prop$predicted * 100
predict_df_prop$conf.low <- predict_df_prop$conf.low * 100
predict_df_prop$conf.high <- predict_df_prop$conf.high * 100
predict_df_prop$sig <- NA
predict_df_prop$sig[c(2,4)] <- "Sig"
predict_df_prop$sig[c(3,6)] <- "Insig"

predict_df_prop$metrics <- factor(predict_df_prop$metrics,levels=c("Taxonomic","Phylogenetic"))
p_prop <- ggplot(data=predict_df_prop)+
  geom_point(aes(x=metrics,y=predicted,colour=x,shape=sig),position=position_dodge(width=0.5))+
  geom_errorbar(aes(x=metrics,ymin=conf.low,ymax=conf.high,colour=x),position=position_dodge(width=0.5),width=0.1)+
  ylab(expression("%"~of~beta[turn]))+
  xlab("")+
  labs(colour="Habitat")+
  guides(shape="none")+
  scale_colour_manual(values=c("green4","darkgoldenrod2","grey25"),labels=c("Natural","Agricultural","Urban"))+
  scale_fill_manual(values=c("green4","darkgoldenrod2","grey25"),labels=c("Natural","Agricultural","Urban"))+
  scale_shape_manual(values= c(4,19),na.value=19)+
  theme_classic()+
  theme(legend.position="bottom")
plot(p_prop)

df_prop_T <- round(data.frame(summary(m_prop_T)$coefficients,R2c=r2T$R2_conditional,R2m=r2T$R2_marginal),3)
df_prop_P <- round(data.frame(summary(m_prop_P)$coefficients,R2c=r2P$R2_conditional,R2m=r2P$R2_marginal),3)
prop_table <- rbind(df_prop_T,df_prop_P)
write.csv(prop_table,"Table/prop_table.csv")
ggsave("Figure/beta_prop.tiff",dpi=600,units="cm",width=12.8,height=12.8,compression="lzw")
