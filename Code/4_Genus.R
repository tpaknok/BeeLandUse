library(tidyverse)
library(BBmisc)
library(reshape2)
library(ape)
source("Code/get_R2.R")

genustree <- read.tree("Data/Tree.tre")
sp.long.df <- read.csv("Data/Cleaned_Genus.csv")

library(spaMM)
response_df_full <- plot_df <- NULL
for (j in 1:length(unique(sp.long.df$genus))) {
  newdata = NULL
  message(j,"/",length(unique(sp.long.df$genus)),unique(sp.long.df$genus)[[j]])
  
  sp.long.df.subset <- subset(sp.long.df,genus == unique(sp.long.df$genus)[[j]])
  remove_study <- names(which(by(paste0(sp.long.df.subset$Long,"_",sp.long.df.subset$Lat),sp.long.df.subset$studyID,function(x) length(unique(x))) == 1))
  const <- 1

  sp.long.df.subset <- sp.long.df.subset[!sp.long.df.subset$studyID %in% remove_study,]
  tl_abundance <- aggregate(orig_value~studyID,data=sp.long.df.subset,FUN=sum)
  tl_abundance <- tl_abundance[tl_abundance$orig_value >= 10, ]
  
  if(nrow(tl_abundance) < 5) {
    message("no study available")
    next
  }
  unique.study <- unique(sp.long.df.subset$studyID)
  unique.study.agr <- unique(sp.long.df.subset$studyID[sp.long.df.subset$LC.author.coarse == "Agriculture"])
  unique.study.urb <- unique(sp.long.df.subset$studyID[sp.long.df.subset$LC.author.coarse == "Urban"])
  unique.study.nat <- unique(sp.long.df.subset$studyID[sp.long.df.subset$LC.author.coarse == "Natural"])
  
  LC.available <- c(length(unique.study.nat),length(unique.study.agr),length(unique.study.urb))
  names(LC.available) <- c("Natural","Agriculture","Urban")
  LC.available <- LC.available[LC.available >= 5]
  
  if (length(LC.available) <= 1 | (!("Natural" %in% names(LC.available))) | length(unique.study) < 5) {
    message("Genus ",unique(sp.long.df$genus)[[j]]," not analyzed")
    response <- data.frame(Taxa = unique(sp.long.df$genus)[[j]],
                           Mean = mean(sp.long.df.subset$value),
                           SD = sd(sp.long.df.subset$value),
                           Int= NA,
                           Agr = NA, 
                           Urb = NA,
                           Agr_p = NA,
                           Urb_p = NA,
                           R2m = NA,
                           R2c = NA,
                           Time = NA,
                           form = 4)
    form=4
  } else {
    Natural.pos <- which(names(LC.available) == "Natural")
    
    if (is_empty(Natural.pos)) {
      LC.available <- LC.available
    } else {
      LC.available <- LC.available[-Natural.pos]
    }
    
    LC.available <- names(LC.available)
    RHS <- paste(LC.available,collapse="+")
    
    if (length(unique(sp.long.df.subset$abundanceMethod)) > 1) {
      formula <- paste0("value~LC.author.coarse+abundanceMethod+PC1+PC2+offset(log(SamplingEffort))+(",RHS,"||studyID)+Matern(1|Long+Lat %in% studyID)")
      form = 1
    } 
    
    if (length(unique(sp.long.df.subset$abundanceMethod)) == 1) {
    formula <- paste0("value~LC.author.coarse+PC1+PC2+offset(log(SamplingEffort))+(",RHS,"||studyID)+Matern(1|Long+Lat %in% studyID)")
    form = 2
    }

        sp.long.df.subset$PC1 <- scale(sp.long.df.subset$PC1)
    sp.long.df.subset$PC2 <- scale(sp.long.df.subset$PC2)
    
    start <- Sys.time()
    mod_genus <- try(fitme(as.formula(formula),data=sp.long.df.subset,family=negbin,method="ML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth")))
    end <- Sys.time()
    b_table <- summary(mod_genus,details=T)$beta_table
    
    aov_table <- NA
    if (length(unique(LC.available)) == 2) {
    aov_table <- anova(mod_genus,method="t-Chisq")
    
    p_df <- c(NA,NA)
    names(p_df) <- c("LC.author.coarseAgriculture","LC.author.coarseUrban")
    if (aov_table$`Pr(>Chisq.)`[2] < 0.05) {
      sp.long.df.subset$LC.author.coarse1 <- relevel(sp.long.df.subset$LC.author.coarse,"Agriculture")
      mod_genus1 <- try(fitme(value~LC.author.coarse1+PC1+PC2+offset(log(SamplingEffort))+(Natural+Urban||studyID)+Matern(1|Long+Lat %in% studyID),
                              data=sp.long.df.subset,
                              family=negbin,
                              method="ML",
                              verbose=c(TRACE=TRUE),
                              control.dist=list(dist.method="Earth")))
      
      p_df <- p.adjust(c(b_table[2:3,4],summary(mod_genus1,details=T,verbose=F)$beta_table[3,4]),"fdr")
      }
    } else {
      p_df <- b_table[2,4]
      names(p_df) <- rownames(b_table)[2]
    }
    
    newdata = expand.grid(LC.author.coarse=c("Natural","Agriculture","Urban"),abundanceMethod="Abundance",PC1=mean(sp.long.df.subset$PC1),PC2=mean(sp.long.df.subset$PC2),SamplingEffort=1)

    if (!"Urban" %in% LC.available) 
      newdata = subset(newdata,LC.author.coarse != "Urban")
    if (!"Agriculture" %in% LC.available) 
      newdata = subset(newdata,LC.author.coarse != "Agriculture")
    
    newdata$predicted_y <- as.numeric(predict(mod_genus,newdata=newdata,re.form=NA,type="response"))
    newdata$predicted_y_link <- as.numeric(predict(mod_genus,newdata=newdata,re.form=NA,type="link"))
  
    newdata <- cbind(newdata,get_intervals(object=mod_genus,newdata=newdata,level=0.95,re.form=NA,type="link"),get_intervals(object=mod_genus,newdata=newdata,level=0.68,re.form=NA,type="link"))
    newdata$indicator <- 0
    newdata$LC <- c("Natural",unique(LC.available))
    natural_ref <- c(newdata[newdata$LC == "Natural","predicted_y"])
    newdata$ratio_mean <- newdata$predicted_y/natural_ref
    newdata$ratio_low <- newdata$predVar_0.025/natural_ref
    newdata$ratio_high <- newdata$predVar_0.975/natural_ref
    
    newdata$Taxa = unique(sp.long.df$genus)[[j]]
    
    plot_df <- rbind(plot_df,newdata)

        if (is.error(mod_genus)){
      message("error model")
      response <- data.frame(Taxa = unique(sp.long.df.subset$genus)[[j]],
                             Mean = mean(sp.long.df.subset$value),
                             SD = sd(sp.long.df.subset$value),
                             Int= NA,
                             Agr = NA, 
                             Urb = NA,
                             Agr_p = NA,
                             Urb_p = NA,
                             R2m = NA,
                             R2c = NA,
                             Time = NA,
                             form = 5)
    } else {
      
      R2 <-  get_R2(mod_genus)
      response <- data.frame(Taxa = unique(sp.long.df$genus)[[j]],
                             Mean = mean(sp.long.df.subset$value),
                             SD = sd(sp.long.df.subset$value),
                             Int = fixef(mod_genus)[[1]],
                             Agr=ifelse("Agriculture" %in% LC.available,fixef(mod_genus)[names(fixef(mod_genus)) == "LC.author.coarseAgriculture"],NA),
                             Urb=ifelse("Urban" %in% LC.available,fixef(mod_genus)[names(fixef(mod_genus)) == "LC.author.coarseUrban"],NA),
                             Agr_se=ifelse("Agriculture" %in% LC.available,b_table[rownames(b_table) == "LC.author.coarseAgriculture","Cond. SE"],NA),
                             Urb_se=ifelse("Urban" %in% LC.available,b_table[rownames(b_table) == "LC.author.coarseUrban","Cond. SE"],NA),
                             Global_LC_p = ifelse(length(unique(LC.available)) == 2,aov_table$`Pr(>Chisq.)`[2],NA),                  
                             Agr_p=ifelse("Agriculture" %in% LC.available,p_df[names(p_df) == "LC.author.coarseAgriculture"],NA),
                             Urb_p=ifelse("Urban" %in% LC.available,p_df[names(p_df) == "LC.author.coarseUrban"],NA),
                             R2m = R2[[1]],
                             R2c = R2[[2]],
                             Time = as.numeric(end-start),
                             form = form)
    }
    
    response_df_full <- rbind(response,response_df_full)
    print(response)
  }
}

response_df <- response_df_full[!is.na(response_df_full$Agr),]
response_df$adjusted_p_agr <- p.adjust(response_df$Agr_p,"fdr")
response_df$adjusted_p_urb <- p.adjust(response_df$Urb_p,"fdr")

response_plot_df_agr <- subset(plot_df,LC.author.coarse == "Agriculture")
response_df$ratio_mean_agr <- response_plot_df_agr[match(response_df$Taxa,response_plot_df_agr$Taxa),"ratio_mean"]
response_df$log_ratio_mean_agr <- log10(response_df$ratio_mean_agr)
response_df$ratio_propoagated_se_agr <- response_plot_df_agr[match(response_df$Taxa,response_plot_df_agr$Taxa),"ratio_propoagated_se2"]
  
response_plot_df_urb <- subset(plot_df,LC.author.coarse == "Urban")
response_df$ratio_mean_urb <- response_plot_df_urb[match(response_df$Taxa,response_plot_df_urb$Taxa),"ratio_mean"]
response_df$log_ratio_mean_urb <- log10(response_df$ratio_mean_urb)
response_df$ratio_propoagated_se_urb <- response_plot_df_urb[match(response_df$Taxa,response_plot_df_urb$Taxa),"ratio_propoagated_se2"]

agr_response_df <- response_df[!is.na(response_df$log_ratio_mean_agr),]

genus_tree <- genustree
genus_tree$tip.label <- word(genus_tree$tip.label,1,sep="_")
tree <- genus_tree

pruned.tree.agr <- drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% agr_response_df$Taxa])
agr <- agr_response_df$log_ratio_mean_agr
names(agr) <- agr_response_df$Taxa
phylosig(pruned.tree.agr,agr,method="lambda",test=T)

urb_response_df <- response_df[!is.na(response_df$Urb),]

pruned.tree.urb <- drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% urb_response_df$Taxa])
urb <- urb_response_df$log_ratio_mean_urb

names(urb) <- urb_response_df$Taxa

phylosig(pruned.tree.urb,urb,method="lambda",test=T)

genus_agr <- agr_response_df
genus_urb <- urb_response_df

library(ggtree)
library(ggnewscale)

agr_response_df$sig <- ifelse(agr_response_df$adjusted_p_agr< 0.05,"sig","insig")
urb_response_df$sig <- ifelse(urb_response_df$adjusted_p_Urb < 0.05,"sig","insig")

pruned.tree.agr <- drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% agr_response_df$Taxa])
sig_gen <- agr_response_df$Taxa[agr_response_df$sig == "sig"]
agr_response_df_tree <- agr_response_df

log_ratio_mean_c <- agr_response_df_tree$ratio_mean_agr
names(log_ratio_mean_c) <- agr_response_df_tree$Taxa
agr_an <- phytools::fastAnc(pruned.tree.agr, log_ratio_mean_c, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(pruned.tree.agr,agr_response_df_tree$Taxa),trait = agr_response_df_tree$ratio_mean_agr)
nd <- data.frame(node = names(agr_an$ace),trait = agr_an$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
pruned.tree.agr <- full_join(pruned.tree.agr, d, by = 'node')

sig <- rep("insig",nrow(pruned.tree.agr@data))
sig[which(pruned.tree.agr@phylo$tip.label %in% subset(agr_response_df,sig=="sig" & Agr > 0)$Taxa)] <- "Positive"
sig[which(pruned.tree.agr@phylo$tip.label %in% subset(agr_response_df,sig=="sig" & Agr < 0)$Taxa)] <- "Negative"

pruned.tree.agr@data$sig <- sig

p_tree_agr <- ggtree(pruned.tree.agr,layout="circular",ladderize=T,size=1.2,hang=0.1) %<+% agr_response_df_tree+
  geom_tree(aes(color=trait),size=1)+
  scale_color_gradient2(high="#ca0020",mid="#ffffbf",low="#0571b0",trans="log10",breaks=c(0.1,0.3,1,3),midpoint=0,limits=c(0.06,6.6),name=expression(paste("Ratio (Agricultural or Urban / Natural)")))+
  new_scale_color()+
  geom_tiplab(aes(colour=sig),fontface=3,align=T,linetype=NA, hjust = -0.05 ,show.legend=F,size=2.8)+
  scale_color_manual(values=c("black","#0571b0","red"))+
  theme(legend.text = element_text(size = 6), 
        legend.title = element_text(size =8), 
        legend.key.size = unit(0.35, 'cm'),legend.position="none",
        legend.box.spacing= unit(-1, 'cm'),
        plot.margin = unit(c(1,1,1,1), "cm"))

plot(p_tree_agr)  
ggsave("C:/Users/pakno/OneDrive - University of Toronto/GRF Bee/Figure/agr_tree.tiff",height=6,width=6,compression="lzw",bg="white",dpi=600)

pruned.tree.urb <- drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% urb_response_df$Taxa])
urb_response_df$sig <- ifelse(urb_response_df$adjusted_p_urb < 0.05,"sig","insig")
sig_gen <- urb_response_df$Taxa[urb_response_df$sig == "sig"]
urb_response_df_tree <- urb_response_df

Urb_change_c <- urb_response_df_tree$ratio_mean_urb
names(Urb_change_c) <- urb_response_df_tree$Taxa
Urb_an <- phytools::fastAnc(pruned.tree.urb, Urb_change_c, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(pruned.tree.urb,urb_response_df_tree$Taxa),trait = urb_response_df_tree$ratio_mean_urb)
nd <- data.frame(node = names(Urb_an$ace),trait = Urb_an$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
pruned.tree.urb <- full_join(pruned.tree.urb, d, by = 'node')

sig <- rep("insig",nrow(pruned.tree.urb@data))
sig[which(pruned.tree.urb@phylo$tip.label %in% subset(urb_response_df,sig=="sig" & Urb < 0)$Taxa)] <- "Negative"
sig[which(pruned.tree.urb@phylo$tip.label %in% subset(urb_response_df,sig=="sig" & Urb > 0)$Taxa)] <- "Positive"

pruned.tree.urb@data$sig <- sig

p_tree_urb <- ggtree(pruned.tree.urb,layout="circular",ladderize=T,size=1.2) %<+% urb_response_df_tree+
  geom_tree(aes(colour=trait),size=1)+
  scale_color_gradient2(high="#ca0020",mid="#ffffbf",low="#0571b0",trans="log10",breaks=c(0.1,0.3,1,3),limits=c(0.07,6.5),name=expression(paste("Ratio (Agricultural or Urban / Natural)")))+
  new_scale_color()+
  geom_tiplab(aes(colour=sig),fontface=3,align=T,linetype=NA, hjust = -0.05 ,show.legend=F,size=2.8)+
  scale_color_manual(values=c("black","#0571b0","red"))+
  theme(legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        legend.key.size = unit(0.5, 'cm'),
        legend.position="bottom",
        legend.box.spacing= unit(0, 'cm'),
        plot.margin = unit(c(1,1,1,1), "cm"))
plot(p_tree_urb)  
ggsave("C:/Users/pakno/OneDrive - University of Toronto/GRF Bee/Figure/urb_tree.tiff",height=6,width=6,compression="lzw",bg="white",dpi=600)

library(ggpubr)

write.csv(response_df,"C:/Users/pakno/OneDrive - University of Toronto/GRF Bee/Table/response_df.csv")
########correlation

cor_df <- agr_response_df[!is.na(agr_response_df$Urb),]
pruned_genus_tree <- keep.tip(genus_tree,cor_df$Taxa)

library(phyr)
cor_phylo(variates=~log_ratio_mean_urb+log_ratio_mean_agr,
          data=cor_df,
          phy=pruned_genus_tree,
          genus=cor_df$Taxa)

cor(cor_df$log_ratio_mean_agr,cor_df$log_ratio_mean_urb)

###bar plot for figure S2

response_df$Family <- sapply(1:nrow(response_df), function(x) getmode(bee[which(grepl(response_df$Taxa[[x]],bee$animalID)),"family"]))

barplot_df <- rbind(response_df[,c("Taxa","Family")],response_df[,c("Taxa","Family")])
barplot_df$LC <- c(rep("Agricultural",nrow(response_df)),rep("Urban",nrow(response_df)))
barplot_df$log_ratio <- c(response_df$log_ratio_mean_agr,response_df$log_ratio_mean_urb)
barplot_df$Simp_Taxa <- substr(barplot_df$Taxa,1,3)
barplot_df$Simp_Taxa[barplot_df$Taxa == "Augochloropsis"] <- "Aug3"
barplot_df$Simp_Taxa[barplot_df$Taxa == "Augochlorella"] <- "Aug2"
barplot_df$Simp_Taxa[barplot_df$Taxa == "Augochlora"] <- "Aug1"

barplot_df$Simp_Taxa[barplot_df$Taxa == "Anthidium"] <- "Ant1"
barplot_df$Simp_Taxa[barplot_df$Taxa == "Anthophora"] <- "Ant2"

barplot_df$Simp_Taxa[barplot_df$Taxa == "Melitta"] <- "Mel3"
  barplot_df$Simp_Taxa[barplot_df$Taxa == "Melissodes"] <- "Mel2"
barplot_df$Simp_Taxa[barplot_df$Taxa == "Melipona"] <- "Mel1"

barplot_df$Simp_Taxa[barplot_df$Taxa == "Nomia"] <- "Nom2"
barplot_df$Simp_Taxa[barplot_df$Taxa == "Nomada"] <- "Nom1"

barplot_df$Simp_Taxa[barplot_df$Taxa == "Panurginus"] <- "Pan1"
barplot_df$Simp_Taxa[barplot_df$Taxa == "Panurgus"] <- "Pan2"

barplot_df$Simp_Taxa[barplot_df$Taxa == "Partamona"] <- "Par2"
barplot_df$Simp_Taxa[barplot_df$Taxa == "Paratrigona"] <- "Par1"

barplot_df$Simp_Taxa[barplot_df$Taxa == "Triepeolus"] <- "Tri2"
barplot_df$Simp_Taxa[barplot_df$Taxa == "Trigona"] <- "Tri1"

urb_na_genus <- barplot_df[is.na(barplot_df$log_ratio),"Simp_Taxa"]
barplot_df$label_colour <- barplot_df$Simp_Taxa %in% urb_na_genus 

family_plot1 <- ggplot(data=subset(barplot_df,Family=="Apidae" | Family == "Andrenidae"),aes(x=Simp_Taxa,y=log_ratio))+
  expand_limits(y=-1.35)+
  geom_bar(aes(fill=LC),stat="identity",position = "dodge2")+
  geom_text(aes(colour=label_colour,label=Simp_Taxa),y=-1.25,vjust=0.7,angle = 20,size=3.5)+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=-1.2)+
  geom_segment(aes(x=0, y=0.75, xend=0, yend=-1.20))+
  scale_colour_manual(values=c("black","blue"),guide=F)+
  scale_fill_manual(values=c("darkgoldenrod2","grey25"),name = "Land use")+
  facet_grid(.~Family,scales="free_x",space="free")+
  ylab("Log ratio")+
  xlab("")+
  theme_classic2()+
  theme(legend.position="bottom",
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=10))


plot(family_plot1)

family_plot2 <- ggplot(data=subset(barplot_df,Family!="Apidae" & Family != "Andrenidae"),aes(x=Simp_Taxa,y=log_ratio))+
  expand_limits(y=-1.35)+
  geom_bar(aes(fill=LC),stat="identity",position = "dodge2")+
  geom_text(aes(colour=label_colour,label=Simp_Taxa),y=-1.25,vjust=0.7,angle = 20,size=3.5)+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=-1.2)+
  geom_segment(aes(x=0, y=0.75, xend=0, yend=-1.20))+
  scale_colour_manual(values=c("black","blue"),guide=F)+
  scale_fill_manual(values=c("darkgoldenrod2","grey25"),name = "Land use")+
  facet_grid(.~Family,scales="free_x",space="free")+
  ylab("Log ratio")+
  xlab("Genus")+
  theme_classic2()+
  theme(legend.position="bottom",
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=10))


plot(family_plot2)

combined_plot <- ggarrange(family_plot1,family_plot2,common.legend=T,nrow=2,legend="bottom")
plot(combined_plot)
ggsave("C:/Users/pakno/OneDrive - University of Toronto/GRF Bee/Figure/bar_plot.tiff",height=10,width=10,compression="lzw",bg="white",dpi=600)
