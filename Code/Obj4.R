library(phytools)
library(sjmisc)
library(BBmisc)
response_df_full <- plot_df <- NULL
sp.long.df <- read.csv("Data/cleaned_genus.csv")
genustree <- read.tree("Data/Tree.tre")

for (j in 1:length(unique(sp.long.df$genus))) {
  newdata = NULL
  message(j,"/",length(unique(sp.long.df$genus)),unique(sp.long.df$genus)[[j]])
  
  sp.long.df.subset <- subset(sp.long.df,genus == unique(sp.long.df$genus)[[j]])
  remove_study <- names(which(by(paste0(sp.long.df.subset$Long,"_",sp.long.df.subset$Lat),sp.long.df.subset$studyID,function(x) length(unique(x))) == 1))
  const <- 1

  sp.long.df.subset <- sp.long.df.subset[!sp.long.df.subset$studyID %in% remove_study,]
  tl_abundance <- aggregate(abundance~studyID,data=sp.long.df.subset,FUN=sum)
  tl_abundance <- tl_abundance[tl_abundance$abundance >= 10, ] #first filter: remove studies with total abundance < 10 across habitats
  
  if(nrow(tl_abundance) < 5) {
    message("no study available") #fewer than five studies = stop
    next
  }
  unique.study <- unique(sp.long.df.subset$studyID)
  unique.study.agr <- unique(sp.long.df.subset$studyID[sp.long.df.subset$LC.author.coarse == "Agriculture"])
  unique.study.urb <- unique(sp.long.df.subset$studyID[sp.long.df.subset$LC.author.coarse == "Urban"])
  unique.study.nat <- unique(sp.long.df.subset$studyID[sp.long.df.subset$LC.author.coarse == "Natural"])
  
  LC.available <- c(length(unique.study.nat),length(unique.study.agr),length(unique.study.urb))
  names(LC.available) <- c("Natural","Agriculture","Urban")
  LC.available <- LC.available[LC.available >= 5] #find land use with more than 5 studies
  
  if (length(LC.available) <= 1 | (!("Natural" %in% names(LC.available))) | length(unique.study) < 5) {
    message("Genus ",unique(sp.long.df$genus)[[j]]," not analyzed")
    response <- data.frame(Taxa = unique(sp.long.df$genus)[[j]],
                           Int= NA,
                           Agr = NA, 
                           Urb = NA,
                           Agr_p = NA,
                           Urb_p = NA,
                           Time = NA,
                           form = 4) #stop if one of the three conditions were met (only one LC, natural is not available, number of unique study < 5)
    form <- 4
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
      formula <- paste0("abundance~",RHS,"+abundanceMethod+PC1+PC2+offset(log(SamplingEffort))+(",RHS,"||studyID)+Matern(1|Long+Lat %in% studyID)")
      null_formula_r2m <- paste0(".~1+offset(log(SamplingEffort))+(",RHS,"||studyID)+Matern(1|Long+Lat %in% studyID)")
      null_formula_r2c <- paste0(".~1+offset(log(SamplingEffort))")
      form = 1
    } 
    
    if (length(unique(sp.long.df.subset$abundanceMethod)) == 1) {
      formula <- paste0("abundance~",RHS,"+PC1+PC2+offset(log(SamplingEffort))+(",RHS,"||studyID)+Matern(1|Long+Lat %in% studyID)") #note the use of offset here. so original abundance is the response
      null_formula_r2m <- paste0(".~1+offset(log(SamplingEffort))+(",RHS,"||studyID)+Matern(1|Long+Lat %in% studyID)")
      null_formula_r2c <- paste0(".~1+offset(log(SamplingEffort))")
      form = 2
    }
    
    sp.long.df.subset$PC1 <- scale(sp.long.df.subset$PC1)
    sp.long.df.subset$PC2 <- scale(sp.long.df.subset$PC2)
    
    start <- Sys.time()
    mod_genus <- try(fitme(as.formula(formula),data=sp.long.df.subset,family=negbin,method="ML",verbose=c(TRACE=TRUE),control.dist=list(dist.method="Earth")))
    end <- Sys.time()
    b_table <- summary(mod_genus,details=T)$beta_table
    
    newdata = expand.grid(Agriculture=c(0,1),Urban=c(0,1),abundanceMethod="Abundance",PC1=mean(sp.long.df.subset$PC1),PC2=mean(sp.long.df.subset$PC2),SamplingEffort=1)
    newdata = newdata[!(newdata$Agriculture == 1 & newdata$Urban == 1),]
    
    if (!"Urban" %in% LC.available) 
      newdata = newdata[!newdata$Urban == 1,]
    if (!"Agriculture" %in% LC.available) 
      newdata = newdata[!newdata$Agriculture == 1,]
    
    newdata$predicted_y <- as.numeric(predict(mod_genus,newdata=newdata,re.form=NA,type="response"))
    newdata$predicted_y_link <- as.numeric(predict(mod_genus,newdata=newdata,re.form=NA,type="link"))
    
    newdata <- cbind(newdata,get_intervals(object=mod_genus,newdata=newdata,level=0.95,re.form=NA,type="link"),get_intervals(object=mod_genus,newdata=newdata,level=0.68,re.form=NA,type="link"))
    newdata$indicator <- 0
    newdata$LC <- c("Natural",unique(LC.available))
    natural_ref <- c(newdata[newdata$LC == "Natural","predicted_y"])
    newdata$ratio_mean <- newdata$predicted_y/natural_ref
    newdata$ratio_se <- newdata$predicted_y_link-newdata$predVar_0.16 #mean-low 84CI = SE?
    newdata$ratio_propoagated_se1 <- newdata$ratio_se+newdata$ratio_se[1]
    newdata$ratio_propoagated_se2 <- sqrt(newdata$ratio_se^2+newdata$ratio_se[1]^2)
    newdata$ratio_low <- newdata$predVar_0.025/natural_ref
    newdata$ratio_high <- newdata$predVar_0.975/natural_ref
    
    newdata$Taxa = unique(sp.long.df$genus)[[j]]
    
    plot_df <- rbind(plot_df,newdata)
    
    mgenus_R2c <- pseudoR2(mod_genus,as.formula(null_formula_r2c))
    mgenus_R2m <- pseudoR2(mod_genus,as.formula(null_formula_r2m))
    
    if (is.error(mod_genus)){
      message("error model")
      response <- data.frame(Taxa = unique(sp.long.df.subset$genus)[[j]],
                             Int= NA,
                             Agr = NA, 
                             Urb = NA,
                             Agr_p = NA,
                             Urb_p = NA,
                             Time = NA,
                             R2m = NA,
                             R2c = NA,
                             form = 5)
    } else {
      response <- data.frame(Taxa = unique(sp.long.df$genus)[[j]],
                             Mean = mean(sp.long.df.subset$abundance),
                             SD = sd(sp.long.df.subset$abundance),
                             Int = fixef(mod_genus)[[1]],
                             Agr=ifelse("Agriculture" %in% LC.available,fixef(mod_genus)[names(fixef(mod_genus)) == "Agriculture"],NA),
                             Urb=ifelse("Urban" %in% LC.available,fixef(mod_genus)[names(fixef(mod_genus)) == "Urban"],NA),
                             Agr_se=ifelse("Agriculture" %in% LC.available,b_table[rownames(b_table) == "Agriculture","Cond. SE"],NA),
                             Urb_se=ifelse("Urban" %in% LC.available,b_table[rownames(b_table) == "Urban","Cond. SE"],NA),
                             Agr_p=ifelse("Agriculture" %in% LC.available,b_table[rownames(b_table) == "Agriculture","p-value"],NA),
                             Urb_p=ifelse("Urban" %in% LC.available,b_table[rownames(b_table) == "Urban","p-value"],NA),
                             Time = as.numeric(end-start),
                             R2m = mgenus_R2m,
                             R2c = mgenus_R2c,
                             form = form)
    }
    
    response_df_full <- rbind(response,response_df_full)
    print(response)
  }
}

response_df <- response_df_full[!is.na(response_df_full$Agr),]
response_df$adjusted_p_agr <- p.adjust(response_df$Agr_p,"fdr")
response_df$adjusted_p_urb <- p.adjust(response_df$Urb_p,"fdr")

response_plot_df_agr <- subset(plot_df,Agriculture==1)
response_df$ratio_mean_agr <- response_plot_df_agr[match(response_df$Taxa,response_plot_df_agr$Taxa),"ratio_mean"]
response_df$log_ratio_mean_agr <- log10(response_df$ratio_mean_agr)
response_df$ratio_propoagated_se_agr <- response_plot_df_agr[match(response_df$Taxa,response_plot_df_agr$Taxa),"ratio_propoagated_se1"]

response_plot_df_urb <- subset(plot_df,Urban==1)
response_df$ratio_mean_urb <- response_plot_df_urb[match(response_df$Taxa,response_plot_df_urb$Taxa),"ratio_mean"]
response_df$log_ratio_mean_urb <- log10(response_df$ratio_mean_urb)
response_df$ratio_propoagated_se_urb <- response_plot_df_urb[match(response_df$Taxa,response_plot_df_urb$Taxa),"ratio_propoagated_se1"]

agr_response_df <- response_df[!is.na(response_df$log_ratio_mean_agr),]

genus_tree <- genustree
genus_tree$tip.label <- word(genus_tree$tip.label,1,sep="_")
tree <- genus_tree

pruned.tree.agr <- drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% agr_response_df$Taxa])
agr <- agr_response_df$log_ratio_mean_agr
names(agr) <-  agr_response_df$Taxa
phylosig(pruned.tree.agr,agr,method="lambda",test=T) #doesn't converge with se...very different results across runs

urb_response_df <- response_df[!is.na(response_df$Urb),]

pruned.tree.urb <- drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% urb_response_df$Taxa])
urb <- urb_response_df$log_ratio_mean_urb

names(urb) <- urb_response_df$Taxa

phylosig(pruned.tree.urb,urb,method="lambda",test=T) #doesn't converge with se...very different results across runs

agr_response_df$adjusted_p <- p.adjust(agr_response_df$Agr_p,"fdr")
urb_response_df$adjusted_p <- p.adjust(urb_response_df$Urb_p,"fdr")

genus_agr <- agr_response_df
genus_urb <- urb_response_df

library(ggtree)
library(ggnewscale)
#agr_response_df$Agr_change <- agr_response_df$Agr_change*100
#urb_response_df$Urb_change <- urb_response_df$Urb_change*100

agr_response_df$sig <- ifelse(agr_response_df$adjusted_p < 0.05,"sig","insig")
urb_response_df$sig <- ifelse(urb_response_df$adjusted_p < 0.05,"sig","insig")

#getMRCA(pruned.tree,agr_response_df$Taxa[agr_response_df$Family == "Halictidae"])
#getMRCA(pruned.tree,agr_response_df$Taxa[agr_response_df$Family == "Melittidae"])
#getMRCA(pruned.tree,agr_response_df$Taxa[agr_response_df$Family == "Andrenidae"])
#getMRCA(pruned.tree,agr_response_df$Taxa[agr_response_df$Family == "Megachilidae"])
#getMRCA(pruned.tree,agr_response_df$Taxa[agr_response_df$Family == "Colletidae"])
#getMRCA(pruned.tree,agr_response_df$Taxa[agr_response_df$Family == "Apidae"])

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
sig[which(pruned.tree.agr@phylo$tip.label %in% subset(agr_response_df,sig=="sig")$Taxa)] <- "sig"
sig[sig=="sig" & pruned.tree.agr@data$trait <1] <- "neg"
sig[sig=="sig" & pruned.tree.agr@data$trait >1] <- "pos"

pruned.tree.agr@data$sig <- sig

p_tree_agr <- ggtree(pruned.tree.agr,layout="circular",ladderize=T,size=1.2,hang=0.1) %<+% agr_response_df_tree+
  geom_tree(aes(color=trait),size=1)+
  scale_color_gradient2(high="#ca0020",mid="#ffffbf",low="#0571b0",trans="log10",breaks=c(0.1,0.3,1,3),midpoint=0,limits=c(0.06,6.6),name=expression(paste("Ratio (Agricultural or Urban / Natural)")))+
  new_scale_color()+
  geom_tiplab(aes(colour=sig),fontface=3,align=T,linetype=NA, hjust = -0.05 ,show.legend=F,size=2.8)+
  scale_color_manual(values=c("black","#0571b0","#ca0020"))+
  theme(legend.text = element_text(size = 6), 
        legend.title = element_text(size =8), 
        legend.key.size = unit(0.35, 'cm'),legend.position="none",
        legend.box.spacing= unit(-1, 'cm'),
        plot.margin = unit(c(1,1,1,1), "cm"))

plot(p_tree_agr)  
ggsave("Figure/agr_tree.tiff",height=6,width=6,compression="lzw",bg="white",dpi=600)

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
sig[which(pruned.tree.urb@phylo$tip.label %in% subset(urb_response_df,sig=="sig")$Taxa)] <- "sig"
sig[sig=="sig" & pruned.tree.urb@data$trait <1] <- "neg"
sig[sig=="sig" & pruned.tree.urb@data$trait >1] <- "pos"

pruned.tree.urb@data$sig <- sig
#aes(label=paste0("italic(",label,")")),parse=T,

p_tree_urb <- ggtree(pruned.tree.urb,layout="circular",ladderize=T,size=1.2) %<+% urb_response_df_tree+
  geom_tree(aes(colour=trait),size=1)+
  scale_color_gradient2(high="#ca0020",mid="#ffffbf",low="#0571b0",trans="log10",breaks=c(0.1,0.3,1,3),limits=c(0.07,6.5),name=expression(paste("Ratio (Agricultural or Urban / Natural)")))+
  new_scale_color()+
  geom_tiplab(aes(colour=sig),fontface=3,align=T,linetype=NA, hjust = -0.05 ,show.legend=F,size=2.8)+
  scale_color_manual(values=c("black","#0571b0","#ca0020"))+
  theme(legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        legend.key.size = unit(0.5, 'cm'),
        legend.position="bottom",
        legend.box.spacing= unit(0, 'cm'),
        plot.margin = unit(c(1,1,1,1), "cm"))
plot(p_tree_urb)  
ggsave("Figure/urb_tree.tiff",height=6,width=6,compression="lzw",bg="white",dpi=600)

library(ggpubr)

#p_tree_combined <- ggarrange(p_tree_agr,p_tree_urb,
#labels=c("(a) Agricultural","(b) Urban"),vjust=1,hjust=0,
#nrow=2,common.legend=F,legend="bottom") #seems that ggarrange can't work with ggnewscale

#ggsave("G:/My Drive/GRF Bee/Figure/tree.tiff",height=24,width=11,units="cm",compression="lzw",bg="white")

write.csv(response_df,"Table/response_df.csv")
########correlation

cor_df <- na.omit(agr_response_df)
pruned_genus_tree <- keep.tip(genus_tree,cor_df$Taxa)

library(phyr)
cor_phylo(variates=~log_ratio_mean_urb+log_ratio_mean_agr,
          data=cor_df,
          phy=pruned_genus_tree,
          genus=cor_df$Taxa)