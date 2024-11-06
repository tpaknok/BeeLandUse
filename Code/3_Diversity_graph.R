library(ggplot2)

### Alpha diversity prediction ##
newdataTD = expand.grid(LC.author.coarse=c("Natural","Agriculture","Urban"),abundanceMethod="Abundance",PC1=mean(analysis.df.PD$PC1),PC2=mean(analysis.df.PD$PC2))

newdataTD$predicted_y <- predict(mTD2,newdata=newdataTD,re.form=NA,type="response")
newdataTD <- cbind(newdataTD,get_intervals(object=mTD2,newdata=newdataTD,levels=0.95,re.form=NA,type="response"))
newdataTD$LC <- c("Natural","Agricultural","Urban")
newdataTD$LC <- factor(newdataTD$LC, levels = c("Natural","Agricultural","Urban"))
newdataTD[,c("predicted_y","predVar_0.025","predVar_0.975")] <- exp(newdataTD[,c("predicted_y","predVar_0.025","predVar_0.975")])
natural_ref <- c(newdataTD[newdataTD$LC == "Natural","predicted_y"])
newdataTD$percent_mean <- (c(newdataTD$predicted_y)-natural_ref)/natural_ref*100
newdataTD$percent_low <- (c(newdataTD$predVar_0.025)-natural_ref)/natural_ref*100
newdataTD$percent_high <- (c(newdataTD$predVar_0.975)-natural_ref)/natural_ref*100
newdataTD$x <- newdataTD$LC

newdataPD = expand.grid(LC.author.coarse=c("Natural","Agriculture","Urban"),abundanceMethod="Abundance",PC1=mean(analysis.df.PD$PC1),PC2=mean(analysis.df.PD$PC2))
newdataPD$predicted_y <- predict(mPD.S2,newdata=newdataPD,re.form=NA,type="response")
newdataPD <- cbind(newdataPD,get_intervals(object=mPD.S2,newdata=newdataPD,levels=0.95,re.form=NA,type="response"))
newdataPD$LC <- c("Natural","Agricultural","Urban")
newdataPD$LC <- factor(newdataPD$LC, levels = c("Natural","Agricultural","Urban"))
newdataPD[,c("predicted_y","predVar_0.025","predVar_0.975")] <- exp(newdataPD[,c("predicted_y","predVar_0.025","predVar_0.975")])
natural_ref <- c(newdataPD[newdataPD$LC == "Natural","predicted_y"])
newdataPD$percent_mean <- (c(newdataPD$predicted_y)-natural_ref)/natural_ref*100
newdataPD$percent_low <- (c(newdataPD$predVar_0.025)-natural_ref)/natural_ref*100
newdataPD$percent_high <- (c(newdataPD$predVar_0.975)-natural_ref)/natural_ref*100
newdataPD$x <- newdataPD$LC

newdataPD_C = expand.grid(LC.author.coarse=c("Natural","Agriculture","Urban"),abundanceMethod="Abundance",PC1=mean(analysis.df.PD$PC1),PC2=mean(analysis.df.PD$PC2),Covariate.Est.Shannon=mean(analysis.df.PD$Covariate.Est.Shannon))
newdataPD_C$predicted_y <- predict(mPD.S.controlled2,newdata=newdataPD_C,re.form=NA,type="response")
newdataPD_C <- cbind(newdataPD_C,get_intervals(object=mPD.S.controlled2,newdata=newdataPD_C,levels=0.95,re.form=NA,type="response"))
newdataPD_C$LC <- c("Natural","Agricultural","Urban")
newdataPD_C$LC <- factor(newdataPD_C$LC, levels = c("Natural","Agricultural","Urban"))
newdataPD_C[,c("predicted_y","predVar_0.025","predVar_0.975")] <- exp(newdataPD_C[,c("predicted_y","predVar_0.025","predVar_0.975")])
natural_ref <- c(newdataPD_C[newdataPD_C$LC == "Natural","predicted_y"])
newdataPD_C$percent_mean <- (c(newdataPD_C$predicted_y)-natural_ref)/natural_ref*100
newdataPD_C$percent_low <- (c(newdataPD_C$predVar_0.025)-natural_ref)/natural_ref*100
newdataPD_C$percent_high <- (c(newdataPD_C$predVar_0.975)-natural_ref)/natural_ref*100
newdataPD_C$x <- newdataPD_C$LC

### Beta #SC_pca = -1.33943234 = when all assemblages have 100% completeness for the complete dataset
library(ggeffects)

newdata_Tbeta <- ggemmeans(m4,"LC",condition=list(abundaceMethod="Abundance",SC_pca=-1.33943234,SamplingEffort.sd=0))
newdata_Tbeta$percent_mean <- ((newdata_Tbeta$predicted)-(newdata_Tbeta$predicted[newdata_Tbeta$x == "Natural"]))/(newdata_Tbeta$predicted[newdata_Tbeta$x == "Natural"])*100
newdata_Tbeta$percent_low <- ((newdata_Tbeta$conf.low)-(newdata_Tbeta$predicted[newdata_Tbeta$x == "Natural"]))/(newdata_Tbeta$predicted[newdata_Tbeta$x == "Natural"])*100
newdata_Tbeta$percent_high <- ((newdata_Tbeta$conf.high)-(newdata_Tbeta$predicted[newdata_Tbeta$x == "Natural"]))/(newdata_Tbeta$predicted[newdata_Tbeta$x == "Natural"])*100

newdata_Pbeta <- ggemmeans(m5,"LC",condition=list(abundaceMethod="Abundance",SC_pca=-1.33943234,SamplingEffort.sd=0))
newdata_Pbeta$percent_mean <- ((newdata_Pbeta$predicted)-(newdata_Pbeta$predicted[newdata_Pbeta$x == "Natural"]))/(newdata_Pbeta$predicted[newdata_Pbeta$x == "Natural"])*100
newdata_Pbeta$percent_low <- ((newdata_Pbeta$conf.low)-(newdata_Pbeta$predicted[newdata_Pbeta$x == "Natural"]))/(newdata_Pbeta$predicted[newdata_Pbeta$x == "Natural"])*100
newdata_Pbeta$percent_high <- ((newdata_Pbeta$conf.high)-(newdata_Pbeta$predicted[newdata_Pbeta$x == "Natural"]))/(newdata_Pbeta$predicted[newdata_Pbeta$x == "Natural"])*100

newdata_Pbeta_C <- ggemmeans(m6,"LC",condition=list(abundaceMethod="Abundance",SC_pca=-1.33943234,SamplingEffort.sd=0))
newdata_Pbeta_C$percent_mean <- ((newdata_Pbeta_C$predicted)-(newdata_Pbeta_C$predicted[newdata_Pbeta_C$x == "Natural"]))/(newdata_Pbeta_C$predicted[newdata_Pbeta_C$x == "Natural"])*100
newdata_Pbeta_C$percent_low <- ((newdata_Pbeta_C$conf.low)-(newdata_Pbeta_C$predicted[newdata_Pbeta_C$x == "Natural"]))/(newdata_Pbeta_C$predicted[newdata_Pbeta_C$x == "Natural"])*100
newdata_Pbeta_C$percent_high <- ((newdata_Pbeta_C$conf.high)-(newdata_Pbeta_C$predicted[newdata_Pbeta_C$x == "Natural"]))/(newdata_Pbeta_C$predicted[newdata_Pbeta_C$x == "Natural"])*100


### Gamma
newdata_Tgamma <- ggemmeans(m1,"LC",condition=list(abundaceMethod="Abundance",SC_pca=-1.33943234,SamplingEffort.sd=0))
newdata_Tgamma$percent_mean <- ((newdata_Tgamma$predicted)-(newdata_Tgamma$predicted[newdata_Tgamma$x == "Natural"]))/(newdata_Tgamma$predicted[newdata_Tgamma$x == "Natural"])*100
newdata_Tgamma$percent_low <- ((newdata_Tgamma$conf.low)-(newdata_Tgamma$predicted[newdata_Tgamma$x == "Natural"]))/(newdata_Tgamma$predicted[newdata_Tgamma$x == "Natural"])*100
newdata_Tgamma$percent_high <- ((newdata_Tgamma$conf.high)-(newdata_Tgamma$predicted[newdata_Tgamma$x == "Natural"]))/(newdata_Tgamma$predicted[newdata_Tgamma$x == "Natural"])*100

newdata_Pgamma <- ggemmeans(m2,"LC",condition=list(abundaceMethod="Abundance",SC_pca=-1.33943234,SamplingEffort.sd=0))
newdata_Pgamma$percent_mean <- ((newdata_Pgamma$predicted)-(newdata_Pgamma$predicted[newdata_Pgamma$x == "Natural"]))/(newdata_Pgamma$predicted[newdata_Pgamma$x == "Natural"])*100
newdata_Pgamma$percent_low <- ((newdata_Pgamma$conf.low)-(newdata_Pgamma$predicted[newdata_Pgamma$x == "Natural"]))/(newdata_Pgamma$predicted[newdata_Pgamma$x == "Natural"])*100
newdata_Pgamma$percent_high <- ((newdata_Pgamma$conf.high)-(newdata_Pgamma$predicted[newdata_Pgamma$x == "Natural"]))/(newdata_Pgamma$predicted[newdata_Pgamma$x == "Natural"])*100

newdata_Pgamma_C <- ggemmeans(m3,"LC",condition=list(abundaceMethod="Abundance",SC_pca=-1.33943234,SamplingEffort.sd=0))
newdata_Pgamma_C$percent_mean <- ((newdata_Pgamma_C$predicted)-(newdata_Pgamma_C$predicted[newdata_Pgamma_C$x == "Natural"]))/(newdata_Pgamma_C$predicted[newdata_Pgamma_C$x == "Natural"])*100
newdata_Pgamma_C$percent_low <- ((newdata_Pgamma_C$conf.low)-(newdata_Pgamma_C$predicted[newdata_Pgamma_C$x == "Natural"]))/(newdata_Pgamma_C$predicted[newdata_Pgamma_C$x == "Natural"])*100
newdata_Pgamma_C$percent_high <- ((newdata_Pgamma_C$conf.high)-(newdata_Pgamma_C$predicted[newdata_Pgamma_C$x == "Natural"]))/(newdata_Pgamma_C$predicted[newdata_Pgamma_C$x == "Natural"])*100

### form df for graphs
column_names <- c("percent_mean","percent_low","percent_high","x")
graph_df <- rbind(newdataTD[,column_names],
                  newdataPD[,column_names],
                  newdataPD_C[,column_names],
                  as.data.frame(newdata_Tbeta)[,column_names],
                  as.data.frame(newdata_Pbeta)[,column_names],
                  as.data.frame(newdata_Pbeta_C)[,column_names],
                  as.data.frame(newdata_Tgamma)[,column_names],
                  as.data.frame(newdata_Pgamma)[,column_names],
                  as.data.frame(newdata_Pgamma_C)[,column_names])

graph_df$metrics <- c(rep(paste("(b)~alpha~diversity"),9),rep(paste("(c)~beta~diversity"),9),rep(paste("(d)~gamma~diversity"),9))
graph_df$facet <- rep(c(rep("Taxonomic",3),rep("Phylogenetic",3),rep("Phylogenetic (Taxonomic controlled)",3)),3)

graph_df$label <- c("a","b","b",
                    "a","b","b",
                    "a","a","a",
                    "a*","b*","c*",
                    "a","a","b",
                    "a","a","b",
                    "a","b","ab",
                    "a","b","ab",
                    "a","a","a")
graph_df$label_y <- ifelse(graph_df$percent_high >-5, graph_df$percent_high+5,graph_df$percent_high+3.5)
### Result graph for fig.1

graph_df <- as.data.frame(graph_df)

library(ggplot2)
graph_df_p1 <- subset(graph_df,facet != "Phylogenetic (Taxonomic controlled)")
graph_df_p1$facet <- relevel(factor(graph_df_p1$facet,ordered=F),ref="Taxonomic")

p1 <- ggplot(data=graph_df_p1,aes(x=facet,colour=x,y=percent_mean,group=x))+
  geom_hline(yintercept=0,linetype=2)+
  geom_point(position=position_dodge(width = 0.9),size=2)+
  geom_errorbar(aes(ymin=percent_low,ymax=percent_high),position=position_dodge(width = 0.9),width=0.5)+
  geom_text(aes(y=label_y,label=label),position=position_dodge(width=0.9),show.legend  = FALSE)+
  ylab(expression(paste("Differences (%)")))+
  xlab("")+
  scale_y_continuous(limits = c(-40, 25), breaks = seq(-40, 25, by = 10))+
  scale_colour_manual(values=c("green4","darkgoldenrod2","grey25"),name = "Land use")+
  theme(axis.text = element_text(size=10))+
  facet_wrap(~metrics,labeller=label_parsed)+
  theme_classic()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        strip.text = element_text(hjust = 0,size=12),
        strip.background = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))
plot(p1)

ggsave("Figure/p1_obs.tiff",width=8,height=4,compression="lzw",dpi=600)

####
graph_df_p2 <- subset(graph_df,facet == "Phylogenetic (Taxonomic controlled)")
graph_df_p2$label_y <- graph_df_p2$percent_high+2.5
graph_df_p2$metrics <- c(rep(paste("(a)~Phylogenetic~alpha~diversity"),3),
                       rep(paste("(b)~Phylogenetic~beta~diversity"),3),
                       rep(paste("(c)~Phylogenetic~gamma~diversity"),3))

p2 <- ggplot(data=graph_df_p2,aes(x=x,colour=x,y=percent_mean,group=x))+
  geom_hline(yintercept=0,linetype=2)+
  geom_point(position=position_dodge(width = 0.9),size=2)+
  geom_errorbar(aes(ymin=percent_low,ymax=percent_high),position=position_dodge(width = 0.9),width=0.5)+
  geom_text(aes(y=label_y,label=label),position=position_dodge(width=0.9))+
  ylab(expression(paste("Differences after controlling for taxonomic diversity (%)")))+
  xlab("")+
  scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, by = 5))+
  scale_colour_manual(values=c("green4","darkgoldenrod2","grey25"))+
  theme(axis.text = element_text(size=10))+
  facet_wrap(~metrics,labeller=label_parsed)+
  theme_classic()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        strip.text = element_text(hjust = 0,size=12),
        strip.background = element_blank(),
        legend.position="none")
plot(p2)
ggsave("Figure/p2_p_controlled.tiff",width=8,height=5,compression="lzw",dpi=600)

####

Lat_df <- analysis.df.PD[,c("studyID","Lat")]
Lat_df$type <- "Assemblage"

study_df_Lat <- analysis.df.PD %>% group_by(studyID) %>% summarize(Lat = mean(Lat))
study_df_Lat$type <- "Study"

Lat_df_plot <- rbind(Lat_df,study_df_Lat)

plot_density <- ggplot(Lat_df_plot,aes(x=Lat,linetype=type))+
  geom_density()+
  xlab("Latitude")+
  ylab("Density")+
  theme_classic()+
  ylim(0,0.04)+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        strip.text = element_text(hjust = 0,size=12),
        strip.background = element_blank(),
        legend.position="none") # a plot with 1) assemblage and 2) study density. not included in the manuscript


plot(plot_density)

plot_density <- ggplot(analysis.df.PD,aes(x=Lat))+
  geom_density()+
  xlab("Latitude")+
  ylab("Assemblage density")+
  theme_classic()+
  ylim(0,0.04)+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        strip.text = element_text(hjust = 0,size=12),
        strip.background = element_blank(),
        legend.position="none")


plot(plot_density) #we used this in Fig. 1

ggsave("Figure/density.tiff",width=2,height=2,compression="lzw",dpi=600)
