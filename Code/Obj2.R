m_pp <- lm(log(qPD.S)~log(Exclude.AM.Est.Shannon),data=analysis.df.PD)
summary(m_pp)

predict_df <- ggpredict(m_pp,terms="Exclude.AM.Est.Shannon[1:75]")

r2 = round(summary(m_pp)$r.squared,2)
p_alpha_corr <- ggplot(data=analysis.df.PD,aes(y=qPD.S,x=Exclude.AM.Est.Shannon))+
  geom_point()+ 
  geom_line(data=predict_df,aes(y=predicted,x=x),colour="blue")+
  geom_ribbon(data=predict_df,aes(y=predicted,ymin=conf.low,ymax=conf.high,x=x),alpha=0.25,fill="blue")+
  scale_y_continuous(limits=c(0.9,8),breaks = seq(1,8,by=2))+
  ylab(expression(paste("Phylogenetic ",alpha,"-diversity")))+
  xlab(expression(paste("Taxonomic ",alpha,"-diversity")))+
  annotate('text', label=bquote(~(a)~R^2*" = "*.(r2)), x=-Inf, y=Inf, hjust=0, vjust=1,size=4)+
  #annotate('text', label=expression(paste(" (a) ",R^2,"=",)), x=-Inf, y=Inf, hjust=0, vjust=1,size=4)+
  theme_classic()

plot(p_alpha_corr)

m_total_b <- lm(log(PBtotal)~log(TBtotal),data=cor_df_Btotal)
summary(m_total_b)
r2 = round(summary(m_total_b)$r.squared,2)

predict_total_b <- ggpredict(m_total_b,"TBtotal[0.01:1,by=0.01]")
predict_total_b <- subset(predict_total_b,x >= min(cor_df_Btotal$TBtotal) & x <= max(cor_df_Btotal$TBtotal))

p_corr_Btotal <- ggplot(data=cor_df_Btotal,aes(y=PBtotal,x=TBtotal))+
  geom_point()+
  geom_line(data=predict_total_b,aes(y=predicted,x=x),colour="blue")+
  geom_ribbon(data=predict_total_b,aes(y=predicted,ymin=conf.low,ymax=conf.high,x=x),alpha=0.25,fill="blue")+
  annotate('text', label=bquote(~(b)~R^2*" = "*.(r2)), x=-Inf, y=Inf, hjust=0, vjust=1,size=4)+
  ylab(expression(Phylogenetic ~ beta[total]))+
  xlab(expression(Taxonomic ~ beta[total]))+
  ylim(0,1)+
  xlim(0,1)+
  theme_classic()
plot(p_corr_Btotal)

cor_df_Bturn <- data.frame(TBturn=P_beta$T_beta[P_beta$metrics == "Bturn"],PBturn=P_beta$mean.dist[P_beta$metrics == "Bturn"])
m_turn_b <- lm(log(PBturn)~log(TBturn),data=cor_df_Bturn)
summary(m_turn_b)
r2 = round(summary(m_turn_b)$r.squared,2)

predict_turn_b <- ggpredict(m_turn_b,"TBturn[0.01:1,by=0.01]")
predict_turn_b <- subset(predict_turn_b,x >= min(cor_df_Bturn$TBturn) & x <= max(cor_df_Bturn$TBturn))

p_corr_Bturn <- ggplot(data=cor_df_Bturn,aes(y=PBturn,x=TBturn))+geom_point()+
  geom_line(data=predict_turn_b,aes(y=predicted,x=x),colour="blue")+
  geom_ribbon(data=predict_turn_b,aes(y=predicted,ymin=conf.low,ymax=conf.high,x=x),alpha=0.25,fill="blue")+
  annotate('text', label=bquote(~(c)~R^2*" = "*.(r2)), x=-Inf, y=Inf, hjust=0, vjust=1,size=4)+
  ylab(expression(Phylogenetic ~ beta[turn]))+
  xlab(expression(Taxonomic ~ beta[turn]))+
  ylim(0,1)+
  xlim(0,1)+
  theme_classic()
plot(p_corr_Bturn)

cor_df_Bab <- data.frame(TBab=P_beta$T_beta[P_beta$metrics == "Bab"],PBab=P_beta$mean.dist[P_beta$metrics == "Bab"])
m_rich_b <- lm(log(PBab)~log(TBab),data=cor_df_Bab)
summary(m_rich_b)
r2 = round(summary(m_rich_b)$r.squared,2)

predict_rich_b <- ggpredict(m_rich_b,"TBab[0.01:1,by=0.01]")
predict_rich_b <- subset(predict_rich_b,x >= min(cor_df_Bab$TBab) & x <= max(cor_df_Bab$TBab))

p_corr_Bab <- ggplot(data=cor_df_Bab,aes(y=PBab,x=TBab))+geom_point()+
  geom_line(data=predict_rich_b,aes(y=predicted,x=x),colour="blue")+
  geom_ribbon(data=predict_rich_b,aes(y=predicted,ymin=conf.low,ymax=conf.high,x=x),alpha=0.25,fill="blue")+
  annotate('text', label=bquote(~(d)~R^2*" = "*.(r2)), x=-Inf, y=Inf, hjust=0, vjust=1,size=4)+
  ylab(expression(Phylogenetic ~ beta[ab]))+
  xlab(expression(Taxonomic ~ beta[ab]))+
  ylim(0,1)+
  xlim(0,1)+
  theme_classic()
plot(p_corr_Bab)

p3 <- ggarrange(p_alpha_corr,p_corr_Btotal,p_corr_Bturn,p_corr_Bab,
                ncol=2,nrow=2,
                common.legend=T,legend="bottom")
plot(p3)

ggsave("Figure/p3.tiff",dpi=600,compression="lzw",width=16.8,height=16.8,units="cm")

loglogtable <- rbind(round(summary(m_pp)$coefficients,3),round(summary(m_total_b)$coefficients,3),round(summary(m_turn_b)$coefficients,3),round(summary(m_rich_b)$coefficients,3))
write.csv(loglogtable,"Table/loglogtable.csv")
