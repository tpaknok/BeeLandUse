m_pp <- lm(log10(qPD.S)~log10(Exclude.AM.Est.Shannon),data=analysis.df.PD)
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

m_total_b <- lm(log10(PBtotal)~log10(TBtotal),data=cor_df_Btotal)
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

cor_df_Brepl <- data.frame(TBrepl=P_beta$T_beta[P_beta$metrics == "Brepl"],PBrepl=P_beta$mean.dist[P_beta$metrics == "Brepl"])
m_repl_b <- lm(log10(PBrepl)~log10(TBrepl),data=cor_df_Brepl)
summary(m_repl_b)
r2 = round(summary(m_repl_b)$r.squared,2)

predict_repl_b <- ggpredict(m_repl_b,"TBrepl[0.01:1,by=0.01]")
predict_repl_b <- subset(predict_repl_b,x >= min(cor_df_Brepl$TBrepl) & x <= max(cor_df_Brepl$TBrepl))

p_corr_Brepl <- ggplot(data=cor_df_Brepl,aes(y=PBrepl,x=TBrepl))+geom_point()+
  geom_line(data=predict_repl_b,aes(y=predicted,x=x),colour="blue")+
  geom_ribbon(data=predict_repl_b,aes(y=predicted,ymin=conf.low,ymax=conf.high,x=x),alpha=0.25,fill="blue")+
  annotate('text', label=bquote(~(c)~R^2*" = "*.(r2)), x=-Inf, y=Inf, hjust=0, vjust=1,size=4)+
  ylab(expression(Phylogenetic ~ beta[repl]))+
  xlab(expression(Taxonomic ~ beta[repl]))+
  ylim(0,1)+
  xlim(0,1)+
  theme_classic()
plot(p_corr_Brepl)

cor_df_Brich <- data.frame(TBrich=P_beta$T_beta[P_beta$metrics == "Brich"],PBrich=P_beta$mean.dist[P_beta$metrics == "Brich"])
m_rich_b <- lm(log10(PBrich)~log10(TBrich),data=cor_df_Brich)
summary(m_rich_b)
r2 = round(summary(m_rich_b)$r.squared,2)

predict_rich_b <- ggpredict(m_rich_b,"TBrich[0.01:1,by=0.01]")
predict_rich_b <- subset(predict_rich_b,x >= min(cor_df_Brich$TBrich) & x <= max(cor_df_Brich$TBrich))

p_corr_Brich <- ggplot(data=cor_df_Brich,aes(y=PBrich,x=TBrich))+geom_point()+
  geom_line(data=predict_rich_b,aes(y=predicted,x=x),colour="blue")+
  geom_ribbon(data=predict_rich_b,aes(y=predicted,ymin=conf.low,ymax=conf.high,x=x),alpha=0.25,fill="blue")+
  annotate('text', label=bquote(~(d)~R^2*" = "*.(r2)), x=-Inf, y=Inf, hjust=0, vjust=1,size=4)+
  ylab(expression(Phylogenetic ~ beta[ab]))+
  xlab(expression(Taxonomic ~ beta[ab]))+
  ylim(0,1)+
  xlim(0,1)+
  theme_classic()
plot(p_corr_Brich)

p3 <- ggarrange(p_alpha_corr,p_corr_Btotal,p_corr_Brepl,p_corr_Brich,
                ncol=2,nrow=2,
                common.legend=T,legend="bottom")
plot(p3)

ggsave("Figure/p3.tiff",dpi=600,compression="lzw",width=16.8,height=16.8,units="cm")

loglogtable <- rbind(round(summary(m_pp)$coefficients,3),round(summary(m_total_b)$coefficients,3),round(summary(m_repl_b)$coefficients,3),round(summary(m_rich_b)$coefficients,3))
write.csv(loglogtable,"Table/loglogtable.csv")
