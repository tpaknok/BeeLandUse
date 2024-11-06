library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)

### Gamma
Beta_Gamma_df <- read.csv("Data/Cleaned_beta.csv")
Beta_Gamma_df$LC <- relevel(factor(Beta_Gamma_df$LC,ordered=F),ref="Natural")
RobustPCA <- rrcov::PcaCov(scale(data.frame(Beta_Gamma_df$SC.sd,log(Beta_Gamma_df$mean_SC)))) #log-transformed to handle non-linearities between SC and diversity
Beta_Gamma_df$SC_PCA <- RobustPCA$scores[,1]
Beta_Gamma_df$SC_PCA2 <- RobustPCA$scores[,2]

m1full <- lmer(T.gamma~LC*log(extent)+abundanceMethod+mean_PC1+mean_PC2+SamplingEffort.sd+log(site)+SC_PCA+(1|study),
             data=Beta_Gamma_df,
             weights=T.gamma.SC,
             REML=T)
summary(m1full)
anova(m1full)
m1 <- lmer(T.gamma~LC+log(extent)+abundanceMethod+mean_PC1+mean_PC2+SamplingEffort.sd+log(site)+SC_PCA+(1|study),
           data=Beta_Gamma_df,
           weights=T.gamma.SC,
           REML=T)
summary(m1)
plot(simulateResiduals(m1))

table_Tgamma <- anova(m1)
pair_comp_TD_gamma <- as.data.frame(pairs(emmeans(m1,"LC",lmer.df = "satterthwaite"),adjust="fdr"))

rTgamma <-performance::r2(m1)
table_Tgamma <- cbind(table_Tgamma,R2m=rTgamma$R2_marginal,R2c=rTgamma$R2_conditional)

new_data <- data.frame(LC=c("Natural","Agriculture","Urban"),
                       abundanceMethod = "Abundance",
                       extent=mean(Beta_Gamma_df$extent),
                       mean_PC1=mean(Beta_Gamma_df$mean_PC1),
                       mean_PC2=mean(Beta_Gamma_df$mean_PC2),
                       SamplingEffort.sd = 0,
                       SC_PCA = mean(Beta_Gamma_df$SC_PCA),
                       site = mean(Beta_Gamma_df$site),
                       study=NA)

m2full <- lmer(P.gamma~LC*log(extent)+abundanceMethod+mean_PC1+mean_PC2++SamplingEffort.sd+log(site)+SC_PCA+(1|study),
             data=Beta_Gamma_df,
             weights=T.gamma.SC,
             REML=T)
summary(m2full)
anova(m2full)
m2 <- lmer(P.gamma~LC+log(extent)+abundanceMethod+mean_PC1+mean_PC2++SamplingEffort.sd+log(site)+SC_PCA+(1|study),
               data=Beta_Gamma_df,
               weights=T.gamma.SC,
               REML=T)
summary(m2)
anova(m2)
plot(simulateResiduals(m2))
pair_comp_PD_gamma <- as.data.frame(pairs(emmeans(m2,"LC",lmer.df = "satterthwaite"),adjust="fdr"))

table_Pgamma <- anova(m2)

rPgamma <-performance::r2(m2)
table_Pgamma <- cbind(table_Pgamma,R2m=rPgamma$R2_marginal,R2c=rPgamma$R2_conditional)

m3 <- lmer(log(P.gamma)~LC+log(T.gamma)+LC+log(extent)+abundanceMethod+mean_PC1+mean_PC2+SamplingEffort.sd+log(site)+SC_PCA+(1|study),
               data=Beta_Gamma_df,
               weights=T.gamma.SC,
               REML=T)
summary(m3)
anova(m3)
plot(simulateResiduals(m3))
table_PgammaControlled <- anova(m3)
rPgammaControlled <-performance::r2(m3)
table_PgammaControlled <- cbind(table_PgammaControlled,R2m=rPgammaControlled$R2_marginal,R2c=rPgammaControlled$R2_conditional)
gamma_table <- rbind(table_Tgamma,
                    table_Pgamma,
                    table_PgammaControlled)
write.csv(gamma_table ,"Table/Gamma_table.csv")
write.csv(rbind(pair_comp_TD_gamma,pair_comp_PD_gamma),"Table/Gamma_table_pair_comp.csv")

################### 

Beta_Gamma_df$T_beta_AGP <- (Beta_Gamma_df$T.gamma)/(Beta_Gamma_df$mean_T_alpha)
Beta_Gamma_df$Beta_Gamma_df_AGP <- (Beta_Gamma_df$P.gamma)/(Beta_Gamma_df$mean_P_alpha)

m4full <- lmer(log(T_beta_AGP)~LC*log(extent)+abundanceMethod+mean_PC1+mean_PC2+SamplingEffort.sd+log(site)+SC_PCA+(1|study),
           data=Beta_Gamma_df,
           weights=T.gamma.SC,
           REML=T)
summary(m4full)
anova(m4full)
m4 <- lmer(log(T_beta_AGP)~LC+log(extent)+abundanceMethod+mean_PC1+mean_PC2+SamplingEffort.sd+log(site)+SC_PCA+(1|study),
           data=Beta_Gamma_df,
           weights=T.gamma.SC,
           REML=T)
summary(m4)
anova(m4)
plot(simulateResiduals(m4))
pair_comp_TD_beta <- pairs(emmeans(m4,"LC",lmer.df = "satterthwaite"),adjust="fdr")
table_Tbeta <- anova(m4)
rTbeta <-performance::r2(m4)
table_Tbeta <- cbind(table_Tbeta,R2m=rTbeta$R2_marginal,R2c=rTbeta$R2_conditional)

m5full <- lmer(log(Beta_Gamma_df_AGP)~LC*log(extent)+abundanceMethod+mean_PC1+mean_PC2+SamplingEffort.sd+log(site)+SC_PCA+(1|study),
               data=Beta_Gamma_df,
               weights=T.gamma.SC,
               REML=T)
summary(m5full)
plot(simulateResiduals(m5full))
anova(m5full)

m5 <- lmer(log(Beta_Gamma_df_AGP)~LC+log(extent)+abundanceMethod+mean_PC1+mean_PC2+SamplingEffort.sd+log(site)+SC_PCA+(1|study),
           data=Beta_Gamma_df,
           weights=T.gamma.SC,
           REML=T)
plot(simulateResiduals(m5))
pair_comp_PD_beta <- pairs(emmeans(m5,"LC",lmer.df = "satterthwaite"),adjust="fdr") 

table_Pbeta <- anova(m5)
rPbeta <-performance::r2(m5)
table_Pbeta <- cbind(table_Pbeta,R2m=rPbeta$R2_marginal,R2c=rPbeta$R2_conditional)

m6 <- lmer(log(Beta_Gamma_df_AGP)~LC+log(T_beta_AGP)+log(extent)+abundanceMethod+mean_PC1+mean_PC2+SamplingEffort.sd+log(site)+SC_PCA+(1|study),
           data=Beta_Gamma_df,
           weights=T.gamma.SC,
           REML=T)
summary(m6)
plot(simulateResiduals(m6))

table_PbetaControlled <- anova(m6)
pair_comp_PD_beta_C <- pairs(emmeans(m6,"LC",lmer.df = "satterthwaite"),adjust="fdr") #no change in the conclusion

rPbetaControlled <-performance::r2(m6)
table_PbetaControlled <- cbind(table_PbetaControlled,R2m=rPbetaControlled$R2_marginal,R2c=rPbetaControlled$R2_conditional)

beta_table <- rbind(table_Tbeta,
                    table_Pbeta,
                    table_PbetaControlled)
write.csv(beta_table,"Table/Beta_table.csv")
write.csv(rbind(as.data.frame(pair_comp_TD_beta),as.data.frame(pair_comp_PD_beta),as.data.frame(pair_comp_PD_beta_C)),"Table/Beta_table_pair_comp.csv")

###