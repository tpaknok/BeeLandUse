library(lme4)
library(lmerTest)
set.seed(123)
xi <- c(rep(0,10),rep(1,10))
ci <<- c(1:154)

x <- expand.grid(xi,ci)

b0=10
b1=0
missing_n <- 111

result <- list()
for (i in 1:1000) {
  message(i)
  z0 <- rep(rnorm(length(ci)),each=20)
  z1 <- rep(rnorm(length(ci)),each=20)
  x$z0 <- z0
  x$zi <- z1
  
  e <- rnorm(n=length(ci)*20)
  x$y <- (b0+z0)+(b1+z1)*x[,1]+e
  
  m <- lmer(y~Var1+(Var1||Var2),data=x)
  summary(m)
  
  missing_df <- data.frame(Var2=sample(1:154,missing_n),Var2_miss = sample(1,missing_n,replace=T))
  missing_df$label <- paste0(missing_df$Var2,"_",missing_df$Var2_miss)
  
  x$label <- paste0(x$Var2,"_",x$Var1)
  x_reduced <- x[!x$label %in% missing_df$label,]
  
  m_reduced <- lmer(y~Var1+(Var1||Var2),data=x_reduced)
  summary(m_reduced)
  
  x_reduced2 <- x[!x$Var2 %in% missing_df$Var2,]
  m_reduced2 <- lmer(y~Var1+(Var1||Var2),data=x_reduced2)
  summary(m_reduced2)
  
  result[[i]]<- data.frame(p_full=coefficients(summary(m))[2,5],
                           slope_full = coefficients(summary(m))[2,1],
                           slope_full_se = coefficients(summary(m))[2,2],
                           intercept_full = coefficients(summary(m))[1,1],
                           p_reduced=coefficients(summary(m_reduced))[2,5],
                           slope_reduced = coefficients(summary(m_reduced))[2,1],
                           slope_full_se_reduced = coefficients(summary(m_reduced))[2,2],
                           intercept_reduced = coefficients(summary(m_reduced))[1,1],
                           p_reduced2=coefficients(summary(m_reduced2))[2,5],
                           slope_reduced2 = coefficients(summary(m_reduced2))[2,1],
                           slope_full_se_reduced_2 = coefficients(summary(m_reduced2))[2,2],
                           intercept_reduced2 = coefficients(summary(m_reduced2))[1,1])
}

result <- do.call(rbind,result)
type_I <- apply(result[,c("p_full","p_reduced","p_reduced2")],2, function(x) sum(x<0.05))/1000
slope <- apply(result[,c("slope_full","slope_reduced","slope_reduced2")],2, range)
apply(result[,c("slope_full","slope_reduced","slope_reduced2")],2, mean)
intercept <- apply(result[,c("intercept_full","intercept_reduced","intercept_reduced2")],2, range)
apply(result[,c("intercept_full","intercept_reduced","intercept_reduced2")],2, mean)
