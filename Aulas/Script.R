require(flexsurv)
require(tidyverse)
require(ggplot2)
require(qqplotr)

#########################################################################
###################### SEM COVARIÁVEIS ##################################
#########################################################################

#### AJUSTE DOS MODELOS

ajuste_exp<-flexsurvreg(Surv(tempos, cens) ~ 1, 
                        dist = "exponential", 
                        data=Dados)
ajuste_wei<-flexsurvreg(Surv(tempos, cens) ~ 1, 
                        dist = "weibull", 
                        data=Dados)
ajuste_lognorm<-flexsurvreg(Surv(tempos, cens) ~ 1, 
                        dist = "lognormal", 
                        data=Dados)
ajuste_gama<-flexsurvreg(Surv(tempos, cens) ~ 1, 
                            dist = "gamma", 
                            data=Dados)
ajuste_gengama<-flexsurvreg(Surv(tempos, cens) ~ 1, 
                            dist = "gengamma", 
                            data=Dados)

##### KAPLAN-MEIER

ekm <- survfit(Surv(tempos,cens)~1, data=Dados)

##### MÉTODO GRÁFICO KAPLAN-MEIER

Surv_ekm <- data.frame(ekm$time,ekm$surv)
colnames(Surv_ekm) <- c("Tempo", "Surv_EKM")

Surv_pam <- data.frame(summary(ajuste_exp)[[1]][,c(1,2)],
                       summary(ajuste_wei)[[1]][,2],
                       summary(ajuste_lognorm)[[1]][,2],
                       summary(ajuste_gama)[[1]][,2],
                       summary(ajuste_gengama)[[1]][,2])
colnames(Surv_pam) <- c("Tempo", "Surv_exp", "Surv_wei", 
                        "Surv_lognorm", "Surv_gama", "Surv_gengama")

Surv_Final <- merge(Surv_pam,Surv_ekm,by="Tempo",all.x=T)

par(mfrow=c(2,3))
plot(Surv_Final$Surv_EKM, Surv_Final$Surv_exp, 
     xlab= "S(t): Kaplan-Meier", ylab = "S(t): Exponencial")
lines(c(0,1), c(0,1), type="l", lty=1)
plot(Surv_Final$Surv_EKM, Surv_Final$Surv_wei, 
     xlab= "S(t): Kaplan-Meier", ylab = "S(t): Weibull" )
lines(c(0,1), c(0,1), type="l", lty=1)
plot(Surv_Final$Surv_EKM, Surv_Final$Surv_lognorm, 
     xlab= "S(t): Kaplan-Meier", ylab = "S(t): LogNormal")
lines(c(0,1), c(0,1), type="l", lty=1)
plot(Surv_Final$Surv_EKM, Surv_Final$Surv_gama, 
     xlab= "S(t): Kaplan-Meier", ylab = "S(t): Gama" )
lines(c(0,1), c(0,1), type="l", lty=1)
plot(Surv_Final$Surv_EKM, Surv_Final$Surv_gengama, 
     xlab= "S(t): Kaplan-Meier", ylab = "S(t): Gama Generalizada" )
lines(c(0,1), c(0,1), type="l", lty=1)


##### MÉTODO GRÁFICO LINEARIZAÇÃO

EixoY_exp = -log(Surv_ekm$Surv_EKM)
EixoX_exp = Surv_ekm$Tempo

EixoY_wei = log(-log(Surv_ekm$Surv_EKM))
EixoX_wei = log(Surv_ekm$Tempo)

EixoY_lognorm = qnorm(Surv_ekm$Surv_EKM)
EixoX_lognorm = log(Surv_ekm$Tempo)

par(mfrow=c(1,3))
plot(EixoX_exp, EixoY_exp, 
     xlab = "t",
     ylab = expression(log(S(t))),
     main = "Exponencial")
plot(EixoX_wei, EixoY_wei,
     xlab = expression(log(t)),
     ylab = expression(log(-log(S(t)))),
     main = "Weibull")
plot(EixoX_lognorm, EixoY_lognorm,
     xlab = expression(log(t)),
     ylab = expression(Phi^{-1}(S(t))),
     main = "Log-Normal")


##### AIC E BIC

AIC<-c(ajuste_exp$AIC,
       ajuste_wei$AIC,
       ajuste_lognorm$AIC,
       ajuste_gama$AIC,
       ajuste_gengama$AIC)
BIC<-c(ajuste_exp$BIC,
       ajuste_wei$BIC,
       ajuste_lognorm$BIC,
       ajuste_gama$BIC,
       ajuste_gengama$BIC)
Modelos<-c("Modelo Exponencial", "Modelo Weibull",
           "Modelo Log-Normal",
           "Modelo Gama", "Modelo Gama Generalizado")
Comparacao<-data.frame(Modelos,AIC,BIC)
Comparacao


##### TRV

TRV_exp = 2*(ajuste_gengama$loglik-ajuste_exp$loglik)
P_valor_exp = pchisq(TRV_exp,2,lower.tail=F)
TRV_wei = 2*(ajuste_gengama$loglik-ajuste_wei$loglik)
P_valor_wei = pchisq(TRV_wei,1,lower.tail=F)
TRV_lognorm = 2*(ajuste_gengama$loglik-ajuste_lognorm$loglik)
P_valor_lognorm = pchisq(TRV_lognorm,1,lower.tail=F)
TRV_gama = 2*(ajuste_gengama$loglik-ajuste_gama$loglik)
P_valor_gama = pchisq(TRV_gama,1,lower.tail=F)
Comparacao<-c("Exponencial x GG", "Weibull x GG",
              "LogNormal x GG",
              "Gama x GG")
P_valor = c(P_valor_exp,P_valor_wei,
            P_valor_lognorm,P_valor_gama)
TRV<-data.frame(Comparacao,P_valor)

###### ESTIMATIVAS DE FUNÇÕES BÁSICAS E MEDIDAS-RESUMO PARA MODELO FINAL

t <- seq(min(Dados$tempos), max(Dados$tempos), length.out = 100)
S <- summary(ajuste, t = t, type = "survival")[[1]]
h <- summary(ajuste, t = t, type = "hazard")[[1]]
H <- summary(ajuste, t = t, type = "cumhaz")[[1]]

df <- rbind(
  data.frame(t = t, valor = S$est, li = S$lcl, ls = S$ucl, funcao = "S(t)"),
  data.frame(t = t, valor = h$est, li = h$lcl, ls = h$ucl, funcao = "h(t)"),
  data.frame(t = t, valor = H$est, li = H$lcl, ls = H$ucl, funcao = "H(t)")
)

df$funcao <- factor(
  df$funcao,
  levels = c("S(t)", "h(t)", "H(t)"),
  labels = c("S(t)", "lambda(t)", "Lambda(t)")
)

ggplot(df, aes(t, valor, color = funcao, fill = funcao)) +
  geom_ribbon(aes(ymin = li, ymax = ls), alpha = 0.2, color = NA) +
  geom_line(size = 1.1) +
  facet_wrap(~funcao, scales = "free_y", labeller = label_parsed) +
  labs(x = "Tempo", y = NULL,
       title = "Funções do Modelo Log-normal") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

Mean<-summary(ajuste,type="mean",cl=0.95,se=T)
Median<-summary(ajuste,type="median",cl=0.95,se=T)

#########################################################################
###################### UMA COVARIÁVEL ###################################
#########################################################################

#### AJUSTE DOS MODELOS

ajuste_exp<-flexsurvreg(Surv(tempos, cens) ~ covariavel, 
                        dist = "exponential", 
                        data=Dados)
ajuste_wei<-flexsurvreg(Surv(tempos, cens) ~ covariavel, 
                        dist = "weibull", 
                        data=Dados)
ajuste_lognorm<-flexsurvreg(Surv(tempos, cens) ~ covariavel, 
                            dist = "lognormal", 
                            data=Dados)
ajuste_gama<-flexsurvreg(Surv(tempos, cens) ~ covariavel, 
                         dist = "gamma", 
                         data=Dados)
ajuste_gengama<-flexsurvreg(Surv(tempos, cens) ~ covariavel, 
                            dist = "gengamma", 
                            data=Dados)


tidy(ajuste)

#### AIC E BIC

AIC<-c(ajuste_exp$AIC,
       ajuste_wei$AIC,
       ajuste_lognorm$AIC,
       ajuste_gama$AIC,
       ajuste_gengama$AIC)
BIC<-c(ajuste_exp$BIC,
       ajuste_wei$BIC,
       ajuste_lognorm$BIC,
       ajuste_gama$BIC,
       ajuste_gengama$BIC)
Modelos<-c("Modelo Exponencial", "Modelo Weibull",
           "Modelo Log-Normal",
           "Modelo Gama", "Modelo Gama Generalizado")
Comparacao<-data.frame(Modelos,AIC,BIC)
Comparacao

#### TRV

TRV_exp = 2*(ajuste_gengama$loglik-ajuste_exp$loglik)
P_valor_exp = pchisq(TRV_exp,2,lower.tail=F)
TRV_wei = 2*(ajuste_gengama$loglik-ajuste_wei$loglik)
P_valor_wei = pchisq(TRV_wei,1,lower.tail=F)
TRV_lognorm = 2*(ajuste_gengama$loglik-ajuste_lognorm$loglik)
P_valor_lognorm = pchisq(TRV_lognorm,1,lower.tail=F)
TRV_gama = 2*(ajuste_gengama$loglik-ajuste_gama$loglik)
P_valor_gama = pchisq(TRV_gama,1,lower.tail=F)
Comparacao<-c("Exponencial x GG", "Weibull x GG",
              "LogNormal x GG",
              "Gama x GG")
P_valor = c(P_valor_exp,P_valor_wei,
            P_valor_lognorm,P_valor_gama)
TRV<-data.frame(Comparacao,P_valor)

#### RESIDUO COX-SNELL

ResiduoCS = coxsnell_flexsurvreg(ajuste)$est
par(mfrow=c(1,2))
ekm<- survfit (Surv(ResiduoCS,Dados$cens)~1,type="kaplan-meier")
plot(ekm, fun="cumhaz",xlab="Resíduos de Cox-Snell",ylab="Função de risco acumulado")
abline(0, 1, col="red")

plot(ekm, xlab="Resíduos de Cox-Snell",ylab="Função de sobrevivência")
ResiduoCS<-sort(ResiduoCS)
exp1<-exp(-ResiduoCS)
points(ResiduoCS,exp1, col="red",type="l")

#### RESIDUO QUANTILICO ALEATORIZADO

survival = NULL
quantilic = NULL
cdf = NULL
Residuo = NULL
for(i in 1:nrow(Dados)){
  survival[i] <- summary(
    ajuste,
    type = "survival",
    newdata = Dados[i, , drop = FALSE],
    t = Dados$tempos[i]
  )[[1]]$est
  cdf[i] = 1-survival[i]
  if(Dados$cens[i] == 1){
    quantilic[i] <- cdf[i]
  } else {
    quantilic[i] <- runif(1, min = cdf[i], max = 1)
  }
  Residuo[i] <- qnorm(quantilic[i])
}

ggplot(data.frame(Residuo), aes(sample=Residuo)) +
  stat_qq_band(fill="lightblue") +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Quantis Teóricos")+
  ylab("Quantis Amostrais")+
  theme_bw() +
  theme(text = element_text(size = 14))


###### INTERPRETAÇÃO DOS COEFICIENTES, 
###### ESTIMATIVAS DE FUNÇÕES BÁSICAS E 
###### MEDIDAS-RESUMO PARA MODELO FINAL

## Se o modelo for gamma ou exponencial:
Est = exp(-ajuste_gengama$coefficient)
IC = exp(-confint(ajuste_gengama)) 

## Para os demais modelos
Est = exp(ajuste_gengama$coefficient)
IC = exp(confint(ajuste_gengama)) 


par(mfrow=c(1,2))
plot(ajuste,newdata = list(covariavel=covariavel_individuo_1),type="survival",
     col="red",col.obs="white", ci=F, xlab="Semanas",ylab="Função de sobrevivência estimada",t=c(0:30))
lines(ajuste,newdata = list(covariavel=covariavel_individuo_2),type="survival",
     col="green",col.obs="white", ci=F, t=c(0:30))
legend(2,0.3,col=c("red","green"),legend=c("Indivíduo 1", "Indivíduo 2"),lty=1)

plot(ajuste,newdata = list(covariavel=covariavel_individuo_1),type="hazard",
     col="red",col.obs="white", ci=F, xlab="Semanas",ylab="Função de sobrevivência estimada",t=c(0:30))
lines(ajuste,newdata = list(covariavel=covariavel_individuo_2),type="hazard",
     col="green",col.obs="white", ci=F, t=c(0:30))
legend(2,0.3,col=c("red","green"),legend=c("Indivíduo 1", "Indivíduo 2"),lty=1)


media<-summary(ajuste,newdata = list(covariavel=covariavel_individuo),type="mean")
median<-summary(ajuste,newdata = list(covariavel=covariavel_individuo),type="median")
