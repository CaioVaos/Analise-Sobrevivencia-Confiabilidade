require(flexsurv)
require(reslife)
library(ggplot2)



tempos<-c(3,5,6,7,8,9,10,10,12,15,15,18,19,20,22,25,28,30,40,45)
cens<-c(1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0)
dados1<-data.frame(tempos,cens)

ajuste_exp<-flexsurvreg(Surv(tempos, cens) ~ 1, dist = "exp", data=dados1)
ajuste_wei<-flexsurvreg(Surv(tempos, cens) ~ 1, dist = "weibull", data=dados1)
ajuste_lognorm<-flexsurvreg(Surv(tempos, cens) ~ 1, dist = "lnorm", data=dados1)
ajuste_gama<-flexsurvreg(Surv(tempos, cens) ~ 1, dist = "gamma", data=dados1)
ajuste_gengama<-flexsurvreg(Surv(tempos, cens) ~ 1, dist = "gengamma", data=dados1)

# 1- Metodo grafico ----
# Aberto a subjetividade

ekm <- survfit(Surv(tempos,cens)~1, data=dados1)
Surv_ekm <- data.frame(ekm$time,ekm$surv)
colnames(Surv_ekm) <- c("Tempo", "Surv_EKM")

Surv_pam <- data.frame(summary(ajuste_exp)[[1]][,c(1,2)],
                       summary(ajuste_wei)[[1]][,2],
                       summary(ajuste_lognorm)[[1]][,2],
                       summary(ajuste_gama)[[1]][,2],
                       summary(ajuste_gengama)[[1]][,2])
colnames(Surv_pam) <- c("Tempo", "Surv_exp", "Surv_wei", "Surv_lognorm", "Surv_gama", "Surv_gengama")

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

# 2- Linearização ----
# Verificar se ao passsar a função de ligação vira um preditor linear

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

# 3- AIC BIC ----
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

# 4 - TRV ----

TRV_exp = 2*(ajuste_gengama$loglik-ajuste_exp$loglik)
P_valor_exp = pchisq(TRV_exp,2,lower.tail=F)

TRV_wei = 2*(ajuste_gengama$loglik-ajuste_wei$loglik)
P_valor_wei = pchisq(TRV_wei,1,lower.tail=F)

TRV_lognorm = 2*(ajuste_gengama$loglik-ajuste_lognorm$loglik)
P_valor_lognorm = pchisq(TRV_lognorm,1,lower.tail=F)

TRV_gama = 2*(ajuste_gengama$loglik-ajuste_gama$loglik)
P_valor_gama = pchisq(TRV_gama,1,lower.tail=F)

# Exercicio ----

dados1 <-  read.csv("Aulas/Data/DadosAIDS.csv")

ajuste_exp<-flexsurvreg(Surv(tempo, status) ~ 1, dist = "exp", data=dados1)
ajuste_wei<-flexsurvreg(Surv(tempo, status) ~ 1, dist = "weibull", data=dados1)
ajuste_lognorm<-flexsurvreg(Surv(tempo, status) ~ 1, dist = "lnorm", data=dados1)
ajuste_gama<-flexsurvreg(Surv(tempo, status) ~ 1, dist = "gamma", data=dados1)
ajuste_gengama<-flexsurvreg(Surv(tempo, status) ~ 1, dist = "gengamma", data=dados1)

# 1- Metodo grafico ----
# Aberto a subjetividade

ekm <- survfit(Surv(tempo,status)~1, data=dados1)
Surv_ekm <- data.frame(ekm$time,ekm$surv)
colnames(Surv_ekm) <- c("Tempo", "Surv_EKM")

Surv_pam <- data.frame(summary(ajuste_exp)[[1]][,c(1,2)],
                       summary(ajuste_wei)[[1]][,2],
                       summary(ajuste_lognorm)[[1]][,2],
                       summary(ajuste_gama)[[1]][,2],
                       summary(ajuste_gengama)[[1]][,2])
colnames(Surv_pam) <- c("Tempo", "Surv_exp", "Surv_wei", "Surv_lognorm", "Surv_gama", "Surv_gengama")

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

# 2- Linearização ----
# Verificar se ao passsar a função de ligação vira um preditor linear

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

# 3- AIC BIC ----
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

# 4 - TRV ----

TRV_exp = 2*(ajuste_gengama$loglik-ajuste_exp$loglik)
P_valor_exp = pchisq(TRV_exp,2,lower.tail=F)

TRV_wei = 2*(ajuste_gengama$loglik-ajuste_wei$loglik)
P_valor_wei = pchisq(TRV_wei,1,lower.tail=F)

TRV_lognorm = 2*(ajuste_gengama$loglik-ajuste_lognorm$loglik)
P_valor_lognorm = pchisq(TRV_lognorm,1,lower.tail=F)

TRV_gama = 2*(ajuste_gengama$loglik-ajuste_gama$loglik)
P_valor_gama = pchisq(TRV_gama,1,lower.tail=F)
