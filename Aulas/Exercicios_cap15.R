require(flexsurv)
require(reslife)
require(dplyr)
require(ggplot2)
library(qqplotr) 
library(coin)


temp<-c(65,156,100,134,16,108,121,4,39,143,56,26,22,1,1,5,65)
cens<-c(rep(1,16),0)
log10_wbc<-c(3.36,2.88,3.63,3.41,3.78,4.02,
             4.00,4.23,3.73,3.85,3.97,
             4.51,4.54,5.00,5.00,4.72,5.00)
dados1<-read.csv("Downloads/AIDS.csv")

ajuste_exp<-flexsurvreg(data=dados1, 
                        Surv(tempo,status)~tratam,
                        dist="exponential")
ajuste_wei<-flexsurvreg(data=dados1,  
                        Surv(tempo,status)~tratam,
                        dist="weibull")
ajuste_lognorm<-flexsurvreg(data=dados1,
                            Surv(tempo,status)~tratam,
                            dist="lognormal")
ajuste_gama<-flexsurvreg(data=dados1,
                         Surv(tempo,status)~tratam,
                         dist="gamma")
ajuste_gengama<-flexsurvreg(data=dados1,
                            Surv(tempo,status)~tratam,
                            dist="gengamma")


tidy(ajuste_exp)
Est = exp(-ajuste_exp$coefficient)
IC = exp(-confint(ajuste_exp)) 
Est
IC


tidy(ajuste_wei)
Est = exp(-ajuste_wei$coefficient)
IC = exp(-confint(ajuste_wei)) 
Est
IC


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
TRV


ResiduoCS = coxsnell_flexsurvreg(ajuste_lognorm)$est
par(mfrow=c(1,2))
ekm<- survfit (Surv(ResiduoCS,dados1$status)~1,type="kaplan-meier")
plot(ekm, fun="cumhaz",xlab="Resíduos de Cox-Snell",ylab="Função de risco acumulado")
abline(0, 1, col="red")

plot(ekm, xlab="Resíduos de Cox-Snell",ylab="Função de sobrevivência")
ResiduoCS<-sort(ResiduoCS)
exp1<-exp(-ResiduoCS)
points(ResiduoCS,exp1, col="red",type="l")



survival = NULL
quantilic = NULL
cdf = NULL
Residuo = NULL
for(i in 1:nrow(dados1)){
  survival[i] <- summary(
    ajuste_lognorm,
    type = "survival",
    newdata = dados1[i, , drop = FALSE],
    t = dados1$tempo[i]
  )[[1]]$est
  cdf[i] = 1-survival[i]
  if(dados1$status[i] == 1){
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

