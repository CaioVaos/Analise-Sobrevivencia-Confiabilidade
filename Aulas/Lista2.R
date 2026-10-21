# Lista 2 - Analise de Sobrevivencia e Confiabilidade

# Setup ----

# install.packages("tidyverse")
# install.packages("survival")
# install.packages("flexsurv")
# install.packages("reslife")
# install.packages("qqplotr")
# install.packages("coin")

library(tidyverse)
library(survival)
library(flexsurv)
library(reslife)
library(ggplot2)
library(qqplotr) 
library(coin)

# 1 ----
#
## EXPONENCIAL:
# Risco constante: sem memória. Não depende do tempo já decorrido.
#
## WEIBULL:
# Risco monótono:
# gamma < 1: decrescente
# gamma = 1: constante (Exponencial)
# gamma > 1: crescente
#
## LOG-NORMAL:
# Risco não monótono (unimodal):
# cresce até um máximo e depois decresce, tendendo a zero.
#
## GAMA:
# Risco monótono, tendendo a um limite constante:
# alpha < 1: decrescente
# alpha = 1: constante (Exponencial)
# alpha > 1: crescente

# 2 ----
#
## GAMA GENERALIZADO:
# Família guarda-chuva: inclui Exponencial, Weibull, Gama e Log-Normal.
#
# Serve como modelo de referência para comparar modelos mais simples
# pelo Teste da Razão de Verossimilhanças (TRV).
# Se o modelo simples não é rejeitado, prefere-se o mais parcimonioso.
#
# Função de risco flexível: pode ser crescente, decrescente,
# unimodal ou apresentar formato de banheira.
#
# Possui mais parâmetros (3): maior flexibilidade, mas maior
# risco de sobreajuste e instabilidade de estimação.

#  3 ----

alpha <- 100
gam   <- 2

S_weibull <- function(t) exp(-(t / alpha)^gam)

## a) ----

S_weibull(30)
S_weibull(45)
# 91,4% dos ratos estão livres de tumor aos 30 dias 
# 81,7% dos ratos estão livres de tumor aos 45 dias.

## b) ----

alpha * gamma(1 + 1 / gam)
# 88.62 dias

## c) ----

alpha * log(2)^(1 / gam)
# 83.26 dias
qweibull(0.5, shape = gam, scale = alpha)

## d) ----
h_weibull <- function(t) (gam / alpha) * (t / alpha)^(gam - 1)
h_weibull(c(30, 45, 60))          # 0.006 0.009 0.012

# Interpretação: 
# o risco é crescente (gamma = 2 > 1)
# Entre os ratos que ainda não têm tumor no instante t, a taxa instantânea de aparecimento
# do tumor é de 0,6% por dia aos 30 dias, 0,9% aos 45 dias e 1,2% aos 60 dias
# quanto mais tempo de exposição à substância, maior a propensão ao tumor.

# 4 ----
#
# No papel

# 5 ----
#
# Nos modelos Weibull, LogNormal e Gama, as equações
# de máxima verossimilhança são não lineares e não permitem
# isolar os parâmetros por fórmulas explícitas.
#
# Utilizamos métodos de otimização numérica iterativos,
# como Newton-Raphson, BFGS ou Nelder-Mead.
# O algoritmo parte de valores iniciais e maximiza a
# log-verossimilhança até atingir a convergência.
#
# No R, esse processo é realizado automaticamente por
# funções de ajuste dos modelos de sobrevivência.

# 6 ----
#
## 1) KAPLAN-MEIER x MODELO PARAMÉTRICO:
# Compara S_KM(t) com S_modelo(t).
# Quanto mais próximas as curvas, melhor o ajuste gráfico.
#
## 2) LINEARIZAÇÃO GRÁFICA:
# Transforma S(t) para obter uma relação linear sob o modelo:
# Exponencial:    -log(S(t))       x t
# Weibull:        log(-log(S(t)))  x log(t)
# Log-Normal:     Phi^(-1)(S(t))   x log(t)
# Log-Logística:  log((1-S(t))/S(t)) x log(t)
# Quanto mais próximos os pontos de uma reta, melhor o ajuste.
#
## 3) AIC E BIC:
# AIC = -2*l(theta_hat) + 2*k
# BIC = -2*l(theta_hat) + k*log(n)
# Comparam a qualidade do ajuste penalizando a complexidade.
# Menor valor = modelo preferível.
# Podem comparar modelos aninhados ou não aninhados.
#
## 4) TESTE DA RAZÃO DE VEROSSIMILHANÇAS (TRV):
# Compara apenas modelos aninhados, como modelos simples
# com a Gama Generalizada.
# H0: o modelo mais simples é suficiente.
# p-valor > alpha: não rejeita H0; prefere o modelo simples.
# p-valor <= alpha: rejeita H0; prefere o modelo completo.

# 7 ----

tempos<-c(
    151, 164, 336, 365, 403, 454, 455, 473, 538, 577, 592, 628,
    632, 647, 675, 727, 785, 801, 811, 816, 867, 893, 930, 937,
    976, 1008, 1040, 1051, 1060, 1183, 1329, 1334, 1379, 1380,
    1633, 1769, 1827, 1831, 1849, 2016, 2282, 2415, 2430, 2686,
    2729, 2729, 2729, 2729, 2729, 2729, 2729, 2729, 2729,
    2729, 2729, 2729, 2729, 2729, 2729, 2729
)

cens<-c(rep(1,45),rep(0,15))
dados<-data.frame(tempos,cens)


## a) ----

ajuste_wei <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "weibull", data = dados)

### estimativas

ajuste_wei$res
# shape = 1.28131
# scale = 1993.21512

### Por máxima verossimilhança 

### erros-padrão
ajuste_wei$res
# shape = 0.1614857
# scale = 231.9023378


### Intervalo de Wald na escala log
ajuste_wei$res
# est   ;     L95%    ;    U95%
# shape = 1.000865    ; 1.640335
# scale = 1586.793431 ; 2503.732655

### Tempo mediano de falha
summary(ajuste_wei, type = "median", cl = 0.95, se = TRUE)
# [1181.497 ; 1896.046]

## b) ----

ajuste_exp     <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "exp",      data = dados)
ajuste_lognorm <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "lnorm",    data = dados)
ajuste_llogis  <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "llogis",   data = dados)
ajuste_gama    <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "gamma",    data = dados)
ajuste_gengama <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "gengamma", data = dados)

ajustes7 <- list(
  "Exponencial"       = ajuste_exp,
  "Weibull"           = ajuste_wei,
  "Log-Normal"        = ajuste_lognorm,
  "Log-Logística"     = ajuste_llogis,
  "Gama"              = ajuste_gama,
  "Gama Generalizada" = ajuste_gengama
)

### AIC/BIC 
data.frame(
  Modelo = names(ajustes7),
  logLik = sapply(ajustes7, function(a) a$loglik),
  AIC    = sapply(ajustes7, function(a) a$AIC),
  BIC    = sapply(ajustes7, function(a) a$BIC)
)

### TRV  ----
TRV_exp    <- 2 * (ajuste_gengama$loglik - ajuste_exp$loglik)
TRV_wei    <- 2 * (ajuste_gengama$loglik - ajuste_wei$loglik)
TRV_lnorm  <- 2 * (ajuste_gengama$loglik - ajuste_lognorm$loglik)
TRV_gama   <- 2 * (ajuste_gengama$loglik - ajuste_gama$loglik)

data.frame(
  Comparacao = c("Exponencial x GG", "Weibull x GG", "Log-Normal x GG", "Gama x GG"),
  TRV        = c(TRV_exp, TRV_wei, TRV_lnorm, TRV_gama),
  gl         = c(2, 1, 1, 1),
  p_valor    = c(pchisq(TRV_exp,   2, lower.tail = FALSE),
                 pchisq(TRV_wei,   1, lower.tail = FALSE),
                 pchisq(TRV_lnorm, 1, lower.tail = FALSE),
                 pchisq(TRV_gama,  1, lower.tail = FALSE))
)

# 8 ----

dados1 <- read.csv("Aulas/Data/DadosAIDS.csv")

# Ajuste dos modelos ----
ajuste_exp <- flexsurvreg(Surv(tempo, status) ~ 1, dist = "exp", data = dados1)
ajuste_wei <- flexsurvreg(Surv(tempo, status) ~ 1, dist = "weibull", data = dados1)
ajuste_lognorm <- flexsurvreg(Surv(tempo, status) ~ 1, dist = "lnorm", data = dados1)
ajuste_gama <- flexsurvreg(Surv(tempo, status) ~ 1, dist = "gamma", data = dados1)
ajuste_gengama <- flexsurvreg(Surv(tempo, status) ~ 1, dist = "gengamma", data = dados1)

# 1 - Graficos: sobrevivencia, risco e risco acumulado ----

ekm <- survfit(Surv(tempo, status) ~ 1, data = dados1)

# Sobrevivencia
plot(ekm, xlab = "Tempo", ylab = "S(t)", main = "Sobrevivencia",
     conf.int = FALSE, col = "black")
lines(ajuste_exp, type = "survival", col = 2, ci = FALSE)
lines(ajuste_wei, type = "survival", col = 3, ci = FALSE)
lines(ajuste_lognorm, type = "survival", col = 4, ci = FALSE)
lines(ajuste_gama, type = "survival", col = 5, ci = FALSE)
lines(ajuste_gengama, type = "survival", col = 6, ci = FALSE)

legend("topright",
       legend = c("Kaplan-Meier", "Exponencial", "Weibull",
                  "LogNormal", "Gama", "Gama Generalizada"),
       col = 1:6, lty = 1, cex = 0.7)

# Risco acumulado: Nelson-Aalen
ena <- survfit(Surv(tempo, status) ~ 1, data = dados1, type = "fh")

plot(ena, fun = "cumhaz", xlab = "Tempo",
     ylab = expression(Lambda(t)), main = "Risco acumulado",
     conf.int = FALSE, col = "black")
lines(ajuste_exp, type = "cumhaz", col = 2, ci = FALSE)
lines(ajuste_wei, type = "cumhaz", col = 3, ci = FALSE)
lines(ajuste_lognorm, type = "cumhaz", col = 4, ci = FALSE)
lines(ajuste_gama, type = "cumhaz", col = 5, ci = FALSE)
lines(ajuste_gengama, type = "cumhaz", col = 6, ci = FALSE)

legend("topleft",
       legend = c("Nelson-Aalen", "Exponencial", "Weibull",
                  "LogNormal", "Gama", "Gama Generalizada"),
       col = 1:6, lty = 1, cex = 0.7)

# Taxa de risco
plot(ajuste_exp, type = "hazard", col = 2, ci = FALSE,
     xlab = "Tempo", ylab = expression(lambda(t)),
     main = "Taxa de risco")
lines(ajuste_wei, type = "hazard", col = 3, ci = FALSE)
lines(ajuste_lognorm, type = "hazard", col = 4, ci = FALSE)
lines(ajuste_gama, type = "hazard", col = 5, ci = FALSE)
lines(ajuste_gengama, type = "hazard", col = 6, ci = FALSE)

legend("topright",
       legend = c("Exponencial", "Weibull", "LogNormal",
                  "Gama", "Gama Generalizada"),
       col = 2:6, lty = 1, cex = 0.7)

# 2 - Probabilidade de obito antes de 800 dias ----

t <- 800

prob_800 <- data.frame(
  Modelo = c("Kaplan-Meier", "Exponencial", "Weibull",
             "LogNormal", "Gama", "Gama Generalizada"),
  Probabilidade = c(
    1 - summary(ekm, times = t)$surv,
    1 - summary(ajuste_exp, t = t)[[1]]$est,
    1 - summary(ajuste_wei, t = t)[[1]]$est,
    1 - summary(ajuste_lognorm, t = t)[[1]]$est,
    1 - summary(ajuste_gama, t = t)[[1]]$est,
    1 - summary(ajuste_gengama, t = t)[[1]]$est
  )
)

prob_800$Porcentagem <- round(100 * prob_800$Probabilidade, 2)
prob_800

# 3 - AIC e BIC ----

AIC <- c(AIC(ajuste_exp), AIC(ajuste_wei),
         AIC(ajuste_lognorm), AIC(ajuste_gama),
         AIC(ajuste_gengama))

BIC <- c(BIC(ajuste_exp), BIC(ajuste_wei),
         BIC(ajuste_lognorm), BIC(ajuste_gama),
         BIC(ajuste_gengama))

Modelos <- c("Exponencial", "Weibull", "LogNormal",
             "Gama", "Gama Generalizada")

Comparacao <- data.frame(Modelos, AIC, BIC)
Comparacao[order(Comparacao$AIC), ]
Comparacao[order(Comparacao$BIC), ]

# 4 - Teste da Razao de Verossimilhancas (TRV) ----

TRV_exp <- 2 * (ajuste_gengama$loglik - ajuste_exp$loglik)
P_valor_exp <- pchisq(TRV_exp, df = 2, lower.tail = FALSE)

TRV_wei <- 2 * (ajuste_gengama$loglik - ajuste_wei$loglik)
P_valor_wei <- pchisq(TRV_wei, df = 1, lower.tail = FALSE)

TRV_lognorm <- 2 * (ajuste_gengama$loglik - ajuste_lognorm$loglik)
P_valor_lognorm <- pchisq(TRV_lognorm, df = 1, lower.tail = FALSE)

TRV_gama <- 2 * (ajuste_gengama$loglik - ajuste_gama$loglik)
P_valor_gama <- pchisq(TRV_gama, df = 1, lower.tail = FALSE)

Resultados_TRV <- data.frame(
  Modelo = c("Exponencial", "Weibull", "LogNormal", "Gama"),
  TRV = c(TRV_exp, TRV_wei, TRV_lognorm, TRV_gama),
  P_valor = c(P_valor_exp, P_valor_wei,
              P_valor_lognorm, P_valor_gama)
)

Resultados_TRV

# 5 - Graficos de comparacao Kaplan-Meier x modelos ----

tempo_ekm <- ekm$time
surv_ekm <- ekm$surv

s_exp <- summary(ajuste_exp, t = tempo_ekm)[[1]]$est
s_wei <- summary(ajuste_wei, t = tempo_ekm)[[1]]$est
s_lognorm <- summary(ajuste_lognorm, t = tempo_ekm)[[1]]$est
s_gama <- summary(ajuste_gama, t = tempo_ekm)[[1]]$est
s_gengama <- summary(ajuste_gengama, t = tempo_ekm)[[1]]$est

par(mfrow = c(2, 3))

plot(surv_ekm, s_exp, xlab = "Kaplan-Meier", ylab = "Exponencial")
abline(0, 1, lty = 2)

plot(surv_ekm, s_wei, xlab = "Kaplan-Meier", ylab = "Weibull")
abline(0, 1, lty = 2)

plot(surv_ekm, s_lognorm, xlab = "Kaplan-Meier", ylab = "LogNormal")
abline(0, 1, lty = 2)

plot(surv_ekm, s_gama, xlab = "Kaplan-Meier", ylab = "Gama")
abline(0, 1, lty = 2)

plot(surv_ekm, s_gengama, xlab = "Kaplan-Meier", ylab = "Gama Generalizada")
abline(0, 1, lty = 2)

par(mfrow = c(1, 1))

# 9 ----

require(survival)
require(flexsurv)
require(dplyr)
require(ggplot2)
require(qqplotr)

dados1 <- read.csv("Aulas/Data/DadosAIDS2.csv")
dados1$tratam <- factor(dados1$tratam)

# 1 - Ajuste dos modelos ----

ajuste_exp <- flexsurvreg(Surv(tempo, status) ~ tratam,
                          dist = "exp", data = dados1)

ajuste_wei <- flexsurvreg(Surv(tempo, status) ~ tratam,
                          dist = "weibull", data = dados1)

ajuste_lognorm <- flexsurvreg(Surv(tempo, status) ~ tratam,
                              dist = "lognormal", data = dados1)

ajuste_gama <- flexsurvreg(Surv(tempo, status) ~ tratam,
                           dist = "gamma", data = dados1)

ajuste_gengama <- flexsurvreg(Surv(tempo, status) ~ tratam,
                              dist = "gengamma", data = dados1)

# 2 - Comparacao dos modelos: AIC e BIC ----

Modelos <- c("Exponencial", "Weibull", "LogNormal",
             "Gama", "Gama Generalizada")

Comparacao <- data.frame(
  Modelo = Modelos,
  AIC = c(AIC(ajuste_exp), AIC(ajuste_wei),
          AIC(ajuste_lognorm), AIC(ajuste_gama),
          AIC(ajuste_gengama)),
  BIC = c(BIC(ajuste_exp), BIC(ajuste_wei),
          BIC(ajuste_lognorm), BIC(ajuste_gama),
          BIC(ajuste_gengama))
)

Comparacao[order(Comparacao$AIC), ]
Comparacao[order(Comparacao$BIC), ]

# 3 - Teste da Razao de Verossimilhancas (TRV) ----

TRV_exp <- 2 * (ajuste_gengama$loglik - ajuste_exp$loglik)
P_valor_exp <- pchisq(TRV_exp, df = 2, lower.tail = FALSE)

TRV_wei <- 2 * (ajuste_gengama$loglik - ajuste_wei$loglik)
P_valor_wei <- pchisq(TRV_wei, df = 1, lower.tail = FALSE)

TRV_lognorm <- 2 * (ajuste_gengama$loglik - ajuste_lognorm$loglik)
P_valor_lognorm <- pchisq(TRV_lognorm, df = 1, lower.tail = FALSE)

TRV_gama <- 2 * (ajuste_gengama$loglik - ajuste_gama$loglik)
P_valor_gama <- pchisq(TRV_gama, df = 1, lower.tail = FALSE)

TRV <- data.frame(
  Comparacao = c("Exponencial x GG", "Weibull x GG",
                 "LogNormal x GG", "Gama x GG"),
  Estatistica = c(TRV_exp, TRV_wei, TRV_lognorm, TRV_gama),
  P_valor = c(P_valor_exp, P_valor_wei,
              P_valor_lognorm, P_valor_gama)
)

TRV

# 4 - Definir o modelo final ----
# Escolha o modelo com base em AIC, BIC, TRV e diagnosticos.

modelo_final <- ajuste_lognorm

# 5 - Residuos de Cox-Snell ----

ResiduoCS <- coxsnell_flexsurvreg(modelo_final)$est

ekm <- survfit(Surv(ResiduoCS, dados1$status) ~ 1)

par(mfrow = c(1, 2))

plot(ekm, fun = "cumhaz",
     xlab = "Residuos de Cox-Snell",
     ylab = "Risco acumulado",
     main = "Cox-Snell: risco acumulado")

abline(0, 1, col = "red", lty = 2)

plot(ekm,
     xlab = "Residuos de Cox-Snell",
     ylab = "Sobrevivencia",
     main = "Cox-Snell: sobrevivencia")

ResiduoCS <- sort(ResiduoCS)
lines(ResiduoCS, exp(-ResiduoCS), col = "red", lty = 2)

par(mfrow = c(1, 1))

# 6 - Residuos quantilicos aleatorizados ----

set.seed(123)

survival <- cdf <- quantilic <- Residuo <- numeric(nrow(dados1))

for (i in 1:nrow(dados1)) {

  survival[i] <- summary(
    modelo_final,
    type = "survival",
    newdata = dados1[i, , drop = FALSE],
    t = dados1$tempo[i]
  )[[1]]$est

  cdf[i] <- 1 - survival[i]

  if (dados1$status[i] == 1) {
    quantilic[i] <- cdf[i]
  } else {
    quantilic[i] <- runif(1, min = cdf[i], max = 1)
  }

  Residuo[i] <- qnorm(quantilic[i])
}

ggplot(data.frame(Residuo), aes(sample = Residuo)) +
  stat_qq_band(fill = "lightblue") +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Quantis teoricos") +
  ylab("Quantis amostrais") +
  theme_bw()

# 7 - Razao de tempos medianos ----
# Para modelo AFT LogNormal

coef_tratam <- coef(modelo_final)["tratam1"]

RTM <- exp(coef_tratam)

IC_RTM <- exp(
  coef_tratam +
    c(-1, 1) * 1.96 *
    sqrt(vcov(modelo_final)["tratam1", "tratam1"])
)

RTM
IC_RTM


# 8 - Graficos de sobrevivencia e risco por tratamento ----

grupos <- data.frame(
  tratam = factor(levels(dados1$tratam),
                  levels = levels(dados1$tratam))
)

par(mfrow = c(1, 2))

plot(modelo_final, newdata = grupos,
     type = "survival", ci = FALSE,
     xlab = "Tempo", ylab = "S(t)",
     main = "Sobrevivencia por tratamento",
     col = seq_len(nrow(grupos)), lwd = 2)

legend("topright", legend = levels(dados1$tratam),
       col = seq_len(nrow(grupos)), lty = 1, lwd = 2,
       title = "Tratamento")

plot(modelo_final, newdata = grupos,
     type = "hazard", ci = FALSE,
     xlab = "Tempo", ylab = expression(lambda(t)),
     main = "Risco por tratamento",
     col = seq_len(nrow(grupos)), lwd = 2)

legend("topright", legend = levels(dados1$tratam),
       col = seq_len(nrow(grupos)), lty = 1, lwd = 2,
       title = "Tratamento")

par(mfrow = c(1, 1))

# 9 - Tempos medios e medianos por grupo ----
# Para o modelo LogNormal

mu0 <- coef(modelo_final)["meanlog"]
beta <- coef(modelo_final)["tratam1"]
sigma <- exp(coef(modelo_final)["sdlog"])

# Grupo de referencia
mu <- c(mu0, mu0 + beta)

Tabela_Tempos <- data.frame(
  Tratamento = levels(dados1$tratam),
  Mediana = exp(mu),
  Media = exp(mu + sigma^2 / 2)
)

Tabela_Tempos