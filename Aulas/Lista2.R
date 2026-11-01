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
dados8 <- read.csv("Aulas/Data/DadosAIDS_cap13.csv")

## Ajuste dos modelos ----

ajuste_exp     <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "exp",      data = dados8)
ajuste_wei     <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "weibull",  data = dados8)
ajuste_lognorm <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "lnorm",    data = dados8)
ajuste_gama    <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "gamma",    data = dados8)
ajuste_gengama <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "gengamma", data = dados8)

ajustes8 <- list(
  "Exponencial"       = ajuste_exp,
  "Weibull"           = ajuste_wei,
  "Log-Normal"        = ajuste_lognorm,
  "Gama"              = ajuste_gama,
  "Gama Generalizada" = ajuste_gengama
)
cores <- 2:6

## 1) S(t), risco e risco acumulado x não paramétrico ----
# (o plot do flexsurv já traz Kaplan-Meier / risco suavizado / Nelson-Aalen)

tipos <- c("survival", "hazard", "cumhaz")
ylabs <- c(expression(hat(S)(t)), expression(hat(lambda)(t)), expression(hat(Lambda)(t)))

par(mfrow = c(1, 3))
for (j in 1:3) {
  plot(ajustes8[[1]], type = tipos[j], xlab = "Tempo (dias)", ylab = ylabs[j],
       lwd.ci = 0, lty.ci = 0, col = cores[1])
  for (i in 2:5)
    lines(ajustes8[[i]], type = tipos[j], lwd.ci = 0, lty.ci = 0, col = cores[i])
  legend("topright", legend = names(ajustes8), fill = cores, cex = 0.7)
}
par(mfrow = c(1, 1))

## 2) P(morrer antes de 800 dias) = 1 - S(800) ----

ekm8 <- survfit(Surv(tempos, cens) ~ 1, data = dados8)

prob800 <- data.frame(
  Modelo = c("Kaplan-Meier", names(ajustes8)),
  Prob   = c(1 - summary(ekm8, times = 800)$surv,
             sapply(ajustes8, function(a) 1 - summary(a, t = 800)[[1]]$est)),
  row.names = NULL
)
prob800$Porcentagem <- round(100 * prob800$Prob, 2)
prob800

## 3) Escolha do modelo ----

### AIC e BIC
Comparacao <- data.frame(
  Modelos = names(ajustes8),
  AIC = sapply(ajustes8, function(a) a$AIC),
  BIC = sapply(ajustes8, function(a) a$BIC),
  row.names = NULL
)
Comparacao[order(Comparacao$AIC), ]
Comparacao[order(Comparacao$BIC), ]

### TRV (modelos x Gama Generalizada)
TRV_exp    <- 2 * (ajuste_gengama$loglik - ajuste_exp$loglik)
TRV_wei    <- 2 * (ajuste_gengama$loglik - ajuste_wei$loglik)
TRV_lognorm <- 2 * (ajuste_gengama$loglik - ajuste_lognorm$loglik)
TRV_gama   <- 2 * (ajuste_gengama$loglik - ajuste_gama$loglik)

data.frame(
  Comparacao = c("Exponencial x GG", "Weibull x GG", "LogNormal x GG", "Gama x GG"),
  TRV     = c(TRV_exp, TRV_wei, TRV_lognorm, TRV_gama),
  P_valor = c(pchisq(TRV_exp,    2, lower.tail = FALSE),
              pchisq(TRV_wei,    1, lower.tail = FALSE),
              pchisq(TRV_lognorm, 1, lower.tail = FALSE),
              pchisq(TRV_gama,   1, lower.tail = FALSE))
)

### Gráfico: Kaplan-Meier x modelo (pontos próximos da reta = bom ajuste)
Surv_ekm <- data.frame(Tempo = ekm8$time, Surv_EKM = ekm8$surv)

par(mfrow = c(2, 3))
for (nm in names(ajustes8)) {
  s_par <- summary(ajustes8[[nm]], t = Surv_ekm$Tempo)[[1]]$est
  plot(Surv_ekm$Surv_EKM, s_par, xlab = "S(t): Kaplan-Meier", ylab = paste("S(t):", nm))
  abline(0, 1)
}
par(mfrow = c(1, 1))

### Gráfico: linearização (Exponencial, Weibull e Log-Normal)
ok <- Surv_ekm$Surv_EKM > 0 & Surv_ekm$Surv_EKM < 1
S  <- Surv_ekm$Surv_EKM[ok]
t  <- Surv_ekm$Tempo[ok]

par(mfrow = c(1, 3))
plot(t,      -log(S),      xlab = "t",      ylab = "-log(S(t))",      main = "Exponencial")
plot(log(t), log(-log(S)), xlab = "log(t)", ylab = "log(-log(S(t)))", main = "Weibull")
plot(log(t), qnorm(S),     xlab = "log(t)", ylab = "Phi^-1(S(t))",    main = "Log-Normal")
par(mfrow = c(1, 1))

# Conclusão:
# (modelo apoiado por gráficos + AIC/BIC + TRV; em empate, o mais simples)


# 9 ----
# dados9 <- read.csv("Aulas/Data/AIDS_cap15.csv")        # colunas: tempo, status, tratam

dados9$tratam <- factor(dados9$tratam)
table(dados9$tratam)

## Teste prévio: tratam discrimina as curvas? ----

survdiff(Surv(tempo, status) ~ tratam, data = dados9)                         # LogRank
logrank_test(Surv(tempo, status) ~ tratam, type = "Peto-Peto", data = dados9) # Peto-Peto

## 1) Ajuste dos 5 modelos ----

ajuste_exp     <- flexsurvreg(Surv(tempo, status) ~ tratam, dist = "exponential", data = dados9)
ajuste_wei     <- flexsurvreg(Surv(tempo, status) ~ tratam, dist = "weibull",     data = dados9)
ajuste_lognorm <- flexsurvreg(Surv(tempo, status) ~ tratam, dist = "lognormal",   data = dados9)
ajuste_gama    <- flexsurvreg(Surv(tempo, status) ~ tratam, dist = "gamma",       data = dados9)
ajuste_gengama <- flexsurvreg(Surv(tempo, status) ~ tratam, dist = "gengamma",    data = dados9)

ajustes9 <- list(
  "Exponencial"       = ajuste_exp,
  "Weibull"           = ajuste_wei,
  "Log-Normal"        = ajuste_lognorm,
  "Gama"              = ajuste_gama,
  "Gama Generalizada" = ajuste_gengama
)

tidy(ajuste_exp)
tidy(ajuste_wei)
tidy(ajuste_lognorm)
tidy(ajuste_gama)
tidy(ajuste_gengama)

## 2) Escolha do modelo: AIC, BIC e TRV ----

Comparacao <- data.frame(
  Modelos = names(ajustes9),
  AIC = sapply(ajustes9, function(a) a$AIC),
  BIC = sapply(ajustes9, function(a) a$BIC),
  row.names = NULL
)
Comparacao[order(Comparacao$AIC), ]
Comparacao[order(Comparacao$BIC), ]

TRV_exp     <- 2 * (ajuste_gengama$loglik - ajuste_exp$loglik)
TRV_wei     <- 2 * (ajuste_gengama$loglik - ajuste_wei$loglik)
TRV_lognorm <- 2 * (ajuste_gengama$loglik - ajuste_lognorm$loglik)
TRV_gama    <- 2 * (ajuste_gengama$loglik - ajuste_gama$loglik)

data.frame(
  Comparacao = c("Exponencial x GG", "Weibull x GG", "LogNormal x GG", "Gama x GG"),
  TRV     = c(TRV_exp, TRV_wei, TRV_lognorm, TRV_gama),
  P_valor = c(pchisq(TRV_exp,     2, lower.tail = FALSE),
              pchisq(TRV_wei,     1, lower.tail = FALSE),
              pchisq(TRV_lognorm, 1, lower.tail = FALSE),
              pchisq(TRV_gama,    1, lower.tail = FALSE))
)

# Modelo final (troque o nome conforme AIC/BIC/TRV e resíduos)
nome_final   <- "Log-Normal"
modelo_final <- ajustes9[[nome_final]]
modelo_final

## 3) Resíduos ----

### Cox-Snell
ResiduoCS <- coxsnell_flexsurvreg(modelo_final)$est

par(mfrow = c(1, 2))
ekm_cs <- survfit(Surv(ResiduoCS, dados9$status) ~ 1, type = "kaplan-meier")
plot(ekm_cs, fun = "cumhaz", xlab = "Resíduos de Cox-Snell", ylab = "Função de risco acumulado")
abline(0, 1, col = "red")

plot(ekm_cs, xlab = "Resíduos de Cox-Snell", ylab = "Função de sobrevivência")
ResiduoCS <- sort(ResiduoCS)
lines(ResiduoCS, exp(-ResiduoCS), col = "red")
par(mfrow = c(1, 1))

### Quantílico aleatorizado
set.seed(123)
Residuo <- numeric(nrow(dados9))

for (i in 1:nrow(dados9)) {
  surv_i <- summary(modelo_final, type = "survival",
                    newdata = dados9[i, , drop = FALSE],
                    t = dados9$tempo[i])[[1]]$est
  cdf_i <- 1 - surv_i
  u <- if (dados9$status[i] == 1) cdf_i else runif(1, min = cdf_i, max = 1)
  Residuo[i] <- qnorm(u)
}

ggplot(data.frame(Residuo), aes(sample = Residuo)) +
  stat_qq_band(fill = "lightblue") +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Quantis Teóricos") +
  ylab("Quantis Amostrais") +
  theme_bw()

## 4) Distribuição aproximada do tempo de vida ----
# Log-Normal: T ~ LogNormal(mu_i = b0 + b1*B (+ b2*C + b3*D), sigma)
# Weibull:    T ~ Weibull(kappa_i = exp(b0 + b1*B + ...), gamma)
# Exp./Gama:  rate_i = exp(b0 + ...) | Gama Gen.: mu_i = b0 + ...
modelo_final$res
coef(modelo_final)

## 5) Razão de tempos medianos (RTM) ----
# Exponencial e Gama: exp(-coef), com intervalo invertido. Demais: exp(coef).

sinal <- ifelse(nome_final %in% c("Exponencial", "Gama"), -1, 1)
termos <- grep("^tratam", names(coef(modelo_final)), value = TRUE)

Est <- exp(sinal * coef(modelo_final)[termos])
IC  <- exp(sinal * confint(modelo_final)[termos, , drop = FALSE])
if (sinal == -1) IC <- IC[, 2:1, drop = FALSE]
colnames(IC) <- c("2.5 %", "97.5 %")
cbind(RTM = Est, IC)

# Interpretação: RTM > 1 -> tempo mediano maior que o do grupo de referência (A);
# RTM < 1 -> menor. Significativo se o IC não contém 1.

## 6) Sobrevivência e risco por tratamento ----

grupos <- levels(dados9$tratam)
tgraf  <- seq(1, max(dados9$tempo), length.out = 100)

par(mfrow = c(1, 2))
plot(modelo_final, newdata = list(tratam = grupos[1]), type = "survival", t = tgraf,
     col = 1, col.obs = "white", ci = FALSE,
     xlab = "Dias", ylab = "Função de sobrevivência estimada")
for (i in 2:length(grupos))
  lines(modelo_final, newdata = list(tratam = grupos[i]), type = "survival",
        t = tgraf, col = i, ci = FALSE)
legend("topright", legend = grupos, col = seq_along(grupos), lty = 1, title = "Tratamento")

plot(modelo_final, newdata = list(tratam = grupos[1]), type = "hazard", t = tgraf,
     col = 1, col.obs = "white", ci = FALSE,
     xlab = "Dias", ylab = "Função de risco estimada")
for (i in 2:length(grupos))
  lines(modelo_final, newdata = list(tratam = grupos[i]), type = "hazard",
        t = tgraf, col = i, ci = FALSE)
legend("topright", legend = grupos, col = seq_along(grupos), lty = 1, title = "Tratamento")
par(mfrow = c(1, 1))

## 7) Tempos médio e mediano por grupo ----

resumo <- do.call(rbind, lapply(grupos, function(g) {
  nd  <- list(tratam = g)
  med <- summary(modelo_final, newdata = nd, type = "median")[[1]]
  mea <- tryCatch(summary(modelo_final, newdata = nd, type = "mean")[[1]],
                  error = function(e) data.frame(est = NA, lcl = NA, ucl = NA))
  data.frame(Tratamento = g,
             Media = mea$est, Media_L95 = mea$lcl, Media_U95 = mea$ucl,
             Mediana = med$est, Mediana_L95 = med$lcl, Mediana_U95 = med$ucl)
}))
resumo