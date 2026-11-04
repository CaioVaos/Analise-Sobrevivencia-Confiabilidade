# Lista 2 - Analise de Sobrevivencia e Confiabilidade
# (rodar com a raiz do projeto como diretório de trabalho: dados em Aulas/Data/)

# Setup ----

# install.packages("survival")
# install.packages("flexsurv")
# install.packages("ggplot2")
# install.packages("qqplotr")
# install.packages("coin")

library(survival)
library(flexsurv)
library(ggplot2)
library(qqplotr)
library(coin)

# 1 ----
#
## EXPONENCIAL:
# Risco constante: h(t) = alpha. Sem memória: não depende do tempo já decorrido.
#
## WEIBULL (forma gamma):
# Risco monótono:
# gamma < 1: decrescente
# gamma = 1: constante (Exponencial)
# gamma > 1: crescente
#
## LOG-NORMAL:
# Risco não monótono (unimodal):
# parte de 0, cresce até um máximo e depois decresce, tendendo a 0.
#
## GAMA (forma alpha, taxa beta):
# Risco monótono, convergindo para beta quando t -> infinito:
# alpha < 1: decrescente
# alpha = 1: constante (Exponencial)
# alpha > 1: crescente

# 2 ----
#
## GAMA GENERALIZADO (3 parâmetros):
# Família guarda-chuva: contém Exponencial, Weibull e Gama como casos
# particulares e a Log-Normal como caso limite
# (flexsurv: Weibull Q = 1 | Gama Q = sigma | Exp. Q = sigma = 1 | Log-Normal Q -> 0).
#
# Papel: modelo de referência. Ajusta-se o GG e testa-se, pelo TRV, se um
# modelo mais simples é suficiente (se não rejeita H0, prefere-se o mais simples).
#
# Risco muito flexível: crescente, decrescente, constante, unimodal ou banheira.
# Custo: mais parâmetros -> maior risco de sobreajuste e instabilidade numérica.

# 3 ----

alpha <- 100
gam   <- 2

S_weibull <- function(t) exp(-(t / alpha)^gam)

## a) ----

S_weibull(c(30, 45))
pweibull(c(30, 45), shape = gam, scale = alpha, lower.tail = FALSE)  # conferência
# 91,4% dos ratos estão livres de tumor aos 30 dias
# 81,7% dos ratos estão livres de tumor aos 45 dias

## b) ----

alpha * gamma(1 + 1 / gam)
# 88,62 dias

## c) ----

alpha * log(2)^(1 / gam)
qweibull(0.5, shape = gam, scale = alpha)  # conferência
# 83,26 dias

## d) ----

h_weibull <- function(t) (gam / alpha) * (t / alpha)^(gam - 1)
h_weibull(c(30, 45, 60))
# 0.006 0.009 0.012

# Interpretação:
# gamma = 2 > 1: risco crescente (h(t) = 0,0002*t).
# Entre os ratos ainda sem tumor no instante t, a taxa instantânea de
# aparecimento do tumor é de cerca de 0,6% por dia aos 30 dias,
# 0,9% aos 45 e 1,2% aos 60: quanto maior a exposição, maior a propensão ao tumor.

# 4 ----
# Papel

# 5 ----
#
# Nos modelos Weibull, Log-Normal e Gama, a log-verossimilhança é não linear
# nos parâmetros (ex.: (t/kappa)^gamma na Weibull; Phi(.) na Log-Normal;
# Gamma(alpha) e gama incompleta na Gama). O sistema de equações de
# verossimilhança U(theta) = 0 não permite isolar os parâmetros por fórmulas
# explícitas. Na Exponencial há um só parâmetro e a equação se resolve direto.
#
# Solução: métodos numéricos iterativos (Newton-Raphson, BFGS...), que partem de
# valores iniciais e maximizam l(theta) até convergir. O flexsurv usa
# optim (BFGS) e faz isso automaticamente.

# 6 ----
#
## 1) KAPLAN-MEIER x MODELO PARAMÉTRICO (gráfico):
# Plota S_KM(t) contra S_modelo(t) nos mesmos tempos.
# Pontos sobre a reta y = x: bom ajuste.
# Descarta bem modelos ruins, mas é subjetivo.
#
## 2) LINEARIZAÇÃO (gráfico):
# Transforma S(t) de modo que fique linear se o modelo for adequado:
# Exponencial:    -log(S(t))          x t
# Weibull:        log(-log(S(t)))     x log(t)
# Log-Normal:     Phi^(-1)(S(t))      x log(t)
# Log-Logística:  log((1-S(t))/S(t))  x log(t)
# Quanto mais próximo de uma reta, melhor. Não existe para Gama e GG;
# também é subjetivo.
#
## 3) AIC E BIC:
# AIC = -2*l(theta_hat) + 2*k
# BIC = -2*l(theta_hat) + k*log(n)
# Premiam o ajuste e penalizam o nº de parâmetros k (BIC penaliza mais se n > 7).
# Escolhe-se o MENOR valor. Servem para modelos aninhados ou não, mas são
# apenas relativos: se todos os modelos forem ruins, não avisam.
#
## 4) TESTE DA RAZÃO DE VEROSSIMILHANÇAS (TRV):
# Só para modelos aninhados (modelo simples x Gama Generalizada).
# Lambda = 2*[l(complexo) - l(simples)] ~ qui-quadrado(gl = p - q)
# H0: o modelo simples é suficiente.
# p-valor > alpha: não rejeita H0, prefere o simples.
# p-valor <= alpha: rejeita H0, prefere o completo.
# Não serve, por ex., para Weibull x Log-Normal ou para a Log-Logística.

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
ajuste_wei

### Estimativas
ajuste_wei$res
# shape (gamma) = 1.281
# scale (kappa) = 1993.2
#
# Obtidas por máxima verossimilhança: L = prod f(t_i)^d_i * S(t_i)^(1-d_i)
# (a censura entra via S), maximizada numericamente (optim/BFGS) na escala log
# dos parâmetros, a partir de valores iniciais automáticos.

### Erros-padrão
# Matriz de covariâncias = inversa da (-Hessiana) da log-verossimilhança no EMV
# (normalidade assintótica do EMV). Calculada na escala log (res.t e cov);
# na escala original, método delta: se(theta) = theta_hat * se(log(theta_hat)).
ep_log <- sqrt(diag(ajuste_wei$cov))
ep_log
ajuste_wei$res[, "est"] * ep_log
# shape = 0.161 ; scale = 231.9  (iguais à coluna se da saída)

### Intervalos de confiança
# Wald na escala log (garante parâmetros > 0), depois exponencial:
# exp( log(theta_hat) +- z_0.975 * se_log )
cbind(inf = exp(ajuste_wei$res.t[, "est"] - qnorm(0.975) * ep_log),
      sup = exp(ajuste_wei$res.t[, "est"] + qnorm(0.975) * ep_log))
# shape: [1.001 ; 1.640]   scale: [1586.8 ; 2503.7]
# (limite inferior do shape muito próximo de 1: a Exponencial fica no limiar)

### Tempo mediano de falha
# Princípio da invariância: mediana = kappa_hat * log(2)^(1/gamma_hat)
ajuste_wei$res["scale", "est"] * log(2)^(1 / ajuste_wei$res["shape", "est"])
# 1497.4 horas

set.seed(123)
summary(ajuste_wei, type = "median", cl = 0.95, se = TRUE)
# IC: o flexsurv simula B = 1000 vetores de parâmetros de N(theta_hat, cov),
# calcula a mediana em cada um e toma os percentis 2,5% e 97,5%.
# Por ser simulação, o IC muda um pouco entre execuções (por isso a semente).
# IC ~ [1170 ; 1900] horas

## b) ----

ajuste_exp     <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "exp",      data = dados)
ajuste_lognorm <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "lnorm",    data = dados)
ajuste_llogis  <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "llogis",   data = dados)
ajuste_gama    <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "gamma",    data = dados)
ajuste_gengama <- flexsurvreg(Surv(tempos, cens) ~ 1, dist = "gengamma", data = dados)

sum(dados$cens) / sum(dados$tempos)
ajuste_exp$res["rate", "est"]
# 0.0004956 nos dois

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
  BIC    = sapply(ajustes7, function(a) a$BIC),
  row.names = NULL
)
# Menores AIC e BIC: Log-Normal (769.4 / 773.5),
# seguida da Log-Logística (770.5 / 774.7).

### TRV (x Gama Generalizada; a Log-Logística não é aninhada)
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
# Rejeitados: Exponencial (p = 0.006), Weibull (0.010), Gama (0.024).
# Não rejeitado: Log-Normal (p = 0.449).

# Conclusão: Log-Normal.
# Critério principal: AIC/BIC (menor valor), confirmado pelo TRV.

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
# Descritivamente: a Exponencial é a que mais se afasta da estimativa não
# paramétrica; Log-Normal e Gama Generalizada são as que mais se aproximam.

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
# KM 34,3% | Exp 32,8% | Weibull 35,5% | Log-Normal 37,9% | Gama 34,8% | GG 37,9%

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
# Menores AIC e BIC: Log-Normal (1539.8 / 1546.3)

### TRV (modelos x Gama Generalizada)
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
# Rejeitados: Exponencial (p = 0.001), Weibull (0.026), Gama (0.012).
# Não rejeitado: Log-Normal (p = 0.936).

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
ok   <- Surv_ekm$Surv_EKM > 0 & Surv_ekm$Surv_EKM < 1
s_km <- Surv_ekm$Surv_EKM[ok]
t_km <- Surv_ekm$Tempo[ok]

par(mfrow = c(1, 3))
plot(t_km,      -log(s_km),      xlab = "t",      ylab = "-log(S(t))",      main = "Exponencial")
plot(log(t_km), log(-log(s_km)), xlab = "log(t)", ylab = "log(-log(S(t)))", main = "Weibull")
plot(log(t_km), qnorm(s_km),     xlab = "log(t)", ylab = "Phi^-1(S(t))",    main = "Log-Normal")
par(mfrow = c(1, 1))

# Conclusão: Log-Normal.
# Apoiado por AIC/BIC (menores) e TRV (não rejeitado contra a GG);
# nos gráficos a Exponencial é a que pior se ajusta.
# Obs.: em t = 800 a Log-Normal superestima um pouco a probabilidade de óbito
# (37,9% contra 34,3% da KM), embora seja o melhor ajuste global.

# 9 ----
dados9 <- read.csv("Aulas/Data/AIDS_cap15.csv")   # colunas: tempo, status, tratam

dados9$tratam <- factor(dados9$tratam)            # referência: A (nenhum)
table(dados9$tratam)

## Teste prévio: tratam discrimina as curvas? ----

survdiff(Surv(tempo, status) ~ tratam, data = dados9)                         # LogRank
logrank_test(Surv(tempo, status) ~ tratam, type = "Peto-Peto", data = dados9) # Peto-Peto
# LogRank p = 7e-15 | Peto-Peto p = 2.4e-10: rejeitam igualdade das curvas
# -> tratam entra na regressão.

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
# Menores AIC e BIC: Log-Normal (1489.2 / 1505.5)

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
# Rejeitados: Exponencial (p = 0.019), Weibull (0.016), Gama (0.008).
# Não rejeitado: Log-Normal (p = 0.234).

# Modelo final: Log-Normal (menor AIC e BIC, e não rejeitado pelo TRV)
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

modelo_final$res
# T_i ~ LogNormal(mu_i = 5.518 + 1.825*B + 3.124*C + 2.626*D ; sigma = 1.537)
# (A = referência; B, C, D = indicadoras do tratamento)

## 5) Razão de tempos medianos (RTM) ----
# Exponencial e Gama: exp(-coef), com intervalo invertido. Demais: exp(coef).

sinal <- ifelse(nome_final %in% c("Exponencial", "Gama"), -1, 1)
termos <- grep("^tratam", names(coef(modelo_final)), value = TRUE)

Est <- exp(sinal * coef(modelo_final)[termos])
IC  <- exp(sinal * confint(modelo_final)[termos, , drop = FALSE])
if (sinal == -1) IC <- IC[, 2:1, drop = FALSE]
colnames(IC) <- c("2.5 %", "97.5 %")
cbind(RTM = Est, IC)
# B: 6.21 [3.43 ; 11.21]
# C: 22.75 [9.04 ; 57.23]
# D: 13.82 [3.38 ; 56.57]

# Interpretação: em relação ao grupo A (sem terapia), o tempo mediano de
# sobrevivência é ~6,2 vezes maior com B (mono), ~22,7 vezes com C (combinada)
# e ~13,8 vezes com D (potente). Todas significativas (IC não contém 1).
# C e D têm IC muito amplos (D tem só 1 óbito em 14 pacientes).

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
# A sobrevivência é menor no A e maior nos grupos com terapia.
# Risco unimodal em todos; o do grupo A é o mais alto e tem pico mais cedo.

## 7) Tempos médio e mediano por grupo ----

set.seed(123)
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
# Mediana (dias): A 249 | B 1546 | C 5666 | D 3443
# Média (dias):   A 812 | B 5035 | C 18460 | D 11218
# Cuidado: a média da Log-Normal é muito sensível à cauda e as medianas de
# C e D passam do tempo máximo de acompanhamento (3228 dias): extrapolação.
# Preferir a mediana na interpretação.