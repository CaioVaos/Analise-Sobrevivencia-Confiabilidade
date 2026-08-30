# Cap. 7 Exercicios

# Setup ----

require(dplyr)
require(ggplot2)
require(survival)
require(ggsurvfit)
require(survRM2)


# 1 ----

tempos <- c(
  0.19, 0.78, 0.96, 1.31, 2.78, 3.16, 4.67, 4.85,
  6.50, 7.35, 8.27, 12.07, 32.52, 33.91, 36.71,
  rep(36.71, 10)
)

cens <- c(rep(1, 15), rep(0, 10))

Dados2 <- data.frame(tempos, cens)


## Kaplan-Meier ----

ekm_Dados2 <- survfit2(
  Surv(tempos, cens) ~ 1,
  data = Dados2,
  conf.type = "log-log",
  conf.int = 0.95
)

summary(ekm_Dados2)

ggsurvfit(ekm_Dados2, linewidth = 1.2) +
  add_confidence_interval() +
  add_risktable() +
  scale_ggsurvfit() +
  labs(
    x = "Tempo (minutos)",
    y = "Sobrevivência estimada"
  ) +
  theme_classic(base_size = 16)

# 2 ----

## a) ----

Dados3 <- readr::read_csv("Aulas/Data/DadosAIDS.csv")


## b) ----

ekm_Dados3_sem_estrat <- survfit2(
  Surv(tempo, status) ~ 1,
  data = Dados3,
  conf.type = "log-log",
  conf.int = 0.95
)

summary(ekm_Dados3_sem_estrat)

ggsurvfit(ekm_Dados3_sem_estrat, linewidth = 1.2) +
  add_confidence_interval() +
  add_risktable() +
  scale_ggsurvfit() +
  labs(
    x = "Tempo (dias)",
    y = "Sobrevivência estimada"
  ) +
  theme_classic(base_size = 16)


# Sobrevivência em 1000 e 1750 dias

summary(
  ekm_Dados3_sem_estrat,
  times = c(1000, 1750)
)

# Resultados:
# P(T > 1000) = 0.561
# P(T > 1750) = 0.443

# Portanto:
# P(T < 1750) = 1 - 0.443 = 0.557

# Interpretação:
# Aproximadamente 56,1% dos pacientes sobrevivem por mais
# de 1000 dias.
#
# Aproximadamente 55,7% dos pacientes morrem antes de 1750 dias.


## c) ----

ekm_Dados3_com_estrat <- survfit2(
  Surv(tempo, status) ~ sexo,
  data = Dados3,
  conf.type = "log-log",
  conf.int = 0.95
)

summary(ekm_Dados3_com_estrat)

ggsurvfit(ekm_Dados3_com_estrat, linewidth = 1.2) +
  add_confidence_interval() +
  add_risktable() +
  scale_ggsurvfit() +
  labs(
    x = "Tempo (dias)",
    y = "Sobrevivência estimada"
  ) +
  theme_classic(base_size = 16)


# Sobrevivência em 1000 e 1750 dias

summary(
  ekm_Dados3_com_estrat,
  times = c(1000, 1750)
)

# Feminino:
# P(T > 1000) = 0.674
# P(T < 1750) = 1 - 0.626 = 0.374

# Masculino:
# P(T > 1000) = 0.522
# P(T < 1750) = 1 - 0.390 = 0.610

# Interpretação:
# Em 1000 dias, a sobrevivência estimada é maior entre as mulheres
# (67,4%) do que entre os homens (52,2%).
#
# Até 1750 dias, a probabilidade estimada de morte é menor
# entre as mulheres (37,4%) do que entre os homens (61,0%).
#
# Descritivamente, os resultados sugerem maior sobrevivência
# para as mulheres.

## d) ----


## e) ----

# Tempo até que 20% dos pacientes tenham apresentado o evento

quantile(
  ekm_Dados3_com_estrat,
  probs = 0.2
)

p <- 0.2

ggsurvfit(ekm_Dados3_com_estrat, linewidth = 1.2) +
  add_confidence_interval() +
  add_quantile(1 - p) +
  add_risktable() +
  scale_ggsurvfit() +
  labs(
    x = "Tempo (dias)",
    y = "Sobrevivência estimada"
  ) +
  theme_classic(base_size = 16)

# Resultados:
# Feminino: 418 dias (IC95%: 145–944)
# Masculino: 285 dias (IC95%: 152–419)

# Interpretação:
# O tempo em que aproximadamente 20% dos pacientes já morreram
# é menor entre os homens (285 dias) do que entre as mulheres
# (418 dias), sugerindo menor sobrevivência inicial entre os homens.


## f) ----

# Mediana da sobrevivência

quantile(
  ekm_Dados3_com_estrat,
  probs = 0.5
)

p <- 0.5

ggsurvfit(ekm_Dados3_com_estrat, linewidth = 1.2) +
  add_confidence_interval() +
  add_quantile(1 - p) +
  add_risktable() +
  scale_ggsurvfit() +
  labs(
    x = "Tempo (dias)",
    y = "Sobrevivência estimada"
  ) +
  theme_classic(base_size = 16)

# Resultados:
# Feminino: mediana não estimável (NA)
# Masculino: 1116 dias (IC95%: 855–1506)

# Interpretação:
# Para os homens, a mediana estimada é de 1116 dias, ou seja,
# aproximadamente metade dos homens já apresentou o evento
# até esse momento.
#
# Para as mulheres, a mediana não é estimável porque a curva
# de sobrevivência não caiu abaixo de 50% durante o período
# observado.


## g) ----

# Nelson-Aalen – risco acumulado

fit_na <- survfit2(
  Surv(tempo, status) ~ sexo,
  data = Dados3,
  conf.type = "log",
  type = "fh"
)

summary(fit_na, times = 1500)


## h) ----

ggsurvfit(
  fit_na,
  linewidth = 1.2,
  type = "cumhaz"
) +
  add_confidence_interval() +
  labs(
    x = "Tempo (dias)",
    y = "Risco acumulado estimado"
  ) +
  theme_classic(base_size = 16)

# Interpretação:
# O risco acumulado cresce ao longo do tempo, pois os eventos
# vão sendo acumulados durante o acompanhamento.
#
# Como a curva masculina esta acima da feminina, isso indica
# maior risco acumulado de morte entre os homens.


## i) ----

Tempo <- Dados3$tempo
Cens <- Dados3$status
Grupo <- as.numeric(as.factor(Dados3$sexo)) - 1

rmst_resultado <- rmst2(Tempo, Cens, Grupo, tau = 3000)

rmst_resultado

# Interpretação:
# O RMST representa o tempo médio de sobrevivência restrito
# aos primeiros 3000 dias de acompanhamento.