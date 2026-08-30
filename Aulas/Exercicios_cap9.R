# Cap. 9 Exercicios

# Setup ----

require(dplyr)
require(ggplot2)
require(survival)
require(ggsurvfit)
require(survRM2)
require(coin)
require(survminer)
require(multcompView)

# 1 ----

tempos <- c(
  28, 89, 175, 195, 309, 377, 393, 421,
  447, 462, 709, 744, 770, 1106, 1206,
  34, 88, 137, 199, 280, 291, 299, 300, 309,
  351, 358, 369, 369, 370, 375, 382, 392,
  429, 451, 1119
)

cens <- c(
  1,1,1,1,1,0,0,0,
  0,1,0,0,0,0,0,
  1,1,1,1,1,1,0,0,1,
  1,1,1,1,1,1,1,1,
  0,1,0
)

grupos <- c(
  rep("Tumor Grande", 15),
  rep("Tumor Pequeno", 20)
)

DadosO <- data.frame(tempos, cens, grupos)


## 1.1 Kaplan-Meier ----

ekm_O <- survfit2(
  Surv(tempos, cens) ~ grupos,
  data = DadosO,
  conf.type = "log-log",
  conf.int = 0.95
)

summary(ekm_O)

ggsurvfit(ekm_O, linewidth = 1.2) +
  add_confidence_interval() +
  add_risktable() +
  add_censor_mark(size = 2, alpha = 1) +
  scale_ggsurvfit() +
  labs(
    x = "Tempo (dias)",
    y = "Sobrevivência estimada"
  ) +
  theme_classic(base_size = 16)

## 1.2 Teste Log-Rank ----

survdiff(
  Surv(tempos, cens) ~ grupos,
  data = DadosO
)

## 1.3 Teste Peto-Peto ----

logrank_test(
  Surv(tempos, cens) ~ factor(grupos),
  data = DadosO,
  type = "Peto-Peto"
)

# 2 ----

tempos <- c(
  1, 4, 5, 6, 7, 7,
  3, 5, 5, 5, 6,
  1, 3, 4, 7, 7, 7,
  3, 5, 7, 7, 7,
  3, 5, 5, 7, 7, 7
)

cens <- c(
  1,1,1,1,0,0,
  1,1,0,0,1,
  1,1,0,1,1,0,
  1,0,1,0,0,
  1,1,1,0,0,0
)

grupos <- c(
  rep("A", 6),
  rep("B", 5),
  rep("C", 6),
  rep("D", 5),
  rep("Controle", 6)
)

DadosL <- data.frame(tempos, cens, grupos)


## 2.1 Kaplan-Meier ----

ekm_L <- survfit2(
  Surv(tempos, cens) ~ grupos,
  data = DadosL,
  conf.type = "log-log",
  conf.int = 0.95
)

summary(ekm_L)

ggsurvfit(ekm_L, linewidth = 1.2) +
  add_confidence_interval() +
  add_risktable() +
  add_censor_mark(size = 2, alpha = 1) +
  scale_ggsurvfit() +
  labs(
    x = "Tempo (meses)",
    y = "Sobrevivência estimada"
  ) +
  theme_classic(base_size = 16)

## 2.2 Teste Log-Rank ----

survdiff(
  Surv(tempos, cens) ~ grupos,
  data = DadosL
)

## 2.3 Teste Peto-Peto ----

logrank_test(
  Surv(tempos, cens) ~ factor(grupos),
  data = DadosL,
  type = "Peto-Peto"
)

# 3 ----

DadosAIDS <- readr::read_csv("Aulas/Data/DadosAIDS2.csv")

## a) ----

ekm_AIDS_sexo <- survfit2(
  Surv(tempo, status) ~ sexo,
  data = DadosAIDS,
  conf.type = "log-log",
  conf.int = 0.95
)

summary(ekm_AIDS_sexo)

ggsurvfit(ekm_AIDS_sexo, linewidth = 1.2) +
  add_confidence_interval() +
  add_risktable() +
  add_censor_mark(size = 2, alpha = 1) +
  scale_ggsurvfit() +
  labs(
    x = "Tempo (dias)",
    y = "Sobrevivência estimada"
  ) +
  theme_classic(base_size = 16)

# Interpretação:
# A curva feminina apresenta maior sobrevivência estimada
# ao longo do acompanhamento, embora os intervalos de confiança
# sejam amplos. Descritivamente, as curvas parecem diferentes.


## b)  ----

# Teste Log-Rank

survdiff(
  Surv(tempo, status) ~ sexo,
  data = DadosAIDS
)

logrank_test(
  Surv(tempo, status) ~ factor(sexo),
  data = DadosAIDS,
  type = "Peto-Peto"
)

# Pelo Log-Rank, rejeitamos H0 (p = 0,04), indicando diferença significativa
# entre as curvas de sobrevivência por sexo.
# Pelo Peto-Peto, não rejeitamos H0 (p = 0,0547), não havendo diferença
# significativa ao nível de 5%. Os testes podem apresentar resultados diferentes
# por atribuírem pesos distintos aos eventos ao longo do tempo

## c) ----


ekm_AIDS_trat <- survfit2(
  Surv(tempo, status) ~ tratam,
  data = DadosAIDS,
  conf.type = "log-log",
  conf.int = 0.95
)

summary(ekm_AIDS_trat)

ggsurvfit(ekm_AIDS_trat, linewidth = 1.2) +
  add_confidence_interval() +
  add_risktable() +
  add_censor_mark(size = 2, alpha = 1) +
  scale_ggsurvfit() +
  labs(
    x = "Tempo (dias)",
    y = "Sobrevivência estimada"
  ) +
  theme_classic(base_size = 16)

## d) ----

# Teste Log-Rank

survdiff(
  Surv(tempo, status) ~ tratam,
  data = DadosAIDS
)

logrank_test(
  Surv(tempo, status) ~ factor(tratam),
  data = DadosAIDS,
  type = "Peto-Peto"
)

# Interpretação:
# Tanto o Log-Rank (p = 7e-15) quanto o Peto-Peto (p = 2,42e-10)
# rejeitam H0, indicando diferença significativa entre as curvas de
# sobrevivência dos tratamentos.

## e)  ----

# Bonferroni

pairwise_survdiff(
  Surv(tempo, status) ~ tratam,
  data = DadosAIDS,
  p.adjust.method = "bonferroni",
  rho = 0
)

pairwise_survdiff(
  Surv(tempo, status) ~ tratam,
  data = DadosAIDS,
  p.adjust.method = "holm",
  rho = 0
)

# Interpretação:
# Com Bonferroni e Holm, houve diferença significativa entre os tratamentos
# 0×1, 0×2, 0×3 e 1×2 (p < 0,05).
# Não houve diferença significativa entre 1×3 e 2×3.