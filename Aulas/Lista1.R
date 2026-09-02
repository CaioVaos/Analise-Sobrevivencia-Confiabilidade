# Lista 1 - Analise de Sobrevivencia e Confiabilidade

# Setup ----

# install.packages("survival")
# install.packages("ggsurvfit")
# install.packages("survRM2")
# install.packages("coin")
# install.packages("survminer")
# install.packages("multcompView")

library(dplyr)
library(ggplot2)
library(survival)
library(ggsurvfit)
library(survRM2)
library(coin)
library(survminer)
library(multcompView)

# Q1 ----
#
# Análise de Sobrevivência: estuda o tempo até a ocorrência de um evento de 
# interesse, normalmente em indivíduos/pacientes.
#
# Confiabilidade: estuda o tempo até a ocorrência de uma falha em componentes, 
# sistemas ou equipamentos.
#
# Semelhanças:
# - as duas estudam um tempo até um evento;
# - ambas podem possuir observações censuradas;
# - utilizam funções como S(t), f(t), lambda(t) e Lambda(t);
# - utilizam métodos como Kaplan-Meier e Nelson-Aalen.
#
# Diferença principal:
# - na Sobrevivência, o evento geralmente está relacionado a indivíduos, 
# por exemplo morte ou recaída;
# - na Confiabilidade, o evento geralmente é uma falha de um equipamento, 
# componente ou sistema.

# Q2 ----
#
# As técnicas usuais não lidam adequadamente com censuras.
# Quando existe censura, não observamos o verdadeiro tempo do evento 
# para todos os indivíduos.
#
# Por exemplo:
# Se um paciente ainda estiver vivo no fim do estudo, sabemos apenas que seu tempo 
# de vida é maior que o tempo observado.
# Ignorar isso e aplicar técnicas usuais trataria esse tempo censurado como se 
# fosse um tempo completo, produzindo estimativas viesadas.

# Q3 ----
#
# Censura ocorre quando o evento de interesse não é observado para determinado 
# indivíduo, mas sabemos alguma informação parcial sobre seu tempo até o evento.
#
# O dado censurado não pode simplesmente ser removido porque ele ainda fornece 
# informação sobre a sobrevivência.
# Um indivíduo censurado em t contribuiu para o estudo enquanto esteve sob risco, 
# até o instante t.

# Q4 ----
#
# Exemplos em Sobrevivência:
# - paciente termina o estudo vivo;
# - paciente abandona o estudo antes de apresentar o evento;
# - paciente é perdido no acompanhamento.
#
# Exemplos em Confiabilidade:
# - equipamento ainda funcionando no fim do teste;
# - componente retirado do teste antes de falhar;
# - teste interrompido depois que um número pré-estabelecido de componentes falhou.

# Q5 ----

## a) ----
#
# Estudos Observacionais:
# O pesquisador só observa os indivíduos e registra as informações, 
# sem alterar a evolução do estudo.

## b) ----
#
# Estudos Experimentais:
# O pesquisador intervém na evolução do extudo, experimentando 
# tratamento no indivíduo.

## c) ----
#
# Prospectivo:
# os indivíduos são acompanhados para o futuro, a partir do início do estudo.
#
# Retrospectivo:
# o pesquisador utiliza informações já ocorridas e registradas anteriormente.

## d) ----
#
# Descritivo:
# Apenas descreve características de uma população, 
# sem necessariamente investigar relações causais.
#
# Caso-controle:
# Parte-se dos indivíduos com e sem o desfecho e investiga-se 
# retrospectivamente a exposição.
#
# Coorte:
# Acompanha-se um grupo de indivíduos ao longo do tempo, 
# observando a ocorrência do evento.
#
# Clínico aleatorizado:
# Indivíduos são aleatoriamente distribuídos entre tratamentos e 
# posteriormente acompanhados para comparação dos desfechos.

# Q6 ----

## a) ----
#
# T = tempo, em semanas, desde a exposição ao material cancerígeno 
# até o desenvolvimento do tumor de determinado tamanho.

## b) ----
#
# Rato A: evento em 10 semanas
# Rato B: evento em 15 semanas
# Rato C: evento em 25 semanas

## c) ----
#
# Rato D: censurado em 20 semanas.
#
# Rato E: censurado em 30 semanas.
# Rato F: censurado em 30 semanas.

# Q7 ----
#
# I. Tempo inicial: instante a partir do qual o acompanhamento começa.
#
# II. Escala de medida: unidade em que o tempo é medido.
#
# III. Evento de interesse: acontecimento cujo tempo até sua ocorrência 
# está sendo estudado.

# Q8 ----
#
# Censura à Direita: sabemos que o evento não ocorreu até determinado instante.
# Exemplo: paciente termina o estudo sem morrer.
#
# Censura à Esquerda: sabemos que o evento já ocorreu antes de determinado instante
# Exemplo: no primeiro acompanhamento o paciente já apresentava o sintoma.
#
# Censura Intervalar: sabemos apenas que o evento ocorreu entre dois instantes.
# Exemplo: num acompanhamento períodico, o 1º é negativo e o 2º é positivo.

# Q9 ----
#
# Tipo I: todos são acompanhados até um tempo fixo de encerramento.
# Quem não apresentou o evento até esse instante é censurado.
#
# Tipo II: o estudo termina quando ocorre um número fixado de eventos.
# Os indivíduos que ainda não falharam nesse instante são censurados.
#
# Aleatória: o tempo de censura varia de indivíduo para indivíduo e é
# determinado por algum mecanismo aleatório, como perda de acompanhamento.

# Q10 ----
#
# Não Informativa: o mecanismo de censura é independente do tempo de ocorrência 
# do evento, condicionado às informações observadas.
# Exemplo: estudo termina em uma data pré-determinada para todos.
#
# Censura informativa: a probabilidade de censura está relacionada ao tempo de 
# ocorrência do evento.
# Exemplo: pacientes com maior risco de morte abandonam o estudo com maior 
# frequência por causa da própria condição clínica.

# Q11 ----

## a) ----
#
# Tempo inicial: momento da infecção pela malária;
# Escala: dias;
# Evento de interesse: morte do camundongo.

## b) ----
#
# Tempo inicial: nascimento da criança;
# Escala: meses;
# Evento de interesse: desmame completo.

# Q12 ----

## a) ----
#
# Censura à direita aleatória.

## b) ----
#
# Censura à esquerda.

## c) ----
#
# Censura intervalar.

## d) ----
#
# Censura à direita tipo I

# Q13 ----

## a) ----
#
# Truncamento ocorre quando a própria inclusão do indivíduo na amostra depende 
# do seu tempo de ocorrência do evento. 
# Assim, certos indivíduos não aparecem na base.

## b) ----
#
# Truncamento à esquerda:
# indivíduos somente entram na amostra depois de determinado instante/critério 
# inicial, fazendo com que parte do histórico anterior não seja observada.
#
# Truncamento à direita:
# indivíduos entram na amostra somente se o evento de interesse
# tiver ocorrido antes de determinado limite/critério final.

## c) ----
#
# Porque os indivíduos que não apresentaram o evento antes do limite simplesmente 
# não entram na amostra. Eles não aparecem como censurados; 
# são excluídos pelo próprio mecanismo de seleção.

# Q14 ----

## a) ----
# f(t): função densidade de probabilidade de T.
# Indica a concentração de probabilidade do tempo do evento.

## b) ----
# F(t): função distribuição acumulada de T
# F(t) = P(T <= t).
# É a probabilidade de o evento ocorrer até t.

## c) ----
# S(t): função de sobrevivência de T
# S(t) = P(T > t).
# É a probabilidade de o indivíduo permanecer sem o evento até t.

## d) ----
# lambda(t): função de risco
# Representa a taxa instantânea de ocorrência do evento no tempo t,
# condicionada a o indivíduo ter sobrevivido até t.

## e) ----
# Lambda(t): função de risco acumulado
# Lambda(t) = integral de lambda(u) du de 0 até t.
# Acumula o risco ao longo do tempo.

## f) ----
# tm: tempo médio até a ocorrência do evento:
# tm = E(T) = integral de S(t) dt de 0 até infinito, quando a integral é finita.

## g) ----
# vmr(t): tempo médio residual:
# vmr(t) = E(T-t | T>t).
# É o tempo médio que ainda resta após t entre os indivíduos que 
# chegaram vivos/sem falha até t.

# Q15 -----

# Q16 ----

# Q17 ----

# Q18 ----

# Q19 ----

# Q20 ----

## a) ----

S_weibull <- function(t, rho, kappa) {
  exp(-rho * t^kappa)
}

S_weibull(30, 0.0001, 2) #0.9139312
S_weibull(45, 0.0001, 2) #0.8166865

# 91.39% não apresentaram tumor até 30 dias.
# 81.67% não apresentaram tumor até 45 dias.

## b) ----
# Media

# Mediana
S_weibull(83.25, 0.0001, 2) #83.25 dias

## c) ----
# lambda = -deriv(log((t)))

lambda_weibull <- function(t, rho, kappa) {
  rho * kappa * t^(kappa - 1)
}

lambda_weibull(30, 0.0001, 2) #0.006
lambda_weibull(45, 0.0001, 2) #0.009
lambda_weibull(60, 0.0001, 2) #0.012

# Interpretação:
# A taxa instantânea de ocorrência do tumor aumenta com o tempo.
# Assim, entre os indivíduos que ainda não desenvolveram o tumor, a taxa 
# instantânea de desenvolvimento do tumor é maior em tempos mais avançados.

# Q21 ----

## a) ----
# Aproximadamente 10%

## b) ----
# Aproximadamente 5 dias.

## c) ----
# Aproximadamente 4 dias.

## Q22 ----
#
# O estimador de Kaplan-Meier constrói a função de sobrevivência de forma 
# acumulativa.
# Calcula-se, em cada tempo de ocorrência do evento, a proporção 
# de indivíduos que permanecem sem o evento entre aqueles que ainda estavam sob 
# risco.
# As censuras não são consideradas como eventos, elas apenas indicam que o 
# indivíduo deixa de ser acompanhado a partir daquele momento.
# Assim, o estimador utiliza as informações disponíveis em cada instante para 
# atualizar progressivamente a estimativa de sobrevivência.

# Q23 ----
#
# Os intervalos de confiança precisam respeitar a natureza da função de 
# sobrevivência, que deve estar entre 0 e 1.
# A transformação log-log é conveniente pois coloca em temos do intervalo 
# (0,1) sem resvalar em valores impossiveis (>1).

# Q24 ----
#
# O estimador de Nelson-Aalen estima diretamente o risco acumulado,
# somando as taxas de eventos ao longo do tempo.
#
# A transformação log é apropriada para intervalos de confiança, pois mantém os 
# limites positivos, já que o risco acumulado não pode ser negativo.

# Q25 ----
#
# Percentis são obtidos identificando, na curva de sobrevivência, o tempo 
# correspondente à probabilidade de sobrevivência desejada.
# Por exemplo, a mediana, é o tempo em que S(t)=0,5, indicando que metade da 
# população já apresentou o evento.

# Q26 ----
#
# 1. estima-se a curva de sobrevivência;
# 2. constrói-se o IC da curva;
# 3. transforma-se a posição horizontal correspondente ao percentil desejado 
# em limites de tempo;
# 4. os pontos de interseção com os limites inferior e superior da curva geram o
# intervalo de confiança para o percentil.

# Q27 ----
#
# O RMST é o tempo médio de sobrevivência até um determinado limite de tempo.
# Ele corresponde à área sob a curva de sobrevivência até esse limite.
#
# Ele pode ser preferido à média tradicional porque não exige que todos os 
# indivíduos apresentem o evento e continua sendo interpretável mesmo quando 
# a distribuição tem forte censura ou média não estimável.

# Q28 (cap.7 - 2)----

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


## i) ----

Tempo <- Dados3$tempo
Cens <- Dados3$status
Grupo <- as.numeric(as.factor(Dados3$sexo)) - 1

rmst_resultado <- rmst2(Tempo, Cens, Grupo, tau = 3000)

rmst_resultado

# 29.
#
# O teste Log-Rank compara as curvas de sobrevivência dos grupos.
#
# H0:
# S1(t) = S2(t) = ... = Sk(t), para todo t.
#
# H1:
# pelo menos uma das funções difere em algum tempo.

# Q30 ----
#
# Na Log-Rank, os tempos de ocorrência recebem o mesmo peso 
# na construção da estatística.
# Em alguns problemas, as curvas podem diferir 
# principalmente no início do acompanhamento.
# Nesse caso, testes ponderados podem atribuir maior peso 
# aos tempos iniciais, quando há mais indivíduos sob risco.
#
# - Gehan-Breslow: peso pelo número de indivíduos em risco
# - Tarone-Ware: peso pela raiz quadrada do número em risco
# - Peto-Peto: peso pela curva de sobrevivência

# Q31 ----
#
# Ao comparar mais de dois grupos, vários testes pareados 
# aumentam o número de hipóteses testadas, inflacionando o 
# erro do tipo I. Por isso, são necessárias correções para 
# múltiplas comparações.

# Q32 ----
#
# Bonferroni: controla o erro dividindo alfa pelo 
# número de comparações,mas é conservador.
#
# Holm-Bonferroni: controla o erro dividindo pelo 
# ordem o p-valor(1ª, 2ª,...), sendo menos conservador.

# Q33 (cap.9 - 3) ----

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