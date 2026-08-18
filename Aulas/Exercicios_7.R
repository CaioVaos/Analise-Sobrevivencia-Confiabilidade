# Setup ----
# install.packages("ggsurvfit")
# install.packages("survRM2")

require(dplyr)
require(ggplot2)
require(survival)    
require(ggsurvfit)
require(survRM2)

# 1 ----

tempos <- c(0.19,0.78,0.96,1.31,2.78,3.16,4.67,4.85,6.50,7.35,8.27,12.07,32.52,33.91,36.71,rep(36.71,10))
cens <- c(rep(1,15),rep(0,10)) 
Dados2 <- data.frame(tempos,cens)

ekm_Dados2<- survfit2(data=Dados2,
               formula = Surv(tempos,cens)~1,
               conf.type = "log-log", 
               conf.int=0.95)
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

# 1 ----
Dados3 <- readr::read_csv("Aulas/Data/DadosAIDS.csv")

## Sem Estrato ----

ekm_Dados3_sem_estrat<- survfit2(data=Dados3,
               formula = Surv(tempo,status)~1,
               conf.type = "log-log", 
               conf.int=0.95)
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

summary(ekm_Dados3_sem_estrat,times = c(1000,1750))

# > summary(ekm_Dados3_sem_estrat,times = c(1000,1750))
# Call: survfit(formula = Surv(tempo, status) ~ 1, data = Dados3, conf.type = "log-log", 
#     conf.int = 0.95)

#  time n.risk n.event survival std.err lower 95% CI upper 95% CI
#  1000     75      75    0.561  0.0389        0.481        0.633
#  1750     29      13    0.443  0.0428        0.358        0.524

# P(T>1000) = 0.561
# P(T<1759) = 1 - P(T<1759) = 1 - 0.443 = 0.557

## Com Estrato ----

ekm_Dados3_com_estrat<- survfit2(data=Dados3,
               formula = Surv(tempo,status)~sexo,
               conf.type = "log-log", 
               conf.int=0.95)
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

summary(ekm_Dados3_com_estrat,times = c(1000,1750))

# > summary(ekm_Dados3_com_estrat,times = c(1000,1750))
# Call: survfit(formula = Surv(tempo, status) ~ sexo, data = Dados3, 
#     conf.type = "log-log", conf.int = 0.95)

#                 sexo=F 
#  time n.risk n.event survival std.err lower 95% CI upper 95% CI
#  1000     21      14    0.674  0.0737        0.507        0.795
#  1750      7       1    0.626  0.0826        0.443        0.763

#                 sexo=M 
#  time n.risk n.event survival std.err lower 95% CI upper 95% CI
#  1000     54      61    0.522  0.0452        0.430        0.606
#  1750     22      12    0.390  0.0477        0.297        0.482

# Feminino
# P(T>1000) = 0.674
# P(T<1759) = 1 - P(T<1759) = 1 - 0.626 = 0.374

# Masculino
# P(T>1000) = 0.522
# P(T<1759) = 1 - P(T<1759) = 1 - 0.390 = 0.61