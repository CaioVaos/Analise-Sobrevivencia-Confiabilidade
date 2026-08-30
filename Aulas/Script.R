require(dplyr)
require(ggplot2)
require(survival)    
require(ggsurvfit)
require(survRM2)
require(coin)
require(survminer)
require(multcompView)

ekm<- survfit2(data=Dados,
               formula = Surv(tempos,cens)~grupos,
               conf.type = "plain" ou "log" ou "log-log", 
               conf.int=0.95)
summary(ekm)

summary(ekm,times = c(???))

ggsurvfit(ekm, linewidth = 1.2) +
  add_confidence_interval() +
  add_censor_mark(size = 2, alpha = 1)+
  add_risktable() +
  scale_ggsurvfit() +
  labs(
    x = "Tempo (semanas)",
    y = "Sobrevivência estimada"
  ) +
  theme_classic(base_size = 16)

quantile(ekm,p)

fit_na <- survfit2(
  Surv(tempos, cens) ~ grupos,
  data = Dados,
  conf.type = "plain" ou "log" ou "log-log", 
  type = "fh"  
)

tabela_nelson_aalen <- function(fit, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level)/2)  # valor crítico da normal
  resultados <- list()
  grupos_nomes <- names(fit$strata)
  
  for(i in seq_along(grupos_nomes)) {
    if(i == 1) {
      idx <- 1:fit$strata[i]
    } else {
      idx <- (sum(fit$strata[1:(i-1)]) + 1):sum(fit$strata[1:i])
    }
    
    tabela_grupo <- data.frame(
      tj = fit$time[idx],
      dj = fit$n.event[idx],
      nj = fit$n.risk[idx],
      Lambda = fit$cumhaz[idx],
      Std.err = fit$std.chaz[idx]
    ) %>%
      dplyr::filter(dj > 0) %>%
      dplyr::mutate(
        IC_lower = Lambda / exp(z * Std.err / Lambda),
        IC_upper = Lambda * exp(z * Std.err / Lambda)
      ) %>%
      dplyr::select(tj, dj, nj, Lambda, Std.err, IC_lower, IC_upper) %>%
      dplyr::mutate(
        dplyr::across(where(is.numeric), ~ round(., 6))
      )
    
    resultados[[grupos_nomes[i]]] <- tabela_grupo
  }
  
  return(resultados)
}

tabelas_por_grupo <- tabela_nelson_aalen(fit_na, conf_level = 0.95) 

for(grupo in names(tabelas_por_grupo)) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("GRUPO:", grupo, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  print(tabelas_por_grupo[[grupo]])
}

ggsurvfit(fit_na, 
          linewidth = 1.2,
          type = "cumhaz") +
  add_censor_mark(size = 2, alpha = 1)+
  add_confidence_interval() +
  labs(
    x = "Tempo (semanas)",
    y = "Função de risco acumulada estimada"
  ) +
  theme_classic(base_size = 16)

Tempo = Dados$tempos
Cens = Dados$cens
Grupo <- as.numeric(as.factor(Dados$grupos)) - 1
resultado = rmst2(Tempo, Cens, Grupo,tau=?)

survdiff(Surv(tempos,cens)~grupos, data = Dados)
logrank_test(Surv(tempos,cens)~factor(grupos),type="Peto-Peto",data=Dados)

res = pairwise_survdiff(Surv(tempos,cens)~grupos,
                        data=Dados,
                        p.adjust.method = "bonferroni" ou "holm",
                        rho=0 ou 1)
p = res$p.value
