# =============================================================================
# PARTE II · ETAPA I — MODELOS PARAMÉTRICOS SEM COVARIÁVEIS (Cap. 13 do livro)
# -----------------------------------------------------------------------------
# Execute PRIMEIRO este script, a partir da raiz do projeto (pasta que contém
# "Trabalho/"). Ele centraliza o pré-processamento e gera:
#   Trabalho/data/df_parte2.RDS            -> base analítica (usada na Etapa II)
#   Trabalho/parte2/data/base_analitica.RDS-> números p/ slide "Base analítica"
#   Trabalho/parte2/data/tabelas_parte1.RDS-> tabelas (brutas e formatadas)
#   Trabalho/parte2/plots/*                -> .RDS e .PNG com os nomes do index.qmd
# Metodologia (Cap. 13): 5 distribuições; escolha conjunta por KM x modelo,
# linearização, AIC/BIC e TRV contra a Gama Generalizada.
# (Não há seleção de variáveis nesta etapa — script inalterado.)
# =============================================================================

# 1. Setup --------------------------------------------------------------------
library(dplyr)
library(survival)
library(flexsurv)
library(coin)
library(ggplot2)
library(qqplotr)
library(broom)

dir_dados <- "Trabalho/data"
dir_saida <- "Trabalho/parte2/data"
dir_plots <- "Trabalho/parte2/plots"
for (d in c(dir_saida, dir_plots)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

dists <- c(
  "Exponencial"       = "exponential",
  "Weibull"           = "weibull",
  "Log-normal"        = "lnorm",
  "Gama"              = "gamma",
  "Gama Generalizada" = "gengamma"
)

# 2. Decisões de base (única vez; a Etapa II lê o resultado salvo) -------------
df <- as.data.frame(readRDS(file.path(dir_dados, "df_final.RDS")))

n_zero  <- sum(df$tempo <= 0)
ev_zero <- sum(df$cens[df$tempo <= 0] == 1)
cat("Pacientes com tempo = 0 removidos:", n_zero, "(", ev_zero, "óbitos )\n")
df <- df %>% filter(tempo > 0)

ausentes <- "categoria"
if (ausentes == "completos") {
  df <- df %>%
    filter(ESCOLARI != "Sem informação",
           !ECGRUP %in% c("Sem informação", "Não se aplica"))
}

df <- df %>%
  mutate(
    ESCOLARI = droplevels(ESCOLARI),
    ECGRUP   = droplevels(ECGRUP),
    IDADE10  = IDADE / 10
  )
if ("Analfabeto" %in% levels(df$ESCOLARI)) df$ESCOLARI <- relevel(df$ESCOLARI, ref = "Analfabeto")
if ("Estágio I"  %in% levels(df$ECGRUP))   df$ECGRUP   <- relevel(df$ECGRUP,   ref = "Estágio I")

saveRDS(df, file.path(dir_dados, "df_parte2.RDS"))

# 3. BLOCO COMUM (idêntico nos dois scripts) -----------------------------------
cor_navy    <- "#10241d"
cor_verde   <- "#3f6b52"
cor_dourado <- "#b8862a"
cor_cinza   <- "#6b7d73"
cor_media   <- "#c0392b"

pal_modelos <- c(
  "Exponencial"       = cor_cinza,
  "Weibull"           = cor_dourado,
  "Log-normal"        = cor_verde,
  "Gama"              = "#7a5a12",
  "Gama Generalizada" = cor_media
)

tema_km <- theme_minimal(base_size = 13) +
  theme(
    plot.title.position  = "plot",
    plot.title           = element_text(color = cor_navy, face = "bold",
                                        size = 17, margin = margin(b = 4)),
    plot.subtitle        = element_text(color = cor_cinza, size = 10.5,
                                        margin = margin(b = 10)),
    axis.title           = element_text(color = cor_navy, size = 11),
    axis.title.x         = element_text(margin = margin(t = 8)),
    axis.title.y         = element_text(margin = margin(r = 8)),
    axis.text            = element_text(color = cor_cinza, size = 10),
    axis.line.x          = element_line(color = cor_cinza, linewidth = 0.4),
    axis.ticks.x         = element_line(color = cor_cinza, linewidth = 0.4),
    strip.text           = element_text(color = cor_navy, face = "bold", size = 10.5),
    legend.position      = "top",
    legend.justification = "left",
    legend.title         = element_text(color = cor_navy, face = "bold", size = 10.5),
    legend.text          = element_text(color = cor_navy, size = 10),
    legend.key.width     = unit(1.3, "cm"),
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(),
    panel.grid.major.y   = element_line(color = "grey92", linewidth = 0.5),
    plot.background      = element_rect(fill = "white", color = NA),
    plot.margin          = margin(14, 18, 10, 14)
  )

n_total  <- format(nrow(df), big.mark = ".", decimal.mark = ",")
n_obitos <- format(sum(df$cens == 1), big.mark = ".", decimal.mark = ",")
subtitulo <- paste0("Câncer de fígado (C22)  ·  FOSP 2014–2019  ·  n = ", n_total,
                    "  ·  ", n_obitos, " óbitos")

formatar_p <- function(p) {
  ifelse(is.na(p), "—",
         ifelse(p < 0.001, "< 0,001",
                formatC(p, format = "f", digits = 3, decimal.mark = ",")))
}
fmt_num <- function(x, d = 1) formatC(x, format = "f", digits = d,
                                      big.mark = ".", decimal.mark = ",")

fml <- function(vars) {
  as.formula(paste("Surv(tempo, cens) ~",
                   if (length(vars) == 0) "1" else paste(vars, collapse = " + ")))
}
ajusta <- function(vars, dist = "gengamma", data = df) {
  flexsurvreg(fml(vars), data = data, dist = dist)
}
trv <- function(m1, m0) {
  est <- 2 * (m1$loglik - m0$loglik)
  gl  <- m1$npars - m0$npars
  data.frame(TRV = est, gl = gl, P_valor = pchisq(est, gl, lower.tail = FALSE))
}
tabela_ic <- function(ajustes) {
  data.frame(
    Modelo = names(ajustes),
    logLik = sapply(ajustes, function(a) a$loglik),
    AIC    = sapply(ajustes, function(a) a$AIC),
    BIC    = sapply(ajustes, function(a) a$BIC),
    row.names = NULL
  )
}
tabela_trv_gg <- function(ajustes) {
  gg <- ajustes[["Gama Generalizada"]]
  do.call(rbind, lapply(setdiff(names(ajustes), "Gama Generalizada"), function(nm) {
    cbind(Comparacao = paste(nm, "x GG"), Modelo = nm, trv(gg, ajustes[[nm]]))
  }))
}
forma_risco <- function(nome, fit) {
  r <- fit$res
  switch(nome,
    "Exponencial" = "risco constante ao longo do tempo",
    "Weibull" = sprintf("risco %s ao longo do tempo (forma = %s)",
      if (r["shape", "est"] > 1) "crescente" else
        if (r["shape", "est"] < 1) "decrescente" else "constante",
      fmt_num(r["shape", "est"], 2)),
    "Log-normal" = "risco unimodal: cresce até um pico e depois decresce",
    "Gama" = sprintf("risco %s ao longo do tempo (forma = %s)",
      if (r["shape", "est"] > 1) "crescente" else "decrescente",
      fmt_num(r["shape", "est"], 2)),
    "Gama Generalizada" = sprintf(
      "risco flexível (Q = %s); engloba as demais como casos particulares",
      fmt_num(r["Q", "est"], 2)),
    "forma de risco não mapeada"
  )
}
salva_plot <- function(p, nome, w = 10, h = 6.5) {
  saveRDS(p, file.path(dir_plots, paste0(nome, ".RDS")))
  ggsave(file.path(dir_plots, paste0(nome, ".png")), p,
         width = w, height = h, dpi = 200, bg = "white")
}

# 4. Etapa I — Modelos sem covariáveis -----------------------------------------
ajustes1 <- lapply(dists, function(d) ajusta(character(0), dist = d))
ajustes1

tab_ic1  <- tabela_ic(ajustes1)
tab_trv1 <- tabela_trv_gg(ajustes1)
tab_trv1$Decisao <- ifelse(tab_trv1$P_valor > 0.05, "compatível com GG", "rejeitado")
tab_ic1[order(tab_ic1$AIC), ]
tab_ic1[order(tab_ic1$BIC), ]
tab_trv1

ekm  <- survfit(Surv(tempo, cens) ~ 1, data = df)
t_km <- ekm$time
s_km <- ekm$surv

df_km <- do.call(rbind, lapply(names(ajustes1), function(nm) {
  data.frame(
    Modelo   = nm,
    S_KM     = s_km,
    S_Modelo = summary(ajustes1[[nm]], t = t_km, type = "survival",
                       ci = FALSE)[[1]]$est
  )
}))
df_km$Modelo <- factor(df_km$Modelo, levels = names(ajustes1))

p_km_modelos <- ggplot(df_km, aes(x = S_KM, y = S_Modelo)) +
  geom_abline(slope = 1, intercept = 0, color = cor_media, linetype = "dashed") +
  geom_point(color = cor_verde, alpha = 0.55, size = 1.3) +
  facet_wrap(~ Modelo, nrow = 2) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Kaplan-Meier × Modelos Paramétricos", subtitle = subtitulo,
       x = "S(t) — Kaplan-Meier", y = "S(t) — modelo paramétrico") +
  tema_km
salva_plot(p_km_modelos, "km_vs_modelos", w = 11, h = 6.5)

ok <- s_km > 0 & s_km < 1
df_lin <- rbind(
  data.frame(Painel = "Exponencial: −log S(t) × t",
             x = t_km[ok],      y = -log(s_km[ok])),
  data.frame(Painel = "Weibull: log(−log S(t)) × log t",
             x = log(t_km[ok]), y = log(-log(s_km[ok]))),
  data.frame(Painel = "Log-normal: Φ⁻¹(S(t)) × log t",
             x = log(t_km[ok]), y = qnorm(s_km[ok]))
)
df_lin$Painel <- factor(df_lin$Painel, levels = unique(df_lin$Painel))

p_linearizacao <- ggplot(df_lin, aes(x = x, y = y)) +
  geom_point(color = cor_verde, alpha = 0.5, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = cor_media,
              linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ Painel, nrow = 1, scales = "free") +
  labs(title = "Método da Linearização",
       subtitle = "Pontos alinhados à reta indicam adequação do modelo",
       x = "Transformação do tempo", y = "Transformação de S(t)") +
  tema_km
salva_plot(p_linearizacao, "linearizacao", w = 12, h = 4.5)

t_grid <- seq(0.05, max(df$tempo), length.out = 200)
df_curvas <- do.call(rbind, lapply(names(ajustes1), function(nm) {
  s <- summary(ajustes1[[nm]], t = t_grid, type = "survival", ci = FALSE)[[1]]
  data.frame(Modelo = nm, tempo = s$time, S = s$est)
}))
df_curvas$Modelo <- factor(df_curvas$Modelo, levels = names(ajustes1))
df_ekm <- data.frame(tempo = c(0, ekm$time), S = c(1, ekm$surv))

p_curvas_modelos <- ggplot() +
  geom_step(data = df_ekm, aes(x = tempo, y = S), color = "black", linewidth = 0.9) +
  geom_line(data = df_curvas, aes(x = tempo, y = S, color = Modelo), linewidth = 0.9) +
  scale_color_manual(values = pal_modelos) +
  scale_x_continuous(breaks = seq(0, 144, 24)) +
  labs(title = "Sobrevida: Kaplan-Meier (preto) e modelos paramétricos",
       subtitle = subtitulo, x = "Tempo desde o diagnóstico (meses)",
       y = "S(t)", color = "Modelo") +
  tema_km
salva_plot(p_curvas_modelos, "curvas_modelos", w = 10, h = 6.5)

compat1   <- c("Gama Generalizada", tab_trv1$Modelo[tab_trv1$P_valor > 0.05])
sugestao1 <- tab_ic1$Modelo[tab_ic1$Modelo %in% compat1][
  which.min(tab_ic1$AIC[tab_ic1$Modelo %in% compat1])]
sugestao1

decisoes1 <- list(
  sugerida        = sugestao1,
  compativeis_trv = compat1,
  menor_aic       = tab_ic1$Modelo[which.min(tab_ic1$AIC)],
  menor_bic       = tab_ic1$Modelo[which.min(tab_ic1$BIC)],
  forma_risco     = forma_risco(sugestao1, ajustes1[[sugestao1]])
)

vars_cat <- intersect(c("SEXO", "ESCOLARI", "ECGRUP", "CATEATEND",
                        "CIRURGIA", "RADIO", "QUIMIO"), names(df))
base_analitica <- list(
  n = nrow(df), obitos = sum(df$cens == 1), censuras = sum(df$cens == 0),
  pct_obitos = 100 * mean(df$cens == 1),
  mediana_km = unname(summary(ekm)$table["median"]),
  removidos_tempo0 = n_zero, obitos_tempo0 = ev_zero,
  ausentes = ausentes,
  niveis = lapply(df[, vars_cat, drop = FALSE], table),
  texto = sprintf(
    paste0("Base analítica: n = %s (%s óbitos por câncer; %s censuras). ",
           "Removidos %s pacientes com tempo = 0 (%s óbitos). Ausentes: %s."),
    fmt_num(nrow(df), 0), fmt_num(sum(df$cens == 1), 0), fmt_num(sum(df$cens == 0), 0),
    fmt_num(n_zero, 0), fmt_num(ev_zero, 0),
    if (ausentes == "categoria")
      "'Sem informação'/'Não se aplica' mantidos como categorias próprias (mesma amostra em todos os modelos)"
    else "análise de casos completos"
  )
)
saveRDS(base_analitica, file.path(dir_saida, "base_analitica.RDS"))

tab_ic1_fmt <- data.frame(Modelo = tab_ic1$Modelo,
                          logLik = fmt_num(tab_ic1$logLik, 1),
                          AIC = fmt_num(tab_ic1$AIC, 1),
                          BIC = fmt_num(tab_ic1$BIC, 1))
tab_trv1_fmt <- data.frame(Comparacao = tab_trv1$Comparacao,
                           TRV = fmt_num(tab_trv1$TRV, 2), gl = tab_trv1$gl,
                           p_valor = formatar_p(tab_trv1$P_valor),
                           Decisao = tab_trv1$Decisao)

tabelas_parte1 <- list(ic = tab_ic1, trv = tab_trv1,
                       ic_fmt = tab_ic1_fmt, trv_fmt = tab_trv1_fmt,
                       decisoes = decisoes1, base = base_analitica)
saveRDS(tabelas_parte1, file.path(dir_saida, "tabelas_parte1.RDS"))

cat("\nEtapa I concluída.",
    "\nDistribuição sugerida:", sugestao1, "|", decisoes1$forma_risco,
    "\nMenor AIC:", decisoes1$menor_aic, "| Menor BIC:", decisoes1$menor_bic,
    "\nArquivos para os slides em:", dir_saida, "e", dir_plots, "\n")