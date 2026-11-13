# =============================================================================
# PARTE II · ETAPA II — MODELOS PARAMÉTRICOS COM COVARIÁVEIS (Cap. 16 do livro)
# -----------------------------------------------------------------------------
# REQUER: rode antes "analise_sem_covariavel.R" (gera df_parte2.RDS,
# base_analitica.RDS e tabelas_parte1.RDS). Execute da raiz do projeto.
# Metodologia (Cap. 16): filtro p < 0,25 (Log-Rank/Peto-Peto); seleção em 5
# passos por TRV a 10% na Gama Generalizada; troca da distribuição ao final;
# diagnóstico (Cox-Snell e quantílicos aleatorizados); RTM = exp(±coef).
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

# 2. Base analítica (gerada na Etapa I; não reprocessar aqui) ------------------
arq_df <- file.path(dir_dados, "df_parte2.RDS")
if (!file.exists(arq_df))
  stop("Rode primeiro 'analise_sem_covariavel.R' (gera df_parte2.RDS e base_analitica.RDS).")
df             <- readRDS(arq_df)
base_analitica <- readRDS(file.path(dir_saida, "base_analitica.RDS"))
ausentes       <- base_analitica$ausentes

# 3. BLOCO COMUM (idêntico ao do script da Etapa I) -----------------------------
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

# Funções de seleção (Cap. 16) -------------------------------------------------
# Passos 1, 3, 4 e 5: acrescenta, uma a uma, cada variável de `novas` ao `base`
trv_adicionar <- function(base, novas, dist = "gengamma", alfa = 0.10) {
  if (length(novas) == 0) return(NULL)
  m0  <- ajusta(base, dist)
  out <- do.call(rbind, lapply(novas, function(v) {
    cbind(Variavel = v, trv(ajusta(c(base, v), dist), m0))
  }))
  out$Decisao <- ifelse(out$P_valor < alfa, "entra", "não entra")
  out
}
# Passo 2: retira, uma a uma, cada variável do modelo múltiplo
trv_retirar <- function(vars, dist = "gengamma", alfa = 0.10) {
  m1  <- ajusta(vars, dist)
  out <- do.call(rbind, lapply(vars, function(v) {
    cbind(Variavel = v, trv(m1, ajusta(setdiff(vars, v), dist)))
  }))
  out$Decisao <- ifelse(out$P_valor < alfa, "fica", "sai")
  out
}

# 5. Etapa II — Modelos com covariáveis ----------------------------------------
alfa_sel <- 0.10   # nível do TRV na seleção (Cap. 16)

# 5.1 Filtro inicial (Log-Rank e Peto-Peto, p < 0,25) ----
qualitativas <- c("SEXO", "ESCOLARI", "ECGRUP", "CATEATEND",
                  "CIRURGIA", "RADIO", "QUIMIO")
filtro <- do.call(rbind, lapply(qualitativas, function(v) {
  f_lr <- as.formula(paste("Surv(tempo, cens) ~", v))
  f_pt <- as.formula(paste0("Surv(tempo, cens) ~ factor(", v, ")"))
  lr   <- survdiff(f_lr, data = df, rho = 0)
  pt   <- logrank_test(f_pt, data = df, type = "Peto-Peto")
  data.frame(
    Variavel       = v,
    Logrank_pvalor = pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE),
    Peto_pvalor    = as.numeric(pvalue(pt))
  )
}))
filtro$Passa_filtro <- pmin(filtro$Logrank_pvalor, filtro$Peto_pvalor) < 0.25
filtro

# 5.2 Covariáveis candidatas (justificativa clínica) ----
# Fora: RADIO (poucos tratados) e CATEATEND (Particular muito pequeno).
candidatas <- c("ECGRUP", "CIRURGIA", "QUIMIO", "ESCOLARI", "SEXO", "IDADE10")
justificativas <- c(
  ECGRUP    = "Extensão da doença no diagnóstico (principal fator prognóstico)",
  CIRURGIA  = "Único tratamento potencialmente curativo",
  QUIMIO    = "Tratamento sistêmico (curvas se cruzam na Parte I)",
  ESCOLARI  = "Marcador de acesso a diagnóstico e cuidado",
  SEXO      = "Diferença modesta na Parte I",
  RADIO     = "Poucos tratados (fora por justificativa clínica)",
  CATEATEND = "Categoria Particular muito pequena (fora)"
)

# 5.3 Passo 1: uma covariável por vez contra o modelo nulo (Gama Generalizada) -
passo1   <- trv_adicionar(character(0), candidatas, alfa = alfa_sel)
passo1
retidas1 <- passo1$Variavel[passo1$P_valor < alfa_sel]
stopifnot(length(retidas1) > 0)

# 5.4 Passo 2: modelo múltiplo, retirando uma variável por vez ----
passo2   <- trv_retirar(retidas1, alfa = alfa_sel)
passo2
retidas2 <- passo2$Variavel[passo2$P_valor < alfa_sel]
stopifnot(length(retidas2) > 0)

# 5.5 Passo 3: reentrada das descartadas no passo 2 ----
descartadas2 <- setdiff(retidas1, retidas2)
passo3   <- trv_adicionar(retidas2, descartadas2, alfa = alfa_sel)
passo3
retidas3 <- c(retidas2, passo3$Variavel[passo3$P_valor < alfa_sel])

# 5.6 Passo 4: reavaliação das excluídas no passo 1 ----
excluidas1 <- setdiff(candidatas, retidas1)
passo4   <- trv_adicionar(retidas3, excluidas1, alfa = alfa_sel)
passo4
retidas4 <- c(retidas3, passo4$Variavel[passo4$P_valor < alfa_sel])

# 5.7 Passo 5: interações duplas, uma a uma, com justificativa clínica ----
# ECGRUP × CIRURGIA: acesso à cirurgia depende do estágio
# ECGRUP × QUIMIO  : quimioterapia mais frequente nos estágios avançados
interacoes <- list(c("ECGRUP", "CIRURGIA"), c("ECGRUP", "QUIMIO"))
elegiveis  <- Filter(function(p) all(p %in% retidas4), interacoes)
passo5 <- NULL
if (length(elegiveis) > 0) {
  termos_inter <- vapply(elegiveis, function(p) paste(p, collapse = "*"), character(1))
  passo5 <- trv_adicionar(retidas4, termos_inter, alfa = alfa_sel)
}
passo5
inter_ret  <- if (!is.null(passo5)) passo5$Variavel[passo5$P_valor < alfa_sel] else character(0)
vars_final <- c(retidas4, inter_ret)
vars_final

# 5.8 Escolha da distribuição com as covariáveis selecionadas (Cap. 16) ----
ajustes2 <- lapply(dists, function(d) ajusta(vars_final, dist = d))
tab_ic2  <- tabela_ic(ajustes2)
tab_trv2 <- tabela_trv_gg(ajustes2)
tab_trv2$Decisao <- ifelse(tab_trv2$P_valor > 0.05, "compatível com GG", "rejeitado")
tab_ic2[order(tab_ic2$AIC), ]
tab_ic2[order(tab_ic2$BIC), ]
tab_trv2

# Menor AIC entre os não rejeitados pelo TRV (p > 0,05); confira também o BIC.
compat2 <- c("Gama Generalizada", tab_trv2$Modelo[tab_trv2$P_valor > 0.05])
dist_final_nome <- tab_ic2$Modelo[tab_ic2$Modelo %in% compat2][
  which.min(tab_ic2$AIC[tab_ic2$Modelo %in% compat2])]
# Override manual (justificar nos slides), se necessário:
# dist_final_nome <- "Log-normal"
dist_final_nome
dist_final   <- unname(dists[dist_final_nome])
modelo_final <- ajustes2[[dist_final_nome]]
modelo_final
saveRDS(modelo_final, file.path(dir_saida, "modelo_final.RDS"))

# 5.9 Reavaliação dos termos na distribuição escolhida ----
# Efeitos principais envolvidos em interações retidas não saem isoladamente
# (hierarquia do modelo). dAIC/dBIC > 0: retirar o termo piora o ajuste.
protegidas <- unique(unlist(strsplit(inter_ret, "*", fixed = TRUE)))
termos_av  <- c(setdiff(retidas4, protegidas), inter_ret)
tab_reavalia <- do.call(rbind, lapply(termos_av, function(v) {
  m0 <- ajusta(setdiff(vars_final, v), dist = dist_final)
  cbind(Variavel = v, trv(modelo_final, m0),
        dAIC = m0$AIC - modelo_final$AIC,
        dBIC = m0$BIC - modelo_final$BIC)
}))
tab_reavalia

# 5.10 Diagnóstico: resíduos de Cox-Snell (Cap. 16.6.1) ----
ResiduoCS <- coxsnell_flexsurvreg(modelo_final)$est
stopifnot(length(ResiduoCS) == nrow(df))   # mesma amostra do ajuste
ekm_cs <- survfit(Surv(ResiduoCS, df$cens) ~ 1, type = "kaplan-meier")

png(file.path(dir_plots, "coxsnell.png"), width = 1800, height = 800, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
plot(ekm_cs, fun = "cumhaz", col = cor_verde, conf.int = FALSE, lwd = 2,
     xlab = "Resíduos de Cox-Snell", ylab = "Função de risco acumulado")
abline(0, 1, col = cor_media, lwd = 2)
legend("topleft", c("KM dos resíduos", "Exp(1) esperada"),
       col = c(cor_verde, cor_media), lwd = 2, bty = "n")
plot(ekm_cs, col = cor_verde, conf.int = FALSE, lwd = 2,
     xlab = "Resíduos de Cox-Snell", ylab = "Função de sobrevivência")
cs_ord <- sort(ResiduoCS)
lines(cs_ord, exp(-cs_ord), col = cor_media, lwd = 2)
legend("topright", c("KM dos resíduos", "Exp(1) esperada"),
       col = c(cor_verde, cor_media), lwd = 2, bty = "n")
dev.off()

# 5.11 Diagnóstico: resíduos quantílicos aleatorizados (Cap. 16.6.2) ----
# Semente fixa: o QQ-plot dos slides não muda a cada renderização.
set.seed(2026)

# Caminho rápido: S(t_i | x_i) calculado direto da distribuição ajustada,
# para todos os pacientes de uma vez (equivale ao loop do livro).
S_vec <- tryCatch({
  pars <- modelo_final$dlist$pars
  loc  <- modelo_final$dlist$location
  res  <- modelo_final$res
  X <- model.matrix(fml(vars_final), df)[, -1, drop = FALSE]
  termos <- setdiff(rownames(res), pars)
  b <- res[termos, "est"]
  stopifnot(all(colnames(X) %in% names(b)))
  eta <- if (length(termos) > 0) as.numeric(X %*% b[colnames(X)]) else rep(0, nrow(df))
  tr  <- modelo_final$dlist$transforms[[loc]]
  inv <- modelo_final$dlist$inv.transforms[[loc]]
  par_loc <- inv(tr(res[loc, "est"]) + eta)   # ex.: meanlog_i, scale_i, rate_i...
  args <- setNames(lapply(pars, function(p) {
    if (p == loc) par_loc else rep(res[p, "est"], nrow(df))
  }), pars)
  1 - do.call(modelo_final$dfns$p, c(list(q = df$tempo), args))
}, error = function(e) NULL)

# Verificação contra o cálculo do Cap. 16 (summary por paciente) em 10 indivíduos
verifica <- function(S) {
  idx <- seq_len(min(10, nrow(df)))
  s_ref <- vapply(idx, function(i) {
    summary(modelo_final, type = "survival", newdata = df[i, , drop = FALSE],
            t = df$tempo[i], ci = FALSE)[[1]]$est
  }, numeric(1))
  isTRUE(all.equal(as.numeric(S[idx]), s_ref, tolerance = 1e-7))
}

if (!is.null(S_vec) && verifica(S_vec)) {
  cdf <- 1 - S_vec
  message("Resíduos quantílicos: cálculo vetorizado (verificado contra o loop do Cap. 16).")
} else {
  message("Resíduos quantílicos: loop do Cap. 16 (pode levar alguns minutos).")
  cdf <- numeric(nrow(df))
  for (i in seq_len(nrow(df))) {
    s_i <- summary(modelo_final, type = "survival", newdata = df[i, , drop = FALSE],
                   t = df$tempo[i], ci = FALSE)[[1]]$est
    cdf[i] <- 1 - s_i
  }
}

# Evento: u = F(t); censura: u ~ U(F(t), 1). Mesma ordem de sorteios do loop.
u <- cdf
idx_cens <- which(df$cens == 0)
u[idx_cens] <- runif(length(idx_cens), min = cdf[idx_cens], max = 1)
u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
Residuo <- qnorm(u)

p_qq_quantilico <- ggplot(data.frame(Residuo), aes(sample = Residuo)) +
  stat_qq_band(fill = "lightblue") +
  stat_qq_line(color = cor_media) +
  stat_qq_point(color = cor_verde, size = 1) +
  xlab("Quantis Teóricos") + ylab("Quantis Amostrais") +
  theme_bw() + theme(text = element_text(size = 14))
salva_plot(p_qq_quantilico, "qq_quantilico", w = 9, h = 7)

# 5.12 Interpretação: razões de tempos medianos (RTM) ----
# Exponencial e Gama: efeito na TAXA -> RTM = exp(-coef) (IC invertido).
# Demais (efeito na escala/locação): RTM = exp(coef). IDADE10: por 10 anos.
tidy(modelo_final)   # estimativas na escala dos coeficientes (como no livro)

res_mod    <- modelo_final$res
termos_cov <- setdiff(rownames(res_mod), modelo_final$dlist$pars)
sinal      <- if (dist_final %in% c("exponential", "gamma")) -1 else 1
est <- exp(sinal * res_mod[termos_cov, "est"])
lo  <- exp(sinal * res_mod[termos_cov, "L95%"])
up  <- exp(sinal * res_mod[termos_cov, "U95%"])
if (sinal == -1) { tmp <- lo; lo <- up; up <- tmp }
# p de Wald (mesmo teste exibido pelo tidy/print do flexsurv)
p_wald <- 2 * pnorm(-abs(res_mod[termos_cov, "est"] / res_mod[termos_cov, "se"]))

tab_rtm <- data.frame(
  Termo    = termos_cov,
  RTM      = as.numeric(est),
  IC_inf   = as.numeric(lo),
  IC_sup   = as.numeric(up),
  Variacao = (as.numeric(est) - 1) * 100,
  p_valor  = as.numeric(p_wald),
  row.names = NULL
)
tab_rtm$Significativo <- tab_rtm$IC_inf > 1 | tab_rtm$IC_sup < 1

# Rótulos legíveis (trata interações "A:B" -> "A × B")
rotulos_var <- c(ECGRUP = "Estadiamento", CIRURGIA = "Cirurgia",
                 QUIMIO = "Quimioterapia", ESCOLARI = "Escolaridade",
                 SEXO = "Sexo", IDADE10 = "Idade (+10 anos)")
rotulo_simples <- function(t) {
  if (t %in% names(rotulos_var)) return(unname(rotulos_var[t]))
  for (v in setdiff(names(rotulos_var), "IDADE10")) {
    if (startsWith(t, v))
      return(paste0(rotulos_var[[v]], ": ", substring(t, nchar(v) + 1)))
  }
  t
}
rotulo_termo <- function(t) {
  partes <- strsplit(t, ":", fixed = TRUE)[[1]]
  paste(vapply(partes, rotulo_simples, character(1)), collapse = " × ")
}
tab_rtm$Rotulo <- vapply(tab_rtm$Termo, rotulo_termo, character(1))
tab_rtm

# Leitura: RTM < 1 -> tempo mediano menor que o da referência; RTM > 1 -> maior.
# Com interações retidas, efeitos principais valem na categoria de referência
# da outra variável: interprete em conjunto.
p_forest_rtm <- ggplot(tab_rtm, aes(x = RTM, y = reorder(Rotulo, RTM),
                                    color = Significativo)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = cor_media) +
  geom_pointrange(aes(xmin = IC_inf, xmax = IC_sup), linewidth = 0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("TRUE" = cor_verde, "FALSE" = cor_cinza),
                     labels = c("TRUE" = "IC95% não inclui 1",
                                "FALSE" = "IC95% inclui 1"),
                     name = NULL) +
  labs(title = "Razão de Tempos Medianos (RTM) e IC95%",
       subtitle = paste0("Modelo ", dist_final_nome,
                         "  ·  RTM < 1: tempo mediano menor que o da referência"),
       x = "RTM (escala logarítmica)", y = NULL) +
  tema_km +
  theme(panel.grid.major.x = element_line(color = "grey92", linewidth = 0.5),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9.5))
salva_plot(p_forest_rtm, "forest_rtm", w = 9.5, h = 7)

# 5.13 Perfis previstos: favorável x desfavorável ----
# Perfis ilustrativos: nível de maior/menor RTM principal (idade: quartis).
# Com interações retidas, revise e defina manualmente abaixo.
vars_principais <- vars_final[!grepl("*", vars_final, fixed = TRUE)]
rtm_de <- function(termo) {
  i <- match(termo, tab_rtm$Termo)
  if (is.na(i)) NA_real_ else tab_rtm$RTM[i]
}
escolhe <- function(v, favoravel = TRUE) {
  if (is.factor(df[[v]])) {
    niv <- levels(df[[v]])
    ef <- vapply(niv, function(n) {
      if (n == niv[1]) 1 else {
        r <- rtm_de(paste0(v, n)); if (is.na(r)) 1 else r
      }
    }, numeric(1))
    factor(niv[if (favoravel) which.max(ef) else which.min(ef)], levels = niv)
  } else {
    r <- rtm_de(v)
    q <- if ((r >= 1) == favoravel) 0.75 else 0.25
    as.numeric(quantile(df[[v]], q))
  }
}
if (length(inter_ret) > 0)
  warning("Interações retidas: perfis automáticos são ilustrativos; revise manualmente.")

perfil_fav <- as.data.frame(setNames(lapply(vars_principais, escolhe, TRUE),  vars_principais))
perfil_des <- as.data.frame(setNames(lapply(vars_principais, escolhe, FALSE), vars_principais))
# Overrides manuais (exemplo):
# perfil_fav$ECGRUP <- factor("Estágio I", levels = levels(df$ECGRUP))
perfil_fav
perfil_des

t_grid <- seq(0.05, max(df$tempo), length.out = 200)   # definido aqui (não depender da Etapa I)
s_fav <- summary(modelo_final, newdata = perfil_fav, type = "survival",
                 t = t_grid, ci = FALSE)[[1]]
s_des <- summary(modelo_final, newdata = perfil_des, type = "survival",
                 t = t_grid, ci = FALSE)[[1]]
df_perfis <- rbind(
  data.frame(Perfil = "Favorável",    tempo = s_fav$time, S = s_fav$est),
  data.frame(Perfil = "Desfavorável", tempo = s_des$time, S = s_des$est)
)
p_perfis <- ggplot(df_perfis, aes(x = tempo, y = S, color = Perfil)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c("Favorável" = cor_verde, "Desfavorável" = cor_media)) +
  scale_x_continuous(breaks = seq(0, 144, 24)) +
  labs(title = "Sobrevida prevista: perfil favorável × desfavorável",
       subtitle = paste0("Modelo ", dist_final_nome, "  ·  ", subtitulo),
       x = "Tempo desde o diagnóstico (meses)", y = "S(t)", color = NULL) +
  tema_km
salva_plot(p_perfis, "curvas_perfis", w = 10, h = 6.5)

# Tempos medianos previstos por perfil (com IC) -> tabela para o slide
med_fav <- summary(modelo_final, newdata = perfil_fav, type = "median")[[1]]
med_des <- summary(modelo_final, newdata = perfil_des, type = "median")[[1]]
tab_medianas <- data.frame(
  Perfil  = c("Favorável", "Desfavorável"),
  Mediana = c(med_fav$est, med_des$est),
  IC_inf  = c(med_fav$lcl, med_des$lcl),
  IC_sup  = c(med_fav$ucl, med_des$ucl)
)
tab_medianas

# 6. Análise de sensibilidade (casos completos) --------------------------------
tab_sens <- NULL
if (ausentes == "categoria") {
  df_cc <- df %>%
    filter(ESCOLARI != "Sem informação",
           !ECGRUP %in% c("Sem informação", "Não se aplica")) %>%
    mutate(ESCOLARI = droplevels(ESCOLARI), ECGRUP = droplevels(ECGRUP))
  m_cc <- tryCatch(ajusta(vars_final, dist = dist_final, data = df_cc),
                   error = function(e) {
                     message("Sensibilidade não ajustada: ", conditionMessage(e)); NULL
                   })
  if (!is.null(m_cc)) {
    res_cc  <- m_cc$res
    term_cc <- setdiff(rownames(res_cc), m_cc$dlist$pars)
    tab_sens <- data.frame(Termo = term_cc,
                           RTM_cc = exp(sinal * res_cc[term_cc, "est"]),
                           row.names = NULL)
    tab_sens$RTM_principal <- tab_rtm$RTM[match(tab_sens$Termo, tab_rtm$Termo)]
    tab_sens$Mesma_direcao <- with(tab_sens, (RTM_cc > 1) == (RTM_principal > 1))
    tab_sens
    # Conclusões estáveis: mesma direção e magnitude parecida nos dois ajustes.
  }
}

# 7. Tabelas para os slides -----------------------------------------------------
# Etapa I: lida do arquivo salvo (fallback: reajusta os modelos nulos)
arq_t1 <- file.path(dir_saida, "tabelas_parte1.RDS")
if (file.exists(arq_t1)) {
  tabelas_parte1 <- readRDS(arq_t1)
  tab_ic1 <- tabelas_parte1$ic; tab_trv1 <- tabelas_parte1$trv
} else {
  a1 <- lapply(dists, function(d) ajusta(character(0), dist = d))
  tab_ic1 <- tabela_ic(a1); tab_trv1 <- tabela_trv_gg(a1)
  tab_trv1$Decisao <- ifelse(tab_trv1$P_valor > 0.05, "compatível com GG", "rejeitado")
  tabelas_parte1 <- NULL
}

tabelas_parte2 <- list(
  ic_etapa1 = tab_ic1, trv_etapa1 = tab_trv1,
  filtro = filtro,
  passo1 = passo1, passo2 = passo2, passo3 = passo3, passo4 = passo4, passo5 = passo5,
  ic_etapa2 = tab_ic2, trv_etapa2 = tab_trv2,
  reavaliacao = tab_reavalia, rtm = tab_rtm,
  sensibilidade = tab_sens, medianas_perfis = tab_medianas,
  base = base_analitica
)
saveRDS(tabelas_parte2, file.path(dir_saida, "tabelas_parte2.RDS"))

# Versões FORMATADAS (vírgula decimal, "< 0,001") para preencher o index.qmd ----
fmt_ic <- function(t) data.frame(Modelo = t$Modelo, logLik = fmt_num(t$logLik, 1),
                                 AIC = fmt_num(t$AIC, 1), BIC = fmt_num(t$BIC, 1))
fmt_trv <- function(t) data.frame(Comparacao = t$Comparacao, TRV = fmt_num(t$TRV, 2),
                                  gl = t$gl, p_valor = formatar_p(t$P_valor),
                                  Decisao = t$Decisao)
fmt_passo <- function(t) {
  if (is.null(t)) return(NULL)
  data.frame(Variavel = gsub("*", " × ", t$Variavel, fixed = TRUE),
             TRV = fmt_num(t$TRV, 2), gl = t$gl,
             p_valor = formatar_p(t$P_valor), Decisao = t$Decisao)
}

filtro_slide <- data.frame(
  Variavel = c(filtro$Variavel, "IDADE10"),
  Justificativa = c(unname(justificativas[filtro$Variavel]),
                    "Quantitativa; entra por justificativa clínica (por 10 anos)"),
  Logrank = c(formatar_p(filtro$Logrank_pvalor), "—"),
  Peto    = c(formatar_p(filtro$Peto_pvalor), "—"),
  Candidata = c(ifelse(filtro$Variavel %in% candidatas, "sim", "não"), "sim")
)

# Tabela combinada do slide "Distribuição com Covariáveis"
p_gg  <- tab_trv2$P_valor[match(tab_ic2$Modelo, tab_trv2$Modelo)]
dec2  <- tab_trv2$Decisao[match(tab_ic2$Modelo, tab_trv2$Modelo)]
dist_etapa2_slide <- data.frame(
  Modelo = tab_ic2$Modelo, AIC = fmt_num(tab_ic2$AIC, 1), BIC = fmt_num(tab_ic2$BIC, 1),
  p_TRV_GG = formatar_p(p_gg),
  Decisao  = ifelse(is.na(p_gg), "referência (GG)", dec2)
)

rtm_slide <- data.frame(
  Efeito = tab_rtm$Rotulo,
  `RTM (IC95%)` = sprintf("%s (%s–%s)", fmt_num(tab_rtm$RTM, 2),
                          fmt_num(tab_rtm$IC_inf, 2), fmt_num(tab_rtm$IC_sup, 2)),
  `Variação (%)` = fmt_num(tab_rtm$Variacao, 1),
  p_valor = formatar_p(tab_rtm$p_valor),
  Significativo = ifelse(tab_rtm$Significativo, "sim", "não"),
  check.names = FALSE
)

decisoes2 <- list(
  dist_etapa1 = if (!is.null(tabelas_parte1)) tabelas_parte1$decisoes$sugerida else NA_character_,
  dist_final  = dist_final_nome,
  vars_final  = vars_final,
  formula     = paste("Surv(tempo, cens) ~", paste(vars_final, collapse = " + ")),
  interacoes  = inter_ret,
  forma_risco = forma_risco(dist_final_nome, modelo_final),
  nota_quimio = if ("QUIMIO" %in% vars_final)
    "QUIMIO: curvas se cruzam na Parte I; a RTM resume o efeito médio e deve ser lida com cautela (possível violação de tempo acelerado)."
  else NULL
)

tabelas_slides <- list(
  base = base_analitica, decisoes = decisoes2,
  ic_etapa1 = if (!is.null(tabelas_parte1)) tabelas_parte1$ic_fmt else fmt_ic(tab_ic1),
  trv_etapa1 = if (!is.null(tabelas_parte1)) tabelas_parte1$trv_fmt else fmt_trv(tab_trv1),
  filtro = filtro_slide,
  passo1 = fmt_passo(passo1), passo2 = fmt_passo(passo2), passo3 = fmt_passo(passo3),
  passo4 = fmt_passo(passo4), passo5 = fmt_passo(passo5),
  ic_etapa2 = fmt_ic(tab_ic2), trv_etapa2 = fmt_trv(tab_trv2),
  dist_etapa2 = dist_etapa2_slide,
  reavaliacao = if (is.null(tab_reavalia)) NULL else data.frame(
    Variavel = gsub("*", " × ", tab_reavalia$Variavel, fixed = TRUE),
    TRV = fmt_num(tab_reavalia$TRV, 2), gl = tab_reavalia$gl,
    p_valor = formatar_p(tab_reavalia$P_valor),
    dAIC = fmt_num(tab_reavalia$dAIC, 1), dBIC = fmt_num(tab_reavalia$dBIC, 1)),
  rtm = rtm_slide,
  sensibilidade = if (is.null(tab_sens)) NULL else data.frame(
    Termo = tab_sens$Termo, RTM_principal = fmt_num(tab_sens$RTM_principal, 2),
    RTM_cc = fmt_num(tab_sens$RTM_cc, 2),
    Mesma_direcao = ifelse(tab_sens$Mesma_direcao, "sim", "não")),
  medianas_perfis = data.frame(
    Perfil = tab_medianas$Perfil,
    `Mediana (IC95%)` = sprintf("%s (%s–%s)", fmt_num(tab_medianas$Mediana, 1),
                                fmt_num(tab_medianas$IC_inf, 1),
                                fmt_num(tab_medianas$IC_sup, 1)),
    check.names = FALSE)
)
saveRDS(tabelas_slides, file.path(dir_saida, "tabelas_slides.RDS"))

cat("\nn analítico:", nrow(df), "| óbitos:", sum(df$cens == 1), "| ausentes:", ausentes,
    "\nCandidatas:", paste(candidatas, collapse = ", "),
    "\nModelo final: ~", paste(vars_final, collapse = " + "),
    "\nDistribuição final:", dist_final_nome, "|", decisoes2$forma_risco,
    "\nSlides: leia 'data/tabelas_slides.RDS' e 'data/base_analitica.RDS' no index.qmd\n")