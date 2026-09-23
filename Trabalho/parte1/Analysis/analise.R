# 1. Setup ---------------------------------------------------------------------

library(arrow)
library(dplyr)
library(survival)
library(survminer)
library(ggsurvfit)
library(ggplot2)
library(gtsummary)
library(multcompView)

Letrinhas <- function(p_mat, alpha = 0.05) {
  p_mat <- as.matrix(p_mat)
  mode(p_mat) <- "numeric"
  rownames(p_mat) <- as.character(rownames(p_mat))
  colnames(p_mat) <- as.character(colnames(p_mat))
  grupos <- union(rownames(p_mat), colnames(p_mat))
  mat_full <- matrix(NA,
                     nrow = length(grupos),
                     ncol = length(grupos),
                     dimnames = list(grupos, grupos))
  mat_full[rownames(p_mat), colnames(p_mat)] <- p_mat
  for (i in seq_len(nrow(mat_full))) {
    for (j in seq_len(ncol(mat_full))) {
      if (is.na(mat_full[i, j]) && !is.na(mat_full[j, i])) {
        mat_full[i, j] <- mat_full[j, i]
      }
    }
  }
  diag(mat_full) <- 1
  rownames(mat_full) <- colnames(mat_full)
  multcompView::multcompLetters(mat_full, threshold = alpha)$Letters
}

# Converte uma letra normal em negrito Unicode (𝐀, 𝐁, 𝐂, ...).
# Não depende do pacote ggtext, então funciona com qualquer versão do ggplot2.
letra_negrito <- function(letra) {
  negrito <- c("𝐀","𝐁","𝐂","𝐃","𝐄","𝐅","𝐆","𝐇","𝐈","𝐉","𝐊","𝐋","𝐌",
               "𝐍","𝐎","𝐏","𝐐","𝐑","𝐒","𝐓","𝐔","𝐕","𝐖","𝐗","𝐘","𝐙")
  idx <- match(toupper(letra), LETTERS)
  negrito[idx]
}

# Constrói os rótulos "Nível (LETRA)": nome do grupo em texto normal,
# só a letra da comparação múltipla em negrito Unicode.
rotulos_com_letras <- function(niveis, letras) {
  paste0(niveis, " (", sapply(letras[niveis], letra_negrito), ")")
}

# 2. Leitura e Filtros ---------------------------------------------------------

df <- read_parquet("Trabalho/parte1/data/FOSP.parquet") %>%
  mutate(
    across(where(is.character), ~ iconv(.x, from = "latin1", to = "UTF-8"))
  ) %>%
  filter(
    TOPOGRUP %in% c("C22"),
    ANODIAG >= 2014,
    ANODIAG <= 2019,
    IDADE >= 18,
    DIAGPREV %in% c(1, 2),
    ERRO == 0,
    !is.na(DTDIAG),
    !is.na(DTULTINFO),
    ULTINFO %in% c(1, 2, 3, 4),
    DTULTINFO >= DTDIAG
  ) %>%
  select(
    TOPOGRUP,
    IDADE, SEXO, ESCOLARI,
    ECGRUP,
    CATEATEND, CIRURGIA, RADIO, QUIMIO,
    DTDIAG, DTULTINFO, ULTINFO
  ) %>%
  mutate(
    tempo = as.numeric(DTULTINFO - DTDIAG) / 30.44,
    cens  = ifelse(ULTINFO == 3, 1, 0)
  )

# 3. Recodificação das Variáveis ----------------------------------------------

df <- df %>%
  mutate(
    SEXO = factor(SEXO,
                  levels = c(1, 2),
                  labels = c("Masculino", "Feminino")),

    ESCOLARI = factor(ESCOLARI,
                      levels = c(1, 2, 3, 4, 5, 9),
                      labels = c("Analfabeto", "Fund. incompleto",
                                 "Fund. completo", "Ensino médio",
                                 "Ensino superior", "Sem informação")),

    ECGRUP = factor(ECGRUP,
                    levels = c("0", "I", "II", "III", "IV", "X", "Y"),
                    labels = c("In situ", "Estágio I", "Estágio II",
                               "Estágio III", "Estágio IV",
                               "Sem informação", "Não se aplica")),

    CATEATEND = factor(CATEATEND,
                       levels = c(1, 2, 3),
                       labels = c("SUS", "Convênio", "Particular")),

    CIRURGIA = factor(CIRURGIA,
                      levels = c(0, 1),
                      labels = c("Não", "Sim")),

    RADIO = factor(RADIO,
                   levels = c(0, 1),
                   labels = c("Não", "Sim")),

    QUIMIO = factor(QUIMIO,
                    levels = c(0, 1),
                    labels = c("Não", "Sim"))
  )

saveRDS(df, "Trabalho/parte1/data/df_final.RDS")

# 4. Análise Descritiva --------------------------------------------------------

nrow(df)
sum(df$cens == 1)
sum(df$cens == 0)

summary(df)

# 5. Análise Sem Covariáveis ---------------------------------------------------

df <- readRDS("Trabalho/parte1/data/df_final.RDS")

# Paleta
cor_navy    <- "#10241d"
cor_verde   <- "#3f6b52"
cor_dourado <- "#b8862a"
cor_cinza   <- "#6b7d73"
cor_media   <- "#c0392b"

# Tema
tema_km <- theme_minimal(base_size = 13) +
  theme(
    plot.title.position = "plot",
    plot.title         = element_text(color = cor_navy, face = "bold",
                                      size = 17, margin = margin(b = 4)),
    plot.subtitle      = element_text(color = cor_cinza, size = 10.5,
                                      margin = margin(b = 10)),
    axis.title         = element_text(color = cor_navy, size = 11),
    axis.title.x       = element_text(margin = margin(t = 8)),
    axis.title.y       = element_text(margin = margin(r = 8)),
    axis.text          = element_text(color = cor_cinza, size = 10),
    axis.line.x        = element_line(color = cor_cinza, linewidth = 0.4),
    axis.ticks.x       = element_line(color = cor_cinza, linewidth = 0.4),
    legend.position    = "top",
    legend.justification = "left",
    legend.title       = element_text(color = cor_navy, face = "bold", size = 10.5),
    legend.text        = element_text(color = cor_navy, size = 10),
    legend.key.width   = unit(1.3, "cm"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.5),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(14, 18, 10, 14)
  )

# Textos de apoio
n_total  <- format(nrow(df), big.mark = ".", decimal.mark = ",")
n_obitos <- format(sum(df$cens == 1), big.mark = ".", decimal.mark = ",")

subtitulo <- paste0(
  "FOSP, diagnósticos de 2014 a 2019  ·  n = ", n_total,
  " pacientes  ·  ", n_obitos, " óbitos por câncer"
)

# Ajuste (IC log-log de 95%)
km_geral <- survfit2(
  Surv(tempo, cens) ~ 1,
  data      = df,
  conf.type = "log-log",
  conf.int  = 0.95
)

mediana <- summary(km_geral)$table["median"]
sobrev  <- summary(km_geral, times = c(12, 60))$surv

texto_resumo <- paste0(
  "Mediana: ", format(round(mediana, 1), decimal.mark = ","), " meses\n",
  "Sobrevida em 12 meses: ", round(100 * sobrev[1]), "%\n",
  "Sobrevida em 60 meses: ", round(100 * sobrev[2]), "%"
)

# Gráfico
curva_sem_covariavel <- ggsurvfit(
  km_geral,
  linewidth = 1.3,
  color     = cor_verde
) +
  add_confidence_interval(fill = cor_verde, alpha = 0.18) +
  add_quantile(
    y_value   = 0.5,
    color     = cor_media,
    linetype  = "dashed",
    linewidth = 0.7
  ) +
  annotate(
    "label",
    x          = 60,
    y          = 0.80,
    label      = texto_resumo,
    hjust      = 0,
    vjust      = 1,
    fill       = "#f7f4ea",
    color      = cor_navy,
    size       = 3.8,
    lineheight = 1.25
  ) +
  add_risktable(
    risktable_stats = "n.risk",
    stats_label     = list(n.risk = "Em risco"),
    size            = 3.6,
    theme           = list(
      theme_risktable_default(),
      theme(plot.title = element_text(face = "bold", color = cor_navy))
    )
  ) +
  scale_ggsurvfit(
    x_scales = list(breaks = seq(0, 144, 24))
  ) +
  labs(
    title    = "Sobrevida global — Câncer de Fígado (C22)",
    subtitle = subtitulo,
    x        = "Tempo desde o diagnóstico (meses)",
    y        = "S(t)"
  ) +
  tema_km

curva_sem_covariavel

saveRDS(
  curva_sem_covariavel,
  "Trabalho/parte1/plots/curva_sem_covariavel.RDS"
)

# 6. Análise por Covariáveis Qualitativas --------------------------------------

df <- readRDS("Trabalho/parte1/data/df_final.RDS")

# Paleta
cor_navy    <- "#10241d"
cor_verde   <- "#3f6b52"
cor_dourado <- "#b8862a"
cor_cinza   <- "#6b7d73"
cor_media   <- "#c0392b"

paleta_2 <- c(cor_verde, cor_dourado)

# Tema
tema_km <- theme_minimal(base_size = 13) +
  theme(
    plot.title.position = "plot",
    plot.title         = element_text(color = cor_navy, face = "bold",
                                      size = 17, margin = margin(b = 4)),
    plot.subtitle      = element_text(color = cor_cinza, size = 10.5,
                                      margin = margin(b = 10)),
    axis.title         = element_text(color = cor_navy, size = 11),
    axis.title.x       = element_text(margin = margin(t = 8)),
    axis.title.y       = element_text(margin = margin(r = 8)),
    axis.text          = element_text(color = cor_cinza, size = 10),
    axis.line.x        = element_line(color = cor_cinza, linewidth = 0.4),
    axis.ticks.x       = element_line(color = cor_cinza, linewidth = 0.4),
    legend.position    = "top",
    legend.justification = "left",
    legend.title       = element_text(color = cor_navy, face = "bold", size = 10.5),
    legend.text        = element_text(color = cor_navy, size = 10),
    legend.key.width   = unit(1.3, "cm"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.5),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(14, 18, 10, 14)
  )

# Textos de apoio
n_total  <- format(nrow(df), big.mark = ".", decimal.mark = ",")
n_obitos <- format(sum(df$cens == 1), big.mark = ".", decimal.mark = ",")

subtitulo <- paste0(
  "Câncer de fígado (C22)  ·  FOSP 2014–2019  ·  n = ", n_total,
  "  ·  ", n_obitos, " óbitos"
)


## SEXO ----

km_sexo <- survfit2(Surv(tempo, cens) ~ SEXO, data = df)

res_sexo    <- pairwise_survdiff(Surv(tempo, cens) ~ SEXO, data = df,
                                  p.adjust.method = "holm", rho = 0)
letras_sexo <- Letrinhas(res_sexo$p.value)

curva_sexo <- ggsurvfit(km_sexo, linewidth = 1.2) +
  add_confidence_interval(alpha = 0.12) +
  add_quantile(
    y_value   = 0.5,
    color     = "grey55",
    linetype  = "dashed",
    linewidth = 0.4
  ) +
  add_risktable(
    risktable_stats = "n.risk",
    stats_label     = list(n.risk = "Em risco"),
    size            = 3.6,
    theme           = list(
      theme_risktable_default(),
      theme(plot.title = element_text(face = "bold", color = cor_navy))
    )
  ) +
  add_risktable_strata_symbol(symbol = "\U25CF", size = 12) +
  add_pvalue(
    location = "annotation",
    caption  = "{p.value}",
    x        = 1,
    y        = 0.1,
    hjust    = 0,
    size     = 4.2,
    fontface = "bold",
    color    = cor_navy
  ) +
  scale_ggsurvfit(x_scales = list(breaks = seq(0, 144, 24))) +
  scale_color_manual(
    values = paleta_2,
    labels = rotulos_com_letras(levels(df$SEXO), letras_sexo)
  ) +
  scale_fill_manual(values = paleta_2) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.6)),
    fill  = "none"
  ) +
  labs(
    title    = "Sobrevida por Sexo",
    subtitle = subtitulo,
    x        = "Tempo desde o diagnóstico (meses)",
    y        = "S(t)",
    color    = "Sexo"
  ) +
  tema_km

curva_sexo

saveRDS(curva_sexo, "Trabalho/parte1/plots/curva_sexo.RDS")


## ESCOLARI ----
# Escala ordinal (dourado -> verde -> navy). "Sem informação" removida do gráfico.

df_escolari <- df %>%
  filter(ESCOLARI != "Sem informação") %>%
  mutate(ESCOLARI = droplevels(ESCOLARI))

# Subtítulo com o n do subconjunto
subtitulo_escolari <- paste0(
  "Câncer de fígado (C22)  ·  FOSP 2014–2019  ·  n = ",
  format(nrow(df_escolari), big.mark = ".", decimal.mark = ","),
  "  ·  ",
  format(sum(df_escolari$cens == 1), big.mark = ".", decimal.mark = ","),
  " óbitos"
)

km_escolari <- survfit2(Surv(tempo, cens) ~ ESCOLARI, data = df_escolari)

res_escolari    <- pairwise_survdiff(Surv(tempo, cens) ~ ESCOLARI, data = df_escolari,
                                      p.adjust.method = "holm", rho = 0)
letras_escolari <- Letrinhas(res_escolari$p.value)

curva_escolari <- ggsurvfit(km_escolari, linewidth = 1.2) +
  add_risktable(
    risktable_stats = "n.risk",
    stats_label     = list(n.risk = "Em risco"),
    size            = 3.4,
    theme           = list(
      theme_risktable_default(),
      theme(plot.title = element_text(face = "bold", color = cor_navy))
    )
  ) +
  add_risktable_strata_symbol(symbol = "\U25CF", size = 12) +
  add_pvalue(
    location = "annotation",
    caption  = "{p.value}",
    x        = 1,
    y        = 0.1,
    hjust    = 0,
    size     = 4.2,
    fontface = "bold",
    color    = cor_navy
  ) +
  scale_ggsurvfit(x_scales = list(breaks = seq(0, 144, 24))) +
  scale_color_manual(
    values = colorRampPalette(c("#d9b054", cor_verde, cor_navy))(5),
    labels = rotulos_com_letras(levels(df_escolari$ESCOLARI), letras_escolari)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(linewidth = 1.6))) +
  labs(
    title    = "Sobrevida por Escolaridade",
    subtitle = subtitulo_escolari,
    x        = "Tempo desde o diagnóstico (meses)",
    y        = "S(t)",
    color    = "Escolaridade"
  ) +
  tema_km

curva_escolari

saveRDS(curva_escolari, "Trabalho/parte1/plots/curva_escolari.RDS")

## ECGRUP ----
# Apenas estágios I a IV ("In situ" não tem observações em C22;
# "Sem informação" e "Não se aplica" foram removidos do gráfico).
# Estágios I a IV: verde -> dourado -> vermelho (gravidade crescente).

df_ecgrup <- df %>%
  filter(ECGRUP %in% c("Estágio I", "Estágio II", "Estágio III", "Estágio IV")) %>%
  mutate(ECGRUP = droplevels(ECGRUP))

# Subtítulo com o n do subconjunto
subtitulo_ecgrup <- paste0(
  "Câncer de fígado (C22)  ·  FOSP 2014–2019  ·  n = ",
  format(nrow(df_ecgrup), big.mark = ".", decimal.mark = ","),
  "  ·  ",
  format(sum(df_ecgrup$cens == 1), big.mark = ".", decimal.mark = ","),
  " óbitos"
)

km_ecgrup <- survfit2(Surv(tempo, cens) ~ ECGRUP, data = df_ecgrup)

res_ecgrup    <- pairwise_survdiff(Surv(tempo, cens) ~ ECGRUP, data = df_ecgrup,
                                    p.adjust.method = "holm", rho = 0)
letras_ecgrup <- Letrinhas(res_ecgrup$p.value)

curva_ecgrup <- ggsurvfit(km_ecgrup, linewidth = 1.2) +
  add_risktable(
    risktable_stats = "n.risk",
    stats_label     = list(n.risk = "Em risco"),
    size            = 3.4,
    theme           = list(
      theme_risktable_default(),
      theme(plot.title = element_text(face = "bold", color = cor_navy))
    )
  ) +
  add_risktable_strata_symbol(symbol = "\U25CF", size = 12) +
  add_pvalue(
    location = "annotation",
    caption  = "{p.value}",
    x        = 1,
    y        = 0.1,
    hjust    = 0,
    size     = 4.2,
    fontface = "bold",
    color    = cor_navy
  ) +
  scale_ggsurvfit(x_scales = list(breaks = seq(0, 144, 24))) +
  scale_color_manual(
    values = colorRampPalette(c(cor_verde, cor_dourado, cor_media))(4),
    labels = rotulos_com_letras(levels(df_ecgrup$ECGRUP), letras_ecgrup)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.6))) +
  labs(
    title    = "Sobrevida por Estadiamento",
    subtitle = subtitulo_ecgrup,
    x        = "Tempo desde o diagnóstico (meses)",
    y        = "S(t)",
    color    = "Estadiamento"
  ) +
  tema_km

curva_ecgrup

saveRDS(curva_ecgrup, "Trabalho/parte1/plots/curva_ecgrup.RDS")

## CATEATEND ----

km_cateatend <- survfit2(Surv(tempo, cens) ~ CATEATEND, data = df)

res_cateatend    <- pairwise_survdiff(Surv(tempo, cens) ~ CATEATEND, data = df,
                                       p.adjust.method = "holm", rho = 0)
letras_cateatend <- Letrinhas(res_cateatend$p.value)

curva_cateatend <- ggsurvfit(km_cateatend, linewidth = 1.2) +
  add_quantile(
    y_value   = 0.5,
    color     = "grey55",
    linetype  = "dashed",
    linewidth = 0.4
  ) +
  add_risktable(
    risktable_stats = "n.risk",
    stats_label     = list(n.risk = "Em risco"),
    size            = 3.6,
    theme           = list(
      theme_risktable_default(),
      theme(plot.title = element_text(face = "bold", color = cor_navy))
    )
  ) +
  add_risktable_strata_symbol(symbol = "\U25CF", size = 12) +
  add_pvalue(
    location = "annotation",
    caption  = "{p.value}",
    x        = 1,
    y        = 0.1,
    hjust    = 0,
    size     = 4.2,
    fontface = "bold",
    color    = cor_navy
  ) +
  scale_ggsurvfit(x_scales = list(breaks = seq(0, 144, 24))) +
  scale_color_manual(
    values = c(cor_verde, cor_dourado, cor_navy),
    labels = rotulos_com_letras(levels(df$CATEATEND), letras_cateatend)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.6))) +
  labs(
    title    = "Sobrevida por Categoria de Atendimento",
    subtitle = subtitulo,
    x        = "Tempo desde o diagnóstico (meses)",
    y        = "S(t)",
    color    = "Atendimento"
  ) +
  tema_km

curva_cateatend

saveRDS(curva_cateatend, "Trabalho/parte1/plots/curva_cateatend.RDS")


## CIRURGIA ----

km_cirurgia <- survfit2(Surv(tempo, cens) ~ CIRURGIA, data = df)

res_cirurgia    <- pairwise_survdiff(Surv(tempo, cens) ~ CIRURGIA, data = df,
                                      p.adjust.method = "holm", rho = 0)
letras_cirurgia <- Letrinhas(res_cirurgia$p.value)

curva_cirurgia <- ggsurvfit(km_cirurgia, linewidth = 1.2) +
  add_confidence_interval(alpha = 0.12) +
  add_quantile(
    y_value   = 0.5,
    color     = "grey55",
    linetype  = "dashed",
    linewidth = 0.4
  ) +
  add_risktable(
    risktable_stats = "n.risk",
    stats_label     = list(n.risk = "Em risco"),
    size            = 3.6,
    theme           = list(
      theme_risktable_default(),
      theme(plot.title = element_text(face = "bold", color = cor_navy))
    )
  ) +
  add_risktable_strata_symbol(symbol = "\U25CF", size = 12) +
  add_pvalue(
    location = "annotation",
    caption  = "{p.value}",
    x        = 1,
    y        = 0.1,
    hjust    = 0,
    size     = 4.2,
    fontface = "bold",
    color    = cor_navy
  ) +
  scale_ggsurvfit(x_scales = list(breaks = seq(0, 144, 24))) +
  scale_color_manual(
    values = paleta_2,
    labels = rotulos_com_letras(levels(df$CIRURGIA), letras_cirurgia)
  ) +
  scale_fill_manual(values = paleta_2) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.6)),
    fill  = "none"
  ) +
  labs(
    title    = "Sobrevida por Cirurgia",
    subtitle = subtitulo,
    x        = "Tempo desde o diagnóstico (meses)",
    y        = "S(t)",
    color    = "Cirurgia"
  ) +
  tema_km

curva_cirurgia

saveRDS(curva_cirurgia, "Trabalho/parte1/plots/curva_cirurgia.RDS")


## RADIO ----

km_radio <- survfit2(Surv(tempo, cens) ~ RADIO, data = df)

res_radio    <- pairwise_survdiff(Surv(tempo, cens) ~ RADIO, data = df,
                                   p.adjust.method = "holm", rho = 0)
letras_radio <- Letrinhas(res_radio$p.value)

curva_radio <- ggsurvfit(km_radio, linewidth = 1.2) +
  add_confidence_interval(alpha = 0.12) +
  add_quantile(
    y_value   = 0.5,
    color     = "grey55",
    linetype  = "dashed",
    linewidth = 0.4
  ) +
  add_risktable(
    risktable_stats = "n.risk",
    stats_label     = list(n.risk = "Em risco"),
    size            = 3.6,
    theme           = list(
      theme_risktable_default(),
      theme(plot.title = element_text(face = "bold", color = cor_navy))
    )
  ) +
  add_risktable_strata_symbol(symbol = "\U25CF", size = 12) +
  add_pvalue(
    location = "annotation",
    caption  = "{p.value}",
    x        = 1,
    y        = 0.1,
    hjust    = 0,
    size     = 4.2,
    fontface = "bold",
    color    = cor_navy
  ) +
  scale_ggsurvfit(x_scales = list(breaks = seq(0, 144, 24))) +
  scale_color_manual(
    values = paleta_2,
    labels = rotulos_com_letras(levels(df$RADIO), letras_radio)
  ) +
  scale_fill_manual(values = paleta_2) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.6)),
    fill  = "none"
  ) +
  labs(
    title    = "Sobrevida por Radioterapia",
    subtitle = subtitulo,
    x        = "Tempo desde o diagnóstico (meses)",
    y        = "S(t)",
    color    = "Radioterapia"
  ) +
  tema_km

curva_radio

saveRDS(curva_radio, "Trabalho/parte1/plots/curva_radio.RDS")


## QUIMIO ----

km_quimio <- survfit2(Surv(tempo, cens) ~ QUIMIO, data = df)

res_quimio    <- pairwise_survdiff(Surv(tempo, cens) ~ QUIMIO, data = df,
                                    p.adjust.method = "holm", rho = 0)
letras_quimio <- Letrinhas(res_quimio$p.value)

curva_quimio <- ggsurvfit(km_quimio, linewidth = 1.2) +
  add_confidence_interval(alpha = 0.12) +
  add_quantile(
    y_value   = 0.5,
    color     = "grey55",
    linetype  = "dashed",
    linewidth = 0.4
  ) +
  add_risktable(
    risktable_stats = "n.risk",
    stats_label     = list(n.risk = "Em risco"),
    size            = 3.6,
    theme           = list(
      theme_risktable_default(),
      theme(plot.title = element_text(face = "bold", color = cor_navy))
    )
  ) +
  add_risktable_strata_symbol(symbol = "\U25CF", size = 12) +
  add_pvalue(
    location = "annotation",
    caption  = "{p.value}",
    x        = 1,
    y        = 0.1,
    hjust    = 0,
    size     = 4.2,
    fontface = "bold",
    color    = cor_navy
  ) +
  scale_ggsurvfit(x_scales = list(breaks = seq(0, 144, 24))) +
  scale_color_manual(
    values = paleta_2,
    labels = rotulos_com_letras(levels(df$QUIMIO), letras_quimio)
  ) +
  scale_fill_manual(values = paleta_2) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.6)),
    fill  = "none"
  ) +
  labs(
    title    = "Sobrevida por Quimioterapia",
    subtitle = subtitulo,
    x        = "Tempo desde o diagnóstico (meses)",
    y        = "S(t)",
    color    = "Quimioterapia"
  ) +
  tema_km

curva_quimio

saveRDS(curva_quimio, "Trabalho/parte1/plots/curva_quimio.RDS")


# 7. Análise por Covariáveis Quantitativas -------------------------------------

df <- readRDS("Trabalho/parte1/data/df_final.RDS")

# Paleta
cor_navy    <- "#10241d"
cor_verde   <- "#3f6b52"
cor_dourado <- "#b8862a"
cor_cinza   <- "#6b7d73"
cor_media   <- "#c0392b"

paleta_2 <- c(cor_verde, cor_dourado)

# Tema (para o boxplot, sem legenda)
tema_km <- theme_minimal(base_size = 13) +
  theme(
    plot.title.position = "plot",
    plot.title         = element_text(color = cor_navy, face = "bold",
                                      size = 17, margin = margin(b = 4)),
    plot.subtitle      = element_text(color = cor_cinza, size = 10.5,
                                      margin = margin(b = 10)),
    axis.title         = element_text(color = cor_navy, size = 11),
    axis.title.y       = element_text(margin = margin(r = 8)),
    axis.text          = element_text(color = cor_cinza, size = 10),
    axis.text.x        = element_text(color = cor_navy, size = 10.5, lineheight = 1.1),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.5),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(14, 18, 10, 14)
  )


## IDADE ----

df_idade_box <- bind_rows(
  df %>%
    mutate(grupo = ifelse(cens == 1, "Óbito por câncer", "Censura")) %>%
    select(grupo, IDADE),

  df %>%
    mutate(grupo = "Total") %>%
    select(grupo, IDADE)
) %>%
  mutate(
    grupo = factor(
      grupo,
      levels = c("Óbito por câncer", "Censura", "Total")
    )
  )

resumo_idade <- df_idade_box %>%
  group_by(grupo) %>%
  summarise(
    n       = n(),
    media   = mean(IDADE, na.rm = TRUE),
    mediana = median(IDADE, na.rm = TRUE),
    q1      = quantile(IDADE, 0.25, na.rm = TRUE),
    q3      = quantile(IDADE, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_num = as.numeric(grupo),
    texto = paste0(
      "Mediana: ", format(round(mediana, 1), decimal.mark = ","), "\n",
      "Q1 – Q3: ", format(round(q1, 1), decimal.mark = ","), " – ",
      format(round(q3, 1), decimal.mark = ","), "\n",
      "Média: ", format(round(media, 1), decimal.mark = ",")
    )
  )

# Rótulos do eixo x com o n de cada grupo
rotulos_x <- setNames(
  paste0(
    as.character(resumo_idade$grupo), "\n(n = ",
    format(resumo_idade$n, big.mark = ".", decimal.mark = ","), ")"
  ),
  as.character(resumo_idade$grupo)
)

paleta_idade <- c(
  "Óbito por câncer" = cor_verde,
  "Censura"          = cor_dourado,
  "Total"            = cor_navy
)

curva_idade_boxplot <- ggplot(
  df_idade_box,
  aes(x = grupo, y = IDADE, fill = grupo)
) +
  geom_violin(
    color = NA,
    alpha = 0.28,
    width = 0.95
  ) +
  geom_boxplot(
    color         = cor_navy,
    alpha         = 0.85,
    width         = 0.20,
    linewidth     = 0.5,
    outlier.shape = NA
  ) +
  geom_point(
    data = resumo_idade,
    aes(x = grupo, y = media),
    inherit.aes = FALSE,
    shape       = 23,
    size        = 3.6,
    stroke      = 1,
    fill        = "white",
    color       = cor_media
  ) +
  geom_text(
    data = resumo_idade,
    aes(x = x_num + 0.30, y = mediana, label = texto),
    inherit.aes = FALSE,
    hjust       = 0,
    size        = 3.3,
    lineheight  = 1.15,
    color       = cor_navy
  ) +
  scale_fill_manual(values = paleta_idade) +
  scale_x_discrete(
    labels = rotulos_x,
    expand = expansion(add = c(0.6, 0.9))
  ) +
  scale_y_continuous(breaks = seq(20, 100, 10)) +
  labs(
    title    = "Distribuição da Idade ao Diagnóstico",
    subtitle = "Losango vermelho = média  ·  caixa = mediana e quartis  ·  formato = densidade",
    x        = NULL,
    y        = "Idade (anos)"
  ) +
  tema_km +
  theme(legend.position = "none")

curva_idade_boxplot

saveRDS(curva_idade_boxplot, "Trabalho/parte1/plots/curva_idade_boxplot.RDS")