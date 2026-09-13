# 1. Setup ---------------------------------------------------------------------

library(arrow)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(gtsummary)
library(cowplot)
library(grid)


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
    cens = ifelse(ULTINFO == 3, 1, 0)
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


# 4. Análise Descritiva --------------------------------------------------------

nrow(df)
sum(df$cens == 1)
sum(df$cens == 0)
mean(df$cens) * 100

summary(df$tempo)


# 5. Paleta, Tema e Borda dos Gráficos ----------------------------------------

cor_navy    <- "#10241d"
cor_verde   <- "#3f6b52"
cor_dourado <- "#b8862a"
cor_cinza   <- "#6b7d73"
cor_media   <- "#c0392b"


tema_km <- theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(color = cor_navy, face = "bold", size = 15),
    axis.title       = element_text(color = cor_navy),
    axis.text        = element_text(color = cor_cinza),
    legend.title     = element_text(color = cor_navy, face = "bold"),
    legend.text      = element_text(color = cor_cinza),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    plot.background  = element_rect(fill = "white", color = NA)
  )


paleta_2 <- c(cor_verde, cor_dourado)

paleta_n <- function(n) {
  colorRampPalette(c(cor_verde, cor_navy, cor_dourado))(n)
}


# 6. Análise Sem Covariáveis ---------------------------------------------------

surv_obj <- Surv(
  time = df$tempo,
  event = df$cens
)

km_geral <- survfit(
  surv_obj ~ 1,
  data = df
)

km_geral

summary(
  km_geral,
  times = c(6, 12, 24, 36, 60)
)


curva_sem_covariavel <- ggsurvplot(
  km_geral,
  data              = df,
  conf.int          = TRUE,
  conf.int.fill     = cor_verde,
  risk.table        = TRUE,
  surv.median.line  = "hv",
  palette           = cor_verde,
  xlab              = "Tempo (meses)",
  ylab              = "S(t)",
  title             = "Kaplan-Meier — Câncer de Fígado (C22)",
  ggtheme           = tema_km,
  risk.table.col    = "black",
  font.title        = c(15, "bold", cor_navy),
  font.x            = c(13, "plain", cor_navy),
  font.y            = c(13, "plain", cor_navy)
)

curva_sem_covariavel

saveRDS(
  curva_sem_covariavel,
  "Trabalho/parte1/plots/curva_sem_covariavel.RDS"
)


# 7. Análise por Covariáveis Qualitativas --------------------------------------

## SEXO ----

curva_sexo <- ggsurvplot(
  survfit(surv_obj ~ SEXO, data = df),
  data         = df,
  conf.int     = FALSE,
  pval         = TRUE,
  pval.coord   = c(1, 0.1),
  risk.table   = TRUE,
  palette      = paleta_2,
  xlab         = "Tempo (meses)",
  ylab         = "S(t)",
  title        = "Kaplan-Meier por Sexo",
  legend.title = "Sexo",
  legend.labs  = levels(df$SEXO),
  ggtheme      = tema_km,
  font.title   = c(15, "bold", cor_navy)
)

curva_sexo

saveRDS(
  curva_sexo,
  "Trabalho/parte1/plots/curva_sexo.RDS"
)


## ESCOLARI ----

curva_escolari <- ggsurvplot(
  survfit(surv_obj ~ ESCOLARI, data = df),
  data         = df,
  conf.int     = FALSE,
  pval         = TRUE,
  pval.coord   = c(1, 0.1),
  risk.table   = TRUE,
  palette      = paleta_n(nlevels(df$ESCOLARI)),
  xlab         = "Tempo (meses)",
  ylab         = "S(t)",
  title        = "Kaplan-Meier por Escolaridade",
  legend.title = "Escolaridade",
  legend.labs  = levels(df$ESCOLARI),
  ggtheme      = tema_km,
  font.title   = c(15, "bold", cor_navy),
  font.legend  = c(10)
)

curva_escolari

saveRDS(
  curva_escolari,
  "Trabalho/parte1/plots/curva_escolari.RDS"
)


## ECGRUP ----

curva_ecgrup <- ggsurvplot(
  survfit(surv_obj ~ ECGRUP, data = df),
  data         = df,
  conf.int     = FALSE,
  pval         = TRUE,
  pval.coord   = c(1, 0.1),
  risk.table   = TRUE,
  palette      = paleta_n(nlevels(df$ECGRUP)),
  xlab         = "Tempo (meses)",
  ylab         = "S(t)",
  title        = "Kaplan-Meier por Estadiamento",
  legend.title = "Estadiamento",
  legend.labs  = levels(df$ECGRUP),
  ggtheme      = tema_km,
  font.title   = c(15, "bold", cor_navy),
  font.legend  = c(10)
)

curva_ecgrup

saveRDS(
  curva_ecgrup,
  "Trabalho/parte1/plots/curva_ecgrup.RDS"
)


## CATEATEND ----

curva_cateatend <- ggsurvplot(
  survfit(surv_obj ~ CATEATEND, data = df),
  data         = df,
  conf.int     = FALSE,
  pval         = TRUE,
  pval.coord   = c(1, 0.1),
  risk.table   = TRUE,
  palette      = paleta_n(nlevels(df$CATEATEND)),
  xlab         = "Tempo (meses)",
  ylab         = "S(t)",
  title        = "Kaplan-Meier por Categoria de Atendimento",
  legend.title = "Atendimento",
  legend.labs  = levels(df$CATEATEND),
  ggtheme      = tema_km,
  font.title   = c(15, "bold", cor_navy)
)

curva_cateatend

saveRDS(
  curva_cateatend,
  "Trabalho/parte1/plots/curva_cateatend.RDS"
)


## CIRURGIA ----

curva_cirurgia <- ggsurvplot(
  survfit(surv_obj ~ CIRURGIA, data = df),
  data         = df,
  conf.int     = FALSE,
  pval         = TRUE,
  pval.coord   = c(1, 0.1),
  risk.table   = TRUE,
  palette      = paleta_2,
  xlab         = "Tempo (meses)",
  ylab         = "S(t)",
  title        = "Kaplan-Meier por Cirurgia",
  legend.title = "Cirurgia",
  legend.labs  = levels(df$CIRURGIA),
  ggtheme      = tema_km,
  font.title   = c(15, "bold", cor_navy)
)

curva_cirurgia

saveRDS(
  curva_cirurgia,
  "Trabalho/parte1/plots/curva_cirurgia.RDS"
)


## RADIO ----

curva_radio <- ggsurvplot(
  survfit(surv_obj ~ RADIO, data = df),
  data         = df,
  conf.int     = FALSE,
  pval         = TRUE,
  pval.coord   = c(1, 0.1),
  risk.table   = TRUE,
  palette      = paleta_2,
  xlab         = "Tempo (meses)",
  ylab         = "S(t)",
  title        = "Kaplan-Meier por Radioterapia",
  legend.title = "Radioterapia",
  legend.labs  = levels(df$RADIO),
  ggtheme      = tema_km,
  font.title   = c(15, "bold", cor_navy)
)

curva_radio

saveRDS(
  curva_radio,
  "Trabalho/parte1/plots/curva_radio.RDS"
)


## QUIMIO ----

curva_quimio <- ggsurvplot(
  survfit(surv_obj ~ QUIMIO, data = df),
  data         = df,
  conf.int     = FALSE,
  pval         = TRUE,
  pval.coord   = c(1, 0.1),
  risk.table   = TRUE,
  palette      = paleta_2,
  xlab         = "Tempo (meses)",
  ylab         = "S(t)",
  title        = "Kaplan-Meier por Quimioterapia",
  legend.title = "Quimioterapia",
  legend.labs  = levels(df$QUIMIO),
  ggtheme      = tema_km,
  font.title   = c(15, "bold", cor_navy)
)

curva_quimio

saveRDS(
  curva_quimio,
  "Trabalho/parte1/plots/curva_quimio.RDS"
)


# 8. Análise por Covariáveis Quantitativas -------------------------------------

## IDADE ----

df_idade_box <- bind_rows(
  df %>%
    mutate(
      grupo = ifelse(cens == 1, "Óbito por câncer", "Censura")
    ) %>%
    select(grupo, IDADE),

  df %>%
    mutate(
      grupo = "Total"
    ) %>%
    select(grupo, IDADE)
) %>%
  mutate(
    grupo = factor(
      grupo,
      levels = c(
        "Óbito por câncer",
        "Censura",
        "Total"
      )
    )
  )


resumo_idade <- df_idade_box %>%
  group_by(grupo) %>%
  summarise(
    media   = mean(IDADE, na.rm = TRUE),
    mediana = median(IDADE, na.rm = TRUE),
    q1      = quantile(IDADE, 0.25, na.rm = TRUE),
    q3      = quantile(IDADE, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_num = as.numeric(grupo),
    x_min = x_num - 0.275,
    x_max = x_num + 0.275
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

  geom_boxplot(
    color         = cor_navy,
    alpha         = 0.75,
    width         = 0.55,
    outlier.color = cor_cinza,
    outlier.alpha = 0.6
  ) +

  geom_segment(
    data = resumo_idade,
    aes(
      x = x_min,
      xend = x_max,
      y = media,
      yend = media
    ),
    inherit.aes = FALSE,
    color = cor_media,
    linetype = "dashed",
    linewidth = 0.7
  ) +

  geom_text(
    data = resumo_idade,
    aes(
      x = grupo,
      y = media,
      label = paste0("Média: ", round(media, 1))
    ),
    inherit.aes = FALSE,
    vjust = -1.0,
    size = 3.4,
    color = cor_media,
    fontface = "bold"
  ) +

  geom_text(
    data = resumo_idade,
    aes(
      x = grupo,
      y = mediana,
      label = paste0("Q2: ", round(mediana, 1))
    ),
    inherit.aes = FALSE,
    vjust = 1.9,
    size = 3.2,
    color = cor_navy
  ) +

  geom_text(
    data = resumo_idade,
    aes(
      x = grupo,
      y = q1,
      label = paste0("Q1: ", round(q1, 1))
    ),
    inherit.aes = FALSE,
    vjust = 1.6,
    hjust = -0.35,
    size = 2.9,
    color = cor_cinza
  ) +

  geom_text(
    data = resumo_idade,
    aes(
      x = grupo,
      y = q3,
      label = paste0("Q3: ", round(q3, 1))
    ),
    inherit.aes = FALSE,
    vjust = -1.0,
    hjust = -0.35,
    size = 2.9,
    color = cor_cinza
  ) +

  scale_fill_manual(values = paleta_idade) +

  labs(
    title = "Distribuição da Idade ao Diagnóstico",
    x     = NULL,
    y     = "Idade (anos)"
  ) +

  tema_km +

  theme(
    legend.position = "none"
  )

curva_idade_boxplot

saveRDS(
  curva_idade_boxplot,
  "Trabalho/parte1/plots/curva_idade_boxplot.RDS"
)
