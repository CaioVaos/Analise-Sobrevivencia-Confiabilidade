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