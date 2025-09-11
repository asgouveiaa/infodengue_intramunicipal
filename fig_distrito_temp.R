  scale_factor <- max(dados_api$receptivo, na.rm = TRUE) / max(dados_api$tempmin, na.rm = TRUE)
  # scale_factor2 <- max(dados_api$receptivo, na.rm = TRUE) / max(distrito_downsc$mean_pred, na.rm = TRUE)
  
  
  dados_api %>% 
    mutate(ano = as.numeric(substr(SE,1,4))) %>% 
    filter(arbo == "dengue" & ano %in% c(ano_atual-1,ano_atual))
  
  dados_api %>% 
    mutate(ano = as.numeric(substr(SE,1,4))) %>% 
    filter(arbo == "dengue" & ano %in% c(ano_atual-1,ano_atual) &
             SE <= as.numeric(as.character(weeks_2_plot[3]))) %>% 
    ggplot() +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1), fill = "grey70") +
    geom_col(aes(x = factor(SE), y = receptivo, group = 1), fill = "steelblue") +
    geom_line(aes(x = factor(SE), y = tempmin * 0.037, group = 1), 
              color = "orange2", linewidth = 1) +
    geom_hline(aes(yintercept = 18 * scale_factor, color = "Limiar Favorável (18°C)"), 
               linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c("red"))+
    scale_x_discrete(breaks = levels(factor(dados_api$SE))[seq(1, nlevels(factor(dados_api$SE)), by = 8)]) +
    scale_y_continuous(
      name = "Receptividade climática",
      sec.axis = sec_axis(~./scale_factor, name = "Temp. Mínima")
    ) +
    labs(x = "Semanas Epidemiológicas", 
         title = "Receptividade climática e Temperatura Mínima por Semana", 
         color = "") +
    theme_minimal() +
    theme(
      axis.title.y.right = element_text(color = "black"),
      axis.text.y.right = element_text(color = "black"),
      axis.title.y.left = element_text(color = "black"),
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}
  
  
ggplot() +
  geom_line(
    data = dados_api %>% 
      mutate(ano = as.numeric(substr(SE,1,4))) %>% 
      filter(arbo == "dengue" & SE %in% ultimas_semanas),
    aes(x = factor(SE), y = tempmin, group = 1, col = "Município")
  ) +
  geom_line(
    data = distrito_downsc %>% 
      mutate(ano = as.numeric(substr(SE,1,4))) %>% 
      filter(SE %in% ultimas_semanas),
  aes(factor(SE), mean_pred, group = Distrito, color = Distrito)
) +
  geom_ribbon(
    data = distrito_downsc %>% 
      mutate(ano = as.numeric(substr(SE,1,4))) %>% 
      filter(SE %in% ultimas_semanas),
    aes(factor(SE),
        ymin = mean_ll,
        ymax = pred_ul, group = Distrito, fill = Distrito),
    alpha = 0.2, show.legend = F
  ) + 
  geom_hline(aes(yintercept = 18, color = "Limiar Favorável (18°C)"), 
             linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("blue","salmon","green3","red","orange2")) +
  scale_fill_manual(values = c("blue","salmon","green3","red","orange2")) +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "black"),
    axis.text.y.right = element_text(color = "black"),
    axis.title.y.left = element_text(color = "black"),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + 
  labs(x = "Semanas Epidemiológicas", 
       y = "Temp. Mínima",
       title = "Receptividade climática e Temperatura Mínima por Semana nos Distritos", 
       color = "") +
  scale_x_discrete(expand = c(0.005,0.005))
  
  
