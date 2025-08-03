## Funções extras boletim submunicipal

# Função para normalizar nomes de bairros
normalizar_bairro <- function(nome) {
  nome %>%
    iconv(to = "ASCII//TRANSLIT") %>%
    toupper() %>%
    str_replace_all("[^A-Z0-9 ]", "") %>%
    str_squish()
}

# Função para ler dados de dengue (caso precise dos anos anteriores)
ler_dados_dengue <- function(caminho) {
  foreign::read.dbf(caminho) %>%
    filter(ID_MUNICIP == "420910", ID_MN_RESI == "420910") %>%  # Joinville
    select(NU_NOTIFIC, NU_ANO, SEM_NOT, SG_UF_NOT, ID_UNIDADE,
           ID_BAIRRO, NM_BAIRRO, ID_RG_RESI, ID_MN_RESI, CLASSI_FIN,
           DT_SIN_PRI, DT_DIGITA)
}

# Função para completar dados (bairros x semanas)
completar_dados_bairros <- function(dados_bairros, bairros_unicos) {
  dados_bairros %>%
    group_by(ano = substr(SEM_NOT, 1, 4)) %>%
    group_modify(~ {
      semanas_ano <- as.numeric(substr(.x$SEM_NOT, 5, 6))
      total_SE <- ifelse(length(semanas_ano) > 0, max(semanas_ano), 0)
      
      if(total_SE > 0) {
        semanas_completas <- sprintf("%s%02d", .y$ano, 1:total_SE)
        
        expand_grid(
          NM_BAIRRO_REF = bairros_unicos,
          SEM_NOT = semanas_completas
        ) %>%
          left_join(.x, by = c("NM_BAIRRO_REF", "SEM_NOT")) %>%
          mutate(notificações = replace_na(notificações, 0))
      } else {
        .x
      }
    }) %>%
    ungroup() %>%
    arrange(NM_BAIRRO_REF, SEM_NOT)
}

# Função para calcular alerta 
calcular_alerta <- function(id, nome, coluna_id, data, gtdist, meangt, sdgt) {
  bdi <- data %>%
    filter(.data[[coluna_id]] == id) %>%  # Filtra usando a coluna passada como string
    arrange(SE) %>%
    AlertTools::Rt(count = "casos", gtdist = gtdist, meangt = meangt, sdgt = sdgt)
  
  # Aplicar critérios de alerta
  crit.x.vector <- structure(as.character(params), names = names(params))
  criteriaU <- setCriteria(rule = params$codmodelo, values = crit.x.vector)
  y <- AlertTools::fouralert(bdi, crit = criteriaU)
  
  bdi$nivel <- y$indices$level
  bdi$receptivo <- y$indices$cytrue
  
  return(bdi)
}

# Função para classificar bairros por distrito
classificar_bairros_distrito <- function(sp_distritos, sp_bairros) {
  # Garantir mesmo CRS
  if(st_crs(sp_distritos) != st_crs(sp_bairros)) {
    sp_bairros <- st_transform(sp_bairros, st_crs(sp_distritos))
  }
  
  # Usar centroides para melhor precisão
  centroides <- st_centroid(sp_bairros)
  
  # Interseção espacial
  intersecao <- st_intersection(centroides, sp_distritos) %>%
    st_drop_geometry() %>%
    select(nome_bairr, id_bairro, Distrito)
  
  # Completar com distritos mais próximos se necessário
  resultado <- sp_bairros %>%
    select(nome_bairr, id_bairro, geometry) %>%
    left_join(intersecao, by = c("nome_bairr", "id_bairro"))
  
  # Para bairros sem classificação, usar distrito mais próximo
  sem_distrito <- which(is.na(resultado$Distrito))
  if(length(sem_distrito) > 0) {
    indices_proximos <- st_nearest_feature(
      st_centroid(resultado[sem_distrito,]), 
      sp_distritos
    )
    resultado$Distrito[sem_distrito] <- sp_distritos$Distrito[indices_proximos]
  }
  
  return(resultado %>% select(-nome_bairr) %>% rename(geometry_distritos = geometry))
}

# Função otimizada para gerar mapas
generate_maps <- function(data, weeks_2_plot, lim_range = NULL,
                          value_col = "p_inc100k", geom_col = "geometry",
                          shape = NULL, area = NULL) {
  
  # Calcular limites automaticamente se não fornecidos
  if (is.null(lim_range)) {
    filtered_data <- data %>% filter(SE %in% weeks_2_plot)
    lim_range <- range(pull(filtered_data, !!sym(value_col)), na.rm = TRUE)
  }
  
  # Definir os padrões de linha desejados
  line_types <- c("dotted", "dashed", "solid")
  
  # Gerar mapas usando map() para maior eficiência
  plots_list <- map(weeks_2_plot, ~ {
    week_data <- data %>% filter(SE == .x)
    
    # Calcular porcentagem de não alocados (se aplicável)
    percent_na <- if(exists("bairros_parecidos") && exists("dg2025")) {
      100 - (bairros_parecidos %>% filter(SE == .x) %>% nrow() / 
               dg2025 %>% filter(SEM_NOT == .x) %>% nrow() * 100)
    } else 0
    
    # Criar plot base
    p <- ggplot() +
      geom_sf(data = week_data %>% st_as_sf(), aes(fill = !!sym(value_col), geometry = !!sym(geom_col))) +
      scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = lim_range)
    
    # Adicionar camada de distritos com os padrões de linha especificados
    if (!is.null(shape)) {
      if ("geom" %in% names(shape) && !("geometry" %in% names(shape))) {
        shape <- shape %>% rename(geometry = geom)
      }
      
      if (!is.null(area)) {
        # Verificar quantos níveis únicos existem na coluna area
        area_levels <- unique(pull(shape, !!sym(area)))
        n_levels <- length(area_levels)
        
        # Criar um mapeamento de linhas para cada nível
        line_scale <- setNames(
          line_types[1:min(n_levels, length(line_types))], 
          area_levels[1:min(n_levels, length(line_types))]
        )
        
        p <- p + 
          geom_sf(data = shape %>% st_as_sf(), 
                  aes(geometry = geometry, linetype = !!sym(area)), 
                  linewidth = 0.5, alpha = 0, color = "black") +
          scale_linetype_manual(values = line_scale)
      } else {
        # Sem cor específica (apenas contorno)
        p <- p + 
          geom_sf(data = shape %>% st_as_sf(), 
                  aes(geometry = geometry), 
                  linewidth = 1, alpha = 0, color = "black", linetype = "solid")
      }
    }
    
    # Aplicar tema e labels
    p <- p +
      theme_pubclean(base_size = 10) +
      theme(axis.text = element_blank(), 
            axis.ticks = element_blank(),
            legend.key.width = unit(1, "cm")) +  # Aumenta espaço para visualizar linhas
      labs(title = paste("Semana de Notificação:", .x),
           subtitle = if(percent_na > 0) paste("Porcentagem de 'Não alocado':", round(percent_na, 2), "%") else NULL)
    
    return(p)
  })
  
  return(plots_list)
}

create_incidence_plot <- function(nowcast_data, observed_data, title_suffix = NULL, 
                                  distrito_nome = NULL, SE = NULL,
                                  ylim_max = NULL, include_zoom = FALSE, pop = pop,
                                  convert_to_incidence = F) {
  
  # Transformar em incidência se solicitado
  if (convert_to_incidence) {
    observed_data <- observed_data %>%
      mutate(total_Y = total_Y / pop * 1e5)
    
    nowcast_data <- nowcast_data %>%
      mutate(across(c(Median, LI, LS), ~ .x / pop * 1e5))
    
    y_label <- "Incidência (casos por 100.000 hab.)"
  } else {
    y_label <- "Incidência (casos por 100.000 hab.)"
  }
  
  # Determinar o título baseado nos parâmetros fornecidos
  if (!is.null(title_suffix)) {
    plot_title <- title_suffix
  } else if (!is.null(distrito_nome) && !is.null(SE)) {
    plot_title <- paste0("Curva de incidência - ", distrito_nome, " até SE", SE)
  } else if (!is.null(distrito_nome)) {
    plot_title <- paste0("Curva de incidência - ", distrito_nome)
  } else {
    plot_title <- "Curva de incidência"
  }
  
  # Criar o gráfico principal
  p_main <- ggplot() +
    geom_segment(
      data = observed_data,
      aes(x = factor(ano_epi), xend = factor(ano_epi), 
          y = 0, yend = total_Y, color = factor(nivel)),
      linewidth = 0.7
    ) +
    # Linha observada
    geom_line(
      data = observed_data,
      aes(x = factor(ano_epi), y = total_Y, group = 1),
      color = "steelblue", linewidth = 1
    ) +
    # Nowcasting
    geom_line(
      data = nowcast_data %>% dplyr::filter(type == "Nowcasting"),
      aes(x = factor(ano_epi), y = Median, color = "Nowcasting", group = type),
      linewidth = 0.5, linetype = 2
    ) +
    # Ribbon nowcasting
    geom_ribbon(
      data = nowcast_data %>% dplyr::filter(type == "Nowcasting"),
      aes(x = factor(ano_epi), ymin = LI, ymax = LS, fill = "Nowcasting", group = type),
      alpha = 0.2
    ) +
    # Escalas
    scale_color_manual(
      name = NULL,
      values = c("Nowcasting" = "red3", "1" = "green", "2" = "yellow", 
                 "3" = "orange", "4" = "red3"),
      labels = c("Nowcasting" = "Nowcasting", "1" = "Baixo Risco", 
                 "2" = "Atenção", "3" = "Risco Moderado", "4" = "Risco Alto"),
      breaks = c("Nowcasting", "1", "2", "3", "4")
    ) +
    scale_fill_manual(values = c("Nowcasting" = "red3"), guide = "none") +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_x_discrete(
      breaks = levels(factor(observed_data$ano_epi))[
        seq(1, nlevels(factor(observed_data$ano_epi)), by = 8)
      ]
    ) +
    labs(
      x = "Semanas Epidemiológicas",
      y = "Incidência (casos por 100.000 hab.)",
      title = title_suffix
    )
    
    # Adicionar ylim se especificado
    if (!is.null(ylim_max)) {
      p_main <- p_main + ylim(0, ylim_max)
    }
  
  # Se não for incluir o zoom, retorne apenas o gráfico principal
  if (!include_zoom) {
    return(p_main)
  }
  
  # Preparar dados para o zoom (últimas 10 semanas)
  ultimas_semanas <- tail(sort(unique(observed_data$ano_epi)), 5)
  zoom_data <- observed_data %>% 
    filter(ano_epi %in% ultimas_semanas)
  
  # Criar gráfico de zoom
  p_zoom <- ggplot(zoom_data) +
    geom_segment(
      aes(x = factor(ano_epi), xend = factor(ano_epi),
          y = 0, yend = total_Y, color = factor(nivel)),
      linewidth = 0.8
    ) +
    geom_line(
      aes(x = factor(ano_epi), y = total_Y, group = 1),
      color = "steelblue", linewidth = 1.2
    ) +
    # Nowcasting
    geom_line(
      data = nowcast_data %>% dplyr::filter(type == "Nowcasting"),
      aes(x = factor(ano_epi), y = Median, color = "Nowcasting", group = type),
      linewidth = 0.5, linetype = 2
    ) +
    # Ribbon nowcasting
    geom_ribbon(
      data = nowcast_data %>% dplyr::filter(type == "Nowcasting"),
      aes(x = factor(ano_epi), ymin = LI, ymax = LS, fill = "Nowcasting", group = type),
      alpha = 0.2
    )+
    # Escalas
    scale_color_manual(
      name = NULL,
      values = c("Nowcasting" = "red3", "1" = "green", "2" = "yellow", 
                 "3" = "orange", "4" = "red"),
      labels = c("Nowcasting" = "Nowcasting", "1" = "Baixo Risco", 
                 "2" = "Atenção", "3" = "Risco Moderado", "4" = "Risco Alto"),
      breaks = c("Nowcasting", "1", "2", "3", "4")
    ) +
    scale_fill_manual(values = c("Nowcasting" = "red3"), guide = "none") +
    theme_minimal() +
    labs(x = NULL, y = NULL, title = "Últimas 5 semanas")+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      title = element_text(size = 6),
      plot.background = element_rect(fill = "white", color = "gray"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_x_discrete(
      breaks = levels(factor(zoom_data$ano_epi))[
        seq(1, nlevels(factor(zoom_data$ano_epi)), by = 3)
      ]
    ) 
  
  # Combinar os gráficos usando patchwork
  final_plot <- p_main + 
    inset_element(p_zoom, left = 0.6, bottom = 0.3, right = 1, top = 1)
  
  return(final_plot)
}

# Função para criar gráfico de receptividade
create_receptivity_plot <- function(api_data, weeks_limit) {
  # Calcular fator de escala
  scale_factor <- max(api_data$receptivo, na.rm = TRUE) / max(api_data$tempmin, na.rm = TRUE)
  
  api_data %>% 
    filter(SE <= as.numeric(as.character(weeks_limit))) %>% 
    ggplot() +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1), fill = "grey70") +
    geom_col(aes(x = factor(SE), y = receptivo, group = 1), fill = "steelblue") +
    geom_line(aes(x = factor(SE), y = tempmin * scale_factor, group = 1), 
              color = "orange2", linewidth = 1) +
    geom_hline(aes(yintercept = 18 * scale_factor, color = "Limiar Favorável (18°C)"), 
               linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c("red"))+
    scale_x_discrete(breaks = levels(factor(api_data$SE))[seq(1, nlevels(factor(api_data$SE)), by = 8)]) +
    scale_y_continuous(
      name = "Receptividade climática",
      sec.axis = sec_axis(~./scale_factor, name = "Temp. Mínima")
    ) +
    labs(x = "Semanas Epidemiológica", 
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

# Função para criar gráfico de Rt (Criar condição para usar facet com distritos)
create_rt_plot <- function(api_data, weeks_limit, title_suffix, facet_by = NULL,
                           lwr = NULL, upr = NULL) {
 p <-  api_data %>% 
    filter(SE <= as.numeric(as.character(weeks_limit))) %>% 
    ggplot() +
    geom_ribbon(
    aes(x = factor(SE), ymin = lwr, ymax = upr, group = facet_by),
     alpha = 0.7, fill = "grey"
   )+
    geom_line(aes(x = factor(SE), y = Rt, group = 1), linetype = 1) +
    geom_hline(yintercept = 1, linetype = 2, col = "red") +
    scale_x_discrete(breaks = levels(factor(api_data$SE))[seq(1, nlevels(factor(api_data$SE)), by = 8)]) +
    labs(x = "Semanas Epidemiológicas", y = "Rt", title = title_suffix) +
    theme_minimal() +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Adicionar facet_wrap condicionalmente
  if (!is.null(facet_by)) {
    # Verificar se a variável existe nos dados
    if (facet_by %in% names(api_data)) {
      p <- p + facet_wrap(as.formula(paste("~", facet_by)), ncol = 1)
    } else {
      warning(paste("Variável", facet_by, "não encontrada nos dados. Facet ignorado."))
    }
  }
  
  return(p)
}
  

# Função para criar painel de mapas
create_map_panel <- function(plots_list, legend_title = "Incidência (100.000 pessoas)") {
  # Extrair legenda
  legend_plot <- plots_list[[length(plots_list) - 2]] + 
    theme(legend.position = "top") + 
    labs(fill = legend_title)
  
  grob_legend <- ggplotGrob(legend_plot)
  legend_grob <- gtable_filter(grob_legend, "guide-box")
  
  # Criar painel
  ggdraw() +
    draw_plot(plots_list[[length(plots_list) - 2]] + theme(legend.position = "none"), 
              x = 0.08, y = 0.1, width = 0.27) +
    draw_plot(plots_list[[length(plots_list) - 1]] + theme(legend.position = "none"), 
              x = 0.65, y = 0.1, width = 0.27) +
    draw_plot(plots_list[[length(plots_list)]] + theme(legend.position = "none"), 
              x = 0.4, y = -0.3, width = 0.27) +
    draw_grob(legend_grob, x = 0.3, y = 0.05,width = 0.4)
}


# Função para nowcasting por distrito
processar_nowcast_distrito <- function(distrito_nome, K = NULL, Dmax = NULL) {
  # Preparar dados
  dados_distrito <- distrito_nowcast %>% 
    left_join(dg_bairros %>%
                select(NM_BAIRRO_REF, Distrito) %>% 
                distinct, 
              by = "NM_BAIRRO_REF") %>% 
    select(NU_NOTIFIC,SEM_NOT,DT_DIGITA,DT_SIN_PRI,Distrito) %>% 
    na.omit() %>% 
    filter(Distrito == distrito_nome,
           year(DT_SIN_PRI) %in% c(ano_atual-2, ano_atual-1, ano_atual)) %>%
    rename_with(tolower) %>%
    select(dt_sin_pri, dt_digita, distrito)
  
  # Gerar nowcasting
  df_nowcast <- nowcasting_inla(
    dataset = dados_distrito,
    data.by.week = TRUE,
    date_onset = dt_sin_pri,
    date_report = dt_digita,
    K = K, Dmax = Dmax
  )
  
  # Processar dados nowcast
  df_nowcast$total <- df_nowcast$total %>%
    mutate(
      ano = year(dt_event),
      epiweek = lubridate::epiweek(dt_event),
      ano_epi = case_when(
        epiweek == 1 & lag(epiweek, default = 52) == 52 ~ 
          paste0(ano_atual, str_pad(epiweek, width = 2, side = "left", pad = "0")),
        TRUE ~ paste0(ano, str_pad(epiweek, width = 2, side = "left", pad = "0"))
      ),
      type = case_when(
        ano_epi > as.numeric(as.character(weeks_2_plot[3])) ~ "Forecast",
        TRUE ~ "Nowcasting"
      )
    )
  
  # Dados observados por semana
  dados_by_week <- dados_distrito %>%
    mutate(
      dt_event = dt_sin_pri,
      epiweek = epiweek(dt_event),
      ano = year(dt_event),
      ano_epi = paste0(ano, str_pad(epiweek, width = 2, pad = "0"))
    ) %>%
    count(ano_epi, dt_event, name = "total_Y") %>%
    filter(total_Y != 0 & ano_epi >= paste0(ano_atual-1, "01")) %>%
    group_by(ano_epi) %>%
    summarise(total_Y = sum(total_Y), .groups = "drop") %>%
    left_join(
      dg_distritos %>% 
        filter(Distrito == distrito_nome) %>%
        mutate(SE = as.character(SE)) %>%
        select(SE, nivel),
      by = c("ano_epi" = "SE")
    )
  
  return(list(nowcast = df_nowcast, observado = dados_by_week))
}



