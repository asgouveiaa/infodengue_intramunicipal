# Analise Joinville

#devtools::install_github("AlertaDengue/AlertTools", dependencies = TRUE)
library(AlertTools)
library(assertthat)
#devtools::install_github("https://github.com/chrismerkord/epical", dependencies = TRUE)
library(epical)

geocode = 4209102

# DADOS ----
# dados casos
load("~/Documents/relatorio/Joinville/JVbair2020_2025.Rdata")

# usando o clima da cidade para todos os bairros
# data clima
#cli <- getClima(geocode) % dados do servidor
#save(cli, file = "clijv.RData")
load("~/Documents/relatorio/Joinville/clijv.RData")
head(cli)
cli %>% filter(substr(SE,1,4) %in% c(ano_atual-1,ano_atual)) %>% ungroup() %>% 
  select(-geocodigo)

# parametros (limiares mem e temp para Joinville)
# foram calculados a nivel de municipio, o ideal
# seria recalcular a nivel mais local
#params <- read.parameters(geocode)
#save(params, file = "paramsjv.RData")
dg_distritos %>% select(!c(geometry,inc,ano)) %>% 
  left_join(cli %>% ungroup() %>% mutate(SE = as.character(SE)) %>% 
              select(-geocodigo),
            by = c("SEM_NOT" = "SE")) -> dg_distritos_nivel

# arrumando os dados (juntando clima e casos) 

## norte
# dg_ditrito_nivel_nort <- dg_distritos_nivel %>%
#   filter(Distrito == "Distrito Norte") %>% 
#   select(SEM_NOT,temp_min,not_total,pop_total) %>% 
#   mutate(SE = as.numeric(SEM_NOT),
#          ano = as.numeric(ano))  %>%# o correto seria usar a SEM_INI
#   rename(casos = not_total) %>%
#   mutate(inc = casos/pop_total*100000) %>% 
#   as.data.frame()


bd <- dg_distritos_nivel %>%
  mutate(id_distrito = case_when(
    Distrito == "Distrito Norte" ~1,
    Distrito == "Distrito Sul" ~2,
    Distrito == "Distrito Centro"~3
  ))%>% 
  select(SEM_NOT,Distrito,id_distrito,temp_min,not_total,pop_total) %>% 
  mutate(SE = as.numeric(SEM_NOT),
         ano = as.numeric(ano))  %>%# o correto seria usar a SEM_INI
  rename(casos = not_total) %>%
  mutate(inc = casos/pop_total*100000) %>% 
  as.data.frame()

str(bd)

# rodando alerta por bairro
ids <- unique(bd$id_distrito)

# dataframe para os resultados
resultado <- data.frame()
  
# calculando alerta por bairro
for(i in 1:length(ids)){
   
   nome <- unique(bd$Distrito)[i]
   bdi <- bd %>%
    filter(id_distrito == ids[i])  %>%  
    arrange(SE)  %>%   # garantir que está em ordem cronologica
    Rt(count = "casos", 
        gtdist = "normal", meangt = 3, sdgt = 1) 
  
   # build rules
   crit.x <- params # parameters
   crit.x.vector <- structure(as.character(crit.x), names = as.character(names(crit.x))) # dataframe -> vector
   criteriaU <- setCriteria(rule = crit.x$codmodelo, values = crit.x.vector) # valued criteria
   
   # Apply alert rules
   y <- fouralert(bdi, crit = criteriaU)  # apply alert 
  
  bdi$nivel <- y$indices$level
  bdi$receptivo <- y$indices$cytrue

  resultado <- rbind(resultado,bdi)
  rm(bdi)
  
}

