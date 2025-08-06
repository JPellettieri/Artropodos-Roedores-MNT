# Cargar tidyverse para manipular los datos
library(dplyr)
library(purrr)
library(ggplot2)
library(MuMIn)
library(tidyr)

load("~/GitHub/Artropodos-Roedores-MNT/montana.small completo.RData") #seteo juli


# Elijo rango de años con datos de MNA y artrop
min.year <- min(montana.small$year[!is.na(montana.small$insect)])
max.year <- max(montana.small$year[!is.na(montana.small$MNA)])

# Defino todas las combinaciones de sitio/a?o/mes y los cruzo
# con las columnas necesarias de montana.small
datos_pc <- expand.grid(month=c(1:12), 
                        year=c(min.year:max.year), 
                        site=c('Cascade','Polson')) %>% 
  select(site,year,month) %>%
  left_join((montana.small %>% select(site,year,month,MNA,insect,arachn,other)),
            by=c('site','year','month'), unmatched='drop')
datos_pc$month <- as.factor(datos_pc$month)

datos_lag <- function(x, lag_meses){
  x %>%
  group_by(site) %>%
    mutate(
      MNA_lag = lead(MNA, lag_meses)
    ) %>%
    ungroup() %>%
    filter(!is.na(insect) & !is.na(arachn) & !is.na(other) & !is.na(MNA_lag))
}

View(datos_pc %>% datos_lag(12))

ajustar_modelo_lag <- function(x, lag_meses) {
  datos_lag <- x %>%
    group_by(site) %>%
    mutate(
      MNA_lag = lead(MNA, lag_meses)
    ) %>%
    ungroup() %>%
    filter(!is.na(insect) & !is.na(arachn) & !is.na(other) & !is.na(MNA_lag))
  full.model <- try(glm.nb(MNA_lag ~ (insect + arachn + other) * site + month,
                           data = datos_lag, na.action = na.fail))
  if(class(full.model)[1] == "try-error"){ 
    return(list(lag=lag_meses,
                n_data = nrow(datos_lag),
                full = NULL,
                nest = NULL,
                aver = NULL))}
  nest.model <- try(dredge(full.model))
  if(class(nest.model)[1] == "try-error") {
    return(list(lag=lag_meses,
                n_data = nrow(datos_lag),
                full = full.model,
                nest = NULL,
                aver = NULL))}
  aver.model <- try(model.avg(nest.model))
  if(class(aver.model)[1] == "try-error"){ 
    return(list(lag=lag_meses,
                n_data = nrow(datos_lag),
                full = full.model,
                nest = nest.model,
                aver = NULL))}
  return(list(lag=lag_meses,
              n_data = nrow(datos_lag),
              full = full.model,
              nest = nest.model,
              aver = aver.model))
}

summary(ajustar_modelo_lag(datos_pc, 1))

lag_meses <- c(0:12)
lista_modelos <- map(lag_meses, \(x) datos_pc %>% ajustar_modelo_lag(x))

######## agrego esto no se si funciona#########
###############################################
importancia_artropodos <- map_dfr(lista_modelos, function(m) {
  if (!is.null(m$nest)) {
    imp <- MuMIn::sw(m$nest)
    tibble(
      lag = m$lag,
      insect = imp["insect"],
      arachn = imp["arachn"],
      other = imp["other"],
      insect_Cascade = imp["insect:siteCascade"],
      arachn_Cascade = imp["arachn:siteCascade"],
      other_Cascade = imp["other:siteCascade"]
    )
  }
}, .id = "modelo_id")

###################################################
#######################################################

modelo_lag <- lista_modelos[[1]]
coeffs_lag <- coefTable(modelo_lag$aver)
rownames(coeffs_lag)
colnames(coeffs_lag)

coeffs_lag[,1]
tabla_modelos <- data.frame(lag = modelo_lag$lag,
                            n_data = modelo_lag$n_data,
                            insect.est = coeffs_lag["insect","Estimate"],
                            arachn.est = coeffs_lag["arachn","Estimate"],
                            other.est = coeffs_lag["other","Estimate"])


plot_results <- function(mod_list, param){
  n_mod <- length(mod_list)
  for (k in c(1:n_mod)) {
    x <- mod_list[[k]]$lag
    if(!is.null(mod_list[[k]]$aver)){
      coef <- coefTable(mod_list[[k]]$aver)
      
    }
  }
}


############aca sigo con lo que falta #####
library(broom)
coef_artropodos <- map_dfr(lista_modelos, function(m) {
  if (!is.null(m$aver)) {
    coefs <- as.data.frame(coefTable(m$aver))
    coefs$term <- rownames(coefs)
    coefs$lag <- m$lag
    coefs$n_data <- m$n_data
    coefs$conf.low <- coefs$Estimate - 1.96 * coefs$`Std. Error`
    coefs$conf.high <- coefs$Estimate + 1.96 * coefs$`Std. Error`
    coefs
  }
}, .id = "modelo_id") %>%
  filter(grepl("insect|arachn|other", term)) %>%
  mutate(
    efecto = case_when(
      term == "insect" ~ "Insectos (Polson)",
      term == "arachn" ~ "Arácnidos (Polson)",
      term == "other" ~ "Otros (Polson)",
      term == "insect:siteCascade" ~ "Insectos (Cascade)",
      term == "arachn:siteCascade" ~ "Arácnidos (Cascade)",
      term == "other:siteCascade" ~ "Otros (Cascade)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(efecto))

coef_por_lag <- map_dfr(lista_modelos, function(m) {
  if (!is.null(m$aver)) {
    coefs <- as.data.frame(coefTable(m$aver))
    coefs$term <- rownames(coefs)
    coefs$lag <- m$lag
    coefs$n_data <- m$n_data
    coefs
  }
}, .id = "modelo_id")


ggplot(coef_artropodos, aes(x = lag, y = Estimate, ymin = conf.low, ymax = conf.high, color = efecto)) +
  geom_line(size = 1) +
  geom_ribbon(aes(fill = efecto), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Efecto estimado de artropodos sobre MNA por lag",
    subtitle = "Modelos promediados por sitio",
    x = "Lag (meses)",
    y = "Coeficiente estimado AIC 95%",
    color = "Efecto",
    fill = "Efecto"
  ) +
  theme_minimal()

#######################grafico de importancia relativa
importancia_artropodos_long <- importancia_artropodos %>%
  pivot_longer(cols = -c(lag, modelo_id), names_to = "efecto", values_to = "importancia")

ggplot(importancia_artropodos_long, aes(x = lag, y = importancia, color = efecto)) +
  geom_line(size = 1.1) +
  geom_point(size = 2) +
  labs(
    title = "Importancia relativa de predictores por lag",
    x = "Lag (meses)",
    y = "Importancia relativa (suma pesos AIC)",
    color = "Efecto"
  ) +
  theme_minimal()


############################################


which(!map_lgl(lista_modelos, ~ is.null(.x$aver)))
summary(lista_modelos[[2]]$aver)  # Porque lag = 1 está en la posición 2 (lag 0 es el primero)

#ver todos los lags
for (i in seq_along(lista_modelos)) {
  cat("\n=========================\n")
  cat("Lag:", lista_modelos[[i]]$lag, "\n")
  if (!is.null(lista_modelos[[i]]$aver)) {
    print(summary(lista_modelos[[i]]$aver))
  } else {
    cat("Modelo promedio no disponible\n")
  }
}







