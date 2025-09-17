# Cargar tidyverse para manipular los datos
library(dplyr)
library(purrr)
library(ggplot2)
library(MuMIn)

load("~/Montana/montana.small completo.RData")


# Elijo rango de a?os con datos de MNA y artrop
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
datos_pc$time <- as.numeric(datos_pc$year)+(as.numeric(datos_pc$month)-1/2)/12

fuzzy.factor <- function(x,x0,w){
  return(pmin(pmax((x-(x0-w/2))/w,0),1)-pmin(pmax((x-(x0+1-w/2))/w,0),1))
}

xp=seq(-2,2,by=1/120)
plot(xp,y=fuzzy.factor(xp,0,3/4),type='l',col='black')
lines(xp,y=fuzzy.factor(xp,-1,3/4),col='red')
lines(xp,y=fuzzy.factor(xp,1,3/4),col='green')


for (yr in seq(min.year,max.year+1)) {
  df <- data.frame(new.col=fuzzy.factor(datos_pc$time, yr,0.5))
  if(yr == min.year)
    plot(x=datos_pc$time, y=df$new.col, col=(1+yr-min.year), type = 'l')
  else
    lines(x=datos_pc$time, y=df$new.col, col=(1+yr-min.year))
  colnames(df) <- sprintf("year%d",yr)
  datos_pc<-cbind(datos_pc,df)
}


datos_pc$insect[datos_pc$month==2] <- 0
insect.lm <- lm(log(1+(insect))~(I(1-cos(2*pi*time))+sin(2*pi*time))*
                  (year2000+year2001+year2002+year2003+year2004+
                     year2005+year2006+year2007+year2008+year2009)*site,
                data=datos_pc)
datos_pc$interp.insect <- exp(predict(insect.lm, datos_pc))-1

plot(insect~time, type='p', data=datos_pc, col=site,
     xlim=c(min.year,max.year+1),log='')
lines(interp.insect~time, data=datos_pc,col=site, subset = site=="Cascade")
lines(interp.insect~time, data=datos_pc,col=site, subset = site=="Polson")
plot(insect.lm)

datos_pc$arachn[datos_pc$month==2] <- 0
arachn.lm <- lm(log(1+(arachn))~(I(1-cos(2*pi*time))+sin(2*pi*time))*
                  (year2000+year2001+year2002+year2003+year2004+
                     year2005+year2006+year2007+year2008+year2009)*site,
                data=datos_pc)
datos_pc$interp.arachn <- exp(predict(arachn.lm, datos_pc))-1

plot(arachn~time, type='p', data=datos_pc, col=site,
     xlim=c(min.year,max.year+1))
lines(interp.arachn~time, data=datos_pc,col=site, subset = site=="Cascade")
lines(interp.arachn~time, data=datos_pc,col=site, subset = site=="Polson")


datos_pc$other[datos_pc$month==2] <- 0
other.lm <- lm(log(1+(other))~(I(1-cos(2*pi*time))+sin(2*pi*time))*
                  (year2000+year2001+year2002+year2003+year2004+
                     year2005+year2006+year2007+year2008+year2009)*site,
                data=datos_pc)
datos_pc$interp.other <- exp(predict(other.lm, datos_pc))-1

plot(other~time, type='p', data=datos_pc, col=site,
     xlim=c(min.year,max.year+1))
lines(interp.other~time, data=datos_pc,col=site, subset = site=="Cascade")
lines(interp.other~time, data=datos_pc,col=site, subset = site=="Polson")


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
  full.model <- try(glm(MNA_lag ~ (insect+arachn+other)*site + month,family = poisson(link = "log"),
                        data = datos_lag,   na.action = na.fail))
  if(class(full.model)[1] == "try-error"){ 
    return(list(lag=lag_meses,
                n_data = nrow(datos_lag),
                disp = NA,
                full = NULL,
                nest = NULL,
                aver = NULL))}
  dispersion <- sum(residuals(full.model, type = "deviance")^2/fitted.values(full.model))/df.residual(full.model)
  nest.model <- try(dredge(full.model))
  if(class(nest.model)[1] == "try-error") {
    return(list(lag=lag_meses,
                n_data = nrow(datos_lag),
                disp = dispersion,
                full = full.model,
                nest = NULL,
                aver = NULL))}
  aver.model <- try(model.avg(nest.model))
  if(class(aver.model)[1] == "try-error"){ 
    return(list(lag=lag_meses,
                n_data = nrow(datos_lag),
                disp = dispersion,
                full = full.model,
                nest = nest.model,
                aver = NULL))}
  return(list(lag=lag_meses,
              n_data = nrow(datos_lag),
              disp = dispersion,
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
      insect_Polson = imp["insect:sitePolson"],
      arachn_Polson = imp["arachn:sitePolson"],
      other_Polson = imp["other:sitePolson"]
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
install.packages("systemfonts")
install.packages("textshaping")
library(tidyverse)

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



ggplot(coef_artropodos, aes(x = lag, y = Estimate, ymin = conf.low, ymax = conf.high)) +
  geom_line(aes(color = efecto), size = 1) +
  geom_ribbon(aes(fill = efecto), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ efecto, ncol = 1, scales = "free_y") +  # Subplots apilados verticalmente
  labs(
    title = "Efecto estimado de artr�podos sobre MNA por lag",
    subtitle = "Modelos promediados por sitio",
    x = "Lag (meses)",
    y = "Coeficiente estimado AIC 95%"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Oculta leyenda redundante
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")  # Espaciado entre subplots
  )
