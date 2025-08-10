library(dplyr)
library(purrr)
library(ggplot2)
library(MuMIn)
library(tidyr)
library(tidyverse)

load("~/GitHub/Artropodos-Roedores-MNT/datos_Pol+Cas+artrop_interp.RData") #seteo juli
montana.small<-datos_pc
ls()
str(montana.small)
#La abundancia neg no tiene sentido biologico asique las reemplazo por 0 
datos <- montana.small %>%
  mutate(across(c(interp.insect, interp.arachn, interp.other),
                ~ ifelse(. < 0, 0, .)))

#Analisis exploratorio
#Grafico en conjunto el interpolado y los datos reales de insectos

datos <- datos%>%
  mutate(fecha = as.Date(paste(year, month, "01", sep = "-")))

# Filtrar un solo sitio (cambiá el nombre si querés otro)
sitio_seleccionado <- "Polson"

datos_largo <- datos %>%
  filter(site == sitio_seleccionado) %>%
  select(fecha, insect, interp.insect) %>%
  pivot_longer(cols = c(insect, interp.insect),
               names_to = "tipo", values_to = "abundancia") %>%
  mutate(tipo = recode(tipo,
                       "insect" = "Observado",
                       "interp.insect" = "Interpolado"))

# Graficar
ggplot(datos_largo, aes(x = fecha, y = abundancia, color = tipo, group = tipo)) +
  geom_line() +
  geom_point(size = 2) +
  labs(title = paste("Abundancia de insectos en", sitio_seleccionado),
       x = "Fecha", y = "Abundancia",
       color = "Tipo de dato") +
  theme_minimal()




# Calcular abundancia de invertebrados
datos_rel <- datos %>%
  mutate(invertebrados = interp.insect + interp.arachn + interp.other) %>%
  select(site, time, MNA, invertebrados) %>%
  pivot_longer(cols = c(MNA, invertebrados),
               names_to = "grupo",
               values_to = "abundancia") %>%
  filter(!is.na(abundancia)) %>%
  group_by(site, grupo) %>%
  mutate(
    abundancia_relativa = abundancia / max(abundancia, na.rm = TRUE)
  ) %>%
  ungroup()

# Gráfico
ggplot(datos_rel, aes(x = time, y = abundancia_relativa, color = grupo)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ site, ncol = 1) +
  labs(
    title = "Abundancia relativa por grupo dentro de cada sitio",
    x = "Año",
    y = "Abundancia relativa (escalada por grupo en cada sitio)",
    color = "Grupo"
  ) +
  theme_minimal()



####### ahora separo para que los bichos esten por separado
datos_rel_sep <- montana.small %>%
  select(site, time, MNA, interp.insect, interp.arachn, interp.other) %>%
  pivot_longer(cols = c(MNA, interp.insect, interp.arachn, interp.other),
               names_to = "grupo",
               values_to = "abundancia") %>%
  filter(!is.na(abundancia)) %>%
  group_by(site, grupo) %>%
  mutate(
    abundancia_relativa = abundancia / max(abundancia, na.rm = TRUE)
  ) %>%
  ungroup()

# Etiquetas más legibles (opcional)
datos_rel_sep <- datos_rel_sep %>%
  mutate(grupo = recode(grupo,
                        MNA = "Roedores",
                        interp.insect = "Insectos",
                        interp.arachn = "Arácnidos",
                        interp.other = "Otros"))

# Graficar
ggplot(datos_rel_sep, aes(x = time, y = abundancia_relativa, color = grupo)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ site, ncol = 1) +
  labs(
    title = "Abundancia relativa por grupo (roedores, insectos, arácnidos y otros) en cada sitio",
    x = "Año",
    y = "Abundancia relativa (escalada por grupo y sitio)",
    color = "Grupo"
  ) +
  theme_minimal()


#Tiene pinta que hay una asociacion clara con insectos y un lag de 3 meses.
# Etiquetas más legibles (opcional)
datos_insectos <- datos_rel_sep %>%
  filter((grupo %in% c("Insectos", "Roedores")))

# Graficar
ggplot(datos_insectos, aes(x = time, y = abundancia_relativa, color = grupo)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ site, ncol = 1) +
  labs(
    title = "Abundancia relativa por grupo (roedores, insectos, arácnidos y otros) en cada sitio",
    x = "Año",
    y = "Abundancia relativa (escalada por grupo y sitio)",
    color = "Grupo"
  ) +
  theme_minimal()

#Objetivo:
#  Encontrar el lag (en meses) que maximiza el ajuste de un modelo donde la abundancia de insectos predice la abundancia de roedores.
library(tidyverse)
library(broom)
library(lme4)
library(dplyr)
library(scales)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(emmeans)

# Crear base limpia
datos_lag <- datos %>%
  mutate(invertebrados = interp.insect + interp.arachn + interp.other) %>%
  select(site, time, month, MNA, invertebrados) %>%
  filter(!is.na(MNA), !is.na(invertebrados)) %>%
  arrange(site, time)

# Función para ajustar modelo con un lag dado
ajustar_glm_lag <- function(lag_meses) {
  datos_lagged <- datos_lag %>%
    group_by(site) %>%
    mutate(insect_lag = lag(invertebrados, n = lag_meses)) %>%
    ungroup() %>%
    filter(!is.na(insect_lag))
  
  modelo <- glmmTMB(MNA ~ insect_lag*month + (1|site), family = nbinom2, data = datos_lagged)
  coef_insect <- fixef(modelo)$cond["insect_lag"]
  
  tibble(
    lag = lag_meses,
    coef = coef_insect,
    aic = AIC(modelo),
    deviance = deviance(modelo),
    n = nrow(datos_lagged)
  )
}

# Probar lags de 0 a 12 meses
resultados_glm <- map_dfr(0:12, ajustar_glm_lag)


# Ordenar por el valor del coeficiente de insect_lag (de mayor a menor)
resultados_glm %>% arrange(desc(coef))


datos_plot <- resultados_glm %>%
  mutate( 
    aic_norm  =  aic / max(aic),  # invertir: AIC más bajo = valor más alto
    n_norm    = n / max(n)              # más observaciones = mejor
  ) %>%
  select(lag, coef, aic_norm, n_norm) %>%
  pivot_longer(cols = c(coef, aic_norm, n_norm),
               names_to = "metrica",
               values_to = "valor")

ggplot(datos_plot, aes(x = lag, y = valor, color = metrica)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = c("coef" = "#1b9e77",
               "aic_norm" = "#d95f02",
               "n_norm" = "#7570b3")) +
  labs(y = "Valor escalado", x = "Lag (meses)", color = "Métrica") +
  theme_minimal()



### Me quedo con el modelo 7
library(lme4)

# Filtrar y crear la variable lag para el modelo
datos_lag7 <- datos_lag %>%
  group_by(site) %>%
  mutate(insect_lag = lag(invertebrados, n = 7)) %>%
  ungroup() %>%
  filter(!is.na(insect_lag)) %>%
  mutate(insect_lag_z = scale(insect_lag))

modelo_lag7 <- glmmTMB (MNA ~ insect_lag_z * month + (1 | site),
                     family = nbinom2,
                       data = datos_lag7)
#####CHEQUEO SUPUESTOS #########
# Función para chequear sobre-dispersión
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  p <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = ratio, rdf = rdf, p = p)
}
overdisp_fun(modelo_lag7) #no hay sobredispersion

# Simulación de residuos 
residuos <- simulateResiduals(modelo_lag7)
plot(residuos) #Todo en orden
testDispersion(residuos) #Todo en orden
testZeroInflation(residuos) #Todo en orden
#########################

#Sigo con el analisis
summary(modelo_lag7)


datos_lag7 <- datos_lag %>%
  group_by(site) %>%
  mutate(insect_lag = lag(invertebrados, n = 7)) %>%
  ungroup() %>%
  filter(!is.na(insect_lag)) 
datos_lag7$MNA_predicha <- predict(modelo_lag7, type = "response")
library(scales)

datos_plot <- datos_lag7 %>%
  select(site, time, month, insect_lag, MNA_predicha) %>%
  mutate(
    insect_lag_scaled = rescale(insect_lag),  # 0 a 1
    MNA_predicha_scaled = rescale(MNA_predicha) # 0 a 1
  ) %>%
  pivot_longer(cols = c(insect_lag_scaled, MNA_predicha_scaled),
               names_to = "variable", values_to = "valor")
ggplot(datos_plot, aes(x = time, y = valor, color = variable)) +
  geom_line(size = 1.2) +
  facet_wrap(~site) +
  scale_color_manual(
    values = c("insect_lag_scaled" = "#1b9e77", "MNA_predicha_scaled" = "#d95f02"),
    labels = c("Invertebrados (lag 7 meses)", "MNA predicha")
  ) +
  labs(
    x = "Tiempo",
    y = "Serie escalada (0-1)",
    color = "Variable"
  ) +
  theme_minimal()


# MNA predicha vs observada
# Crear columna fecha correctamente
datos_lag7 <- datos_lag7 %>%
  mutate(fecha = as.Date(paste(.data$year, .data$month, "01", sep = "-")))

# Asegurarse insect_lag_z esté bien escalada (igual que en modelo)
datos_lag7 <- datos_lag7 %>%
  mutate(insect_lag_z = scale(insect_lag))

# Obtener predicción en escala link con error estándar
pred_link <- predict(modelo_lag7, newdata = datos_lag7, type = "link", se.fit = TRUE)

# Añadir predicción y IC al dataframe
datos_lag7 <- datos_lag7 %>%
  mutate(
    fit_link = pred_link$fit,
    se_link = pred_link$se.fit,
    lower_link = fit_link - 1.96 * se_link,
    upper_link = fit_link + 1.96 * se_link,
    pred = exp(fit_link),
    lower = exp(lower_link),
    upper = exp(upper_link)
  )

# Gráfico con fecha real en eje X
ggplot(datos_lag7, aes(x = time)) +
  geom_point(aes(y = MNA, color = "Observado"), size = 2, alpha = 0.7) +
  geom_line(aes(y = pred, color = "Predicho"), size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
  facet_wrap(~ site) +
  scale_color_manual(
    name = "Abundancia",
    values = c("Observado" = "black", "Predicho" = "red")
  ) +
  labs(
    x = "Fecha",
    y = "Abundancia de roedores (MNA)",
    title = "Abundancia observada y predicha con intervalo de confianza"
  ) +
  theme_minimal()

#############################################################
#Analizo insectos aracnidos y otros por separado, manteniendo el lag de 7 meses##
datos_lag7_sep <- datos %>%
  select(site, year, month, time, MNA, interp.insect, interp.arachn, interp.other) %>%
  arrange(site, time) %>%
  group_by(site) %>%
  mutate(
    lag_insect = lag(interp.insect, 7),
    lag_arachn = lag(interp.arachn, 7),
    lag_other  = lag(interp.other, 7)
  ) %>%
  ungroup() %>%
  filter(!is.na(lag_insect), !is.na(lag_arachn), !is.na(lag_other)) %>%
  mutate(
    lag_insect_z = scale(lag_insect),
    lag_arachn_z = scale(lag_arachn),
    lag_other_z  = scale(lag_other)
  )

modelo_sep <- glmmTMB(
  MNA ~ lag_insect_z + lag_arachn_z + lag_other_z + month + (1 | site),
  family = nbinom2,
  data = datos_lag7_sep
)
#Chequeo supuestos
overdisp_fun(modelo_sep) #no hay sobredispersion

# Simulación de residuos 
residuos <- simulateResiduals(modelo_sep)
plot(residuos) #Todo en orden
testDispersion(residuos) #Todo en orden
testZeroInflation(residuos) #Todo en orden
#Resultados modelo
summary(modelo_sep)

# Predecir en escala link
pred_sep <- predict(modelo_sep, newdata = datos_lag7_sep, type = "link", se.fit = TRUE)

# Agregar predicciones al dataset
datos_lag7_sep <- datos_lag7_sep %>%
  mutate(
    fit_link = pred_sep$fit,
    se_link = pred_sep$se.fit,
    lower_link = fit_link - 1.96 * se_link,
    upper_link = fit_link + 1.96 * se_link,
    pred = exp(fit_link),
    lower = exp(lower_link),
    upper = exp(upper_link)
  )

datos_lag7_sep <- datos_lag7_sep %>%
  mutate(fecha = as.Date(paste(year, month, "01", sep = "-")))

ggplot(datos_lag7_sep, aes(x = fecha)) +
  geom_point(aes(y = MNA, color = "Observado"), size = 2, alpha = 0.7) +
  geom_line(aes(y = pred, color = "Predicho"), size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
  facet_wrap(~ site, scales = "free_y") +
  scale_color_manual(
    name = "Abundancia",
    values = c("Observado" = "black", "Predicho" = "red")
  ) +
  labs(
    x = "Fecha",
    y = "Abundancia de roedores (MNA)",
    title = "Abundancia observada y predicha (con IC 95%) según insectos, arácnidos y otros (lag 7)"
  ) +
  theme_minimal()

