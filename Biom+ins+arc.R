library(dplyr)
library(purrr)
library(ggplot2)
library(MuMIn)
library(tidyr)
library(tidyverse)
library(broom)
library(lme4)
library(scales)
library(glmmTMB)
library(DHARMa)
library(emmeans)


# Cargar el archivo
load("MT_mice+insect+plant.RData")
load("~/GitHub/Artropodos-Roedores-MNT/datos_Pol+Cas+artrop_interp.RData")
str(montana.full)
str(datos_pc)

# Filtrado del dataset
montana_filtrado <- montana.full %>%
  mutate( 
    anio = as.numeric(substr(session, 1, 4)),
    mes = as.numeric(substr(session, 5, 6)),
    fecha = as.Date(paste(anio, mes, "01", sep = "-") )) %>%
  # Aplicar filtros por año y sitio
  filter(site %in% c("Cascade", "Polson"),
         anio >= 2000, anio <= 2012) %>%
  # Seleccionar solo columnas necesarias
  select(site, session, time, MNA, biomF, anio, mes)

# Ver un resumen
glimpse(montana_filtrado)

datos_plot <- montana_filtrado %>%
  mutate(
    anio = as.numeric(substr(session, 1, 4)),
    mes = as.numeric(substr(session, 5, 6)),
    fecha = as.Date(paste(anio, mes, "01", sep = "-")),
    MNA_escalado = rescale(MNA),
    biomF_escalado = rescale(biomF)
  ) %>%
  pivot_longer(cols = c(MNA_escalado, biomF_escalado),
               names_to = "variable", values_to = "valor") %>%
  mutate(
    variable = recode(variable,
                      MNA_escalado = "MNA (roedores)",
                      biomF_escalado = "Biomasa F")
  )
# Graficar
ggplot(datos_plot, aes(x = fecha, y = valor, color = variable)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ site, scales = "free_y") +
  scale_color_manual(values = c("MNA (roedores)" = "#d95f02", "Biomasa F" = "#1b9e77")) +
  labs(
    title = "Evolución temporal de MNA y Biomasa F (2000–2012)",
    x = "Fecha",
    y = "Valor escalado (0–1)",
    color = "Variable"
  ) +
  theme_minimal()


#############Hago fusion de las dos bases de datos!! #############

# Renombrar year y month para que coincidan con tu base
datos_pc_fecha <- datos_pc %>%
  rename(anio = year, mes = month) %>%
  select(site, anio, mes,
         insect, arachn, other,
         interp.insect, interp.arachn, interp.other)

# Hacer left_join con tu base original (la de 418 filas)
datos_completos <- montana_filtrado %>%
  left_join(datos_pc_fecha, by = c("site", "anio", "mes")) %>%
  mutate(fecha = as.Date(paste(anio, mes, "01", sep = "-")))

# Verificar estructura
glimpse(datos_completos)

# Asegurate de tener las columnas 'fecha' ya generadas
datos_grafico <- datos_completos %>%
  mutate(
    fecha = as.Date(paste(anio, mes, "01", sep = "-")),
    MNA_escalado = rescale(MNA),
    biomF_escalado = rescale(biomF),
    interp_insect_escalado = rescale(interp.insect)
  ) %>%
  pivot_longer(cols = c(MNA_escalado, biomF_escalado, interp_insect_escalado),
               names_to = "variable", values_to = "valor") %>%
  mutate(
    variable = recode(variable,
                      MNA_escalado = "MNA (roedores)",
                      biomF_escalado = "Biomasa F",
                      interp_insect_escalado = "Insectos interpolados")
  )

# Graficar
ggplot(datos_grafico, aes(x = fecha, y = valor, color = variable)) +
  geom_line(size = 1.2) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ site, scales = "free_y") +
  scale_color_manual(values = c(
    "MNA (roedores)" = "#d95f02",
    "Biomasa F" = "#1b9e77",
    "Insectos interpolados" = "#7570b3")) +
  labs(
    title = "Evolución de MNA, Biomasa F e Insectos interpolados (2000–2012)",
    x = "Fecha",
    y = "Valor escalado (0–1)",
    color = "Variable"
  ) +
  theme_minimal()

########## chequeo lag vegetacion
# Secuencia de lags a probar
lags_veg <- 0:12  
resultados_veg <- data.frame(
  lag = lags_veg,
  coef = NA,
  pvalue = NA
)

for (l in lags_veg) {
  datos_temp <- datos_completos %>%
    arrange(site, fecha) %>%
    group_by(site) %>%
    mutate(lag_biomF = lag(biomF, l)) %>%
    ungroup() %>%
    filter(!is.na(lag_biomF)) %>%
    mutate(lag_biomF_z = scale(lag_biomF))
  
  mod <- glmmTMB(MNA ~ lag_biomF_z + mes + (1|site),
                 family = nbinom2,
                 data = datos_temp)
  
  resumen <- summary(mod)
  
  resultados_veg$coef[resultados_veg$lag == l] <- 
    resumen$coefficients$cond["lag_biomF_z", "Estimate"]
  
  resultados_veg$pvalue[resultados_veg$lag == l] <- 
    resumen$coefficients$cond["lag_biomF_z", "Pr(>|z|)"]
  
  resultados_insect$AIC[resultados_insect$lag == l] <- AIC(mod)
}

# Ordenar por el coeficiente más grande positivo
resultados_veg[order(-resultados_veg$coef), ] # Me COEF max el lag de 10 meses 
resultados_insect
plot(resultados_veg$lag, resultados_veg$AIC, type="b", # Me AIC optimo el lag de 10 meses 
     xlab="Lag veg (meses)", ylab="AIC")
plot(resultados_veg$lag, resultados_veg$coef, type="b", # Me AIC optimo el lag de 10 meses 
     xlab="Lag veg (meses)", ylab="Coef")
# Lag 10  Coef 0.5096642  p-value 3.830999e-02


########## Ahora veo cual es el lag optimo en artropodos con el lag en plantas de 10 meses #############
lags_insect <- 0:12  
resultados_insect <- data.frame(
  lag = lags_insect,
  coef = NA,
  pvalue = NA
)

for (l in lags_insect) {
  datos_temp <- datos_completos %>%
    arrange(site, fecha) %>%
    group_by(site) %>%
    mutate(
      lag_biomF  = lag(biomF, 10),  # fijo en 10
      lag_insect = lag(interp.insect, l)
    ) %>%
    ungroup() %>%
    filter(!is.na(lag_biomF), !is.na(lag_insect)) %>%
    mutate(
      lag_biomF_z  = scale(lag_biomF),
      lag_insect_z = scale(lag_insect)
    )
  
  mod <- glmmTMB(MNA ~ lag_insect_z + lag_biomF_z + mes + (1|site),
                 family = nbinom2,
                 data = datos_temp)
  resumen <- summary(mod)
  
  resultados_insect$coef[resultados_insect$lag == l] <- 
    resumen$coefficients$cond["lag_insect_z", "Estimate"]
  
  resultados_insect$pvalue[resultados_insect$lag == l] <- 
    resumen$coefficients$cond["lag_insect_z", "Pr(>|z|)"]
  
  resultados_insect$AIC[resultados_insect$lag == l] <- AIC(mod)
}

resultados_insect
plot(resultados_insect$lag, resultados_insect$AIC, type="b",
     xlab="Lag insectos (meses)", ylab="AIC")
plot(resultados_insect$lag, resultados_insect$coef, type="b",
     xlab="Lag insectos (meses)", ylab="coef")

##### Lag 2 meses  COEF 0.089899449 p.value 0.4137739 AIC 534.5757 # Ninguno dio significativo


### Ahora veo cual es el LAG optimo en aracnidos

# Definir los lags a evaluar para arácnidos
lags_arachn <- 0:12  
resultados_arachn <- data.frame(
  lag = lags_arachn,
  coef = NA,
  pvalue = NA,
  AIC = NA
)

for (l in lags_arachn) {
  datos_temp <- datos_completos %>%
    arrange(site, fecha) %>%
    group_by(site) %>%
    mutate(
      lag_biomF  = lag(biomF, 10),        # biomasa fija en 10 meses
      lag_insect = lag(interp.insect, 2), # insectos fijos en 2 meses
      lag_arachn = lag(interp.arachn, l)  # arácnidos variable
    ) %>%
    ungroup() %>%
    filter(!is.na(lag_biomF), !is.na(lag_insect), !is.na(lag_arachn)) %>%
    mutate(
      lag_biomF_z  = scale(lag_biomF),
      lag_insect_z = scale(lag_insect),
      lag_arachn_z = scale(lag_arachn)
    )
  
  # Modelo con arácnidos, biomasa y mes
  mod <- glmmTMB(MNA ~ lag_arachn_z + lag_biomF_z + mes + (1|site),
                 family = nbinom2,
                 data = datos_temp)
  
  resumen <- summary(mod)
  
  resultados_arachn$coef[resultados_arachn$lag == l] <- 
    resumen$coefficients$cond["lag_arachn_z", "Estimate"]
  
  resultados_arachn$pvalue[resultados_arachn$lag == l] <- 
    resumen$coefficients$cond["lag_arachn_z", "Pr(>|z|)"]
  
  resultados_arachn$AIC[resultados_arachn$lag == l] <- AIC(mod)
}

# Resultados
resultados_arachn #Me da 7 meses Max Coef neg

# Gráficos
plot(resultados_arachn$lag, resultados_arachn$AIC, type="b",
     xlab="Lag arácnidos (meses)", ylab="AIC")

plot(resultados_arachn$lag, resultados_arachn$coef, type="b",
     xlab="Lag arácnidos (meses)", ylab="Coeficiente")
####################### Modelo veg con lag de 10 meses y insectos con alg de 7 meses.

# Generar lag de 2 meses para insectos y arácnidos, y de 10 meses para biomasa
datos_lags <- datos_completos %>%
  arrange(site, fecha) %>%
  group_by(site) %>%
  mutate(
    lag_insect = lag(interp.insect, 2),
    lag_arachn = lag(interp.arachn, 7),
    lag_biomF  = lag(biomF, 10)
  ) %>%
  ungroup()
datos_modelo <- datos_lags %>%
  filter(!is.na(MNA), !is.na(lag_insect), !is.na(lag_arachn), !is.na(lag_biomF)) %>%
  mutate(
    lag_insect_z = scale(lag_insect),
    lag_arachn_z = scale(lag_arachn),
    lag_biomF_z  = scale(lag_biomF)
  )
library(glmmTMB)

modelo_lag_veg_insect_arachn <- glmmTMB(
  MNA ~ lag_insect_z + lag_arachn_z + lag_biomF_z + mes + (1 | site),
  family = nbinom2,
  data = datos_modelo
)


#Chequeo supuestos todos dan bien!!
# Simulación de residuos
residuos <- simulateResiduals(modelo_lag_veg_insect_arachn)
plot(residuos)
testDispersion(residuos)
testZeroInflation(residuos)

summary(modelo_lag_veg_insect_arachn)


# Predecir en escala link
pred <- predict(modelo_lag_veg_insect_arachn,
                newdata = datos_modelo,
                type = "link",
                se.fit = TRUE)

# Agregar predicciones al dataset
datos_pred <- datos_modelo %>%
  mutate(
    fit_link = pred$fit,
    se_link = pred$se.fit,
    lower_link = fit_link - 1.96 * se_link,
    upper_link = fit_link + 1.96 * se_link,
    pred = exp(fit_link),
    lower = exp(lower_link),
    upper = exp(upper_link)
  )

# Graficar observado vs predicho
ggplot(datos_pred, aes(x = fecha)) +
  geom_point(aes(y = MNA, color = "Observado"), size = 2, alpha = 0.7) +
  geom_line(aes(y = pred, color = "Predicho"), size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "red", alpha = 0.2) +
  facet_wrap(~ site, scales = "free_y") +
  scale_color_manual(
    name = "Abundancia",
    values = c("Observado" = "black", "Predicho" = "red")
  ) +
  labs(
    x = "Fecha",
    y = "Abundancia de roedores (MNA)",
    title = "Abundancia observada y predicha (IC 95%)\nLag: insectos/arácnidos = 2 meses, biomasa = 10 meses"
  ) +
  theme_minimal()

#######y si pruebo con un lugar q no use en el modelo q pasa?
# 1. Armar dataset de GoldCreek con los mismos lags
datos_goldcreek <- datos_completos %>%
  filter(site == "GoldCreek") %>%
  arrange(site, fecha) %>%
  group_by(site) %>%
  mutate(
    lag_insect = lag(interp.insect, 2),   # insectos fijos en 2 meses
    lag_arachn = lag(interp.arachn, 7),   # arácnidos fijos en 7 meses
    lag_biomF  = lag(biomF, 10)           # biomasa fija en 10 meses
  ) %>%
  ungroup() %>%
  filter(!is.na(lag_insect), !is.na(lag_arachn), !is.na(lag_biomF)) %>%
  mutate(
    lag_insect_z = scale(lag_insect),
    lag_arachn_z = scale(lag_arachn),
    lag_biomF_z  = scale(lag_biomF)
  )

# 2. Usar el modelo entrenado (solo Cascade y Polson) para predecir GoldCreek
pred_gold <- predict(modelo_lag_veg_insect_arachn,
                     newdata = datos_goldcreek,
                     type = "link",
                     se.fit = TRUE)

# 3. Guardar las predicciones
datos_gold_pred <- datos_goldcreek %>%
  mutate(
    fit_link = pred_gold$fit,
    se_link = pred_gold$se.fit,
    lower_link = fit_link - 1.96 * se_link,
    upper_link = fit_link + 1.96 * se_link,
    pred = exp(fit_link),
    lower = exp(lower_link),
    upper = exp(upper_link)
  )

# 4. Graficar abundancia observada vs predicha
ggplot(datos_gold_pred, aes(x = fecha)) +
  geom_point(aes(y = MNA, color = "Observado"), size = 2, alpha = 0.7) +
  geom_line(aes(y = pred, color = "Predicho"), size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "red", alpha = 0.2) +
  scale_color_manual(
    name = "Abundancia",
    values = c("Observado" = "black", "Predicho" = "red")
  ) +
  labs(
    x = "Fecha",
    y = "Abundancia de roedores (MNA)",
    title = "Validación externa en GoldCreek\nModelo entrenado con Cascade y Polson"
  ) +
  theme_minimal()


