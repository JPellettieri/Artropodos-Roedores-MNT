library(dplyr)
library(purrr)
library(ggplot2)
library(MuMIn)
library(tidyr)
library(tidyverse)
library(broom)
library(lme4)
library(dplyr)
library(scales)
library(ggplot2)
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
  left_join(datos_pc_fecha, by = c("site", "anio", "mes"))

# Verificar estructura
glimpse(datos_completos)

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

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

#######################Modelo veg con lag de 9 meses y insectos con alg de 7 meses.

library(dplyr)

# Generar lag de 7 meses para insectos y arácnidos, y de 9 meses para biomasa
datos_lags <- datos_completos %>%
  arrange(site, fecha) %>%
  group_by(site) %>%
  mutate(
    lag_insect = lag(interp.insect, 7),
    lag_arachn = lag(interp.arachn, 7),
    lag_biomF  = lag(biomF, 9)
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
  MNA ~ lag_insect_z + lag_arachn_z + lag_biomF_z + month + (1 | site),
  family = nbinom2,
  data = datos_modelo
)
library(DHARMa)

# Simulación de residuos
residuos <- simulateResiduals(modelo_lag_veg_insect_arachn)
plot(residuos)
testDispersion(residuos)
testZeroInflation(residuos)

summary(modelo_lag_veg_insect_arachn)
