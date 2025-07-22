library(dplyr)
library(purrr)
library(ggplot2)
library(MuMIn)
library(tidyr)
library(tidyverse)

load("~/GitHub/Artropodos-Roedores-MNT/montana.small completo.RData") #seteo juli

#Analisis exploratorio
str(montana.small)
# Calcular abundancia de invertebrados
datos_rel <- montana.small %>%
  mutate(invertebrados = insect + arachn + other) %>%
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
  select(site, time, MNA, insect, arachn, other) %>%
  pivot_longer(cols = c(MNA, insect, arachn, other),
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
                        insect = "Insectos",
                        arachn = "Arácnidos",
                        other = "Otros"))

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

# Base con solo site, time, insectos y roedores
datos_lag <- montana.small%>%
  mutate(invertebrados = insect + arachn + other) %>%
  select(site, time, MNA, insect) %>%
  filter(!is.na(MNA), !is.na(insect)) %>%
  arrange(site, time)
# Para cada lag, calcular modelo lineal de MNA ~ insectos_lag
ajustar_modelo_con_lag <- function(lag_meses) {
  datos_lagged <- datos_lag %>%
    group_by(site) %>%
    arrange(time) %>%
    mutate(insect_lag = lag(insect, n = lag_meses)) %>%
    ungroup() %>%
    filter(!is.na(insect_lag))
  
  modelo <- lm(MNA ~ insect_lag, data = datos_lagged)
  
  tibble(
    lag = lag_meses,
    r2 = summary(modelo)$r.squared,
    aic = AIC(modelo)
  )
}
resultados_lags <- map_dfr(0:6, ajustar_modelo_con_lag)

# Ver resultados: Evaluar qué lag maximiza el R² o minimiza el AIC.
print(resultados_lags)
#    lag         r2   aic
#<int>      <dbl> <dbl>
#  1     0 0.00762    1737.
#2     1 0.00000436 1684.
#3     2 0.0000149  1630.
#4     3 0.00200    1575.
#5     4 0.00677    1519.
#6     5 0.00688    1464.    Este tiene la mejor relacion R2 y AIC
#7     6 0.000152   1404.
#
ggplot(resultados_lags, aes(x = lag, y = r2)) +
  geom_line() +
  geom_point() +
  labs(title = "R² del modelo MNA ~ insectos_lag",
       x = "Lag (meses)",
       y = "R²") +
  theme_minimal()

#### Hago lo mismo pero un poco mas complejo planteo poison y comparo modelos GLM
library(tidyverse)

# Crear base limpia
datos_lag <- montana.small%>%
  select(site, time, MNA, insect) %>%
  arrange(site, time) %>%
  filter(!is.na(MNA), !is.na(insect))

# Función para ajustar modelo con un lag dado
ajustar_glm_lag <- function(lag_meses) {
  datos_lagged <- datos_lag %>%
    group_by(site) %>%
    mutate(insect_lag = lag(insect, n = lag_meses)) %>%
    ungroup() %>%
    filter(!is.na(insect_lag))
  
  modelo <- glm(MNA ~ insect_lag, family = poisson, data = datos_lagged)
  
  tibble(
    lag = lag_meses,
    aic = AIC(modelo),
    deviance = deviance(modelo),
    n = nrow(datos_lagged)
  )
}

# Probar lags de 0 a 6 meses
resultados_glm <- map_dfr(0:12, ajustar_glm_lag)

# Ver resultados ordenados
resultados_glm %>% arrange(aic)
#     # A tibble: 7 × 4
#lag   aic deviance     n
#<int> <dbl>    <dbl> <int>
#  1     6 4270.    3563.   146
#2     5 4552.    3824.   152
#3     4 4724.    3973.   158
#4     3 4912.    4138.   164
#5     2 5066.    4268.   170
#6     1 5228.    4410.   176
#7     0 5332.    4493.   182
#
ggplot(resultados_glm, aes(x = lag, y = aic)) +
  geom_line() +
  geom_point() +
  labs(title = "AIC del modelo GLM (Poisson) según lag",
       x = "Lag (meses)",
       y = "AIC") +
  theme_minimal()

#Cae en picada es claro que no esta funcionando simplemente disminuye a medida que aumentan los años. voy a agregar variables aleatorias y ver que pasa 
library(tidyverse)
library(broom)
library(lme4)
library(dplyr)
library(scales)
library(ggplot2)

# Crear base limpia
datos_lag <- montana.small%>%
  select(site, time, MNA, insect,month) %>%
  arrange(site, time) %>%
  filter(!is.na(MNA), !is.na(insect))

# Función para ajustar modelo con un lag dado
ajustar_glm_lag <- function(lag_meses) {
  datos_lagged <- datos_lag %>%
    group_by(site) %>%
    mutate(insect_lag = lag(insect, n = lag_meses)) %>%
    ungroup() %>%
    filter(!is.na(insect_lag))
  
  modelo <- glmer(MNA ~ insect_lag*month + (1|site), family = poisson, data = datos_lagged)
  coef_insect <- fixef(modelo)["insect_lag"]
  
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

##    lag     coef   aic deviance     n
#<int>    <dbl> <dbl>    <dbl> <int>
#  1    12  0.122   1650.    1060.   110
#2     8  0.107   2041.    1339.   134
#3     3  0.0972  2529.    1709.   164
#4     2  0.0766  2608.    1765.   170
#5     4  0.0642  2438.    1642.   158
#6     1  0.0424  2673.    1809.   176
#7     0  0.0243  2772.    1886.   182
#8     7  0.00713 2133.    1404.   140
#9     9 -0.0144  2007.    1333.   128
#10    11 -0.0350  1808.    1188.   116
#11     6 -0.0925  2224.    1472.   146
#12     5 -0.165   2283.    1509.   152
#13    10 -0.345   1864.    1217.   122
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



