Biom+ins+arc.R es el archivo que tiene el modelo los otros archivos son de analisis exploratorio

> modelo_lag_veg_insect_arachn <- glmmTMB(
+   MNA ~ lag_insect_z + lag_arachn_z + lag_biomF_z + mes + (1 | site),
+   family = nbinom2,
+   data = datos_modelo
+ )
> summary(modelo_lag_veg_insect_arachn)
 Family: nbinom2  ( log )
Formula:          
MNA ~ lag_insect_z + lag_arachn_z + lag_biomF_z + mes + (1 |      site)
Data: datos_modelo

      AIC       BIC    logLik -2*log(L)  df.resid 
    521.4     535.9    -253.7     507.4        51 

Random effects:

Conditional model:
 Groups Name        Variance Std.Dev.
 site   (Intercept) 0.3988   0.6315   
Number of obs: 58, groups:  site, 2

Dispersion parameter for nbinom2 family (): 1.97 

Conditional model:
             Estimate Std. Error z value Pr(>|z|)    
(Intercept)   3.56681    0.56482   6.315  2.7e-10 ***
lag_insect_z  0.24624    0.11541   2.134   0.0329 *       # Proteinas importantes para reproduccion
lag_arachn_z -0.55966    0.13306  -4.206  2.6e-05 ***    # Posible explicacion-> efecto indirecto por la competencia con los roedores por los  artropodos 
lag_biomF_z   0.43366    0.17507   2.477   0.0132 *  
mes          -0.02497    0.03928  -0.636   0.5249    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

<img width="1064" height="576" alt="image" src="https://github.com/user-attachments/assets/85c4ce68-6e3f-48ac-82b4-a01a7f95cc3e" />
<img width="1655" height="784" alt="image" src="https://github.com/user-attachments/assets/d3059352-e1c1-4709-888e-92352acc16b6" />
<img width="1647" height="780" alt="image" src="https://github.com/user-attachments/assets/0c849370-78ff-43e2-8fdb-9477377d52ac" />
LAG insectos max en 2 meses
<img width="1618" height="678" alt="image" src="https://github.com/user-attachments/assets/d3bd3223-7a8f-4908-874d-dd69d03b59aa" />
LAG max en 10 meses
<img width="1607" height="721" alt="image" src="https://github.com/user-attachments/assets/9707fdf1-45f0-4ab4-910d-aa0bffecc084" />
LAG min en 7 meses
<img width="1204" height="668" alt="image" src="https://github.com/user-attachments/assets/49dae4e0-9789-401a-9f28-a6b86a521019" />
