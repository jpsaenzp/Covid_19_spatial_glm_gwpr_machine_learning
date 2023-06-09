# http://rstudio-pubs-static.s3.amazonaws.com/483800_14662f9098974b31b20110a0ca3
# 45708.html#025_residuales_espacialmente_correlacionados

rm(list=ls())
gc(reset = TRUE)
# Cargar librerías
library(sp)   
library(spdep)
library(RColorBrewer)
library(classInt)  
library(grid)     
library(gridExtra)
library(ggplot2)
library(rgdal)
library(readxl)
library(dplyr)
library(tmap)
library(viridis)
library(hglm) # HGLM model
library(spatialreg)
library(Matrix)
library(sf)
library(R2BayesX)
library(INLA)
# Cargar archivos
barrios_Cali <- readOGR('databases_and_shapefiles/Barrios Cali.shp')
datos_barrio_cali <- barrios_Cali@data
datos_cali <- read_excel('databases_and_shapefiles/Cali.xlsx')

datos_cali$Estrato_moda <- as.factor(datos_cali$Estrato_moda)

datos_cali <- datos_cali %>% mutate(cat_comorb = as.numeric(cut_number(datos_cali$Comorbilidad,3)),
                                    cat_persondifi = as.numeric(cut_number(datos_cali$Personas_con_dificultades,3)),
                                    cat_ics_bajo = as.numeric(cut_number(datos_cali$ind_cobert_serv_bajo,3)),
                                    cat_adultmayor = as.numeric(cut_number(datos_cali$Adulto_mayor,3)),
                                    cat_prop_trabaj = as.numeric(cut_number(datos_cali$Proporcion_trabajando,3)))

barrios_Cali@data = left_join(barrios_Cali@data, datos_cali, by = 'ID_BARRIO')
# Análisis exploratorio de los datos al igual que selección de variables
corrplot::corrplot.mixed(cor(datos_cali[,c(4:9,11:16)], use = "complete.obs"))
Variables_Cali <- barrios_Cali[,c('ID_BARRIO','Barrio_Urbanizacion_o_Sector',
                                  'Comuna','Poblacion_total','Covid',
                                  'Adulto_mayor','ind_cobert_serv_bajo',
                                  'Comorbilidad', 'razon_letalidad',
                                  'Proporcion_trabajando','Estrato_moda',
                                  'area_ha','perimetro_m',
                                  'Centroide_X','Centroide_Y',
                                  'cat_comorb', 'cat_prop_trabaj',
                                  'cat_adultmayor', 'cat_ics_bajo',
                                  'cat_persondifi')]
Variables_Cali <- Variables_Cali[-c(65,112,304,326),]
View(Variables_Cali@data)

Variables_Cali$cat_comorb <- as.factor(Variables_Cali$cat_comorb)
Variables_Cali$cat_prop_trabaj <- as.factor(Variables_Cali$cat_prop_trabaj)
Variables_Cali$cat_adultmayor <- as.factor(Variables_Cali$cat_adultmayor)
Variables_Cali$cat_ics_bajo <- as.factor(Variables_Cali$cat_ics_bajo)
Variables_Cali$cat_persondifi <- as.factor(Variables_Cali$cat_persondifi)

#Variables_Cali = spTransform(Variables_Cali,CRS("+init=epsg:21897"))
#(l <- length(Variables_Cali))

#xy0 = data.frame(x = Variables_Cali$Centroide_X, y = Variables_Cali$Centroide_Y)
#coordinates(xy0) <- c('x','y')
#proj4string(xy0) <- CRS("+init=epsg:4326")
#xy0 = spTransform(xy0, CRS("+init=epsg:21897"))

# Ver mapa por casos de covid-19 por barrio en Cali
library(mapview)
mapView(Variables_Cali, zcol="Covid", legend=T)

## CAR random effects with hglm
# Matrices de vecindad para el modelo hglm
nb_q_93 <- spdep::poly2nb(Variables_Cali, 
                          row.names=unique(as.character(Variables_Cali$ID_BARRIO)))
W <- as(spdep::nb2listw(nb_q_93, style = "B"), "CsparseMatrix")

form <- formula(Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb + Estrato_moda)

HGLM_car <- hglm(fixed = form, offset = log(perimetro_m), random = ~ 1 | ID_BARRIO,
                 data = Variables_Cali@data,
                 family = poisson(), rand.family = CAR(D=W))
Variables_Cali$HGLM_ss <- HGLM_car$ranef[,1]
HGLM_car

## SAR random effects with hglm
HGLM_sar <- hglm(fixed = form, offset = log(perimetro_m), random = ~ 1 | ID_BARRIO,
                 data = Variables_Cali@data,
                 family = poisson(), rand.family = SAR(D=W))
Variables_Cali$HSAR_ss <- HGLM_sar$ranef[,1]
HGLM_sar

#ICAR random effects with R2BayesX
suppressPackageStartupMessages(library(R2BayesX))

#Comunas_Cali <- aggregate(st_as_sf(Variables_Cali[, "Comuna"]), 
#                          list(ids = Variables_Cali$Comuna), unique)
#nb_q_93 <- spdep::poly2nb(Comunas_Cali, 
#                          row.names=unique(as.character(Comunas_Cali$Comuna)))

#RBX_gra <- nb2gra(nb_q_93)
#all.equal(row.names(RBX_gra), attr(nb_q_93, "region.id"))
#all.equal(unname(diag(RBX_gra)), spdep::card(nb_q_93))

# Más matrices de vecindad
RBX_gra <- nb2gra(nb_q_93)
all.equal(row.names(RBX_gra), attr(nb_q_93, "region.id"))
all.equal(unname(diag(RBX_gra)), spdep::card(nb_q_93))

Variables_Cali$Estrato_moda <- as.numeric(Variables_Cali$Estrato_moda)
# Modelo BayesX
BX_mrf <- bayesx(update(form, . ~ . + sx(ID_BARRIO, bs = "mrf", map = RBX_gra)),
                 offset = log(Variables_Cali$perimetro_m),
                 family = "poisson", data = st_as_sf(Variables_Cali), 
                 method = "MCMC", 
                 iterations = 12000, burnin = 2000, step = 2, seed = 123)
bayesx_logfile(BX_mrf)
BX_mrf$effects

Variables_Cali$BX_ss <- BX_mrf$effects["sx(ID_BARRIO):mrf"][[1]]$Mean

# Un hallazgo es que la variable Covid presenta mucha dispersión.
#Hay que Agregarle offset

##ICAR and Leroux random effects with INLA
suppressPackageStartupMessages(library(INLA))

Variables_Cali$Estrato_moda <- as.factor(Variables_Cali$Estrato_moda)
ID2 <- as.integer(as.factor(Variables_Cali$ID_BARRIO))

INLA_ss <- inla(formula = update(form, . ~ . + f(ID2, model = "besag", graph = W)), 
                family = "poisson", offset = log(Variables_Cali$perimetro_m),
                data = data.frame(Variables_Cali@data), verbose = TRUE)

Variables_Cali$INLA_ss <- INLA_ss$summary.random$ID2$mean

#The sparse Leroux representation 
M <- Diagonal(nrow(W), rowSums(W)) - W
Cmatrix <- Diagonal(nrow(M), 1) -  M

INLA_lr <- inla(formula = update(form, . ~ . + f(ID2, model = "generic1",
                                                 Cmatrix = Cmatrix)),
                family = "poisson", offset = log(Variables_Cali$perimetro_m),
                data = data.frame(Variables_Cali@data))
Variables_Cali$INLA_lr <- INLA_lr$summary.random$ID2$mean

### ICAR random effects with mgcv::gam()
library(mgcv)
Comunas_Cali <- aggregate(st_as_sf(Variables_Cali[, "Comuna"]), 
                          list(ids = Variables_Cali$Comuna), unique)
nb_q_93 <- spdep::poly2nb(Comunas_Cali, 
                          row.names=unique(as.character(Comunas_Cali$Comuna)))

names(nb_q_93) <- attr(nb_q_93, "region.id")
#Variables_Cali$ID_BARRIO <- as.factor(Variables_Cali$ID_BARRIO)
Variables_Cali$Comuna <- as.factor(Variables_Cali$Comuna)

#The `"REML"` method of `bayesx()` gives the same results as `gam()` using `"REML"`
GAM_MRF <- gam(formula = update(form, . ~ . + s(Comuna, bs = "mrf",
                                                xt = list(nb = nb_q_93))),
               family = poisson, offset = log(Variables_Cali$perimetro_m),
               data = st_as_sf(Variables_Cali), method = "REML")

Variables_Cali$Estrato_moda <- as.numeric(Variables_Cali$Estrato_moda)
GAM_mrf <- bayesx(update(form, . ~ . + sx(ID_BARRIO, bs = "mrf", map = RBX_gra)),
                  offset = log(Variables_Cali$perimetro_m),
                  family = "poisson", data = st_as_sf(Variables_Cali),
                  method = "REML",
                  iterations = 12000, burnin = 2000, step = 2, seed = 123)
bayesx_logfile(GAM_mrf)
# Es muy raro porque si se hiciera con el ID_BARRIO no se puede estimar,
# hay mas coeficientes que datos. El GAM_MRF se estimó con el id de comuna.
Variables_Cali$Estrato_moda <- as.factor(Variables_Cali$Estrato_moda)

ssre <- predict(GAM_MRF, type = "terms", se = FALSE)[, "s(Comuna)"]
all(sapply(tapply(ssre, list(Variables_Cali$Comuna), c), 
           function(x) length(unique(x)) == 1))

Comunas_Cali$GAM_ss <- aggregate(ssre, list(Variables_Cali$Comuna), head, n=1)$x
Variables_Cali$GAM_ss <- ssre

### Upper level random effects: summary
View(Variables_Cali@data)

tm_shape(Variables_Cali) + tm_fill(c("HGLM_ss", "HSAR_ss", "INLA_lr", "INLA_ss", 
                                     "BX_ss", "GAM_ss"), midpoint = 0, 
                                   title = "SSRE")  + 
  tm_facets(free.scales = FALSE) + tm_borders(lwd = 0.3, alpha = 0.4) + 
  tm_layout(panel.labels = c("hglm CAR", "hsar SAR", "inla Leroux", "inla ICAR", 
                             "bayesx ICAR", "gam ICAR"))
