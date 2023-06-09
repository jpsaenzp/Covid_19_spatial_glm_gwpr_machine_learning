############# Definir directorio y cargar datos ----------------
rm(list=ls())
gc(reset = T)
require(rgdal)
require(pscl)
require(sf)
require(spdep)
require(spatialreg)
require(stringr)
require(performance)
require(AER)
require(ggplot2)
require(vcdExtra)
library(readxl)
library(maptools)
library(tmap)
library(dplyr)
library(lme4)
library(raster)
library(caret)
library(CAST)
library(fasterize)
library(randomForest)
library(kernlab)
library(leaflet)
library(maptools)
require(RColorBrewer)
library(texreg)
library(stargazer)
library(R2BayesX)
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

Variables_Cali$cat_comorb <- as.factor(Variables_Cali$cat_comorb)
Variables_Cali$cat_prop_trabaj <- as.factor(Variables_Cali$cat_prop_trabaj)
Variables_Cali$cat_adultmayor <- as.factor(Variables_Cali$cat_adultmayor)
Variables_Cali$cat_ics_bajo <- as.factor(Variables_Cali$cat_ics_bajo)
Variables_Cali$cat_persondifi <- as.factor(Variables_Cali$cat_persondifi)

#Variables_Cali = spTransform(Variables_Cali,CRS("+init=epsg:21897"))
(l <- length(Variables_Cali))

xy0 = data.frame(x = Variables_Cali$Centroide_X, y = Variables_Cali$Centroide_Y)
coordinates(xy0) <- c('x','y')
#proj4string(xy0) <- CRS("+init=epsg:4326")
#xy0 = spTransform(xy0, CRS("+init=epsg:21897"))

tri_nb_b = nb2listw(tri2nb(xy0), style="B",zero.policy = TRUE)


modelo3 <- glm(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb + Estrato_moda,
               offset = log(Variables_Cali$perimetro_m), #offset = log(Poblacion_total)
               family = poisson, data = Variables_Cali@data)
summary(modelo3) # Mejor?
lmtest::coeftest(modelo3, vcov.=vcovHC(modelo3, type="const"))
pR2(modelo3)

#https://rpubs.com/arquez9512/spatial-statistics
#Fit a Bayesian GLM
#You can include a spatially correlated term based on the adjacency structure 
#by adding a term to the formula specifying a spatially correlated model.

nb_q_93 <- spdep::poly2nb(Variables_Cali, 
                          row.names=unique(as.character(Variables_Cali$ID_BARRIO)))
RBX_gra <- nb2gra(nb_q_93)

Variables_Cali$Estrato_moda <- as.numeric(Variables_Cali$Estrato_moda)
mod_spatial <- bayesx(Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb + Estrato_moda +
                        sx(ID_BARRIO, bs = "spatial", map = RBX_gra),
                      offset = log(Variables_Cali$perimetro_m),
                      family = "poisson", data = st_as_sf(Variables_Cali),
                      control = bayesx.control(seed = 17610407))

bayesx_logfile(mod_spatial)
summary(mod_spatial)

# Map the fitted spatial term only
Variables_Cali$spatial <- fitted(mod_spatial, term = "sx(ID_BARRIO):mrf")[, "Mean"]
spplot(Variables_Cali, zcol = "spatial")

# Map the residuals
Variables_Cali$spatial_resid <- residuals(mod_spatial)[, "mu"]
spplot(Variables_Cali, zcol = "spatial_resid")

# Test residuals for spatial correlation
moran.mc(Variables_Cali$spatial_resid, nb2listw(nb_q_93), 999)
