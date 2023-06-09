############# Paquetes ----------------------
#install.packages("sppois")
#install.packages("devtools")
#devtools::install_github("gregmacfarlane/sppois")
# Este script es para versi?n 4.1.2 de R
rm(list=ls())
gc(reset = TRUE)
# cargue de las librerias
require(rgdal);require(pscl);require(sf);require(spdep);require(spatialreg)
require(stringr);require(performance);require(AER);require(ggplot2)
require(vcdExtra);library(readxl);library(maptools);library(tmap)
library(dplyr);library(lme4);library(raster);library(caret);library(xtable)
library(CAST);library(fasterize);library(randomForest);library(kernlab)
library(leaflet);library(maptools);require(RColorBrewer);library(texreg)
library(stargazer);library(xtable)
#library(sppois) # No le han hecho actualización
############# Definir directorio y cargar datos ----------------
barrios_Cali <- readOGR('databases_and_shapefiles/Barrios Cali.shp')
datos_barrio_cali <- barrios_Cali@data
datos_cali <- read_excel('databases_and_shapefiles/Cali.xlsx')
datos_cali$Estrato_moda <- as.factor(datos_cali$Estrato_moda)

datos_cali <- datos_cali %>%
  mutate(cat_comorb = as.numeric(cut_number(datos_cali$Comorbilidad,3)),
         cat_persondifi = as.numeric(cut_number(datos_cali$Personas_con_dificultades,3)),
         cat_ics_bajo = as.numeric(cut_number(datos_cali$ind_cobert_serv_bajo,3)),
         cat_adultmayor = as.numeric(cut_number(datos_cali$Adulto_mayor,3)),
         cat_prop_trabaj = as.numeric(cut_number(datos_cali$Proporcion_trabajando,3)))

table(datos_cali$cat_comorb)
table(datos_cali$cat_prop_trabaj)
table(datos_cali$cat_adultmayor)
table(datos_cali$cat_ics_bajo)
tapply(datos_cali$ind_cobert_serv_bajo, datos_cali$cat_ics_bajo, summary)
tapply(datos_cali$Comorbilidad, datos_cali$cat_comorb, summary)
tapply(datos_cali$Proporcion_trabajando, datos_cali$cat_prop_trabaj, summary)
tapply(datos_cali$Adulto_mayor, datos_cali$cat_adultmayor, summary)

barrios_Cali@data = left_join(barrios_Cali@data, datos_cali, by = 'ID_BARRIO')

corrplot::corrplot.mixed(cor(datos_cali[,c(4:9,11:16)], use = "complete.obs"))
Variables_Cali <- barrios_Cali[,c('ID_BARRIO','Barrio_Urbanizacion_o_Sector',
                                  'Comuna','Poblacion_total','Covid',
                                  'Adulto_mayor','ind_cobert_serv_bajo',
                                  'Comorbilidad', 'razon_letalidad',
                                  'Proporcion_trabajando','Estrato_moda',
                                  'Personas_con_dificultades',
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

# Primera estimacion del GLM
modelo3 <- glm(Covid ~ Adulto_mayor + ind_cobert_serv_bajo + log(Comorbilidad) + 
                 razon_letalidad + Proporcion_trabajando%in%Estrato_moda,
               family = poisson, data = Variables_Cali@data)
summary(modelo3) # este ya no

modelo3 <- glm(Covid ~ Adulto_mayor + ind_cobert_serv_bajo + Comorbilidad +
                 Proporcion_trabajando%in%Estrato_moda,
               offset = log(perimetro_m), #log(Poblacion_total)
               family = poisson, data = Variables_Cali@data)
summary(modelo3) # proporcion trabajando anidado en estrato moda no es significativo

modelo3 <- glm(Covid ~ Adulto_mayor + ind_cobert_serv_bajo + Comorbilidad +
                 Estrato_moda,
               offset = log(perimetro_m), #log(Poblacion_total)
               family = poisson, data = Variables_Cali@data)
summary(modelo3)

#Adulto_mayor + ind_cobert_serv_bajo + cat_comorb + Estrato_moda + cat_adultmayor +
#Comorbilidad
modelo3 <- glm(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb + Estrato_moda,
               offset = log(Variables_Cali$perimetro_m), #offset = log(Poblacion_total)
               family = poisson, data = Variables_Cali@data)
summary(modelo3)
# cat_comorb es importante y significativa.
# los ultimos 2 modelos de arriba presentan mejor significancia. el ultimo es el mejor
# las personas con dificultades no son significativas

#Proporcion_trabajando%in%cat_estrato
# Se ensayo con la proporcion trabajando anidada en la categorica de a dos estratos,
# y no da significativa

# analisis con estrato:
kruskal.test(formula = Covid ~ Estrato_moda, data = Variables_Cali@data)
ggplot(Variables_Cali@data) + geom_boxplot(aes(Estrato_moda, Covid)) 
table(Variables_Cali@data$Estrato_moda)
tapply(Variables_Cali@data$Covid, Variables_Cali@data$Estrato_moda, summary)
#
Variables_Cali$cat_estrato <- ifelse(as.numeric(as.character(Variables_Cali$Estrato_moda)) <= 2,
                                     "Bajo", ifelse(as.numeric(as.character(Variables_Cali$Estrato_moda)) >= 5,
                                                    "Alto", "Medio"))
#Variables_Cali$cat_estrato <- ifelse(as.numeric(as.character(Variables_Cali$Estrato_moda)) <= 3,
#                                     "Bajo", "Medio_alto")
#
Variables_Cali$cat_estrato <- as.factor(Variables_Cali$cat_estrato)
Variables_Cali$cat_estrato <- relevel(Variables_Cali$cat_estrato, ref = "Bajo")
#modelo3 <- glm(formula = Covid ~ cat_estrato,
#               offset = log(perimetro_m), #offset = log(Poblacion_total)
#               family = poisson, data = Variables_Cali@data)
#summary(modelo3)
# El estrato no es suficiente significativo 

# modelo por zonas por ej por estrato
table(Variables_Cali@data$cat_estrato)
modelo4 <- glm(formula = Covid ~ Adulto_mayor + cat_comorb,
               offset = log(perimetro_m), 
               family = poisson, data = Variables_Cali@data,
               subset = (cat_estrato == 'Alto'))
summary(modelo4)
lmtest::coeftest(modelo4, vcov.=vcovHC(modelo4, type="const"))
Devianza = ((modelo4$null.deviance - modelo4$deviance)/modelo4$null.deviance)*100;Devianza
cbind(exp(coef(modelo4)),exp(confint(modelo4)))

#estrato bajo: Adulto_mayor 1.00006939 , cat_comorb2  2.00276299 , cat_comorb3  3.58817204
#estrato medio: Adulto_mayor 1.00003478 , cat_comorb2  2.13387205 , cat_comorb3  3.46338178
#estrato alto: Adulto_mayor 1.00036229 , cat_comorb2  2.12955075 , cat_comorb3  2.40327927

# ahora con proporcion trabajando categorica
#modelo3 <- glm(Covid ~ Adulto_mayor + ind_cobert_serv_bajo + cat_comorb +
#                 Estrato_moda + cat_prop_trabaj,
#               offset = log(Poblacion_total),
#               family = poisson, data = Variables_Cali@data)
#summary(modelo3)
# la proporcion trabajando categorica no es significativa

cbind(exp(coef(modelo3)),exp(confint(modelo3)))


#https://rpubs.com/jaimeisaacp/774744
#cov.modelo3 <- sandwich::vcovHC(modelo3, type = "HC0")
#std.err <- sqrt(diag(cov.modelo3))
#r.est <- cbind(
#  Estimate = coef(modelo3), "Robust SE" = std.err,
#  "Pr(>|z|)" = 2 * pnorm(abs(coef(modelo3) / std.err), lower.tail = FALSE),
#  LL = coef(modelo3) - 1.96 * std.err,
#  UL = coef(modelo3) + 1.96 * std.err
#)
#r.est
## drop the p-value and exp the coefficients
#rexp.est <- exp(r.est[, -3])
## compute the SE for the exp coefficients
#s <- msm::deltamethod(list(
#  ~ exp(x1), ~ exp(x2), ~ exp(x3), ~ exp(x4), ~ exp(x5), ~ exp(x6), ~ exp(x7),
#  ~ exp(x8), ~ exp(x9), ~ exp(x10), ~ exp(x11)
#), coef(modelo3), cov.modelo3)
## replace with the new SE's
#rexp.est[, "Robust SE"] <- s
#rexp.est[, "LL"] <- rexp.est[, 1] - 1.96 * s
#rexp.est[, "UL"] <- rexp.est[, 1] + 1.96 * s
#rexp.est


# Test de sobredispersión
#library(AER)
#dispersiontest(modelo3) # Hay mucha sobredispersion
#modBN <- glm.nb(Covid ~ Adulto_mayor + ind_cobert_serv_bajo + log(Comorbilidad) + 
#                  razon_letalidad + Proporcion_trabajando%in%Estrato_moda,
#                data=Variables_Cali@data)
#summary(modBN)
#cbind(exp(coef(modBN)),exp(confint(modBN)))
#pchisq(modBN$deviance, df=modBN$df.residual, lower.tail=FALSE)

lmtest::coeftest(modelo3, vcov.=vcovHC(modelo3, type="const"))
Devianza = ((modelo3$null.deviance - modelo3$deviance)/modelo3$null.deviance)*100
Devianza

pR2(modelo3)

#coeftest(modelo3, vcov = sandwich::vcovHC(x = modelo3, type = "HC1"))

anova(aov(modelo3))
Anova(modelo3, type="II")
Anova(modelo3, type="III")
xtable(modelo3)
stargazer(modelo3,single.row = T, no.space = T)

#validacion de supuestos 
#distribucion
library(tseries)
library(mvnormtest)
library(MVN)
library(olsrr)
library(lmtest)
library(stats)
par(mfrow =c(1,1))
hist(x = modelo3$residuals, main = "Histograma de los Residuales", axes = TRUE,
     breaks = "Sturges", prob=T, xlab='Residuales del modelo', ylab='Densidad')
lines(density(modelo3$residuals, adjust = 7.7))
par(mfrow=c(2,2))
plot(modelo3)

library(vcd) ## loading vcd package
gf <- goodfit(modelo3$residuals,type= "poisson",method= "MinChisq")
stargazer(summary(gf)) # Si sigue una distribucion poisson

#Homocedasticidad
bartlett.test(rstudent(modelo3)~c(rep(1,67),rep(2,67),rep(3,67),rep(4,67),rep(5,66)))
# No hay homocedasticidad ya que p-valor < 0.05
leveneTest(rstudent(modelo3)~as.factor(c(rep(1,67),rep(2,67),rep(3,67),rep(4,67),
                                         rep(5,66))))
# Si hay homocedasticidad ya que p-valor > 0.05
fligner.test(rstudent(modelo3)~as.factor(c(rep(1,67),rep(2,67),rep(3,67),rep(4,67),
                                           rep(5,66))))
# Si hay homocedasticidad ya que p-valor > 0.05
bptest(modelo3)
# No hay homocedasticidad ya que p-valor < 0.05
#No homogeneidad de varianzas 

#multicolinealidad
vif(modelo3) # No hay multicolinealidad
xtable(vif(modelo3))
#en la variable de interaccion

#Autocorrelacion espacial
Variables_Cali = spTransform(Variables_Cali,CRS("+init=epsg:21897"))
(l <- length(Variables_Cali))
############# Centroides y matrices de vecindad -----------
xy0 = data.frame(x = Variables_Cali$Centroide_X, y = Variables_Cali$Centroide_Y)
coordinates(xy0) <- c('x','y')
proj4string(xy0) <- CRS("+init=epsg:4326")
xy0 = spTransform(xy0, CRS("+init=epsg:21897"))
queen_nb_b = nb2listw(poly2nb(Variables_Cali,queen=TRUE), style="B",
                      zero.policy = TRUE)
queen_nb_w = nb2listw(poly2nb(Variables_Cali,queen=TRUE), style="W",
                      zero.policy = TRUE)
relative_nb_b = nb2listw(graph2nb(relativeneigh(xy0), sym=TRUE), style="B",
                         zero.policy = TRUE)
relative_nb_w = nb2listw(graph2nb(relativeneigh(xy0), sym=TRUE), style="W",
                         zero.policy = TRUE)
knn_nb_b = nb2listw(knn2nb(knearneigh(xy0, k = 5)), style="B",zero.policy = TRUE)
knn_nb_w = nb2listw(knn2nb(knearneigh(xy0, k = 5)), style="W",zero.policy = TRUE)

Variables_Cali$residLR = residuals(modelo3)
moran.mc(x = Variables_Cali$residLR, listw = queen_nb_b, nsim = 999,zero.policy = T,
         alternative = "two.sided")
moran.test(Variables_Cali$residLR, relative_nb_b, alternative="two.sided")
moran.plot(Variables_Cali$residLR, relative_nb_b, zero.policy = T,
           pch=16, col="black", cex=.5, quiet=TRUE,
           labels=Variables_Cali$Barrio_Urbanizacion_o_Sector,
           xlab="Residuales del GLM", 
           ylab="Residuales del GLM (Spatial Lag)")


############# SpatialFiltering --------------
#https://stats.stackexchange.com/questions/193440/spatial-autoregressive-poisson-
#model-in-r
############# spaMM ------------------
library(spaMM)
# Fitting function for fixed- and mixed-effect models with GLM response.
#Covid_spamm <- fitme(Covid ~ Adulto_mayor + ind_cobert_serv_bajo + 
#                       Comorbilidad + razon_letalidad + 
#                       Matern(1 | Centroide_X + Centroide_Y), family = poisson(),
#                     data = Variables_Cali@data) # this take a bit of time
Covid_spamm <- fitme(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                       offset(log(perimetro_m)) +
                       Estrato_moda + Matern(1 | Centroide_X + Centroide_Y),
                     family = poisson(),
                     data = Variables_Cali@data) # this take a bit of time

summary(Covid_spamm)

# CAR by Laplace with 'inner' estimation of rho
W <- as(spdep::nb2listw(tri2nb(coords = xy0,row.names = Variables_Cali$ID_BARRIO), 
                        style="B",zero.policy = TRUE), "CsparseMatrix")
#Fits a (spatially) correlated mixed model, for given correlation parameters
HLCor_Covid <- HLCor(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                       offset(log(perimetro_m)) +
                       Estrato_moda + adjacency(1|ID_BARRIO),
                     adjMatrix=W, family=poisson(), data=Variables_Cali@data, 
                     method="ML")

summary(HLCor_Covid)
# Fits a mixed model, typically a spatial GLMM.
corrHL_Covid <- corrHLfit(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                            offset(log(perimetro_m)) +
                            Estrato_moda + adjacency(1|ID_BARRIO),
                          adjMatrix=W, family=poisson(), 
                          data=Variables_Cali@data, method="ML")
summary(corrHL_Covid)

############# glmmTMB -------------
library(glmmTMB)
#spaMM
#Incluir coordenadas como covariables
#https://datascienceplus.com/spatial-regression-in-r-part-1-spamm-vs-glmmtmb/
library(DHARMa)
library(geoR)
library(viridis)
Variables_Cali$pos <- numFactor(scale(Variables_Cali$Centroide_X), 
                                scale(Variables_Cali$Centroide_Y))
Variables_Cali$ID <- factor(rep(1, nrow(Variables_Cali)))
# Fit a generalized linear mixed model (GLMM) using Template Model Builder (TMB).
glmmTMB_Covid <- glmmTMB(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                           offset(log(perimetro_m)) +
                           Estrato_moda + mat(pos + 0 | ID), 
                         Variables_Cali@data)
summary(glmmTMB_Covid)

#https://aip.scitation.org/doi/pdf/10.1063/5.0041857
#https://rdrr.io/github/gregmacfarlane/sppois/src/R/sarpoisson.R
#https://rpubs.com/corey_sparks/250410
#https://rpubs.com/corey_sparks/111362
#https://stat.ethz.ch/pipermail/r-sig-geo/2016-February/023978.html