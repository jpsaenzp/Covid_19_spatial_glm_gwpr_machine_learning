############# Paquetes ----------------------
# https://rpubs.com/corey_sparks/250410
# Cargar librerias
rm(list=ls())
gc(reset = TRUE)
require(rgdal);require(pscl);require(sf);require(spdep);require(spatialreg);require(stringr)
require(performance);require(AER);require(ggplot2);require(vcdExtra);library(readxl)
library(maptools);library(tmap);library(dplyr);library(lme4);library(raster)
library(caret);library(CAST);library(fasterize);library(randomForest);library(kernlab)
library(leaflet);library(maptools);require(RColorBrewer);library(texreg)
library(stargazer);library(ggspatial)

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

barrios_Cali@data = left_join(barrios_Cali@data, datos_cali, by = 'ID_BARRIO')
# variables
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

mean(na.omit(Variables_Cali$Covid))
sqrt(var(na.omit(Variables_Cali$Covid)))/mean(na.omit(Variables_Cali$Covid))
Variables_Cali = spTransform(Variables_Cali,CRS("+init=epsg:21897"))
(l <- length(Variables_Cali))
############# Centroides y matrices de vecindad -----------
xy0 = data.frame(x = Variables_Cali$Centroide_X, y = Variables_Cali$Centroide_Y)
coordinates(xy0) <- c('x','y')
proj4string(xy0) <- CRS("+init=epsg:4326")
xy0 = spTransform(xy0, CRS("+init=epsg:21897"))
#Create a k=4 nearest neighbor set
nb4 <- knearneigh(coordinates(Variables_Cali), k=4)
nb4 <- knn2nb(nb4)
nb4 <- make.sym.nb(nb4)
us.wt4 <- nb2listw(nb4, style="W")
# evaluar matrices de vecindad para ver cual es mejor
rook_nb_b = nb2listw(poly2nb(Variables_Cali,queen=FALSE), style="B",
                     zero.policy = TRUE)
rook_nb_w = nb2listw(poly2nb(Variables_Cali,queen=FALSE), style="W",
                     zero.policy = TRUE)
queen_nb_b = nb2listw(poly2nb(Variables_Cali,queen=TRUE), style="B",
                      zero.policy = TRUE)
queen_nb_w = nb2listw(poly2nb(Variables_Cali,queen=TRUE), style="W",
                      zero.policy = TRUE)
#Graphs neighbours
trinb = tri2nb(xy0)
options(warn = -1)
tri_nb_b = nb2listw(tri2nb(xy0), style="B",zero.policy = TRUE)
tri_nb_w = nb2listw(tri2nb(xy0), style="W",zero.policy = TRUE)
#soi_nb_b = nb2listw(graph2nb(soi.graph(trinb,xy0)), style="B",zero.policy = TRUE)
#soi_nb_w = nb2listw(graph2nb(soi.graph(trinb,xy0)), style="W",zero.policy = TRUE)
relative_nb_b = nb2listw(graph2nb(relativeneigh(xy0), sym=TRUE), style="B",
                         zero.policy = TRUE)
relative_nb_w = nb2listw(graph2nb(relativeneigh(xy0), sym=TRUE), style="W",
                         zero.policy = TRUE)
gabriel_nb_b = nb2listw(graph2nb(gabrielneigh(xy0), sym=TRUE), style="B",
                        zero.policy = TRUE)
gabriel_nb_w = nb2listw(graph2nb(gabrielneigh(xy0), sym=TRUE), style="W",
                        zero.policy = TRUE)

#Distance neighbours

knn1_nb_b=nb2listw(knn2nb(knearneigh(xy0, k = 5)), style="B",zero.policy = TRUE)
knn1_nb_w=nb2listw(knn2nb(knearneigh(xy0, k = 5)), style="W",zero.policy = TRUE)
knn2_nb_b=nb2listw(knn2nb(knearneigh(xy0, k = 2)), style="B",zero.policy = TRUE)
knn2_nb_w=nb2listw(knn2nb(knearneigh(xy0, k = 2)), style="W",zero.policy = TRUE)
knn3_nb_b=nb2listw(knn2nb(knearneigh(xy0, k = 3)), style="B",zero.policy = TRUE)
knn3_nb_w=nb2listw(knn2nb(knearneigh(xy0, k = 3)), style="W",zero.policy = TRUE)
knn4_nb_b=nb2listw(knn2nb(knearneigh(xy0, k = 4)), style="B",zero.policy = TRUE)
knn4_nb_w=nb2listw(knn2nb(knearneigh(xy0, k = 4)), style="W",zero.policy = TRUE)
knn6_nb_b=nb2listw(knn2nb(knearneigh(xy0, k = 6)), style="B",zero.policy = TRUE)
knn6_nb_w=nb2listw(knn2nb(knearneigh(xy0, k = 6)), style="W",zero.policy = TRUE)

mat=list(us.wt4,rook_nb_b,rook_nb_w,
         queen_nb_b,queen_nb_w,
         tri_nb_b,tri_nb_w,
         gabriel_nb_b,gabriel_nb_w,
         relative_nb_b,relative_nb_w,
         knn1_nb_b,knn1_nb_w,
         knn2_nb_b,knn2_nb_w,
         knn3_nb_b,knn3_nb_w,
         knn4_nb_b,knn4_nb_w,
         knn6_nb_b,knn6_nb_w)

############# Primer modelo --------------
# Modelo lineal generalizado
modelo1 <- glm(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb + Estrato_moda,
               offset = log(perimetro_m), #offset = log(Poblacion_total)
               family = poisson, data = Variables_Cali@data)
summary(modelo1)

lmtest::coeftest(modelo1, vcov.=vcovHC(modelo1, type="const"))
anova(modelo1)
Variables_Cali$residLR = residuals(modelo1)

1-modelo1$deviance/modelo1$null.deviance # R^2

scale <- sqrt(modelo1$deviance/modelo1$df.residual)
round(exp(coef(modelo1)), 3)
hist(fitted(modelo1)/Variables_Cali$Covid,main="Distribution from Poisson Model" )

############# Probando correlación espacial con I de moran -----------
aux=numeric(0)
options(warn = -1)
{
  for(i in 1:length(mat))
    aux[i]=moran.test(Variables_Cali$residLR,mat[[i]],
                      alternative="two.sided")$"p"
}
aux # La sexta
moran.test(Variables_Cali$residLR, mat[[which.min(aux)]], alternative="two.sided")
moran.plot(Variables_Cali$residLR, mat[[which.min(aux)]], zero.policy = T,
           pch=16, col="black", cex=.5, quiet=TRUE,
           labels=Variables_Cali$Barrio_Urbanizacion_o_Sector,
           xlab="Residuales del GLM", 
           ylab="Residuales del GLM (Spatial Lag)")

############# Moran Eigenvector ---------------
lm.morantest(modelo1, listw = knn3_nb_w)
#mejores: knn3_nb_w, knn3_nb_b, queen_nb_b presentaron mayor p-valor
#hay que resaltar que de todas las 17 matrices espaciales con que se probo, TODAS
#no rechazan la hiptesis nula de presencia de autocorrelacion espacial
# relative_nb_b y relative_nb_w presentaron menor p-valor

me.fit <- ME(Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb + Estrato_moda,
             offset = log(perimetro_m), data = Variables_Cali@data,
             family="poisson", listw = relative_nb_w, verbose=T,alpha=0.05)
me.fit

fits <- data.frame(fitted(me.fit))
Variables_Cali$me1 <- fits$vec1
Variables_Cali$me2 <- fits$vec20
Variables_Cali$me3 <- fits$vec31

lm.morantest(modelo1, listw = relative_nb_b)
# Segundo modelo con otra matriz de distancia
me.fit2 <- ME(Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb + Estrato_moda,
              offset = log(perimetro_m), data = Variables_Cali@data,
              family="poisson", listw = relative_nb_b, verbose=T,alpha=0.05)
me.fit2
fits2 <- data.frame(fitted(me.fit2))
Variables_Cali$me4 <- fits2$vec54
Variables_Cali$me5 <- fits2$vec2

############# Mapas --------------------
spplot(Variables_Cali, c("me1","me2","me3"), at=quantile(c(Variables_Cali$me1,
                                                           Variables_Cali$me2,
                                                           Variables_Cali$me3),
                                             p=seq(0,1,length.out = 7)),
       col.regions=brewer.pal(n=7, "Greys"), 
       main="Distribución espacial de los vectores propios de Moran con relative_nb_w")

spplot(Variables_Cali, c("me4", "me5"), at=quantile(c(Variables_Cali$me4,
                                                      Variables_Cali$me5),
                                                    p=seq(0,1,length.out = 7)), 
       col.regions=brewer.pal(n=7, "Greys"), 
       main="Spatial Distribution Moran Eigenvectors with relative_nb_b")

ggplot(data = st_as_sf(Variables_Cali)) + geom_sf(aes(fill=me4), color="black",
                                                  size = 0.2) +
        ggtitle("Spatial Distribution Moran Eigenvectors") + 
        theme(panel.background = element_blank(), plot.title = element_text(size=10),
              panel.border = element_rect(fill = NA, size = 0.2)) +
        annotation_scale(location = "br") + 
        annotation_north_arrow(location = "br", which_north = "true",
                               pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                               style = north_arrow_fancy_orienteering) +
        coord_sf(expand = FALSE)

############# Modelos con componentes espaciales ---------------
modelo2 <- glm(Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb + Estrato_moda +
                 fitted(me.fit), offset = log(Variables_Cali$perimetro_m),
               family = poisson, data = Variables_Cali@data)
summary(modelo2) # Mejor?
stargazer(modelo2, no.space = TRUE, single.row=TRUE)
lmtest::coeftest(modelo2, vcov.=vcovHC(modelo2, type="const"))
# Los vectores propios de moran son significativos
lm.morantest(modelo2, listw = relative_nb_w)
pR2(modelo2)
Devianza = ((modelo2$null.deviance - modelo2$deviance)/modelo2$null.deviance)*100
Devianza

# 5% de significancia, mejoro bastante
#stargazer(lm.morantest(modelo2, listw = relative_nb_w))

modelo3 <- glm(Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb + Estrato_moda +
                 fitted(me.fit2), offset = log(Variables_Cali$perimetro_m),
               family = poisson, data = Variables_Cali@data)
summary(modelo3) # Mejor?
lmtest::coeftest(modelo3, vcov.=vcovHC(modelo3, type="const"))
# Igualmente para el modelo 2, los vectores propios son significativos
lm.morantest(modelo3, listw=relative_nb_b)

pR2(modelo2) # Mejor el modelo 2
pR2(modelo3)

AICs<-c(AIC(modelo1), AIC(modelo2), AIC(modelo3))
plot(AICs, type="l", lwd=1.5, xaxt="n", xlab="")
axis(1, at=1:6,labels=F) 
labels<-c("modelo 1","modelo 2","modelo 3")
text(1:6, par("usr")[3]-.25, srt=45, adj=1, labels=labels, xpd=T)
mtext(side=1, text="Model Specification", line=3)
symbols(x= which.min(AICs), y=AICs[which.min(AICs)], circles=1, fg=2,lwd=2,add=T)

############# Auto-models ---------------
Variables_Cali$lag_rate <- lag.listw(x=relative_nb_w, var=Variables_Cali$Covid)
modelo4 <- glm(Covid ~ lag_rate + cat_ics_bajo + Adulto_mayor + cat_comorb + Estrato_moda,
               offset = log(perimetro_m), family = poisson, data = Variables_Cali)
summary(modelo4)
lmtest::coeftest(modelo4, vcov.=vcovHC(modelo4, type="const"))
# No es significativo
lm.morantest(modelo4, listw=relative_nb_b)

############# Auto-models and eigen vectors -------------
modelo5 <- glm(Covid ~ lag_rate + cat_ics_bajo + Adulto_mayor + cat_comorb + Estrato_moda +
                 fitted(me.fit), offset = log(perimetro_m),
               family = poisson, data = Variables_Cali)
summary(modelo5)
lmtest::coeftest(modelo5, vcov.=vcovHC(modelo5, type="const"))

lm.morantest(modelo5, listw=relative_nb_w)


AICs<-c(AIC(modelo1), AIC(modelo3), AIC(modelo4), AIC(modelo2), AIC(modelo5))
plot(AICs, type="l", lwd=1.5, xaxt="n", xlab="")
axis(1, at=1:6,labels=F) 
labels<-c("Modelo 1", "Modelo 3", "Modelo 4", "Modelo 2", "Modelo 5")
text(1:6, par("usr")[3]-.25, srt=45, adj=1, labels=labels, xpd=T)
mtext(side=1, text="Model Specification", line=3)
symbols(x= which.min(AICs), y=AICs[which.min(AICs)], circles=1, fg=2,lwd=2,add=T)

# El mejos modelo es el segundo.

Variables_Cali$resnb <- rstudent(modelo3)
Variables_Cali$resnb1 <- rstudent(modelo1)
Variables_Cali$resclean <- rstudent(modelo2)
spplot(Variables_Cali,c("resnb1"),at=quantile(c(Variables_Cali$resnb1,
                                                Variables_Cali$resnb, 
                                                Variables_Cali$resclean), 
                                              p=seq(0,1,length.out = 7) ), 
       col.regions=brewer.pal(n=7, "Greys"), 
       main="Spatial Distribution Poisson - GLM", names.attr=c("Poisson Regresssion"))

spplot(Variables_Cali,c("resnb"),at=quantile(c(Variables_Cali$resnb1,
                                               Variables_Cali$resnb, 
                                               Variables_Cali$resclean), 
                                             p=seq(0,1,length.out = 7) ), 
       col.regions=brewer.pal(n=7, "Greys"), 
       main="Spatial Distribution Poisson Residuals - Auto count model", 
       names.attr=c("Auto Poisson Regression"))


spplot(Variables_Cali,c("resclean"),at=quantile(c(Variables_Cali$resnb1,
                                                  Variables_Cali$resnb, 
                                                  Variables_Cali$resclean), 
                                                p=seq(0,1,length.out = 7) ), 
       col.regions=brewer.pal(n=7, "Greys"), 
       main="Distribución espacial de los residuales del GLM con los vectores propios de Moran", 
       names.attr=c( "Filtered Poisson"))

### Es mejor el modelo 2