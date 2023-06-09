# https://zia207.github.io/geospatial-r-github.io/geographically-weighted-poisson-
#regression.html
# https://github.com/hoxo-m/gwpr
# https://rpubs.com/corey_sparks/260983
# http://rstudio-pubs-static.s3.amazonaws.com/483800_14662f9098974b31b20110a0ca345708.
#html025_residuales_espacialmente_correlacionados
rm(list=ls())
gc(reset = T)
# Cargue de librerías
library(GWmodel)      ### GW models
library(sp)           ## Data management
library(spdep)        ## Spatial autocorrelation
library(RColorBrewer) ## Visualization
library(classInt)     ## Class intervals
library(raster)       ## spatial data
library(grid)         # plot
library(gridExtra)    # Multiple plot
library(ggplot2)      # Multiple plot
library(gtable);library(rgdal);library(readxl);library(dplyr);library(xtable)
library(Rcpp);library(tmap);library(viridis);library(ggspatial);library(stargazer)
# Cargue de base de datos
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
# Correlaciones inicialea
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
(l <- length(Variables_Cali))

xy0 = data.frame(x = Variables_Cali$Centroide_X, y = Variables_Cali$Centroide_Y)
coordinates(xy0) <- c('x','y')
#proj4string(xy0) <- CRS("+init=epsg:4326")
#xy0 = spTransform(xy0, CRS("+init=epsg:21897"))

# Mapa de casos de covid
library(mapview)
mapView(Variables_Cali, zcol="Covid", legend=T)

# Bandwidth selection
DM <- gw.dist(dp.locat=coordinates(Variables_Cali))

#Covid ~ Adulto_mayor + ind_cobert_serv_bajo + log(Comorbilidad) + razon_letalidad + 
#Proporcion_trabajando%in%Estrato_moda, family = poisson,data = Variables_Cali@data)

#### Ancho de banda adaptativo
# Primer modelo kernel gaussiano
bw.gwr <- bw.ggwr(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                    Estrato_moda + offset(log(perimetro_m)),
                  data = Variables_Cali,
                  family = "poisson",
                  approach = "AICc",
                  kernel = "gaussian", 
                  adaptive = TRUE,
                  dMat = DM);bw.gwr
#Adaptive bandwidth (number of nearest neighbours): 20 AICc value: 15416.89

# Primer modelo kernel exponencial
bw.gwr <- bw.ggwr(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                    Estrato_moda + offset(log(perimetro_m)),
                  data = Variables_Cali,
                  family = "poisson",
                  approach = "AICc",
                  kernel = "exponential", 
                  adaptive = FALSE,
                  dMat = DM)
bw.gwr #Fixed bandwidth: 0.003830947 AICc value: 6081.316

# SEGUNDO modelo kernel exponencial
bw.gwr <- bw.ggwr(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                    Estrato_moda + offset(log(perimetro_m)),
                  data = Variables_Cali,
                  family = "poisson",
                  approach = "AICc",
                  kernel = "exponential", 
                  adaptive = TRUE,
                  dMat = DM)
bw.gwr #Con el kernel EXPONENCIAL también arrojo 20 vecinos y 14093,18 de AICc

# Fit the model
bgwr.res <- ggwr.basic(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                         Estrato_moda + offset(log(perimetro_m)),
                       data = Variables_Cali,
                       family = "poisson",
                       bw = bw.gwr, 
                       kernel = "gaussian", 
                       adaptive = TRUE,
                       dMat = DM);bgwr.res
#modelo gaussiano bandwidth adaptativo, AICc : 12804.89, Pseudo R-square value:  0.8902721

bgwr.res <- ggwr.basic(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                         Estrato_moda + offset(log(perimetro_m)),
                       data = Variables_Cali,
                       family = "poisson",
                       bw = bw.gwr, 
                       kernel = "exponential", 
                       adaptive = FALSE,
                       dMat = DM);bgwr.res
#modelo exponencial bandwidth fijo, AICc : 3467.999, Pseudo R-square value:  0.9865585

bgwr.res <- ggwr.basic(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                         Estrato_moda + offset(log(perimetro_m)),
                       data = Variables_Cali,
                       family = "poisson",
                       bw = bw.gwr, 
                       kernel = "exponential", 
                       adaptive = TRUE,
                       dMat = DM);bgwr.res
#modelo exponencial bandwidth adaptativo, AICc: 11481.19, Pseudo R-square value:  0.9019524 

bgwr.res ##############################################################################3
stargazer(bgwr.res)
bgwr.res$GW.diagnostic$gw.deviance
###s Save the summary output
#capture.output(print(bgwr.res),file="summary_GWRP.doc")

# Modelos por estrato
mode_1_2 <- ggwr.basic(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                         offset(log(perimetro_m)),
                       data = subset(Variables_Cali, Estrato_moda %in% c("1","2")),
                       family = "poisson",
                       bw = bw.ggwr(formula = Covid ~ cat_ics_bajo + Adulto_mayor +
                                      cat_comorb + offset(log(perimetro_m)),
                                    data = subset(Variables_Cali,
                                                  Estrato_moda %in% c("1","2")),
                                    family = "poisson",
                                    approach = "AICc",
                                    kernel = "exponential", 
                                    adaptive = TRUE,
                                    dMat = gw.dist(dp.locat=coordinates(subset(Variables_Cali,
                                                                               Estrato_moda %in% c("1","2"))))), 
                       kernel = "exponential", 
                       adaptive = TRUE,
                       dMat = gw.dist(dp.locat=coordinates(subset(Variables_Cali,
                                                                  Estrato_moda %in% c("1","2")))));mode_1_2


mode_3_4 <- ggwr.basic(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                         offset(log(perimetro_m)),
                       data = subset(Variables_Cali, Estrato_moda %in% c("3","4")),
                       family = "poisson",
                       bw = bw.ggwr(formula = Covid ~ cat_ics_bajo + Adulto_mayor +
                                      cat_comorb + offset(log(perimetro_m)),
                                    data = subset(Variables_Cali,
                                                  Estrato_moda %in% c("3","4")),
                                    family = "poisson",
                                    approach = "AICc",
                                    kernel = "exponential", 
                                    adaptive = TRUE,
                                    dMat = gw.dist(dp.locat=coordinates(subset(Variables_Cali,
                                                                               Estrato_moda %in% c("3","4"))))), 
                       kernel = "exponential", 
                       adaptive = TRUE,
                       dMat = gw.dist(dp.locat=coordinates(subset(Variables_Cali,
                                                                  Estrato_moda %in% c("3","4")))));mode_3_4


mode_5_6 <- ggwr.basic(formula = Covid ~ cat_ics_bajo + Adulto_mayor + cat_comorb +
                         offset(log(perimetro_m)),
                       data = subset(Variables_Cali, Estrato_moda %in% c("5","6")),
                       family = "poisson",
                       bw = bw.ggwr(formula = Covid ~ cat_ics_bajo + Adulto_mayor +
                                      cat_comorb + offset(log(perimetro_m)),
                                    data = subset(Variables_Cali,
                                                  Estrato_moda %in% c("5","6")),
                                    family = "poisson",
                                    approach = "AICc",
                                    kernel = "exponential", 
                                    adaptive = TRUE,
                                    dMat = gw.dist(dp.locat=coordinates(subset(Variables_Cali,
                                                                               Estrato_moda %in% c("5","6"))))), 
                       kernel = "exponential", 
                       adaptive = TRUE,
                       dMat = gw.dist(dp.locat=coordinates(subset(Variables_Cali,
                                                                  Estrato_moda %in% c("5","6")))));mode_5_6

mode_1_2
mode_3_4
mode_5_6

#Extract GWPR results
### Create spatial data frame
Variables_Cali@data$y <- bgwr.res$SDF$y
Variables_Cali@data$yhat <- bgwr.res$SDF$yhat
Variables_Cali@data$residual <- bgwr.res$SDF$residual

rsd = sd(Variables_Cali@data$residual)
# Agregar Residuales a la base de datos para observarlos
Variables_Cali@data$stdRes <- (Variables_Cali@data$residual)/sd(Variables_Cali@data$residual)
Variables_Cali@data$LLN = Variables_Cali@data$yhat - 1.645*rsd
Variables_Cali@data$ULN = Variables_Cali@data$yhat + 1.645*rsd

# Intercepto Coeficientes estimados
Variables_Cali@data$Intercept <- bgwr.res$SDF$Intercept
Variables_Cali@data$est_Adultomayor <- bgwr.res$SDF$Adulto_mayor
Variables_Cali@data$est_indcobertservbajo2 <- bgwr.res$SDF$cat_ics_bajo2
Variables_Cali@data$est_indcobertservbajo3 <- bgwr.res$SDF$cat_ics_bajo3
Variables_Cali@data$est_comorbilidad2 <- bgwr.res$SDF$cat_comorb2
Variables_Cali@data$est_comorbilidad3 <- bgwr.res$SDF$cat_comorb3
Variables_Cali@data$est_estrato2 <- bgwr.res$SDF$Estrato_moda2
Variables_Cali@data$est_estrato3 <- bgwr.res$SDF$Estrato_moda3
Variables_Cali@data$est_estrato4 <- bgwr.res$SDF$Estrato_moda4
Variables_Cali@data$est_estrato5 <- bgwr.res$SDF$Estrato_moda5
Variables_Cali@data$est_estrato6 <- bgwr.res$SDF$Estrato_moda6

# T-values
Variables_Cali@data$t_Intercept <- bgwr.res$SDF$Intercept_TV
Variables_Cali@data$t_Adultomayor <- bgwr.res$SDF$Adulto_mayor_TV
Variables_Cali@data$t_indcobertservbajo2 <- bgwr.res$SDF$cat_ics_bajo2_TV
Variables_Cali@data$t_indcobertservbajo3 <- bgwr.res$SDF$cat_ics_bajo3_TV
Variables_Cali@data$t_comorbilidad2 <- bgwr.res$SDF$cat_comorb2_TV
Variables_Cali@data$t_comorbilidad3 <- bgwr.res$SDF$cat_comorb3_TV
Variables_Cali@data$t_estrato2 <- bgwr.res$SDF$Estrato_moda2_TV
Variables_Cali@data$t_estrato3 <- bgwr.res$SDF$Estrato_moda3_TV
Variables_Cali@data$t_estrato4 <- bgwr.res$SDF$Estrato_moda4_TV
Variables_Cali@data$t_estrato5 <- bgwr.res$SDF$Estrato_moda5_TV
Variables_Cali@data$t_estrato6 <- bgwr.res$SDF$Estrato_moda6_TV

# Calculate pseudo-t values
Variables_Cali@data$p_Adultomayor <- 2*pt(-abs(bgwr.res$SDF$Adulto_mayor_TV), df=334*5)
Variables_Cali@data$p_indcobertservbajo2 <- 2*pt(-abs(bgwr.res$SDF$cat_ics_bajo2_TV),
                                                 df=334*5)
Variables_Cali@data$p_indcobertservbajo3 <- 2*pt(-abs(bgwr.res$SDF$cat_ics_bajo3_TV),
                                                 df=334*5)
Variables_Cali@data$p_comorbilidad2 <- 2*pt(-abs(bgwr.res$SDF$cat_comorb2_TV),
                                            df=334*5)
Variables_Cali@data$p_comorbilidad3 <- 2*pt(-abs(bgwr.res$SDF$cat_comorb3_TV),
                                            df=334*5)
Variables_Cali@data$p_estrato2 <- 2*pt(-abs(bgwr.res$SDF$Estrato_moda2_TV),
                                       df=334*5)
Variables_Cali@data$p_estrato3 <- 2*pt(-abs(bgwr.res$SDF$Estrato_moda3_TV),
                                       df=334*5)
Variables_Cali@data$p_estrato4 <- 2*pt(-abs(bgwr.res$SDF$Estrato_moda4_TV),
                                       df=334*5)
Variables_Cali@data$p_estrato5 <- 2*pt(-abs(bgwr.res$SDF$Estrato_moda5_TV),
                                       df=334*5)
Variables_Cali@data$p_estrato6 <- 2*pt(-abs(bgwr.res$SDF$Estrato_moda6_TV),
                                       df=334*5)

# Significancia de las variables por barrio
Variables_Cali$sig_Adultomayor <-ifelse(Variables_Cali@data$est_Adultomayor > 0 &
                                          Variables_Cali@data$p_Adultomayor <= 0.05 , 1, 0)
Variables_Cali$sig_indcobertservbajo2 <-ifelse(Variables_Cali@data$est_indcobertservbajo2 > 0 &
                                                 Variables_Cali@data$p_indcobertservbajo2 <= 0.05 , 1, 0)
Variables_Cali$sig_indcobertservbajo3 <-ifelse(Variables_Cali@data$est_indcobertservbajo3 > 0 &
                                                 Variables_Cali@data$p_indcobertservbajo3 <= 0.05 , 1, 0)
Variables_Cali$sig_comorbilidad2 <-ifelse(Variables_Cali@data$est_comorbilidad2 > 0 &
                                            Variables_Cali@data$p_comorbilidad2 <= 0.05 , 1, 0)
Variables_Cali$sig_comorbilidad3 <-ifelse(Variables_Cali@data$est_comorbilidad3 > 0 &
                                            Variables_Cali@data$p_comorbilidad3 <= 0.05 , 1, 0)
Variables_Cali$sig_estrato2 <-ifelse(Variables_Cali@data$est_estrato2 > 0 &
                                                Variables_Cali@data$p_estrato2 <= 0.05 , 1, 0)
Variables_Cali$sig_estrato3 <-ifelse(Variables_Cali@data$est_estrato3 > 0 &
                                                Variables_Cali@data$p_estrato3 <= 0.05 , 1, 0)
Variables_Cali$sig_estrato4 <-ifelse(Variables_Cali@data$est_estrato4 > 0 &
                                                Variables_Cali@data$p_estrato4 <= 0.05 , 1, 0)
Variables_Cali$sig_estrato5 <-ifelse(Variables_Cali@data$est_estrato5 > 0 &
                                                Variables_Cali@data$p_estrato5 <= 0.05 , 1, 0)
Variables_Cali$sig_estrato6 <-ifelse(Variables_Cali@data$est_estrato6 > 0 &
                                                Variables_Cali@data$p_estrato6 <= 0.05 , 1, 0)

# Plot GWRP Statistics
polys <- list("sp.lines", as(Variables_Cali, "SpatialLines"), col="grey", 
              lwd=.8,lty=1)
col.palette <- colorRampPalette(c("red", "orange", "light blue", "dark blue"),space="rgb",
                                interpolate = "linear")

# Plot Local Estimates
# Estimaciones locales para adultos mayores
options(scipen=999)
col.palette <- colorRampPalette(c("dark red", "red", "orange", "beige", "light blue",
                                  "dark blue"),space="rgb",
                                interpolate = "linear") 
Est_Adultomayor <- spplot(Variables_Cali,"est_Adultomayor", main = "Adulto Mayor",
                          sp.layout=list(polys),
                          col="transparent",
                          col.regions=col.palette(100))

Est_Adultomayor <- ggplot(data = st_as_sf(Variables_Cali)) +
  geom_sf(aes(fill=est_Adultomayor), color="black", size = 0.2) +
  ggtitle("Adulto Mayor") +
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        legend.title=element_blank(), axis.text.x = element_text(size = 5),
        panel.border = element_rect(fill = NA, linewidth = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_gradientn(colours = col.palette(100))

# Estimaciones locales para cobertura de servicios públicos en estratos bajos

Est_indcobertservbajo2 <- ggplot(data = st_as_sf(Variables_Cali)) +
  geom_sf(aes(fill=est_indcobertservbajo2), color="black", size = 0.2) +
  ggtitle("Cobertura de servicios estratos bajos 2") +
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        legend.title=element_blank(), axis.text.x = element_text(size = 5),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_gradientn(colours = col.palette(100))

Est_indcobertservbajo3 <- ggplot(data = st_as_sf(Variables_Cali)) +
  geom_sf(aes(fill=est_indcobertservbajo3), color="black", size = 0.2) +
  ggtitle("Cobertura de servicios estratos bajos 3") +
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        legend.title=element_blank(), axis.text.x = element_text(size = 5),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_gradientn(colours = col.palette(100))


# Estimaciones locales para comorbilidades
Est_comorbilidad2 <- ggplot(data = st_as_sf(Variables_Cali)) +
  geom_sf(aes(fill=est_comorbilidad2), color="black", size = 0.2) +
  ggtitle("Comorbilidades 2") +
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        legend.title=element_blank(), axis.text.x = element_text(size = 5),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_gradientn(colours = col.palette(100))

Est_comorbilidad3 <- ggplot(data = st_as_sf(Variables_Cali)) +
  geom_sf(aes(fill=est_comorbilidad3), color="black", size = 0.2) +
  ggtitle("Comorbilidades 3") +
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        legend.title=element_blank(), axis.text.x = element_text(size = 5),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_gradientn(colours = col.palette(100))

# mapas juntos para las variables
grid.arrange(Est_Adultomayor, Est_indcobertservbajo2, Est_indcobertservbajo3,
             Est_comorbilidad2, Est_comorbilidad3,
             ncol = 3, nrow = 2, heights = c(40,6), top = textGrob("Local Estimates",
                                                        gp=gpar(fontsize=25)))

ggpubr::ggarrange(Est_Adultomayor, Est_indcobertservbajo2, Est_indcobertservbajo3,
                  Est_comorbilidad2, Est_comorbilidad3, ncol = 2,nrow = 3)

# Gráfico de las estimaciones para estratos
Est_estrato2 <- ggplot(data = st_as_sf(Variables_Cali)) +
  geom_sf(aes(fill=est_estrato2), color="black", size = 0.2) +
  ggtitle("Estrato 2") +
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        legend.title=element_blank(), axis.text.x = element_text(size = 6),
        panel.border = element_rect(fill = NA, size = 0.2),
        legend.key.size = unit(1, 'cm')) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_gradientn(colours = col.palette(100))

Est_estrato3 <- ggplot(data = st_as_sf(Variables_Cali)) +
  geom_sf(aes(fill=est_estrato3), color="black", size = 0.2) +
  ggtitle("Estrato 3") +
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        legend.title=element_blank(), axis.text.x = element_text(size = 6),
        panel.border = element_rect(fill = NA, size = 0.2),
        legend.key.size = unit(1, 'cm')) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_gradientn(colours = col.palette(100))

Est_estrato4 <- ggplot(data = st_as_sf(Variables_Cali)) +
  geom_sf(aes(fill=est_estrato4), color="black", size = 0.2) +
  ggtitle("Estrato 4") +
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        legend.title=element_blank(), axis.text.x = element_text(size = 6),
        panel.border = element_rect(fill = NA, size = 0.2),
        legend.key.size = unit(1, 'cm')) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_gradientn(colours = col.palette(100))

Est_estrato5 <- ggplot(data = st_as_sf(Variables_Cali)) +
  geom_sf(aes(fill=est_estrato5), color="black", size = 0.2) +
  ggtitle("Estrato 5") +
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        legend.title=element_blank(), axis.text.x = element_text(size = 6),
        panel.border = element_rect(fill = NA, size = 0.2),
        legend.key.size = unit(1, 'cm')) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_gradientn(colours = col.palette(100))

Est_estrato6 <- ggplot(data = st_as_sf(Variables_Cali)) +
  geom_sf(aes(fill=est_estrato6), color="black", size = 0.2) +
  ggtitle("Estrato 6") +
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        legend.title=element_blank(), axis.text.x = element_text(size = 6),
        panel.border = element_rect(fill = NA, size = 0.2),
        legend.key.size = unit(1, 'cm')) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_gradientn(colours = col.palette(100))

# mapas juntos para las variables de estratos
grid.arrange(Est_estrato2, Est_estrato3, Est_estrato4, Est_estrato5, Est_estrato6,
             ncol = 3, nrow = 2, heights = c(30,6),
             top = textGrob("Estratos",
                            gp=gpar(fontsize=25)))

ggpubr::ggarrange(Est_estrato2, Est_estrato3, ncol = 2, nrow = 1)

ggpubr::ggarrange(Est_estrato4, Est_estrato5, Est_estrato6,
                  ncol = 3, nrow = 1)
#Plot Local t-values

col.palette.t<-colorRampPalette(c("dark red", "red", "orange", "beige", "light blue",
                                  "dark blue"),space="rgb", interpolate = "linear") 

T_Adultomayor <- spplot(Variables_Cali,"t_Adultomayor", main = "Adulto mayor",
                        sp.layout=list(polys),
                        col="transparent",
                        col.regions=rev(col.palette.t(100)))

T_indcobertservbajo2 <- spplot(Variables_Cali,"t_indcobertservbajo2", main = "ICS bajo 2", 
                               sp.layout=list(polys),
                               col="transparent",
                               col.regions=rev(col.palette.t(100)))

T_indcobertservbajo3 <- spplot(Variables_Cali,"t_indcobertservbajo3", main = "ICS bajo 3", 
                               sp.layout=list(polys),
                               col="transparent",
                               col.regions=rev(col.palette.t(100)))

T_comorbilidad2 <- spplot(Variables_Cali,"t_comorbilidad2", main = "Comorbilidad 2", 
                          sp.layout=list(polys),
                          col="transparent",
                          col.regions=rev(col.palette.t(100)))

T_comorbilidad3 <- spplot(Variables_Cali,"t_comorbilidad3", main = "Comorbilidad 3", 
                         sp.layout=list(polys),
                         col="transparent",
                         col.regions=rev(col.palette.t(100)))

# mapas juntos para las 4 variables
grid.arrange(T_Adultomayor, T_indcobertservbajo2, T_indcobertservbajo3,
             T_comorbilidad2, T_comorbilidad3,
             ncol=3, nrow=2, heights = c(30,6), top = textGrob("Local t-values",
                                                               gp=gpar(fontsize=25)))

# Para de valores t para los estratos
T_estrato2 <- spplot(Variables_Cali,"t_estrato2", main = "Estrato 2",
                     sp.layout=list(polys), col="transparent",
                     col.regions=rev(col.palette.t(100)))
T_estrato3 <- spplot(Variables_Cali,"t_estrato3", main = "Estrato 3",
                     sp.layout=list(polys), col="transparent",
                     col.regions=rev(col.palette.t(100)))
T_estrato4 <- spplot(Variables_Cali,"t_estrato4", main = "Estrato 4",
                     sp.layout=list(polys), col="transparent",
                     col.regions=rev(col.palette.t(100)))
T_estrato5 <- spplot(Variables_Cali,"t_estrato5", main = "Estrato 5",
                     sp.layout=list(polys), col="transparent",
                     col.regions=rev(col.palette.t(100)))
T_estrato6 <- spplot(Variables_Cali,"t_estrato6", main = "Estrato 6",
                     sp.layout=list(polys), col="transparent",
                     col.regions=rev(col.palette.t(100)))

# mapas juntos para las variables de estratos
grid.arrange(T_estrato2, T_estrato3, T_estrato4, T_estrato5, T_estrato6,
             ncol=3, nrow=2, heights = c(30,6), top = textGrob("Local t-values",
                                                               gp=gpar(fontsize=25)))

#Plot Std-Residuals
myPaletteRes <- colorRampPalette(c("orange", "yellow", "beige",
                                   "light blue", "dark blue"))
std_res <- spplot(Variables_Cali,"stdRes", main = "GWPR Std. Residuals",
                  sp.layout=list(polys),
                  col="transparent",
                  col.regions=myPaletteRes(100))

std_res <- ggplot(data = st_as_sf(Variables_Cali)) +
  geom_sf(aes(fill=stdRes), color="black", size = 0.2) +
  ggtitle("GWPR Std. Residuals") +
  theme(panel.background = element_blank(), plot.title = element_text(size=15),
        legend.title=element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_gradientn(colours = myPaletteRes(100))

#windows(width=4, height=3.5)
#tiff( file="FIG_GWRP_Std_Residuals.tif", 
#      width=4, height=3.5,units = "in", pointsize = 12, res=1600,
#      restoreConsole = T,bg="transparent")
std_res


########################
legend_title = expression("Residual")
map_fbgwr1 = tm_shape(Variables_Cali) +
  tm_fill(col = "residual", title = legend_title, palette = viridis(256), style = "cont") + # add fill
  tm_borders(col = "white", lwd = .1)  + # add borders
  tm_compass(type = "arrow", position = c("left", "bottom") , size = 2) + # add compass
  tm_scale_bar(breaks = c(0,1,2), text.size = 0.7, position =  c("left", "bottom")) + # add scale bar
  tm_layout(bg.color = "white") # change background colour
map_fbgwr1 + tm_shape(barrios_Cali) + # add region boundaries
  tm_borders(col = "white", lwd = .5) # add borders

##### Interpretacion
bgwr.res

legend_title = expression("Adultos mayores")
map_abgwr2 = tm_shape(Variables_Cali) +
  tm_fill(col = "est_Adultomayor", title = legend_title, palette = viridis(256), style = "cont") + # add fill
  tm_borders(col = "white", lwd = .1)  + # add borders
  tm_compass(type = "arrow", position = c("left", "bottom") , size = 2) + # add compass
  tm_scale_bar(breaks = c(0,1,2), text.size = 0.7, position =  c("left", "bottom")) + # add scale bar
  tm_layout(bg.color = "white") # change background colour
#map_abgwr2 = map_abgwr2 + tm_shape(barrios_Cali) + # add region boundaries
#  tm_borders(col = "white", lwd = .5) # add borders
map_abgwr2

#### significancia. Categorizar el t-value
Variables_Cali$t_Adultomayor_cat <- cut(Variables_Cali$t_Adultomayor,
                                        breaks=c(min(Variables_Cali$t_Adultomayor),
                                                 -1.96, 1.96,
                                                 max(Variables_Cali$t_Adultomayor)),
                                        labels=c("Sig","No sig", "Sig"))

Variables_Cali$t_indcobertservbajo2_cat <- cut(Variables_Cali$t_indcobertservbajo2,
                                               breaks=c(min(Variables_Cali$t_indcobertservbajo2),
                                                        -1.96, 1.96,
                                                        max(Variables_Cali$t_indcobertservbajo2)),
                                               labels=c("Sig","No sig", "Sig"))

Variables_Cali$t_indcobertservbajo3_cat <- cut(Variables_Cali$t_indcobertservbajo3,
                                               breaks=c(min(Variables_Cali$t_indcobertservbajo3),
                                                        -1.96, 1.96,
                                                        max(Variables_Cali$t_indcobertservbajo3)),
                                               labels=c("Sig","No sig", "Sig"))

Variables_Cali$t_comorbilidad2_cat <- cut(Variables_Cali$t_comorbilidad2,
                                          breaks=c(min(Variables_Cali$t_comorbilidad2),
                                                   -1.96, 1.96,
                                                   max(Variables_Cali$t_comorbilidad2)),
                                          labels=c("Sig","No sig", "Sig"))

Variables_Cali$t_comorbilidad3_cat <- cut(Variables_Cali$t_comorbilidad3,
                                          breaks=c(min(Variables_Cali$t_comorbilidad3),
                                                   -1.96, 1.96,
                                                   max(Variables_Cali$t_comorbilidad3)),
                                          labels=c("Sig","No sig", "Sig"))
table(Variables_Cali$t_Adultomayor_cat)
table(Variables_Cali$t_indcobertservbajo2_cat)
table(Variables_Cali$t_indcobertservbajo3)

# mapas juntos para la significancia de las palabras más importantes
g1 <- ggplot() + geom_sf(data=sf::st_as_sf(Variables_Cali), aes(fill=t_Adultomayor_cat),
                         color="black", size = 0.2) +
  scale_fill_manual(values = c("beige", "gray")) + 
  ggtitle("Sig. de Adulto mayor") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        #axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 5), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + coord_sf(expand = FALSE)

g2 <- ggplot() + geom_sf(data=sf::st_as_sf(Variables_Cali), aes(fill=t_indcobertservbajo2_cat),
                   color="black", size = 0.2) +
  scale_fill_manual(values = c("beige", "gray")) + 
  ggtitle("Sig. de cobertura de servicios bajo 2") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        #axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 5), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + coord_sf(expand = FALSE)

g3 <- ggplot() + geom_sf(data=sf::st_as_sf(Variables_Cali), aes(fill=t_indcobertservbajo3_cat),
                         color="black", size = 0.2) +
  scale_fill_manual(values = c("beige", "gray")) + 
  ggtitle("Sig. de cobertura de servicios bajo 3") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        #axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 5), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + coord_sf(expand = FALSE)

g4 <- ggplot() + geom_sf(data=sf::st_as_sf(Variables_Cali), aes(fill=t_comorbilidad2_cat),
                         color="black", size = 0.2) +
  scale_fill_manual(values = c("beige", "gray")) + 
  ggtitle("Sig. de Comorbilidad 2") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        #axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 5), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + coord_sf(expand = FALSE)

g5 <- ggplot() + geom_sf(data=sf::st_as_sf(Variables_Cali), aes(fill=t_comorbilidad3_cat),
                         color="black", size = 0.2) +
  scale_fill_manual(values = c("beige", "gray")) + 
  ggtitle("Sig. de Comorbilidad 3") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        #axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 5), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + coord_sf(expand = FALSE)

ggpubr::ggarrange(g1, g2, g3, g4, g5, common.legend = TRUE, legend="bottom", ncol = 3,
                  nrow = 2)
#
# mapas juntos para la significancia del eos estratos anidados con los trabajadores
Variables_Cali$t_estrato2_cat <- cut(Variables_Cali$t_estrato2,
                                     breaks=c(min(Variables_Cali$t_estrato2),
                                              -1.96, 1.96,
                                              max(Variables_Cali$t_estrato2)),
                                     labels=c("Sig","No sig", "Sig"))
table(Variables_Cali$t_estrato2_cat)

Variables_Cali$t_estrato3_cat <- cut(Variables_Cali$t_estrato3,
                                     breaks=c(min(Variables_Cali$t_estrato3),
                                              -1.96, 1.96,
                                              max(Variables_Cali$t_estrato3)),
                                     labels=c("Sig","No sig", "Sig"))
table(Variables_Cali$t_estrato3_cat)

Variables_Cali$t_estrato4_cat <- cut(Variables_Cali$t_estrato4,
                                     breaks=c(min(Variables_Cali$t_estrato4),
                                              -1.96, 1.96,
                                              max(Variables_Cali$t_estrato4)),
                                     labels=c("Sig","No sig", "Sig"))
table(Variables_Cali$t_estrato4_cat)

Variables_Cali$t_estrato5_cat <- cut(Variables_Cali$t_estrato5,
                                     breaks=c(min(Variables_Cali$t_estrato5),
                                              -1.96, 1.96,
                                              max(Variables_Cali$t_estrato5)),
                                     labels=c("Sig","No sig", "Sig"))
table(Variables_Cali$t_estrato5_cat)

Variables_Cali$t_estrato6_cat <- cut(Variables_Cali$t_estrato6,
                                     breaks=c(min(Variables_Cali$t_estrato6),
                                              -1.96, 1.96,
                                              max(Variables_Cali$t_estrato6)),
                                     labels=c("Sig","No sig", "Sig"))
table(Variables_Cali$t_estrato6_cat)

g6 <- ggplot() + geom_sf(data=sf::st_as_sf(Variables_Cali), aes(fill=t_estrato2_cat),
                         color="black", size = 0.2) +
  scale_fill_manual(values = c("beige", "gray")) + 
  ggtitle("Estrato 2") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + coord_sf(expand = FALSE)

g7 <- ggplot() + geom_sf(data=sf::st_as_sf(Variables_Cali), aes(fill=t_estrato3_cat),
                         color="black", size = 0.2) +
  scale_fill_manual(values = c("beige", "gray")) + 
  ggtitle("Estrato 3") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + coord_sf(expand = FALSE)

g8 <- ggplot() + geom_sf(data=sf::st_as_sf(Variables_Cali), aes(fill=t_estrato4_cat),
                         color="black", size = 0.2) +
  scale_fill_manual(values = c("beige", "gray")) + 
  ggtitle("Estrato 4") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + coord_sf(expand = FALSE)

g9 <- ggplot() + geom_sf(data=sf::st_as_sf(Variables_Cali), aes(fill=t_estrato5_cat),
                         color="black", size = 0.2) +
  scale_fill_manual(values = c("beige", "gray")) + 
  ggtitle("Estrato 5") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + coord_sf(expand = FALSE)

g10 <- ggplot() + geom_sf(data=sf::st_as_sf(Variables_Cali), aes(fill=t_estrato6_cat),
                          color="black", size = 0.2) +
  scale_fill_manual(values = c("beige", "gray")) + 
  ggtitle("Estrato 6") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + coord_sf(expand = FALSE)

ggpubr::ggarrange(g6, g7, common.legend = TRUE, legend="bottom")

ggpubr::ggarrange(g8, g9, g10, common.legend = TRUE,
                  legend="bottom", ncol = 3, nrow = 1)
