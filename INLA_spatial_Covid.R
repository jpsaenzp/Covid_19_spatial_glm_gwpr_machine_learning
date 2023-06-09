#https://rpubs.com/corey_sparks/607118
# https://rpubs.com/corey_sparks/439277
rm(list=ls())
gc(reset = TRUE)
# Librerías
library(spdep)
library(RColorBrewer)
library(lattice)
library(INLA)
library(tigris)
library(tidycensus)
library(ggplot2)
library(dplyr)
library(plyr)
library(rgdal)
library(readxl)
# Cargue de la base de datos
barrios_Cali <- readOGR('databases_and_shapefiles/Barrios Cali.shp')
datos_barrio_cali <- barrios_Cali@data
datos_cali <- read_excel('databases_and_shapefiles/Cali.xlsx')
datos_cali$Estrato_moda <- as.factor(datos_cali$Estrato_moda)
barrios_Cali@data = left_join(barrios_Cali@data, datos_cali, by = 'ID_BARRIO')
# gráfico de correlación inicial
corrplot::corrplot.mixed(cor(datos_cali[,4:9], use = "complete.obs"))
Variables_Cali <- barrios_Cali[,c('ID_BARRIO','Barrio_Urbanizacion_o_Sector',
                                  'Comuna','Poblacion_total','Covid',
                                  'Adulto_mayor','ind_cobert_serv_bajo',
                                  'Comorbilidad', 'razon_letalidad',
                                  'Proporcion_trabajando','Estrato_moda',
                                  'Centroide_X','Centroide_Y')]
Variables_Cali <- Variables_Cali[-c(65,112,304,326),]
View(Variables_Cali@data)
# Matriz de vecindades
nbs <- knearneigh(coordinates(Variables_Cali), k = 5, longlat = T) #k=5 nearest neighbors
nbs <- knn2nb(nbs, row.names = Variables_Cali$ID_BARRIO, sym = T) #force symmetry!!
plot(nbs, coords= coordinates(Variables_Cali))

mat <- nb2mat(nbs, style="B",zero.policy=TRUE)
colnames(mat) <- rownames(mat) 
mat <- as.matrix(mat[1:dim(mat)[1], 1:dim(mat)[1]])

nb2INLA("cl_graph",nbs)
am_adj <- paste(getwd(),"/cl_graph",sep="")
H <- inla.read.graph(filename="cl_graph")
image(inla.graph2matrix(H), xlab="", ylab="", main="")

Variables_Cali$struct <- 1:dim(Variables_Cali)[1]
Variables_Cali <- st_as_sf(Variables_Cali)
# Grafico del mapa de cali
Variables_Cali$cofips <- paste(Variables_Cali$Comuna, Variables_Cali$ID_BARRIO, sep="")
Variables_Cali %>%
  ggplot()+geom_sf()+coord_sf(crs = "+proj=aea +lat_1=25 +lat_2=50 +lon_0=-100")

#gráficos importantes para observar la distribución de casos de covid
ggplot(data = Variables_Cali) +geom_histogram(aes(x = Covid, y = 0.5*..density..))+ 
  ggtitle(label = "Distribución de Covid", subtitle = "Barrios de Cali")

ggplot(data = Variables_Cali) +geom_histogram(aes(x = log(Covid), y = 0.5*..density..))+ 
  ggtitle(label = "Distribución de Covid", subtitle = "Barrios de Cali")

# Gráfico de cuartiles de covid
Variables_Cali$E_dCovid <- mean(Variables_Cali$Covid)
Variables_Cali %>%
  mutate(qrr = cut(I(Covid/E_dCovid), breaks = quantile(I(Covid/E_dCovid),
                                                        p=seq(0,1,length.out = 5)),
                   include.lowest = T)) %>%
  ggplot() + geom_sf(aes(fill=qrr)) + scale_colour_brewer(palette = "RdBu" ) +
  scale_fill_brewer(palette = "RdBu", na.value="grey") + 
  guides(fill = guide_legend(title="Covid Quartile")) + 
  ggtitle(label="Covid Quartile - Cali, 2020-2021") + 
  coord_sf(crs = "+proj=laea +lat_0=35 +lon_0=-100")

#We can fit these model using the Bayesian framework with INLA.
#Model specification:
f1 <- Covid ~ Adulto_mayor + ind_cobert_serv_bajo + log(Comorbilidad) +
  razon_letalidad + Proporcion_trabajando%in%Estrato_moda

#Inla Model 1 fit
mod1<-inla(formula = f1,data = Variables_Cali, #linear predictor - fixed effects
           family = "poisson", E = E_dCovid,  #marginal distribution for the outcome, expected count
           control.compute = list(waic=T), # compute DIC or not?
           control.predictor = list(link=1), #estimate predicted values & their marginals or not?
           num.threads = 3, 
           verbose = F)
#model summary
summary(mod1)
# Gráfico de observados vs predichos
plot(x= mod1$summary.fitted.values$mean, y=Variables_Cali$Covid/Variables_Cali$E_dCovid, 
     ylab="Observed", xlab="Estimated" )

# Basic county level random intercept model
f2 <- Covid ~ Adulto_mayor + ind_cobert_serv_bajo + log(Comorbilidad) +
  razon_letalidad + Proporcion_trabajando%in%Estrato_moda + f(ID_BARRIO, model = "iid")
# Inla model 2
mod2<-inla(formula = f2,data = Variables_Cali,
           family = "poisson", E = E_dCovid, 
           control.compute = list(waic=T), 
           control.predictor = list(link=1),
           num.threads = 3, 
           verbose = F)

#total model summary
summary(mod2)

#Marginal Distributions of hyperparameters
m2<- inla.tmarginal(
  function(x) (1/x), #invert the precision to be on variance scale
  mod2$marginals.hyperpar$`Precision for ID_BARRIO`)
#95% credible interval for the variance
inla.hpdmarginal(.95, marginal=m2)

plot(m2, type="l", main=c("Posterior distibution for between ID_barrio variance",
                          "- IID model -"), xlim=c(0.035, 0.06))
Variables_Cali$fitted_m2 <- mod2$summary.fitted.values$mean
# Gráfico de covid cuartile
Variables_Cali%>%
  mutate(qrr=cut(fitted_m2, breaks = quantile(fitted_m2, p=seq(0,1,length.out = 6)),
                 include.lowest = T)) %>%
  ggplot()+geom_sf(aes(fill=qrr))+scale_colour_brewer(palette = "RdBu" )+
  scale_fill_brewer(palette = "RdBu", na.value="grey")+
  guides(fill=guide_legend(title="Covid Quartile"))+
  ggtitle(label="Covid Quartile - IID Model")+
  coord_sf(crs = "+proj=laea +lat_0=35 +lon_0=-100")

library(mapview)

map1 <- Variables_Cali %>% 
  mutate(qrr=cut(fitted_m2, breaks = quantile(fitted_m2, p=seq(0,1,length.out = 8))))
clrs <- colorRampPalette(brewer.pal(8, "RdBu"))
mapView(as(map1, "Spatial"), zcol="qrr", legend=T, col.regions=clrs)

### BYM Model - Besag, York, and Mollie Model
#For the BYM model we must specify the spatial connectivity matrix in the 
#random effect.
nbs <- knearneigh(coordinates(barrios_Cali[-c(65,112,304,326),]), k = 5, longlat = T) #k=5 nearest neighbors
nbs <- knn2nb(nbs, row.names = Variables_Cali$struct, sym = T) #force symmetry!!
plot(nbs, coords= coordinates(barrios_Cali[-c(65,112,304,326),]))

mat <- nb2mat(nbs, style="B",zero.policy=TRUE)
colnames(mat) <- rownames(mat) 
mat <- as.matrix(mat[1:dim(mat)[1], 1:dim(mat)[1]])

nb2INLA("cl_graph",nbs)
am_adj <- paste(getwd(),"/cl_graph",sep="")
H <- inla.read.graph(filename="cl_graph")
image(inla.graph2matrix(H), xlab="", ylab="", main="")
# Tercer modelo con modelo Besag, York, and Mollie
f3 <- Covid ~ Adulto_mayor + ind_cobert_serv_bajo + log(Comorbilidad) +
  razon_letalidad + Proporcion_trabajando%in%Estrato_moda +
  f(struct, model = "bym", constr = T, scale.model = T, graph = H)#+
  #f(Comuna, model="iid")
mod3 <- inla(formula = f3,data = Variables_Cali,
             family = "poisson", E = E_dCovid, 
             control.compute = list(waic=T), 
             num.threads = 3, 
             verbose = F,
             control.predictor = list(link=1))

#total model summary
summary(mod3)

# Distribuciones marginales
m3a<- inla.tmarginal(
  function(x) (1/x),
  mod3$marginals.hyperpar$`Precision for struct (iid component)`)
m3b<- inla.tmarginal(
  function(x) (1/x),
  mod3$marginals.hyperpar$`Precision for struct (spatial component)`)

plot(m3a, type="l", main=c("Posterior distibution for between neigh variance", 
                           "- IID model -"), xlim=c(0, 0.07),  ylim=c(0,300))
lines(m3b, col="red")

inla.hpdmarginal(.95,m3a)
inla.hpdmarginal(.95,m3b)

#Esto indica una varianza correlacionada espacialmente baja en estos datos ya que
# es menos a 5%

#Space mapping of the fitted values
Variables_Cali$fitted_m3 <- mod3$summary.fitted.values$mean

Variables_Cali %>%
  mutate(qrr=cut(fitted_m3, breaks = quantile(fitted_m3, p=seq(0,1,length.out = 6)), 
                 include.lowest = T))%>%
  ggplot()+geom_sf(aes(fill=qrr))+scale_colour_brewer(palette = "RdBu" )+
  scale_fill_brewer(palette = "RdBu", na.value="grey")+
  guides(fill=guide_legend(title="Covid Quartile"))+
  ggtitle(label="Covid Quartile - BYM Model")+
  coord_sf(crs = "+proj=laea +lat_0=35 +lon_0=-100")

map1 <- Variables_Cali %>%
  mutate(qrr=cut(fitted_m3, breaks = quantile(fitted_m3, p=seq(0,1,length.out = 8))))
clrs <- colorRampPalette(brewer.pal(8, "RdBu"))
mapView(as(map1, "Spatial"), zcol="qrr", legend=T, col.regions=clrs)

##########Map of spatial random effects
#Es común mapear los efectos aleatorios del modelo BYM para buscar tendencias 
#espaciales, en este caso, no hay señales espaciales fuertes
Variables_Cali$sp_re <- mod3$summary.random$struct$mean[1:334]

Variables_Cali %>%
  mutate(qse=cut(sp_re, breaks = quantile(sp_re, p=seq(0,1,length.out = 6)), 
                 include.lowest = T))%>%
  ggplot()+geom_sf(aes(fill=qse))+scale_colour_brewer(palette = "RdBu" )+
  scale_fill_brewer(palette = "RdBu", na.value="grey")+
  guides(fill=guide_legend(title="Spatial Excess Risk"))+
  ggtitle(label="Spatial Random Effect - BYM Model")+
  coord_sf(crs = "+proj=laea +lat_0=35 +lon_0=-100")

#Exceedence probabilities
thetastar <- 1.5#theta*
inlaprob <- unlist(lapply(mod3$marginals.fitted.values, function(X){
  1-inla.pmarginal(thetastar, X)
}))
hist(inlaprob)

#So, we see lots of occasions where the exceedence probability is greater than .9. 
#We can visualize these in a map
Variables_Cali$exceedprob <- inlaprob

Variables_Cali %>%
  mutate(qrr=cut(exceedprob, breaks = c(0, .5, .8, .9, .95, 1), 
                 include.lowest = T))%>%
  ggplot()+geom_sf(aes(fill=qrr))+scale_colour_brewer(palette = "Blues" )+
  scale_fill_brewer(palette = "Blues", na.value="grey")+
  guides(fill=guide_legend(title=""))+
  ggtitle(label=expression(paste("Exceedence Probability Covid ",
                                 "Pr( ",theta," >1.5"," )") ))+
  coord_sf(crs = "+proj=laea +lat_0=35 +lon_0=-100")

map1 <- Variables_Cali%>%
  mutate(qrr=cut(exceedprob, breaks = c(0, .5, .8, .9, .95, 1),
                 include.lowest = T))
clrs <- colorRampPalette(brewer.pal(6, "Blues"))
mapView(as(map1, "Spatial"), zcol="qrr", legend=T, col.regions=clrs, 
        map.types="OpenStreetMap")

#Which shows several areas where Covid is higher than the municipal cases.



