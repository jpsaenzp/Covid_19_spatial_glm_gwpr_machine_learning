############# Carga de librer?as ------------------
rm(list=ls())
gc(reset = T)
library(readxl); library(maptools); library(tmap); library(dplyr)
library(rgdal); library(lme4); library(raster); library(caret); #library(mapview)
library(sf); library(CAST); library(fasterize); library(randomForest)
library(leaflet); library(spdep); library(maptools); library(spatialreg)
require(RColorBrewer); library(texreg); library(stargazer);
library(bartMachine)
library(Metrics); library(xgboost); library(Metrics); library(kernlab)
library(ggspatial); library(kernlab); library(ggplot2)
############# Definir directorio y cargar datos ----------------
barrios_Cali <- readOGR('databases_and_shapefiles/Barrios Cali.shp')
datos_barrio_cali <- barrios_Cali@data
datos_cali <- read_excel('databases_and_shapefiles/Covid_19_espacial_machine_learning.xlsx') #Cali_nueva_base.xlsx
barrios_Cali@data = left_join(barrios_Cali@data, datos_cali, by = 'ID_BARRIO')

############# Analisis de mapas --------------------
nas_columna <- apply(X = is.na(barrios_Cali@data), MARGIN = 2, FUN = sum)
View(nas_columna) # Datos faltantes por columna

spplot(obj = barrios_Cali[,77], main = colnames(barrios_Cali@data[77]))
spplot(obj = barrios_Cali[,73], main = colnames(barrios_Cali@data[73]))
spplot(obj = barrios_Cali[,74], main = colnames(barrios_Cali@data[74]))
spplot(obj = barrios_Cali[,70], main = colnames(barrios_Cali@data[70]))
spplot(obj = barrios_Cali[,69], main = colnames(barrios_Cali@data[69]))
spplot(obj = barrios_Cali[,66], main = colnames(barrios_Cali@data[66]))
spplot(obj = barrios_Cali[,65], main = colnames(barrios_Cali@data[65]))
spplot(obj = barrios_Cali[,32], main = colnames(barrios_Cali@data[32]))
spplot(obj = barrios_Cali[,19], main = colnames(barrios_Cali@data[19]))
spplot(obj = barrios_Cali[,29], main = colnames(barrios_Cali@data[29]))
spplot(obj = barrios_Cali[,16], main = colnames(barrios_Cali@data[16]))

### Graficos 
## Grafico de estrato por barrio en Cali
barrios_Cali1 <- sf::st_as_sf(barrios_Cali)

ggplot() + geom_sf(data = barrios_Cali1, aes(fill = as.factor(Estrato_moda)), 
                   color = "grey85") + theme(text = element_text(size = 10))+
  ggtitle ("Estratos en Cali por barrio") + labs(fill = "Estrato")

###### graficar el mapa 
tm_shape(barrios_Cali) + tm_fill("Covid", title=c("Contagios por Covid-19"), n=7, 
                                 palette="viridis") + 
  tm_borders(lwd=0.0) + tm_layout(legend.text.size = 0.6, 
                                  legend.position = c(0.7,0.02), 
                                  panel.labels="Covid-19 por barrios en Cali",
                                  bg="grey90")

library(leaflet)
paleta <- colorRampPalette(c('red', 'orange', 'green'))(100)
q_pal <- colorNumeric(paleta, barrios_Cali$Covid)
#q_pal <- colorNumeric("OrRd", barrios_Cali$Covid)

labels <- sprintf(
  "<strong>%s</strong><br/>%g casos",
  barrios_Cali$NOMBRE.x, barrios_Cali$Covid
) %>% lapply(htmltools::HTML)

leaflet(barrios_Cali) %>% addTiles() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(color = '#444444', weight = 1, smoothFactor = 0.5, opacity = 1.0,
              fillOpacity = 0.5, fillColor = q_pal(barrios_Cali$Covid), 
              highlightOptions = highlightOptions(color = "white", weight = 2, 
                                                  bringToFront = TRUE), 
              label = labels, labelOptions = labelOptions(
                style = list("font-weight" = "normal", padding = "3px 8px"),
                textsize = "15px", direction = "auto")) %>% 
  addLegend(pal = q_pal, values = ~Covid)
#https://rstudio.github.io/leaflet/choropleths.html


############# An?lisis de los datos --------------------
custom_corr_plot <- function(variable1, variable2, df, alpha=0.3){
  p <- df %>%
    mutate(
      # Truco para que se ponga el t?tulo estilo facet
      title = paste(toupper(variable2), "vs", toupper(variable1))
    ) %>%
    ggplot(aes(x = !!sym(variable1), y = !!sym(variable2))) + 
    geom_point(alpha = alpha) +
    # Tendencia no lineal
    geom_smooth(se = FALSE, method = "gam", formula =  y ~ splines::bs(x, 3)) +
    # Tendencia lineal
    geom_smooth(se = FALSE, method = "lm", color = "firebrick") +
    facet_grid(. ~ title) +
    theme_bw() +
    theme(strip.text = element_text(colour = "black", size = 8, face = 2),
          axis.title = element_blank())
  return(p)
}
variables_continuas <- colnames(datos_cali[,6:73])

plots <- purrr::map(
  .x = variables_continuas,
  .f = custom_corr_plot,
  variable2 = "Covid",
  df = datos_cali[,6:73]
)

plots

#ggpubr::ggarrange(plotlist = plots, ncol = 3, nrow = 2) %>%
#  ggpubr::annotate_figure(
#    top = ggpubr::text_grob("Comportamiento variables con Covid", 
#                            face = "bold", size = 16, x = 0.20)
#  )
#

############# Planteamiento de modelos con base a la informaci?n de Python y R, correlaciones --------------
datos_cali$Estrato_moda <- as.factor(datos_cali$Estrato_moda)

modelo2 <- glm(formula = Covid ~ ., family = poisson, data = na.omit(datos_cali[,c(6:12,25:27,44:73)]))
step(modelo2, ~1, data = datos, direction = "both")

summary(glm(formula = Covid ~ `TGP = PEA/PET` + `TGI = PEI/PET` + `TD = Desocupados/PEA` + 
              `TO = Ocupados/PET` + Estrato_moda + Poblacion_total + Hogares + 
              Adulto_mayor + Proporcion_etnia + Proporcion_ngr_afro_pal_raiz + 
              Proporc_personas_dificultades + Asistencia_esc_5 + Asistencia_esc_6a10 + 
              Asistencia_esc_17a21 + Asistencia_esc_5a16 + Asistencia_esc_5a21 + 
              Cobertura_Gas_natural + Cobertura_Recoleccion_basura + `Cobertura_Internet_ fijoomovil` + 
              Proporcion_buscando_trabajo + Proporcion_jubilados + Hacinam_cuarto + 
              Hacinam_dormit + ind_cobert_serv_bajo + Pobre_nbi + Tasa_ind_cobert_serv_bajo + 
              Tasa_Pobre_nbi + Comorbilidad + Muertos + Recuperados + razon_letalidad, 
            family = poisson, data = na.omit(datos_cali[, c(6:12, 25:27, 
                                                            44:73)])))
# Correlaci?n con todas las variables:
source("connected-graphs.R")
datos_cali %>%
  dplyr::select(c(6:10,12,25:27,44:66,70)) %>%
  corr_graph(., use.only = c(6:10,12,25:27,44:66,70),
             vertex.cex = 1, annotate = TRUE, set.seed = 105, height = 8, width = 8,
             ncol = 2, nrow = 1, margin = 0.2, line.col = rep(grey(0.3),2))

# Correlaci?n variables de edad
corrplot::corrplot.mixed(cor(datos_cali[,c(16,20,24,27,28,70)],use = "complete.obs"))
corrplot::corrplot.mixed(cor(datos_cali[,c(13,14,15,17,18,19,21,22,23,27,73)],
                             use = "complete.obs"))
# Correlaci?n asistencia escolar
corrplot::corrplot.mixed(cor(datos_cali[,c(48:53,73)],
                             use = "complete.obs"))
# Hay mucha correlaci?n con las variables de edad, se selecciona "adulto_mayor"
corrplot::corrplot.mixed(cor(datos_cali[,c(27,29:43,73)],
                             use = "complete.obs"))

############# modelo lineal con base en el an?lisis de Python ------------------
# Despues de haber arreglado los datos y agregado variables:
corrplot::corrplot.mixed(cor(datos_cali[,c(25,27,63,39,67:73)],use = "complete.obs"))
Variables_Cali <- barrios_Cali[,c('ID_BARRIO','Covid','Adulto_mayor','ind_cobert_serv_bajo',
                                  'Comorbilidad', 'Muertos', 'razon_letalidad')]

# Se debe transformar la variable de Covid 
# Aunque se usar? el valor que arroja python: 0.1834979059151258
shapiro.test(Variables_Cali$Covid)
Variables_Cali$Covid <- forecast::BoxCox(x = Variables_Cali$Covid, 
                                         lambda = 0.18349790591512582)
shapiro.test(Variables_Cali$Covid)
spplot(Variables_Cali[,2], main = colnames(Variables_Cali@data[2]))
View(Variables_Cali@data)

ggplot(data = Variables_Cali@data, aes(x = Comorbilidad, y = Covid)) + geom_point()
ggplot(data = Variables_Cali@data, aes(x = Muertos, y = Covid)) + geom_point()
ggplot(data = Variables_Cali@data, aes(x = razon_letalidad, y = Covid)) + geom_point()

step(lm(formula = Covid ~ .,data = Variables_Cali@data),
     ~1, data = datos, direction = "both")

modelo_def <- lm(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo +
                   log(Comorbilidad) + Muertos + razon_letalidad, data = Variables_Cali@data)
summary(modelo_def)

modelo_def2 <- lm(formula = Covid ~ Hogares + Adulto_mayor + ind_cobert_serv_bajo + 
             Prop_Edad_M_60a64 + `Cobertura_Internet_ fijoomovil`%in%Estrato_moda,
               data = datos_cali)
summary(modelo_def2)

spplot(Variables_Cali[,1:5], main = colnames(Variables_Cali@data[1:5]))

#validacion de supuestos 
#normalidad
library(tseries)
library(mvnormtest)
library(MVN)
library(olsrr)
library(lmtest)
library(stats)
library(car)
par(mfrow =c(1,2))
hist(modelo_def$residuals, main = "Histograma de los Residuales", axes = TRUE,
     breaks = "Sturges")
qqPlot(modelo_def$residuals, col="black", pch = 1, lwd = 1, cex = 1)

qqnorm(modelo_def$residuals,pch=20)
qqline(modelo_def$residuals)
par(mfrow=c(2,2))
plot(modelo_def)

ks.test(modelo_def$residuals,"pnorm",mean(modelo_def$residuals),
        sd(modelo_def$residuals))  
shapiro.test(rstudent(modelo_def)) 
jarque.bera.test(modelo_def$residuals) 

#Homocedasticidad
bptest(modelo_def)
#multicolinealidad
vif(modelo_def)
ols_vif_tol(modelo_def)
ols_coll_diag(modelo_def) 
ols_eigen_cindex(modelo_def)

############# Modelos de machine learning ---------------
rm(list=ls())
gc(reset = T)
library(readxl); library(maptools); library(tmap); library(dplyr)
library(rgdal); library(lme4); library(raster); library(caret); #library(mapview)
library(sf); library(CAST); library(fasterize); library(randomForest)
library(leaflet); library(spdep); library(maptools); library(spatialreg)
require(RColorBrewer); library(texreg); library(stargazer)
library(Metrics); library(xgboost); library(Metrics); library(kernlab)
library(ggspatial); library(kernlab); library(ggplot2)
############# Definir directorio y cargar datos ----------------
barrios_Cali <- readOGR('databases_and_shapefiles/Barrios Cali.shp')
datos_barrio_cali <- barrios_Cali@data
datos_cali <- read_excel('databases_and_shapefiles/Cali_nueva_base.xlsx')
barrios_Cali@data = left_join(barrios_Cali@data, datos_cali, by = 'ID_BARRIO')
Variables_Cali = barrios_Cali
Variables_Cali$geometry = Variables_Cali@polygons
Variables_Cali = Variables_Cali[,c('Covid','Adulto_mayor','ind_cobert_serv_bajo',
                                   'Comorbilidad', 'razon_letalidad',
                                   'Proporcion_trabajando', 'Estrato_moda',
                                   'Pfizer', 'ss_subsidiado',
                                   'COD_COMUNA')]
colnames(Variables_Cali@data)
View(Variables_Cali@data)
#Variables_Cali@data <- na.omit(Variables_Cali@data)
Variables_Cali = Variables_Cali[-c(65,112,304,326),]
nrow(Variables_Cali@data)
Variables_Cali$Estrato_moda2 <- as.factor(Variables_Cali$Estrato_moda)
str(Variables_Cali@data$Estrato_moda2)

#raster
r <- raster(Variables_Cali, resolution=0.0005) #, resolution=0.00005
et.ras <-rasterize(x = Variables_Cali, y = r)
s <- deratify(et.ras)
plot(s)

# Se realizar? la predicci?n de una comuna
# Predicci?n de la comuna 2
set.seed(100)
#trainids <- createDataPartition(Variables_Cali$Covid,list=FALSE,p=0.8) #Por bloques
trainData <- Variables_Cali[Variables_Cali@data$COD_COMUNA != "02",]
#trainData <- Variables_Cali@data[trainids,]
#testData <- Variables_Cali@data[-trainids,]
testData <- Variables_Cali[Variables_Cali@data$COD_COMUNA == "02",]
plot(trainData)
plot(testData)
corrplot::corrplot.mixed(corr = cor(Variables_Cali@data[,1:9]))
predictors <- c('Adulto_mayor','ind_cobert_serv_bajo',
                'Comorbilidad', 'razon_letalidad',
                'Proporcion_trabajando', 'Estrato_moda',
                'Pfizer', 'ss_subsidiado') #, 'sintomatico_si'
predictors2 <- c('Adulto_mayor','ind_cobert_serv_bajo',
                 'Comorbilidad', 'razon_letalidad',
                 'Proporcion_trabajando', 'Estrato_moda2',
                 'Pfizer', 'ss_subsidiado') #, 'sintomatico_si'
response <- 'Covid'

asd3 = CreateSpacetimeFolds(x = trainData@data,spacevar = "COD_COMUNA",
                            k=21,class = "COD_COMUNA",seed = 25)
#indices <- CreateSpacetimeFolds(trainData,spacevar = "geometry", k=5)
View(asd3)
View(asd3$indexOut)
asd3$indexOut[[21]] # Comuna 22, la ?ltima
#rownames(trainData@data) <- 1:nrow(trainData@data)
rownames(trainData@data[trainData@data$COD_COMUNA == "22",]) # Comuna 22

############# Random Forest -----------------
# Para el modelo de Random Forest se us? el estrato como factor
set.seed(100)
model_rf <- train(x = trainData@data[,predictors2],
                  y = trainData@data[,response],
                  method="rf", metric = 'RMSE', maximize = F,
                  trControl=trainControl(search = "grid",method="cv",number = 21,
                                         indexOut = asd3$indexOut), tuneLength = 10,
                  form = Covid ~ Adulto_mayor + ind_cobert_serv_bajo +
                    log(Comorbilidad) + razon_letalidad +
                    Proporcion_trabajando%in%Estrato_moda + Pfizer + ss_subsidiado,
                  data = trainData@data, subset = asd3$indexOut,
                  tuneGrid = data.frame(mtry = c(2,3,4,5,6,7,8,9,10))) #Este es: 5
print(model_rf)
ggplot(model_rf) + scale_x_log10() ############## Este es el definitivo
plot(model_rf)
print(model_rf$finalModel)
plot(varImp(model_rf))
##### view in raster
prediction_rf <- model_rf %>% predict(testData@data[,predictors2])
testData@data$prediccion_rf <- prediction_rf
comuna_pred <- deratify(rasterize(x = testData, y = raster(testData,
                                                           resolution=0.0005)))
plot(comuna_pred[[c(1,12)]])
spplot(obj = testData[,c(1,12)],
       col.regions=brewer.pal(n = 10, name = "YlGn"), cuts=8)
#prediction <- predict(comuna_pred, model_rf)

mse(testData@data$Covid, testData@data$prediccion_rf)
Metrics::mae(testData@data$Covid, testData@data$prediccion_rf)
caret::postResample(pred = testData@data$prediccion_rf, obs = testData@data$Covid)[2]

rss <- sum((testData@data$prediccion_rf - testData@data$Covid)^2)
tss <- sum((testData@data$Covid - mean(testData@data$Covid))^2)
rsq <- 1 - rss/tss; rsq

cor(testData@data$prediccion_rf, testData@data$Covid)^2 # R_squared
#mse: 11870.85, con pfizer y subsidiados: 11655.77
#mae: 76.48327, con pfizer y subsidiados: 78.27677
#Rsquared: 0.933331 ? 0.8823798, con pfizer y subsidiados: 0.9019047 ? 0.8835527

############# Nerual network -------------

#http://topepo.github.io/caret/train-models-by-tag.html#Neural_Network
#model_nn <- train(x = barrios_Cali@data[,predictors], y = barrios_Cali@data[,response],
#                  method="mlpKerasDropoutCost", metric = 'Rsquared', maximize = T,
#                  trControl=trainControl(search = "grid",method="cv",number = 22,
#                                         indexOut = asd3$indexOut), tuneLength = 10,
#                  form = Covid ~ Adulto_mayor + ind_cobert_serv_bajo + log(Comorbilidad) +
#                    razon_letalidad + Proporcion_trabajando%in%Estrato_moda,
#                  data = Variables_Cali@data, subset = asd3$indexOut,
#                  tuneGrid = data.frame(layer1 = seq(10,100),
#                                        layer2 = seq(10,100),
#                                        layer3 = seq(10,100))) #Este es
#print(model_nn)
#ggplot(model_nn) + scale_x_log10() ############## Este es el definitivo

#asd = createFolds(y = Variables_Cali$Covid%in%Variables_Cali$COD_COMUNA, k = 22,
#                  list = F, returnTrain = F)
#View(cbind(Variables_Cali@data$COD_COMUNA, asd))

#model_nn <- train(x = barrios_Cali@data[,predictors], y = barrios_Cali@data[,response],
#                  method="neuralnet", metric = 'RMSE', maximize = F,
#                  trControl=trainControl(search = "grid",method="cv",number = 22,
#                                         indexOut = asd3$indexOut), tuneLength = 10,
#                  subset = asd3$indexOut,
#                  tuneGrid = data.frame(layer1 = seq(10,100),
#                                        layer2 = seq(10,100),
#                                        layer3 = seq(10,100))) #Este es

#train(x = barrios_Cali@data[0:5,predictors], y = barrios_Cali@data[0:5,response],
#      method="neuralnet", metric = 'RMSE', maximize = F,
#      tuneGrid = data.frame(layer1 = seq(10,100),
#                            layer2 = seq(10,100),
#                            layer3 = seq(10,100)))

############# XGboost --------------
#rownames(trainData@data) <- 1:nrow(trainData@data)
#trainData$Estrato_moda <- as.factor(trainData$Estrato_moda)
# Ejecut? sin la formula, y con el estrato como num?rico
model_XGB <- train(x = trainData@data[,predictors], y = trainData@data[,response],
                   method="xgbDART", importance=TRUE, metric = 'Rsquared', maximize = T,
                   trControl = trainControl(search = "grid",method="cv",number = 21,
                                            indexOut = asd3$indexOut), tuneLength = 10,
                   subset = asd3$indexOut,
                   tuneGrid = expand.grid(nrounds = c(50, 100, 150, 180, 200, 250, 300),
                                          max_depth = c(10, 25, 50, 75, 100),
                                          eta = c(0, 0.01, 0.03, 0.05, 0.1),
                                          gamma = c(0, 3, 5, 10, 20, 30),
                                          subsample = c(0.2, 0.5, 0.8, 1),
                                          colsample_bytree = c(0.2, 0.5, 0.8, 1),
                                          rate_drop = c(0, 0.01, 0.03, 0.05, 0.10),
                                          skip_drop = c(0, 0.01, 0.03, 0.05, 0.10),
                                          min_child_weight = c(1, 10, 15, 20, 30)))
#min_child_weight means something like "stop trying to split once your sample size 
#in a node goes below a given threshold".

model_XGB
print(model_XGB$finalModel)
plot(varImp(model_XGB))

model_XGB <- train(x = trainData@data[,predictors], y = trainData@data[,response],
                   method="xgbDART", importance=TRUE, metric = 'Rsquared', maximize = T,
                   trControl = trainControl(search = "grid",method="cv",number = 21,
                                            indexOut = asd3$indexOut), tuneLength = 10,
                   subset = asd3$indexOut,
                   tuneGrid = expand.grid(nrounds = c(200),
                                          max_depth = c(50),
                                          eta = c(0.05),
                                          gamma = c(10),
                                          subsample = c(1),
                                          colsample_bytree = c(1),
                                          rate_drop = c(0),
                                          skip_drop = c(0.05),
                                          min_child_weight = c(1)))
model_XGB
print(model_XGB$finalModel)
plot(varImp(model_XGB))
##### view in raster
predictionXGB <- model_XGB %>% predict(testData@data[,predictors])
testData@data$prediccion_XGB <- predictionXGB
spplot(obj = testData[,c(1,12,13)],
       col.regions=brewer.pal(n = 10, name = "YlGn"), cuts=8)

comuna_pred <- deratify(rasterize(x = testData, y = raster(testData,
                                                           resolution=0.0005)))
predictionXGB <- predict(comuna_pred,model_XGB)
spplot(predictionXGB)
plot(comuna_pred[[c(1,12,13)]])

mse(testData@data$Covid, testData@data$prediccion_XGB)
Metrics::mae(testData@data$Covid, testData@data$prediccion_XGB)
caret::postResample(pred = testData@data$prediccion_XGB, obs = testData@data$Covid)[2]

rss <- sum((testData@data$prediccion_XGB - testData@data$Covid)^2)
tss <- sum((testData@data$Covid - mean(testData@data$Covid))^2)
rsq <- 1 - rss/tss; rsq

#mse: 17039.17, con pfizer y subsidiados: 13627.73
#mae: 96.16176, con pfizer y subsidiados: 91.99681
#Rsquared: 0.9281048 ? 0.8311704, con pfizer y subsidiados: 0.9182582 ? 0.864972

############# Support Vector Machine --------------
#with Boundrange String Kernel: svmBoundrangeString
#with exponential string Kernel: svmExpoString
#with linear Kernel: svmLinear
#with polynomial Kernel: svmPoly
#with Radial Basis Function Kernel: svmRadial
#with Spectrum String Kernel: svmSpectrumString
model_SVM3 <- train(x = trainData@data[,predictors],
                    y = trainData@data[,response], method="svmLinear",
                    importance=TRUE, metric = 'Rsquared', maximize = T,
                    trControl=trainControl(search = "grid", method="cv", number = 21,
                                           indexOut = asd3$indexOut),
                    form = Covid ~ Adulto_mayor + ind_cobert_serv_bajo + log(Comorbilidad) +
                      razon_letalidad + Proporcion_trabajando%in%Estrato_moda + 
                      Pfizer + ss_subsidiado,
                    data = trainData@data,
                    tuneLength = 10,tuneGrid = data.frame(C=c(10,20,30,40,50,60,70,
                                                              80,90,100,120,140,160))) #,90
model_SVM3
print(model_SVM3$finalModel)
plot(varImp(model_SVM3))

model_SVM3 <- train(x = trainData@data[,predictors],
                    y = trainData@data[,response], method="svmLinear",
                    importance=TRUE, metric = 'Rsquared', maximize = T,
                    trControl=trainControl(search = "grid", method="cv", number = 21,
                                           indexOut = asd3$indexOut),
                    form = Covid ~ Adulto_mayor + ind_cobert_serv_bajo +
                      log(Comorbilidad) + razon_letalidad +
                      Proporcion_trabajando%in%Estrato_moda + Pfizer +
                      ss_subsidiado,
                    data = trainData@data,
                    tuneLength = 10,tuneGrid = data.frame(C=c(40))) #,70 o 90 o 10
model_SVM3
print(model_SVM3$finalModel)
plot(varImp(model_SVM3))
#plot(model_SVM3) # C de 70 es mejor
##### view in raster
predictionSVM3 <- model_SVM3 %>% predict(testData@data[,predictors])
testData@data$prediccion_SVM3 <- predictionSVM3
spplot(obj = testData[,c(1,12,13,14)],
       col.regions=brewer.pal(n = 10, name = "YlGn"), cuts=8)

comuna_pred <- deratify(rasterize(x = testData, y = raster(testData,
                                                           resolution=0.0005)))
predictionSVM3 <- predict(comuna_pred,model_SVM3)
#spplot(predictionSVM3)
plot(comuna_pred[[c(1,12,13,14)]])

mse(testData@data$Covid, testData@data$prediccion_SVM3)
Metrics::mae(testData@data$Covid, testData@data$prediccion_SVM3)
caret::postResample(pred = testData@data$prediccion_SVM3, obs = testData@data$Covid)[2]

rss <- sum((testData@data$prediccion_SVM3 - testData@data$Covid)^2)
tss <- sum((testData@data$Covid - mean(testData@data$Covid))^2)
rsq <- 1 - rss/tss; rsq

# 16848.09
# 91.51734
# 0.8542609, 0.8330637

#mse: 16561.05, con pfizer y subsidiados: 16867.61
#mae: 96.69497, con pfizer y subsidiados: 91.53236
#Rsquared: 0.8621281 ? 0.8359078, con pfizer y subsidiados: 0.8543025 ? 0.8328703

######### SVM 2
model_SVM2 <- train(x = trainData@data[,predictors],
                    y = trainData@data[,response], method="svmPoly",
                    importance=TRUE, metric = 'Rsquared', maximize = T,
                    trControl=trainControl(search = "grid", method="cv", number = 21,
                                           indexOut = asd3$indexOut),
                    #subset = asd3$indexOut,
                    form = Covid ~ Adulto_mayor + ind_cobert_serv_bajo + log(Comorbilidad) +
                      razon_letalidad + Proporcion_trabajando%in%Estrato_moda + Pfizer +
                      ss_subsidiado,
                    data = trainData@data,
                    tuneLength = 10, tuneGrid = expand.grid(C=c(10, 20, 40, 60, 80, 100, 120), #0.1, 1, 10, 20, 
                                                            degree= c(3),
                                                            scale = c(0.01, 0.07, 0.1, 0.3)))
#https://www.quora.com/Support-Vector-Machines-How-do-I-choose-a-kernel-scale-parameter
model_SVM2
print(model_SVM2$finalModel)
plot(varImp(model_SVM2))

model_SVM2 <- train(x = trainData@data[,predictors],
                    y = trainData@data[,response], method="svmPoly",
                    importance=TRUE, metric = 'Rsquared', maximize = T,
                    trControl=trainControl(search = "grid", method="cv", number = 21,
                                           indexOut = asd3$indexOut),
                    #subset = asd3$indexOut,
                    form = Covid ~ Adulto_mayor + ind_cobert_serv_bajo +
                      log(Comorbilidad) + razon_letalidad +
                      Proporcion_trabajando%in%Estrato_moda + Pfizer +
                      ss_subsidiado + sintomatico_si,
                    data = trainData@data,
                    tuneLength = 10,tuneGrid = data.frame(C=c(80), #90
                                                          degree= c(3), #2
                                                          scale = c(0.01))) #2
model_SVM2
print(model_SVM2$finalModel)
plot(varImp(model_SVM2))
#plot(model_SVM3) # C de 70 es mejor
##### view in raster
predictionSVM2 <- model_SVM2 %>% predict(testData@data[,predictors])
testData@data$prediccion_SVM2 <- predictionSVM2
spplot(obj = testData[,c(1,12,13,14,15)],
       col.regions=brewer.pal(n = 10, name = "YlGn"), cuts=8)

comuna_pred <- deratify(rasterize(x = testData, y = raster(testData,
                                                           resolution=0.0005)))
predictionSVM2 <- predict(comuna_pred,model_SVM2)
#spplot(predictionSVM3)
plot(comuna_pred[[c(1,12,13,14,15)]])

mse(testData@data$Covid, testData@data$prediccion_SVM2)
Metrics::mae( testData@data$Covid, testData@data$prediccion_SVM2)
caret::postResample(pred = testData@data$prediccion_SVM2, obs = testData@data$Covid)[2]

rss <- sum((testData@data$prediccion_SVM2 - testData@data$Covid)^2)
tss <- sum((testData@data$Covid - mean(testData@data$Covid))^2)
rsq <- 1 - rss/tss; rsq

#mse: 21030.13, con pfizer y subsidiados 16155.15
#mae: 107.7966, con pfizer y subsidiados 100.9819
#Rsquared: 0.8408321 ? 0.7916268, con pfizer y subsidiados 0.8861105 ? 0.8399296

############ SVM 1
model_SVM1 <- train(x = trainData@data[,predictors],
                    y = trainData@data[,response], method="svmRadial",
                    importance=TRUE, metric = 'Rsquared', maximize = T,
                    trControl=trainControl(search = "grid", method="cv", number = 21,
                                           indexOut = asd3$indexOut),
                    form = Covid ~ Adulto_mayor + ind_cobert_serv_bajo +
                      log(Comorbilidad) + razon_letalidad +
                      Proporcion_trabajando%in%Estrato_moda + Pfizer +
                      ss_subsidiado + sintomatico_si,
                    data = trainData@data,
                    tuneLength = 10,
                    tuneGrid = expand.grid(C=c(1,3,5,10,20,30,40,50,60,70,80,90,100),
                                           sigma = c(0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,
                                                     0.08,0.09,0.1,0.3,0.5,1,2)))

#20, 90
model_SVM1 #  sigma = 0.6 and C = 20.
print(model_SVM1$finalModel)
plot(varImp(model_SVM1))

model_SVM1 <- train(x = trainData@data[,predictors],
                    y = trainData@data[,response], method="svmRadial",
                    importance=TRUE, metric = 'Rsquared', maximize = T,
                    trControl=trainControl(search = "grid", method="cv", number = 21,
                                           indexOut = asd3$indexOut),
                    form = Covid ~ Adulto_mayor + ind_cobert_serv_bajo +
                      log(Comorbilidad) + razon_letalidad +
                      Proporcion_trabajando%in%Estrato_moda + Pfizer +
                      ss_subsidiado + sintomatico_si,
                    data = trainData@data,
                    tuneLength = 10,
                    tuneGrid = expand.grid(C=c(20),
                                           sigma = c(0.01)))

#20, 90
model_SVM1 #  sigma = 0.6 and C = 20.
print(model_SVM1$finalModel)
plot(varImp(object = model_SVM1, scale = F, useModel = F))
varImp(object = model_SVM1, scale = F, nonpara = T)

##### view in raster
predictionSVM1 <- model_SVM1 %>% predict(testData@data[,predictors])
testData@data$prediccion_SVM1 <- predictionSVM1
#spplot(obj = testData[,c(1,12,13,14,15,16)],
#       col.regions=brewer.pal(n = 10, name = "YlGn"), cuts=8)
spplot(obj = testData[,c(1,12,13,14)],
       col.regions=brewer.pal(n = 10, name = "YlGn"), cuts=8)

comuna_pred <- deratify(rasterize(x = testData, y = raster(testData,
                                                           resolution=0.0005)))
predictionSVM1 <- predict(comuna_pred,model_SVM1)
#spplot(predictionSVM1)
#plot(comuna_pred[[c(1,12,13,14,15,16)]])
plot(comuna_pred[[c(1,12,13,14)]])

mse(testData@data$Covid, testData@data$prediccion_SVM1)
Metrics::mae(testData@data$Covid, testData@data$prediccion_SVM1)
caret::postResample(pred = testData@data$prediccion_SVM1, obs = testData@data$Covid)[2]

rss <- sum((testData@data$prediccion_SVM1 - testData@data$Covid)^2)
tss <- sum((testData@data$Covid - mean(testData@data$Covid))^2)
rsq <- 1 - rss/tss; rsq

#mse: 24977.48        con pfizer 15122.47
#mae: 124.0734        con pfizer 92.97329
#Rsquared: 0.760055 ? 0.7525151 con pfizer 0.9004668 ? 0.8501617
testData <- sf::st_as_sf(testData)
g1 <- ggplot(data = testData) +
  geom_sf(aes(fill=Covid), color="black", size = 0.2) +
  ggtitle("Casos de Covid-19 en la comuna 2. Valor real") +
  theme(panel.background = element_blank(), plot.title = element_text(size=14),
        legend.title=element_blank(), panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_distiller(palette ="RdBu")

g1

g2 <- ggplot(data = testData) +
  geom_sf(aes(fill=prediccion_rf), color="black", size = 0.2) +
  ggtitle("Casos de Covid-19 en la comuna 2. Predicci?n RF") +
  theme(panel.background = element_blank(), plot.title = element_text(size=15),
        legend.title=element_blank(), panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_distiller(palette ="RdBu")

g2

g3 <- ggplot(data = testData) +
  geom_sf(aes(fill=prediccion_XGB), color="black", size = 0.2) +
  ggtitle("Casos de Covid-19 en la comuna 2. Predicci?n XGB") +
  theme(panel.background = element_blank(), plot.title = element_text(size=15),
        legend.title=element_blank(), panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_distiller(palette ="RdBu")

g3

g4 <- ggplot(data = testData) +
  geom_sf(aes(fill=prediccion_SVM1), color="black", size = 0.2) +
  ggtitle("Casos de Covid-19 en la comuna 2. Predicci?n SVM1") +
  theme(panel.background = element_blank(), plot.title = element_text(size=15),
        legend.title=element_blank(), panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_distiller(palette ="RdBu")

g4

g5 <- ggplot(data = testData) +
  geom_sf(aes(fill=prediccion_SVM2), color="black", size = 0.2) +
  ggtitle("Casos de Covid-19 en la comuna 2. Predicci?n SVM2") +
  theme(panel.background = element_blank(), plot.title = element_text(size=15),
        legend.title=element_blank(), panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_distiller(palette ="RdBu")

g5

g6 <- ggplot(data = testData) +
  geom_sf(aes(fill=prediccion_SVM3), color="black", size = 0.2) +
  ggtitle("Casos de Covid-19 en la comuna 2. Predicci?n SVM3") +
  theme(panel.background = element_blank(), plot.title = element_text(size=15),
        legend.title=element_blank(), panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_distiller(palette ="RdBu")

g6

prediccion_nn <- c(808.3204346,133.0464783,1184.459595,734.730896,186.657074,1006.644958,
                   1130.816406,460.6136475,220.758667,813.5688477,495.6878662,99.41513824,
                   228.3886719,923.6265259,672.9551392,657.1071167,165.6322479,328.4180603,
                   643.5836792,385.6372986,505.1529846,326.2724609,200.3936768,724.6391602,
                   365.0610962)
testData$prediccion_nn <- prediccion_nn
View(testData[,c(1,2,18)])

g7 <- ggplot(data = testData) +
  geom_sf(aes(fill=prediccion_nn), color="black", size = 0.2) +
  ggtitle("Casos de Covid-19 en la comuna 2. Predicci?n NN") +
  theme(panel.background = element_blank(), plot.title = element_text(size=15),
        legend.title=element_blank(), panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) + scale_fill_distiller(palette ="RdBu")

g7

ggpubr::ggarrange(g1, g2, g3, g4, g5, g7, common.legend = TRUE, legend="bottom", ncol = 3, nrow = 2)
############# Otros modelos de ML con base anterior ------------------
set.seed(100)
model_rf <- train(trainData[,predictors],trainData[,response],
                   method="rf", importance=TRUE,metric = 'RMSE',maximize = F,
                   trControl=trainControl(search = "random",method="cv",
                                          index = indices$index)) #Este es

model_SVM <- train(trainData[,predictors],trainData[,response],
               method="svmLinear", importance=TRUE,metric = 'RMSE',maximize = F,
               trControl=trainControl(search = "random",method="cv",
                                      index = indices$index))

model_GB <- train(trainData[,predictors],trainData[,response],
                   method="xgbDART", importance=TRUE,metric = 'RMSE',maximize = F,
                   trControl=trainControl(search = "random",method="cv",
                                          index = indices$index))

model_knn <- train(trainData[,predictors],trainData[,response],
                  method="kknn", importance=TRUE,metric = 'RMSE',maximize = F,
                  trControl=trainControl(search = "random",method="cv",
                                         index = indices$index))

model_knn2 <- train(trainData[,predictors],trainData[,response],
                   method="knn", importance=TRUE,metric = 'RMSE',maximize = F,
                   trControl=trainControl(search = "random",method="cv",
                                          index = indices$index))

model_brnn <- train(trainData[,predictors],trainData[,response],
                  method="brnn", importance=TRUE,metric = 'RMSE',maximize = F,
                  trControl=trainControl(search = "random",method="cv",
                                         index = indices$index))

model_avnnet <- train(trainData[,predictors],trainData[,response],
                    method="avNNet", importance=TRUE,metric = 'Rsquared',maximize = T,
                    trControl=trainControl(search = "random",method="cv",
                                           index = indices$index))

model_nn <- train(trainData[,predictors],trainData[,response],
                  method="nnet", importance=TRUE,metric = 'Rsquared',maximize = T,
                  trControl=trainControl(search = "random",method="cv",
                                         index = indices$index))

print(model_rf)
print(model_SVM)
print(model_GB)
print(model_knn)
print(model_knn2)
print(model_brnn)
print(model_avnnet)
print(model_nn)

plot(model_rf)
plot(varImp(model_rf))
plot(model_SVM)
plot(varImp(model_SVM))
plot(model_GB)
plot(varImp(model_GB))

plot(Variables_Cali[-c(trainids),])
# Se debe empezar a mirar como hacer particiones espaciales, bien sea por 
# deciles de Covid-19, por estrato, comuna. Probar muchar particiones por variables
# e identificadores de ubicaci?n para mirar cual da mejor ajuste.

print(model_rf$finalModel)
print(model_SVM$finalModel)
print(model_GB$finalModel)
print(model_brnn$finalModel)
print(model_avnnet$finalModel)
prediction_rf <- predict(model_rf, testData)
prediction_GB <- predict(model_GB, testData)
prediction_brnn <- predict(model_brnn, testData)
prediction_knn <- predict(model_knn, testData)
prediction_knn2 <- predict(model_knn2, testData)
prediction_avnnet <- predict(model_avnnet, testData)
prediction_nn <- predict(model_nn)


prediction_rf <- forecast::InvBoxCox(x = prediction_rf, lambda = 0.18349790591512582)
prediction_GB <- forecast::InvBoxCox(x = prediction_GB, lambda = 0.18349790591512582)
prediction_brnn <- forecast::InvBoxCox(x = prediction_brnn, lambda = 0.18349790591512582)
prediction_knn <- forecast::InvBoxCox(x = prediction_knn, lambda = 0.18349790591512582)
prediction_knn2 <- forecast::InvBoxCox(x = prediction_knn2, lambda = 0.18349790591512582)
prediction_avnnet <- forecast::InvBoxCox(x = prediction_avnnet, lambda = 0.18349790591512582)
prediction_nn <- forecast::InvBoxCox(x = prediction_nn, lambda = 0.18349790591512582)
testData$Covid <- forecast::InvBoxCox(x = testData$Covid, 
                                      lambda = 0.18349790591512582)


mse_rf <- (sum((testData$Covid - prediction_rf)^2))/nrow(testData)
mse_GB <- (sum((testData$Covid - prediction_GB)^2))/nrow(testData)
mse_brnn <- (sum((testData$Covid - prediction_brnn)^2))/nrow(testData)
mse_knn <- (sum((testData$Covid - prediction_knn)^2))/nrow(testData)
mse_knn2 <- (sum((testData$Covid - prediction_knn2)^2))/nrow(testData)
mse_avnnet <- (sum((testData$Covid - prediction_avnnet)^2))/nrow(testData)

Comparacion <- cbind(prediction_brnn,prediction_knn,
                     prediction_rf, prediction_GB, testData$Covid)

############# Econometria espacial --------------------
#http://www.econ.uiuc.edu/~lab/workshop/Spatial_in_R.html
Variables_Cali <- barrios_Cali[,c('ID_BARRIO','Covid','Adulto_mayor','ind_cobert_serv_bajo',
                                  'Comorbilidad', 'Muertos', 'razon_letalidad')]
Variables_Cali$Covid <- forecast::BoxCox(x = Variables_Cali$Covid, 
                                         lambda = 0.18349790591512582)

Variables_Cali <- Variables_Cali[-c(65,112,304,326),]
plot(Variables_Cali)
Variables_Cali2 <- Variables_Cali[c(trainids),] # Data de entrenamiento
plot(Variables_Cali2)
class(Variables_Cali)
str(slot(Variables_Cali,"data"))

summary(Variables_Cali@data$Covid) # Tiene transformaci?n Box-cox

leaflet(Variables_Cali) %>%
  addPolygons(stroke = FALSE, fillOpacity = 0.5, smoothFactor = 0.5) %>%
  addTiles() #adds a map tile, the default is OpenStreetMap

qpal<-colorQuantile("OrRd", Variables_Cali@data$Covid, n=9) 
leaflet(Variables_Cali) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, 
              color = ~qpal(Covid)) %>% addTiles()

##### Empezar los modelos espaciales:
modelo_ols <- lm(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo + 
                   + log(Comorbilidad) + Muertos + razon_letalidad, 
                 data = Variables_Cali2@data)
summary(modelo_ols)

list.queen <- poly2nb(Variables_Cali2, queen=TRUE)
W <- nb2listw(list.queen, style="W", zero.policy=TRUE)
W

plot(W,coordinates(Variables_Cali2))

coords <- coordinates(Variables_Cali2)
W_dist <- dnearneigh(coords,0,1,longlat = FALSE)
plot(W_dist,coordinates(Variables_Cali2))

moran.lm <- lm.morantest(modelo_ols, W, alternative="two.sided")
print(moran.lm)
moran.test(Variables_Cali2$Covid, W, alternative="two.sided",randomisation = F)
moran.plot(x = Variables_Cali2$Covid, listw = W, zero.policy = T,
           pch=16, col="black",cex=.5, quiet=TRUE, 
           labels=barrios_Cali[c(trainids),]$Barrio_Urbanizacion_o_Sector,
           xlab="Covid en Cali", 
           ylab="Covid (Spatial Lag)")

moran.test(Variables_Cali2$Comorbilidad, W, alternative="two.sided",randomisation = F)
moran.plot(x = Variables_Cali2$Comorbilidad, listw = W, zero.policy = T,
           pch=16, col="black",cex=.5, quiet=TRUE, 
           labels=barrios_Cali[c(trainids),]$Barrio_Urbanizacion_o_Sector,
           xlab="Comorbilidades en Cali", 
           ylab="Comorbilidades (Spatial Lag)")

moran.test(Variables_Cali2$Adulto_mayor, W, alternative="two.sided",randomisation = F)
moran.plot(x = Variables_Cali2$Adulto_mayor, listw = W, zero.policy = T,
           pch=16, col="black",cex=.5, quiet=TRUE, 
           labels=barrios_Cali[c(trainids),]$Barrio_Urbanizacion_o_Sector,
           xlab="Adultos mayores en Cali", 
           ylab="Adultos mayores (Spatial Lag)")

moran.test(Variables_Cali2$ind_cobert_serv_bajo, W, alternative="two.sided",
           randomisation = F)
moran.plot(x = Variables_Cali2$ind_cobert_serv_bajo, listw = W, zero.policy = T,
           pch=16, col="black",cex=.5, quiet=TRUE, 
           labels=barrios_Cali[c(trainids),]$Barrio_Urbanizacion_o_Sector,
           xlab="ICS en Cali", 
           ylab="ICS (Spatial Lag)")

moran.test(Variables_Cali2$Muertos, W, alternative="two.sided",
           randomisation = F)
moran.plot(x = Variables_Cali2$Muertos, listw = W, zero.policy = T,
           pch=16, col="black",cex=.5, quiet=TRUE, 
           labels=barrios_Cali[c(trainids),]$Barrio_Urbanizacion_o_Sector,
           xlab="Muertos en Cali", 
           ylab="Muertos (Spatial Lag)")

moran.test(Variables_Cali2$razon_letalidad, W, alternative="two.sided",
           randomisation = F)
moran.plot(x = Variables_Cali2$razon_letalidad, listw = W, zero.policy = T,
           pch=16, col="black",cex=.5, quiet=TRUE, 
           labels=barrios_Cali[c(trainids),]$Barrio_Urbanizacion_o_Sector,
           xlab="Raz?n de letalidad en Cali", 
           ylab="Raz?n de letalidad (Spatial Lag)")


LM <- lm.LMtests(modelo_ols, W, test="all")
print(LM) # RLMerr parece ser el mejor

## Error sar lm
modelo_errorsalm <- errorsarlm(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo + 
                                 + log(Comorbilidad) + Muertos + razon_letalidad, 
                               data = Variables_Cali2@data, W)
summary(modelo_errorsalm)

#Spatial Error model by the Feasible Generalized Least Squares (GLS)
modelo_fglserror <- GMerrorsar(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo + 
                                 + log(Comorbilidad) + Muertos + razon_letalidad, 
                               data = Variables_Cali2@data, listw = W)
#listw = nb2listw(poly2nb(Variables_Cali2, queen=TRUE), 
#                 style="W", zero.policy=TRUE), na.action = na.omit)
summary(modelo_fglserror)

## Error SAR Durbin model
modelo_Durbin_error <- errorsarlm(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo + 
                                    + log(Comorbilidad) + Muertos + razon_letalidad, 
                                  data = Variables_Cali2@data, W,Durbin = T)
summary(modelo_Durbin_error)

## Modelo lag-sar
modelo_lagsar <- lagsarlm(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo + 
                            + log(Comorbilidad) + Muertos + razon_letalidad, 
                          data = Variables_Cali2@data,listw = W)
summary(modelo_lagsar)

## Modelo lag-sar Durbin
modelo_Durbin <- lagsarlm(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo + 
                            + log(Comorbilidad) + Muertos + razon_letalidad, 
                          data = Variables_Cali2@data,listw = W,Durbin = T)
summary(modelo_Durbin)

## Modelo SAC 
modelo_SAC <- sacsarlm(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo + 
                         + log(Comorbilidad) + Muertos + razon_letalidad, 
                       data = Variables_Cali2@data,listw = W)
summary(modelo_SAC)

## Modelo SAC Durbin
modelo_cliff_ord <- sacsarlm(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo + 
                         + log(Comorbilidad) + Muertos + razon_letalidad, 
                       data = Variables_Cali2@data,listw = W,Durbin = T)
summary(modelo_cliff_ord)

## Modelo Spatial-lagged x variables
modelo_SLX <- lmSLX(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo + 
                      + log(Comorbilidad) + Muertos + razon_letalidad, 
                    data = Variables_Cali2@data, listw = W, zero.policy = TRUE)
summary(modelo_SLX)

## Spatial Lag model with the Two-Stage Least Squares (2SLS)
modelo_stsls <- stsls(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo + 
                        + log(Comorbilidad) + Muertos + razon_letalidad, 
                      data = Variables_Cali2@data, listw = W, 
                    zero.policy = TRUE)
summary(modelo_stsls)

#https://keen-swartz-3146c4.netlify.app/spatecon.html

#Variables_Cali@data <- na.omit(Variables_Cali@data)
Variables_Cali2@data$ols.res <- resid(modelo_ols) #residuals ols
Variables_Cali2@data$sem.res <- resid(modelo_errorsalm) #residual sar
Variables_Cali2@data$fglserror <- resid(modelo_fglserror)
Variables_Cali2@data$Durbin_error <- resid(modelo_Durbin_error)
Variables_Cali2@data$lagsar <- resid(modelo_lagsar)
Variables_Cali2@data$Durbin <- resid(modelo_Durbin)
Variables_Cali2@data$SAC <- resid(modelo_SAC)
Variables_Cali2@data$cliff_ord <- resid(modelo_cliff_ord)
Variables_Cali2@data$SLX <- resid(modelo_SLX)
Variables_Cali2@data$stsls <- resid(modelo_stsls)
Variables_Cali2@data$brnn <- resid(model_brnn)

spplot(Variables_Cali2[,8], main = colnames(Variables_Cali2@data[8]))
spplot(Variables_Cali2[,9], main = colnames(Variables_Cali2@data[9]))
spplot(Variables_Cali2[,10], main = colnames(Variables_Cali2@data[10]))
spplot(Variables_Cali2[,11], main = colnames(Variables_Cali2@data[11]))
spplot(Variables_Cali2[,12], main = colnames(Variables_Cali2@data[12]))
spplot(Variables_Cali2[,13], main = colnames(Variables_Cali2@data[13]))
spplot(Variables_Cali2[,14], main = colnames(Variables_Cali2@data[14]))
spplot(Variables_Cali2[,15], main = colnames(Variables_Cali2@data[15]))
spplot(Variables_Cali2[,16], main = colnames(Variables_Cali2@data[16]))
spplot(Variables_Cali2[,17], main = colnames(Variables_Cali2@data[17]))
spplot(Variables_Cali2[,18], main = colnames(Variables_Cali2@data[18]))

impacts(modelo_cliff_ord, listw=W) # Toca mirar con se hace la inversa de box-cox
# para interpretar eso

texreg(list(modelo_ols, modelo_lagsar, modelo_errorsalm, modelo_SAC, 
            modelo_Durbin, modelo_Durbin_error, modelo_cliff_ord),
       digits = 3) #### Este es 

##### EFectos marginales con modelo completo:
modelo_clifford_global <- sacsarlm(formula = Covid ~ log(Adulto_mayor) + ind_cobert_serv_bajo + 
                                     + log(Comorbilidad) + Muertos + razon_letalidad, 
                                   data = Variables_Cali@data,
                       listw = nb2listw(poly2nb(Variables_Cali, queen=TRUE), style="W", zero.policy=TRUE),Durbin = T)
summary(modelo_clifford_global)

impacts(modelo_clifford_global, listw=nb2listw(poly2nb(Variables_Cali, queen=TRUE), style="W", zero.policy=TRUE))

Variables_Cali$geometry = Variables_Cali@polygons
set.seed(100)
model_brnn_global <- train(Variables_Cali@data[,predictors],Variables_Cali@data[,response],
                           method="brnn", importance=TRUE,metric = 'RMSE',
                           maximize = F, trControl=trainControl(search = "random",method="cv",
                                           index = CreateSpacetimeFolds(Variables_Cali,spacevar = "geometry",
                                                                                   k=5)$index))
print(model_brnn_global)
Variables_Cali@data$clifford_global <- resid(modelo_clifford_global)
Variables_Cali@data$brnn_global <- resid(model_brnn_global)
spplot(Variables_Cali[,9], main = "Residuales del modelo Cliff-Ord", at=seq(min(Variables_Cali@data$clifford_global),
                                      max(Variables_Cali@data$clifford_global),
                                      length=12),
       col.regions=rev(brewer.pal(11,"RdBu")))


spplot(Variables_Cali[,10], main = "Residuales del modelo Bayesian Regularized Neural Network", 
       at=seq(min(Variables_Cali@data$brnn_global),
              max(Variables_Cali@data$brnn_global),length=12),
       col.regions=rev(brewer.pal(11,"RdBu")))
############# Predicciones econometr?a espacial -------------
#https://github.com/jsay/spatial-pred-R/blob/master/DOC.org

#endogenous predictor
prdKP1 <- function(prdWX, rho= rho, lsw= lsw,
                   power= power, order= order, tol= tol){
  if (power){
    W <- as(as_dgRMatrix_listw(lsw), "CsparseMatrix")
    prdKP1 <- c(as(powerWeights(W, rho= rho, X= as.matrix(prdWX),
                                order= order, tol= tol), "matrix"))
  } else {
    prdKP1 <- c(invIrW(lsw, rho) %*% prdWX)
  }
  prdKP1
}

#exogenous predictor
prdWX <- function(prdX, X= X, Bl= Bl, mod= mod, lsw= lsw){
  if (!mod %in% c("sxm", "sdm", "smc")){
    prdWX <- prdX } else {
      K <- ifelse(colnames(X)[ 1] == "(Intercept)", 2, 1)
      m <- ncol(X) ; WX <- matrix(nrow= length(prdX), ncol= m+ 1- K)
      for (k in K: m){
        WX[, k+ 1- K] <- lag.listw(lsw, X[, k])
      }
      prdWX <- prdX+ (WX %*% Bl)
    } 
  prdWX
}

sppred <- function(object, newdata = NULL, listw = NULL, yobs= object$y,
                   condset= "DEF", blup = NULL, loo = FALSE, rg = NULL, SCrho= FALSE, SClab= FALSE,
                   power = NULL, zero.policy = NULL, legacy = TRUE, order = 250,
                   tol= .Machine$double.eps^(3/5), ..., tst= 1, ctm= F) {
  require(spdep)
  ## USUAL VERIFICATIONS
  if (is.null(zero.policy)) 
    zero.policy <- get("zeroPolicy", envir = spdep:::.spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (is.null(power)) power <- object$method != "eigen"
  stopifnot(is.logical(legacy)) ; stopifnot(is.logical(power))
  ## DETERMINING THE MODEL
  if (object$type== "error"){
    mod <- ifelse(object$etype== "error", "sem", "sxm")
  } else {
    mod <- switch(object$type, "lag"= "sar", "mixed"= "sdm",
                  "sac"= "sac", "sacmixed"= "smc")}
  ## DATA SHAPING
  if (mod %in% c("sem", "sxm")) {lab= object$lambda ; rho= 0         }
  if (mod %in% c("sar", "sdm")) {lab= 0             ; rho= object$rho}
  if (mod %in% c("sac", "smc")) {lab= object$lambda ; rho= object$rho}
  Wlg <- substr(names(object$coefficients), 1, 4)== "lag."
  B <- object$coefficients[ !Wlg] ; Bl <- object$coefficients[ Wlg]
  if (is.null(newdata)){
    X   <- object$X[, !Wlg]
  } else {
    frm <- formula(object$call)
    mt  <- delete.response(terms(frm, data = newdata))
    mf  <- model.frame(mt, newdata)
    X   <- model.matrix(mt, mf)
    if (any(object$aliased)) X <- X[, -which(object$aliased)]
  }
  ## WEIGHT MATRIX
  if (is.null(listw)) lsw <- eval(object$call$listw) else lsw <- listw
  ## PREDICTORS
  if (is.null(blup)){
    pt <- switch(condset, "X"= 1, "XW"= 2, "DEF"= 3, "XWy"= 4)
  } else {
    pt <- switch(blup, "LSP"= 5, "KPG"= 6, "KP2"= 7, "KP3"= 8)
  }
  prdX <- as.vector(X %*% B)
  if (pt> 1) prdWX   <- prdWX(prdX, X, Bl, mod, lsw)
  if (pt> 2 && pt!= 4) prdKP1  <- prdKP1(prdWX, rho, lsw, power, order, tol)
  if (SCrho) rho= 0
  if (SClab) lab= 0
  if (pt> 3 && tst== 1){
    prdWXy <- prdWX+
      rho* lag.listw(lsw, yobs)+ lab* lag.listw(lsw, yobs- prdWX)}
  if (pt> 3 && tst== 2){
    prdWXy <- prdWX+ rho* lag.listw(lsw, yobs)}
  if (pt==5) prdLSP <- prdLSP(prdKP1, rho, lab, lsw, yobs, loo)
  if (pt==6) prdKPG <- prdKPG(prdKP1, yobs, rg, loo)
  if (pt> 6 && !loo) stop("Set loo= TRUE for this blup predictor")
  if (pt==7){
    prdKP2 <- prdKP2(prdKP1, prdWXy,
                     rho, lab, lsw, yobs, power, order, tol)}
  if (pt==8){
    prdKP3 <- prdKP3(prdKP1, prdWXy,
                     rho, lab, lsw, yobs, power, order, tol)}
  prd <- switch(pt, "1"= prdX  , "2"= prdWX , "3"= prdKP1, "4"= prdWXy,
                "5"= prdLSP, "6"= prdKPG, "7"= prdKP2, "8"= prdKP3)
  if (ctm) prd <- prdctm(prdWX, rho, lab, lsw, yobs)
  class(prd) <- "sppred" ; as.vector(prd) 
}


plot(Variables_Cali[-c(trainids),])
nb2listw(poly2nb(Variables_Cali[-c(trainids),], queen=TRUE), style="W", zero.policy=TRUE)
nb2listw(dnearneigh(coordinates(Variables_Cali[-c(trainids),]),0,1,longlat = FALSE), 
         style="W", zero.policy=TRUE)
# Se observa que no hay vecinos cerca, se debe realizar una matriz de distancias.
plot(nb2listw(dnearneigh(coordinates(Variables_Cali[-c(trainids),]),0,1,longlat = FALSE), 
              style="W", zero.policy=TRUE),coordinates(Variables_Cali[-c(trainids),]))

prediccion_ols <- predict(object = modelo_ols, testData[,c(3,4,5,6,7)])
prediccion_ols <- forecast::InvBoxCox(x = prediccion_ols, lambda = 0.18349790591512582)
Comparacion <- cbind(prediccion_ols, Comparacion)
mse_ols <- (sum((testData$Covid - prediccion_ols)^2))/nrow(testData)

prediccion_cliff_ord <- sppred(object = modelo_cliff_ord, newdata = testData[,c(3,4,5,6,7)], 
                            listw = nb2listw(dnearneigh(coordinates(Variables_Cali[-c(trainids),]),0,1,longlat = FALSE), 
                                             style="W", zero.policy=TRUE))

prediccion_cliff_ord <- forecast::InvBoxCox(x = prediccion_cliff_ord, lambda = 0.18349790591512582)
Comparacion <- cbind(prediccion_cliff_ord, Comparacion)
mse_cliff_ord <- (sum((testData$Covid - prediccion_cliff_ord)^2))/nrow(testData)

prediccion_DSEM <- sppred(object = modelo_Durbin_error, newdata = testData[,c(3,4,5,6,7)], 
                            listw = nb2listw(dnearneigh(coordinates(Variables_Cali[-c(trainids),]),0,1,longlat = FALSE), 
                                             style="W", zero.policy=TRUE))
prediccion_DSEM <- forecast::InvBoxCox(x = prediccion_DSEM, lambda = 0.18349790591512582)
Comparacion <- cbind(prediccion_DSEM, Comparacion)
mse_DSEM <- (sum((testData$Covid - prediccion_DSEM)^2))/nrow(testData)

prediccion_errorsalm <- sppred(object = modelo_errorsalm, newdata = testData[,c(3,4,5,6,7)], 
                          listw = nb2listw(dnearneigh(coordinates(Variables_Cali[-c(trainids),]),0,1,longlat = FALSE), 
                                           style="W", zero.policy=TRUE))
prediccion_errorsalm <- forecast::InvBoxCox(x = prediccion_errorsalm, lambda = 0.18349790591512582)
Comparacion <- cbind(prediccion_errorsalm, Comparacion)
mse_errorsalm <- (sum((testData$Covid - prediccion_errorsalm)^2))/nrow(testData)

prediccion_lagsar <- sppred(object = modelo_lagsar, newdata = testData[,c(3,4,5,6,7)], 
                               listw = nb2listw(dnearneigh(coordinates(Variables_Cali[-c(trainids),]),0,1,longlat = FALSE), 
                                                style="W", zero.policy=TRUE))
prediccion_lagsar <- forecast::InvBoxCox(x = prediccion_lagsar, lambda = 0.18349790591512582)
Comparacion <- cbind(prediccion_lagsar, Comparacion)
mse_lagsar <- (sum((testData$Covid - prediccion_lagsar)^2))/nrow(testData)

prediccion_Durbin <- sppred(object = modelo_Durbin, newdata = testData[,c(3,4,5,6,7)], 
                            listw = nb2listw(dnearneigh(coordinates(Variables_Cali[-c(trainids),]),0,1,longlat = FALSE), 
                                             style="W", zero.policy=TRUE))
prediccion_Durbin <- forecast::InvBoxCox(x = prediccion_Durbin, lambda = 0.18349790591512582)
Comparacion <- cbind(prediccion_Durbin, Comparacion)
mse_Durbin <- (sum((testData$Covid - prediccion_Durbin)^2))/nrow(testData)

prediccion_SAC <- sppred(object = modelo_SAC, newdata = testData[,c(3,4,5,6,7)], 
                               listw = nb2listw(dnearneigh(coordinates(Variables_Cali[-c(trainids),]),0,1,longlat = FALSE), 
                                                style="W", zero.policy=TRUE))
prediccion_SAC <- forecast::InvBoxCox(x = prediccion_SAC, lambda = 0.18349790591512582)
Comparacion <- cbind(prediccion_SAC, Comparacion)
mse_SAC <- (sum((testData$Covid - prediccion_SAC)^2))/nrow(testData)


############# Pruebas varias iniciales con modelos de bases nuevas ---------------
## Random forest inicial 
# Habr?a que intentar con todos los datos y adem?s imputando datos faltantes
barrios_Cali = barrios_Cali[-304,]
barrios_Cali = barrios_Cali[,c("Muertos", "Recuperados", "razon_letalidad",
                               'Covid', 'COD_COMUNA')]
View(barrios_Cali@data)
r <- raster(barrios_Cali, resolution=0.0005) #, resolution=0.00005
et.ras <-rasterize(x = barrios_Cali, y = r)
s <- deratify(et.ras)
plot(s)
predictors <- c("Muertos", "Recuperados", "razon_letalidad")
response <- 'Covid'
asd3 = CreateSpacetimeFolds(x = barrios_Cali@data,spacevar = "COD_COMUNA",
                            k=22,class = "COD_COMUNA",seed = 25)
# Random forest

model_rf <- train(x = barrios_Cali@data[,predictors], y = barrios_Cali@data[,response],
                  method="rf", metric = 'RMSE', maximize = F,
                  trControl=trainControl(search = "grid",method="cv",number = 22,
                                         indexOut = asd3$indexOut), tuneLength = 10,
                  subset = asd3$indexOut,
                  tuneGrid = data.frame(mtry = c(2,3,4,5,6,7,8,9,10)))
prediction <- predict(s,model_rf)
spplot(prediction)


# Extreme Gradient Boosting inicial
model_XGB <- train(x = trainData@data[,predictors], y = trainData@data[,response],
                   method="xgbDART", importance=TRUE, metric = 'Rsquared',
                   maximize = T,
                   trControl = trainControl(search = "grid",method="cv",number = 21,
                                            indexOut = asd3$indexOut),
                   tuneLength = 10,
                   form = Covid ~ Adulto_mayor + ind_cobert_serv_bajo +
                     log(Comorbilidad) + razon_letalidad +
                     Proporcion_trabajando%in%Estrato_moda,
                   data = trainData@data, subset = asd3$indexOut,
                   tuneGrid = expand.grid(nrounds = c(10, 20, 30, 40, 50, 100, 200),
                                          max_depth = c(10, 15, 20, 25, 30, 40, 50),
                                          eta = c(0.01, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3),
                                          gamma = c(1,3,5,7,10,12,15),
                                          subsample = c(1,0.9,0.8,0.7,0.6,0.55,0.5),
                                          colsample_bytree = c(1,2,3,4,5,6,7),
                                          rate_drop = c(0,0.005,0.01,0.02,0.05,0.07,0.1),
                                          skip_drop = c(0,0.01,0.03,0.07,0.1,0.15,0.2),
                                          min_child_weight = c(1,3,5,10,16,20,30)))
model_XGB

model_XGB <- train(x = barrios_Cali@data[,predictors], y = barrios_Cali@data[,response],
                   method="xgbDART", importance=TRUE,metric = 'Rsquared',maximize = T,
                   trControl=trainControl(search = "grid",method="cv",number = 22,
                                          indexOut = asd3$indexOut), tuneLength = 10,
                   subset = asd3$indexOut,
                   tuneGrid = data.frame(nrounds = c(200),
                                         max_depth = c(50),
                                         eta = c(0.05),
                                         gamma = c(10),
                                         subsample = c(1),
                                         colsample_bytree = c(1),
                                         rate_drop = c(0),
                                         skip_drop = c(0.05),
                                         min_child_weight = c(1)))
model_XGB
predictionXGB <- predict(s,model_XGB)
spplot(predictionXGB)
"""
Grilla de par?metros xgboost:
,
tuneGrid = data.frame(nrounds = c(10, 20, 30, 40, 50, 100, 200),
                      max_depth = c(10, 15, 20, 25, 30, 40, 50),
                      eta = c(0.01, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3),
                      gamma = c(1,3,5,7,10,12,15),
                      subsample = c(1,1,1,1,1,1,1),
                      colsample_bytree = c(1,2,3,4,5,6,7),
                      rate_drop = c(0,0.005,0.01,0.02,0.05,0.07,0.1),
                      skip_drop = c(0,0,0,0,0,0,0),
                      min_child_weight = c(1,3,5,10,16,20,30))
"""

#### Support vector machine inicial
model_SVM3 <- train(x = barrios_Cali@data[,predictors],
                    y = barrios_Cali@data[,response], method="svmLinear",
                    importance=TRUE, metric = 'Rsquared', maximize = T,
                    trControl=trainControl(search = "grid", method="cv", number = 22,
                                           indexOut = asd3$indexOut),
                    subset = asd3$indexOut,
                    tuneLength = 10,tuneGrid = data.frame(C=c(70,90)))
model_SVM3
plot(model_SVM3) # C de 70 es mejor
predictionSVM3 <- predict(s,model_SVM3)
print(model_SVM3$finalModel)
plot(varImp(model_SVM3))
spplot(predictionSVM3)

predictionSVM1 <- model_SVM3 %>% predict(barrios_Cali@data[,predictors])
barrios_Cali@data$prediccion_SVM3 <- predictionSVM3

r <- raster(barrios_Cali, resolution=0.0005) #, resolution=0.00005
et.ras <-rasterize(x = barrios_Cali, y = r)
s <- deratify(et.ras)
plot(s$prediccion_SVM3)

#predictionSVM3 <- predict(s,model_SVM3)
spplot(predictionSVM3)
