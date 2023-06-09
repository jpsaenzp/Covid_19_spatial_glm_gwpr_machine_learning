##################################################################################
##################### GRÁFICOS DE MORAN ###################################
##################################################################################
rm(list=ls())
gc(reset = T)
# Cargue de librerías
library(readxl)
library(maptools)
library(tmap)
library(dplyr)
library(ggplot2)
library(rgdal)
library(sandwich)
library(stringr)
library(spdep)
library(spatialreg)
library(ggspatial)
# Definir directorio y cargar datos
# Cargue de base de datos
barrios_Cali <- readOGR('databases_and_shapefiles/Barrios Cali.shp')
datos_barrio_cali <- barrios_Cali@data
datos_cali <- read_excel('databases_and_shapefiles/Cali_nueva_base.xlsx')
datos_cali$Estrato_moda <- as.factor(datos_cali$Estrato_moda)
barrios_Cali@data = left_join(barrios_Cali@data, datos_cali, by = 'ID_BARRIO')
#barrios_Cali@data <- barrios_Cali@data[complete.cases(barrios_Cali@data$Adulto_mayor),]
barrios_Cali <- barrios_Cali[-c(65,112,304,326),]
#barrios_Cali <- sf::st_as_sf(barrios_Cali)

################# Indicadores de Moran local bivariado ####
####### Covid - Adulto mayor
x <- barrios_Cali$Adulto_mayor
y <- barrios_Cali$Covid
#==
# Adjacency Matrix (Queen)
nb <- poly2nb(barrios_Cali)
lw <- spdep::nb2listw(nb, style = "B", zero.policy = T)
W  <- as(lw, "symmetricMatrix")
W  <- as.matrix(W/rowSums(W))
W[which(is.na(W))] <- 0

# I de Moran bajo aleatorizacion y por MC
moran.test(barrios_Cali$Covid, listw=lw, zero.policy = TRUE)
set.seed(333)
moran.mc(barrios_Cali$Covid, listw=lw, zero.policy = TRUE,nsim=3000)
# Si hay correlación espacial significativa en global moran para Covid
# local Moran I
locmo = localmoran(barrios_Cali$Covid,lw,zero.policy = TRUE)
#==
# Programming some functions
# Bivariate Moran's I
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  yp <- (y - mean(y, na.rm=T))/sd(y, na.rm=T)
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  list(global = global, local  = as.numeric(local))
}
# Permutations for the Bivariate Moran's I
simula_moran <- function(x, y = NULL, W, nsims = 1000){
  if(is.null(y)) y = x
  n   = nrow(W)
  IDs = 1:n
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  W[which(is.na(W))] <- 0
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  list(global_sims = global_sims,
       local_sims  = local_sims)
}
#==
# Calculating the index and its simulated distribution
# for global and local values
m   <- moran_I(x, y, W)
m[[1]] # global value
m_i <- m[[2]]  # local values
local_sims <- simula_moran(x, y, W)$local_sims
# Identifying the significant values 
alpha <- .05  # for a 95% confidence interval
probs <- c(alpha/2, 1-alpha/2)
intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs)))
sig <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )
#==
# Preparing for plotting
barrios_Cali <- st_as_sf(barrios_Cali)
barrios_Cali$sig <- sig
# Identifying the LISA Temb_Tmatt
xp <- (x-mean(x))/sd(x)
yp <- (y-mean(y))/sd(y)

Covid_Adultomay <- as.character( interaction(xp > 0, W%*%yp > 0) ) 
Covid_Adultomay <- Covid_Adultomay %>% 
  str_replace_all("TRUE","Alto") %>% 
  str_replace_all("FALSE","Bajo")
Covid_Adultomay[barrios_Cali$sig==0] <- "No significativo"
barrios_Cali$Covid_Adultomay <- Covid_Adultomay
# Plotting
g1 <- ggplot() + geom_sf(data=barrios_Cali, aes(fill=Covid_Adultomay), color="black",
                         size = 0.2) +
  scale_fill_manual(values = c("red", "orange", "light blue", "dark blue", "grey95")) + 
  ggtitle("Moran's I Local Covid y Adultos mayores") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)
g1
####### Covid - Indice de cobertura de servicios bajo, moran bivariado
x <- barrios_Cali$ind_cobert_serv_bajo
y <- barrios_Cali$Covid
#==
# Adjacency Matrix (Queen)
nb <- poly2nb(barrios_Cali)
lw <- spdep::nb2listw(nb, style = "B", zero.policy = T)
W  <- as(lw, "symmetricMatrix")
W  <- as.matrix(W/rowSums(W))
W[which(is.na(W))] <- 0
#==
# Programming some functions
# Bivariate Moran's I
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  yp <- (y - mean(y, na.rm=T))/sd(y, na.rm=T)
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  list(global = global, local  = as.numeric(local))
}
# Permutations for the Bivariate Moran's I
simula_moran <- function(x, y = NULL, W, nsims = 1000){
  if(is.null(y)) y = x
  n   = nrow(W)
  IDs = 1:n
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  W[which(is.na(W))] <- 0
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  list(global_sims = global_sims,
       local_sims  = local_sims)
}
#==
# Calculating the index and its simulated distribution
# for global and local values
m   <- moran_I(x, y, W)
m[[1]] # global value
m_i <- m[[2]]  # local values
local_sims <- simula_moran(x, y, W)$local_sims
# Identifying the significant values 
alpha <- .05  # for a 95% confidence interval
probs <- c(alpha/2, 1-alpha/2)
intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs)))
sig <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )
#==
# Preparing for plotting
barrios_Cali <- st_as_sf(barrios_Cali)
barrios_Cali$sig <- sig
# Identifying the LISA Temb_Tmatt
xp <- (x-mean(x))/sd(x)
yp <- (y-mean(y))/sd(y)

Covid_ICSbajo <- as.character( interaction(xp > 0, W%*%yp > 0) ) 
Covid_ICSbajo <- Covid_ICSbajo %>% 
  str_replace_all("TRUE","Alto") %>% 
  str_replace_all("FALSE","Bajo")
Covid_ICSbajo[barrios_Cali$sig==0] <- "No significativo"
barrios_Cali$Covid_ICSbajo <- Covid_ICSbajo
# Plotting
g2 <- ggplot() + geom_sf(data=barrios_Cali, aes(fill=Covid_ICSbajo), color="black",
                         size = 0.2) +
  scale_fill_manual(values = c("red", "orange", "light blue", "dark blue", "grey95")) + 
  ggtitle("Moran's I Local Covid y ICS bajo") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

g2
####### Covid - Comorbilidad
x <- barrios_Cali$Comorbilidad
y <- barrios_Cali$Covid
#==
# Adjacency Matrix (Queen)
nb <- poly2nb(barrios_Cali)
lw <- spdep::nb2listw(nb, style = "B", zero.policy = T)
W  <- as(lw, "symmetricMatrix")
W  <- as.matrix(W/rowSums(W))
W[which(is.na(W))] <- 0
#==
# Programming some functions
# Bivariate Moran's I
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  yp <- (y - mean(y, na.rm=T))/sd(y, na.rm=T)
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  list(global = global, local  = as.numeric(local))
}
# Permutations for the Bivariate Moran's I
simula_moran <- function(x, y = NULL, W, nsims = 1000){
  if(is.null(y)) y = x
  n   = nrow(W)
  IDs = 1:n
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  W[which(is.na(W))] <- 0
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  list(global_sims = global_sims,
       local_sims  = local_sims)
}
#==
# Calculating the index and its simulated distribution
# for global and local values
m   <- moran_I(x, y, W)
m[[1]] # global value
m_i <- m[[2]]  # local values
local_sims <- simula_moran(x, y, W)$local_sims
# Identifying the significant values 
alpha <- .05  # for a 95% confidence interval
probs <- c(alpha/2, 1-alpha/2)
intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs)))
sig <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )
#==
# Preparing for plotting
barrios_Cali <- st_as_sf(barrios_Cali)
barrios_Cali$sig <- sig
# Identifying the LISA
xp <- (x-mean(x))/sd(x)
yp <- (y-mean(y))/sd(y)

Covid_comorb <- as.character( interaction(xp > 0, W%*%yp > 0) ) 
Covid_comorb <- Covid_comorb %>% 
  str_replace_all("TRUE","Alto") %>% 
  str_replace_all("FALSE","Bajo")
Covid_comorb[barrios_Cali$sig==0] <- "No significativo"
barrios_Cali$Covid_comorb <- Covid_comorb
# Plotting
g3 <- ggplot() + geom_sf(data=barrios_Cali, aes(fill=Covid_comorb), color="black",
                         size = 0.2) +
  scale_fill_manual(values = c("red", "orange", "light blue", "dark blue", "grey95")) + 
  ggtitle("Moran's I Local Covid y Comorbilidad") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

g3
####### Covid - razon_letalidad
x <- barrios_Cali$razon_letalidad
y <- barrios_Cali$Covid
#==
# Adjacency Matrix (Queen)
nb <- poly2nb(barrios_Cali)
lw <- spdep::nb2listw(nb, style = "B", zero.policy = T)
W  <- as(lw, "symmetricMatrix")
W  <- as.matrix(W/rowSums(W))
W[which(is.na(W))] <- 0
#==
# Programming some functions
# Bivariate Moran's I
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  yp <- (y - mean(y, na.rm=T))/sd(y, na.rm=T)
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  list(global = global, local  = as.numeric(local))
}
# Permutations for the Bivariate Moran's I
simula_moran <- function(x, y = NULL, W, nsims = 1000){
  if(is.null(y)) y = x
  n   = nrow(W)
  IDs = 1:n
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  W[which(is.na(W))] <- 0
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  list(global_sims = global_sims,
       local_sims  = local_sims)
}
#==
# Calculating the index and its simulated distribution
# for global and local values
m   <- moran_I(x, y, W)
m[[1]] # global value
m_i <- m[[2]]  # local values
local_sims <- simula_moran(x, y, W)$local_sims
# Identifying the significant values 
alpha <- .05  # for a 95% confidence interval
probs <- c(alpha/2, 1-alpha/2)
intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs)))
sig <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )
#==
# Preparing for plotting
barrios_Cali <- st_as_sf(barrios_Cali)
barrios_Cali$sig <- sig
# Identifying the LISA
xp <- (x-mean(x))/sd(x)
yp <- (y-mean(y))/sd(y)

Covid_razonlet <- as.character( interaction(xp > 0, W%*%yp > 0) ) 
Covid_razonlet <- Covid_razonlet %>% 
  str_replace_all("TRUE","Alto") %>% 
  str_replace_all("FALSE","Bajo")
Covid_razonlet[barrios_Cali$sig==0] <- "No significativo"
barrios_Cali$Covid_razonlet <- Covid_razonlet
# Plotting
g4 <- ggplot() + geom_sf(data=barrios_Cali, aes(fill=Covid_razonlet), color="black",
                         size = 0.2) +
  scale_fill_manual(values = c("red", "orange", "light blue", "dark blue", "grey95")) + 
  ggtitle("Moran's I Local Covid y Razón de letalidad del Covid-19") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=7),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) + 
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

g4

ggpubr::ggarrange(g1,g2,g3,g4,common.legend = TRUE, legend="bottom",ncol = 4,
                  nrow = 1)

####### Covid - Pfizer moran bivariado
x <- barrios_Cali$Pfizer
y <- barrios_Cali$Covid
#==
# Adjacency Matrix (Queen)
nb <- poly2nb(barrios_Cali)
lw <- spdep::nb2listw(nb, style = "B", zero.policy = T)
W  <- as(lw, "symmetricMatrix")
W  <- as.matrix(W/rowSums(W))
W[which(is.na(W))] <- 0

# I de Moran bajo aleatorizacion y por MC
moran.test(barrios_Cali$Covid, listw=lw, zero.policy = TRUE)
set.seed(333)
moran.mc(barrios_Cali$Covid, listw=lw, zero.policy = TRUE,nsim=3000)
# Si hay correlación espacial significativa en global moran para Covid
# local Moran I
locmo = localmoran(barrios_Cali$Covid,lw,zero.policy = TRUE)
#==
# Programming some functions
# Bivariate Moran's I
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  yp <- (y - mean(y, na.rm=T))/sd(y, na.rm=T)
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  list(global = global, local  = as.numeric(local))
}
# Permutations for the Bivariate Moran's I
simula_moran <- function(x, y = NULL, W, nsims = 1000){
  if(is.null(y)) y = x
  n   = nrow(W)
  IDs = 1:n
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  W[which(is.na(W))] <- 0
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  list(global_sims = global_sims,
       local_sims  = local_sims)
}
#==
# Calculating the index and its simulated distribution
# for global and local values
m   <- moran_I(x, y, W)
m[[1]] # global value
m_i <- m[[2]]  # local values
local_sims <- simula_moran(x, y, W)$local_sims
# Identifying the significant values 
alpha <- .05  # for a 95% confidence interval
probs <- c(alpha/2, 1-alpha/2)
intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs)))
sig <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )
#==
# Preparing for plotting
barrios_Cali <- st_as_sf(barrios_Cali)
barrios_Cali$sig <- sig
# Identifying the LISA Temb_Tmatt
xp <- (x-mean(x))/sd(x)
yp <- (y-mean(y))/sd(y)

Covid_Pfizer <- as.character( interaction(xp > 0, W%*%yp > 0) ) 
Covid_Pfizer <- Covid_Pfizer %>% 
  str_replace_all("TRUE","Alto") %>% 
  str_replace_all("FALSE","Bajo")
Covid_Pfizer[barrios_Cali$sig==0] <- "No significativo"
barrios_Cali$Covid_Pfizer <- Covid_Pfizer
# Plotting
g5 <- ggplot() + geom_sf(data=barrios_Cali, aes(fill=Covid_Pfizer), color="black",
                         size = 0.2) +
  scale_fill_manual(values = c("red", "orange", "light blue", "dark blue", "grey95")) + 
  ggtitle("Moran's I Local Covid y Vacunados con Pfizer") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)
g5

####### Covid - ss_subsidiado moran bivariado
x <- barrios_Cali$ss_subsidiado
y <- barrios_Cali$Covid
#==
# Adjacency Matrix (Queen)
nb <- poly2nb(barrios_Cali)
lw <- spdep::nb2listw(nb, style = "B", zero.policy = T)
W  <- as(lw, "symmetricMatrix")
W  <- as.matrix(W/rowSums(W))
W[which(is.na(W))] <- 0

# I de Moran bajo aleatorizacion y por MC
moran.test(barrios_Cali$Covid, listw=lw, zero.policy = TRUE)
set.seed(333)
moran.mc(barrios_Cali$Covid, listw=lw, zero.policy = TRUE,nsim=3000)
# Si hay correlación espacial significativa en global moran para Covid
# local Moran I
locmo = localmoran(barrios_Cali$Covid,lw,zero.policy = TRUE)
#==
# Programming some functions
# Bivariate Moran's I
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  yp <- (y - mean(y, na.rm=T))/sd(y, na.rm=T)
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  list(global = global, local  = as.numeric(local))
}
# Permutations for the Bivariate Moran's I
simula_moran <- function(x, y = NULL, W, nsims = 1000){
  if(is.null(y)) y = x
  n   = nrow(W)
  IDs = 1:n
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  W[which(is.na(W))] <- 0
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  list(global_sims = global_sims,
       local_sims  = local_sims)
}
#==
# Calculating the index and its simulated distribution
# for global and local values
m   <- moran_I(x, y, W)
m[[1]] # global value
m_i <- m[[2]]  # local values
local_sims <- simula_moran(x, y, W)$local_sims
# Identifying the significant values 
alpha <- .05  # for a 95% confidence interval
probs <- c(alpha/2, 1-alpha/2)
intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs)))
sig <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )
#==
# Preparing for plotting
barrios_Cali <- st_as_sf(barrios_Cali)
barrios_Cali$sig <- sig
# Identifying the LISA Temb_Tmatt
xp <- (x-mean(x))/sd(x)
yp <- (y-mean(y))/sd(y)

Covid_ss_subsidiado <- as.character( interaction(xp > 0, W%*%yp > 0) ) 
Covid_ss_subsidiado <- Covid_ss_subsidiado %>% 
  str_replace_all("TRUE","Alto") %>% 
  str_replace_all("FALSE","Bajo")
Covid_ss_subsidiado[barrios_Cali$sig==0] <- "No significativo"
barrios_Cali$Covid_ss_subsidiado <- Covid_ss_subsidiado
# Plotting
g6 <- ggplot() + geom_sf(data=barrios_Cali, aes(fill=Covid_ss_subsidiado), color="black",
                         size = 0.2) +
  scale_fill_manual(values = c("red", "orange", "light blue", "dark blue", "grey95")) + 
  ggtitle("Moran's I Local Covid y subsidiados en seguridad social") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)
g6

####### Covid - sintomatico_si moran bivariado
x <- barrios_Cali$sintomatico_si
y <- barrios_Cali$Covid
#==
# Adjacency Matrix (Queen)
nb <- poly2nb(barrios_Cali)
lw <- spdep::nb2listw(nb, style = "B", zero.policy = T)
W  <- as(lw, "symmetricMatrix")
W  <- as.matrix(W/rowSums(W))
W[which(is.na(W))] <- 0

# I de Moran bajo aleatorizacion y por MC
moran.test(barrios_Cali$Covid, listw=lw, zero.policy = TRUE)
set.seed(333)
moran.mc(barrios_Cali$Covid, listw=lw, zero.policy = TRUE,nsim=3000)
# Si hay correlación espacial significativa en global moran para Covid
# local Moran I
locmo = localmoran(barrios_Cali$Covid,lw,zero.policy = TRUE)
#==
# Programming some functions
# Bivariate Moran's I
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  yp <- (y - mean(y, na.rm=T))/sd(y, na.rm=T)
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  list(global = global, local  = as.numeric(local))
}
# Permutations for the Bivariate Moran's I
simula_moran <- function(x, y = NULL, W, nsims = 1000){
  if(is.null(y)) y = x
  n   = nrow(W)
  IDs = 1:n
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  W[which(is.na(W))] <- 0
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  list(global_sims = global_sims,
       local_sims  = local_sims)
}
#==
# Calculating the index and its simulated distribution
# for global and local values
m   <- moran_I(x, y, W)
m[[1]] # global value
m_i <- m[[2]]  # local values
local_sims <- simula_moran(x, y, W)$local_sims
# Identifying the significant values 
alpha <- .05  # for a 95% confidence interval
probs <- c(alpha/2, 1-alpha/2)
intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs)))
sig <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )
#==
# Preparing for plotting
barrios_Cali <- st_as_sf(barrios_Cali)
barrios_Cali$sig <- sig
# Identifying the LISA Temb_Tmatt
xp <- (x-mean(x))/sd(x)
yp <- (y-mean(y))/sd(y)

Covid_sintomatico_si <- as.character( interaction(xp > 0, W%*%yp > 0) ) 
Covid_sintomatico_si <- Covid_sintomatico_si %>% 
  str_replace_all("TRUE","Alto") %>% 
  str_replace_all("FALSE","Bajo")
Covid_sintomatico_si[barrios_Cali$sig==0] <- "No significativo"
barrios_Cali$Covid_sintomatico_si <- Covid_sintomatico_si
# Plotting
g7 <- ggplot() + geom_sf(data=barrios_Cali, aes(fill=Covid_sintomatico_si), color="black",
                         size = 0.2) +
  scale_fill_manual(values = c("red", "orange", "light blue", "dark blue", "grey95")) + 
  ggtitle("Moran's I Local Covid y sintomáticos por barrio") + 
  theme(panel.background = element_blank(), plot.title = element_text(size=10),
        axis.text.x = element_text(size = 6), legend.title=element_blank(),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)
g7

ggpubr::ggarrange(g5,g6,g7,common.legend = TRUE, legend="bottom",ncol = 3,nrow = 1)
################# Graficos de relación de variables
custom_corr_plot <- function(variable1, variable2, df, alpha=0.3){
  p <- df %>%
    mutate(
      # Truco para que se ponga el título estilo facet
      title = paste(toupper(variable2), "vs", toupper(variable1))
    ) %>%
    ggplot(aes(x = !!sym(variable1), y = !!sym(variable2))) + 
    geom_point(alpha = alpha) +
    # Tendencia no lineal
    geom_smooth(se = FALSE, method = "gam", formula =  y ~ splines::bs(x, 3)) +
    # Tendencia lineal
    geom_smooth(se = FALSE, method = "lm", color = "firebrick") +
    labs(x = variable1, y = variable2) +
    facet_grid(. ~ title) +
    theme_bw() +
    theme(strip.text = element_text(colour = "black", size = 8, face = 2))
  return(p)
}
variables_continuas <- colnames(datos_cali[,c(4:11)])

plots <- purrr::map(
  .x = variables_continuas,
  .f = custom_corr_plot,
  variable2 = "Covid",
  df = datos_cali[,c(4:11)]
)
plots

# Gráfico de casos de covid por proporción trabajando, segmentado por estrato
# y población
ggplot(datos_cali, aes(x=Proporcion_trabajando, y=Covid, color = Estrato_moda,
                       size = Poblacion_total)) +
  geom_point(alpha=0.4) +
  scale_size(range = c(1, 7), name="Población") +
  theme(legend.key.size = unit(3, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size


library(ggpubr)
# Grouped Scatter plot with marginal density plots
ggscatterhist(
  datos_cali, x = "Proporcion_trabajando", y = "Covid",
  color = "Estrato_moda", size = 3, alpha = 0.6,
  margin.params = list(fill = "Estrato_moda", color = "black", size = 0.2)
)

##################################################################################
############################## GRÁFICOS INDIVIDUALES##############################
##################################################################################
rm(list=ls())
gc(reset = T)
library(readxl)
library(maptools)
library(tmap)
library(dplyr)
library(ggplot2)
library(rgdal)
library(sandwich)
library(stringr)
library(spdep)
library(spatialreg)
library(ggspatial)
# Definir directorio y cargar datos

barrios_Cali <- readOGR('databases_and_shapefiles/Barrios Cali.shp')
datos_barrio_cali <- barrios_Cali@data
datos_cali <- read_excel('databases_and_shapefiles/Cali_nueva_base.xlsx')
datos_cali$Estrato_moda <- as.factor(datos_cali$Estrato_moda)
barrios_Cali@data = left_join(barrios_Cali@data, datos_cali, by = 'ID_BARRIO')
barrios_Cali <- sf::st_as_sf(barrios_Cali)

# Casos positivos de covid-19 en Cali por barrio
barrios_Cali$Covid[is.na(barrios_Cali$Covid)] <- 0
options(scipen=999)
levels(cut(x = barrios_Cali$Covid, breaks = 10, include.lowest = F,dig.lab = 5))
barrios_Cali %>% mutate(qrr = cut(Covid,
                                  breaks = 10,
                                  labels = c("(0, 375]",
                                             "(375, 751]",
                                             "(751, 1.126]",
                                             "(1.126, 1.502]",
                                             "(1.502, 1.877]",
                                             "(1.877, 2.252]",
                                             "(2.252, 2.628]",
                                             "(2.628, 3.003]",
                                             "(3.003, 3.379]",
                                             "(3.379, 3.758]"), include.lowest = F)) %>%
  ggplot() + geom_sf(aes(fill = qrr), color="black", size = 0.2) +
  scale_colour_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu", na.value="grey") +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  ggtitle ("Casos positivos de covid-19 en Cali por barrio") + labs(fill = "N° de casos") +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

# Casos positivos de covid-19 en Cali por barrio
barrios_Cali$Covid_nueva_base[is.na(barrios_Cali$Covid_nueva_base)] <- 0
options(scipen=999)
levels(cut(x = barrios_Cali$Covid_nueva_base, breaks = 10, include.lowest = F,dig.lab = 5))
barrios_Cali %>% mutate(qrr = cut(Covid_nueva_base,
                                  breaks = 10,
                                  labels = c("(0, 673]",
                                             "(673, 1.347]",
                                             "(1.347, 2.020]",
                                             "(2.020, 2.693]",
                                             "(2.693, 3.366]",
                                             "(3.366, 4.040]",
                                             "(4.040, 4.713]",
                                             "(4.713, 5.386]",
                                             "(5.386, 6.060]",
                                             "(6.060, 6.740]"), include.lowest = F)) %>%
  ggplot() + geom_sf(aes(fill = qrr), color="black", size = 0.2) +
  scale_colour_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu", na.value="grey") +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  ggtitle ("Casos positivos de covid-19 en Cali por barrio") + labs(fill = "N° de casos") +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

# Estratos en Cali por barrio
ggplot() + geom_sf(data = barrios_Cali, aes(fill = Estrato_moda), color="black",
                   size = 0.2) +
  scale_colour_brewer(palette = "RdBu" ) +
  scale_fill_brewer(palette = "RdBu", na.value="grey") +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2))+
  ggtitle ("Estratos en Cali por barrio") + labs(fill = "Estrato") +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

# Población en Cali por barrio
barrios_Cali$Poblacion_total[is.na(barrios_Cali$Poblacion_total)] <- 0
options(scipen=999)
#cut(x = barrios_Cali$Poblacion_total, breaks = 10, include.lowest = F,dig.lab = 5)
barrios_Cali %>% mutate(qrr = cut(Poblacion_total,
                                  breaks = 10,
                                  labels = c("(0, 3.343]",
                                             "(3.343, 6.685]",
                                             "(6.685, 10.028]",
                                             "(10.028, 13.370]",
                                             "(13.370, 16.713]",
                                             "(16.713, 20.056]",
                                             "(20.056, 23.398]",
                                             "(23.398, 26.741]",
                                             "(26.741, 30.083]",
                                             "(30.083, 33.459]"), include.lowest = F)) %>%
  ggplot() + geom_sf(aes(fill = qrr), color="black", size = 0.2) +
  scale_colour_brewer(palette = "RdBu" ) +
  scale_fill_brewer(palette = "RdBu", na.value="grey") +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  ggtitle ("Población en Cali por barrio") + labs(fill = "Población") +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

# Adultos mayores en Cali por barrio
barrios_Cali$Adulto_mayor[is.na(barrios_Cali$Adulto_mayor)] <- 0
levels(cut(x = barrios_Cali$Adulto_mayor, breaks = 10, include.lowest = F,dig.lab = 5))
barrios_Cali %>% mutate(qrr = cut(Adulto_mayor,
                                  breaks = 10,
                                  labels = c("(0, 476]",
                                             "(476, 953]",
                                             "(953, 1.429]",
                                             "(1.429, 1.905]",
                                             "(1.905, 2.381]",
                                             "(2.381, 2.858]",
                                             "(2.858, 3.334]",
                                             "(3.334, 3.810]",
                                             "(3.810, 4.287]",
                                             "(4.287, 4.768]"), include.lowest = F)) %>%
  ggplot() + geom_sf(aes(fill = qrr), color="black", size = 0.2) +
  scale_colour_brewer(palette = "RdBu" ) +
  scale_fill_brewer(palette = "RdBu", na.value="grey") +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  ggtitle ("Adultos mayores en Cali por barrio") + labs(fill = "Adultos mayores") +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

# ersonas de estratos bajos en cobertura de servicios, por barrio
barrios_Cali$ind_cobert_serv_bajo[is.na(barrios_Cali$ind_cobert_serv_bajo)] <- 0
levels(cut(x = barrios_Cali$ind_cobert_serv_bajo, breaks = 10, include.lowest = F,
           dig.lab = 5))
barrios_Cali %>% mutate(qrr = cut(ind_cobert_serv_bajo,
                                  breaks = 10,
                                  labels = c("(0, 834]",
                                             "(834, 1.662]",
                                             "(1.662, 2.490]",
                                             "(2.490, 3.318]",
                                             "(3.318, 4.145]",
                                             "(4.145, 4.973]",
                                             "(4.973, 5.801]",
                                             "(5.801, 6.629]",
                                             "(6.629, 7.457]",
                                             "(7.457, 8.293]"), include.lowest = F)) %>%
  ggplot() + geom_sf(aes(fill = qrr), color="black", size = 0.2) +
  scale_colour_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu", na.value="grey",direction = -1) +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  ggtitle ("Personas de estratos bajos en cobertura de servicios, por barrio") + 
  labs(fill = "Número de personas en CS bajo") +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

# Personas con comorbilidad por barrio
barrios_Cali$Comorbilidad[is.na(barrios_Cali$Comorbilidad)] <- 0
levels(cut(x = barrios_Cali$Comorbilidad, breaks = 10, include.lowest = F, dig.lab = 5))
barrios_Cali %>% mutate(qrr = cut(Comorbilidad,
                                  breaks = 10,
                                  labels = c("(0, 98]",
                                             "(98, 197]",
                                             "(197, 295]",
                                             "(295, 394]",
                                             "(394, 492]",
                                             "(492, 590]",
                                             "(590, 689]",
                                             "(689, 787]",
                                             "(787, 886]",
                                             "(886, 985]"), include.lowest = F)) %>%
  ggplot() + geom_sf(aes(fill = qrr), color="black", size = 0.2) +
  scale_colour_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu", na.value="grey") +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  ggtitle("Personas con comorbilidad por barrio") + labs(fill = "Número de personas con comorbilidad") +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

# Razón de letalidad del Covid por barrio
barrios_Cali$razon_letalidad[is.na(barrios_Cali$razon_letalidad)] <- 0
levels(cut(x = barrios_Cali$razon_letalidad, breaks = 10, include.lowest = F, dig.lab = 5))
barrios_Cali %>% mutate(qrr = cut(razon_letalidad,
                                  breaks = 10,
                                  labels = c("(0, 0.011]",
                                             "(0.011, 0.022]",
                                             "(0.022, 0.033]",
                                             "(0.033, 0.044]",
                                             "(0.044, 0.056]",
                                             "(0.056, 0.067]",
                                             "(0.067, 0.078]",
                                             "(0.078, 0.089]",
                                             "(0.089, 0.1]",
                                             "(0.1, 0.112]"), include.lowest = F)) %>%
  ggplot() + geom_sf(aes(fill = qrr), color="black", size = 0.2) +
  scale_colour_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu", na.value="grey") +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  ggtitle("Razón de letalidad del Covid por barrio") + labs(fill = "Razón de letalidad") +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

#Proporción de ocupados por barrio
#barrios_Cali$Proporcion_trabajando[is.na(barrios_Cali$Proporcion_trabajando)] <- 0
levels(cut(x = barrios_Cali$Proporcion_trabajando, breaks = 10, include.lowest = F, 
           dig.lab = 5))
barrios_Cali %>% mutate(qrr = cut(Proporcion_trabajando,
                                  breaks = 10,
                                  labels = c("(0.358, 0.392]",
                                             "(0.392, 0.426]",
                                             "(0.426, 0.460]",
                                             "(0.460, 0.494]",
                                             "(0.494, 0.528]",
                                             "(0.528, 0.562]",
                                             "(0.562, 0.596]",
                                             "(0.596, 0.629]",
                                             "(0.629, 0.663]",
                                             "(0.663, 0.698]"), include.lowest = F)) %>%
  ggplot() + geom_sf(aes(fill = qrr), color="black", size = 0.2) +
  scale_colour_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu", na.value="grey") +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  ggtitle("Proporción de ocupados por barrio") + labs(fill = "Proporción de ocupados") +
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

# N° de vacunados con Pfizer en Cali por barrio
barrios_Cali$Pfizer[is.na(barrios_Cali$Pfizer)] <- 0
options(scipen=999)
levels(cut(x = barrios_Cali$Pfizer, breaks = 10, include.lowest = F, dig.lab = 5))
barrios_Cali %>% mutate(qrr = cut(Pfizer,
                                  breaks = 10,
                                  labels = c("(0, 209]",
                                             "(209, 418]",
                                             "(418, 627]",
                                             "(627, 836]",
                                             "(836, 1.045]",
                                             "(1.045, 1.254]",
                                             "(1.254, 1.463]",
                                             "(1.463, 1.672]",
                                             "(1.672, 1.881]",
                                             "(1.881, 2.092]"), include.lowest = F)) %>%
  ggplot() + geom_sf(aes(fill = qrr), color="black", size = 0.2) +
  scale_colour_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu", na.value="grey") +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  ggtitle ("N° de vacunados con Pfizer en Cali por barrio") + labs(fill = "N° de vacunas") +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

# Subsidiados seguridad social en Cali por barrio
barrios_Cali$ss_subsidiado[is.na(barrios_Cali$ss_subsidiado)] <- 0
options(scipen=999)
levels(cut(x = barrios_Cali$ss_subsidiado, breaks = 10, include.lowest = F, dig.lab = 5))
barrios_Cali %>% mutate(qrr = cut(ss_subsidiado,
                                  breaks = 10,
                                  labels = c("(0, 81]",
                                             "(81, 161]",
                                             "(161, 242]",
                                             "(242, 323]",
                                             "(323, 403]",
                                             "(403, 484]",
                                             "(484, 565]",
                                             "(565, 646]",
                                             "(646, 726]",
                                             "(726, 808]"), include.lowest = F)) %>%
  ggplot() + geom_sf(aes(fill = qrr), color="black", size = 0.2) +
  scale_colour_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu", na.value="grey") +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  ggtitle ("Subsidiados seguridad social en Cali por barrio") + labs(fill = "N° de subsidiados") +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)

# Sintomáticos en Cali por barrio
barrios_Cali$sintomatico_si[is.na(barrios_Cali$sintomatico_si)] <- 0
options(scipen=999)
levels(cut(x = barrios_Cali$sintomatico_si, breaks = 10, include.lowest = F, dig.lab = 5))
barrios_Cali %>% mutate(qrr = cut(sintomatico_si,
                                  breaks = 10,
                                  labels = c("(0, 347]",
                                             "(347, 695]",
                                             "(695, 1.042]",
                                             "(1.042, 1.389]",
                                             "(1.389, 1.736]",
                                             "(1.736, 2.084]",
                                             "(2.084, 2.431]",
                                             "(2.431, 2.778]",
                                             "(2.778, 3.126]",
                                             "(3.126, 3.476]"), include.lowest = F)) %>%
  ggplot() + geom_sf(aes(fill = qrr), color="black", size = 0.2) +
  scale_colour_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu", na.value="grey") +
  theme(panel.background = element_blank(), text = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.2)) +
  ggtitle ("Sintomáticos en Cali por barrio") + labs(fill = "N° de sintomáticos") +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE)
######################################################################################
######################################################################################
rm(list=ls())
gc(reset = T)
library(readxl)
library(maptools)
library(tmap)
library(dplyr)
library(ggplot2)
library(rgdal)
library(sandwich)
# Definir directorio y cargar datos
barrios_Cali <- readOGR('databases_and_shapefiles/Barrios Cali.shp')
datos_barrio_cali <- barrios_Cali@data
datos_cali <- read_excel('databases_and_shapefiles/Base_final.xlsx')

# hacer el cruce de codigos de barrios del DataFrame con el mapa
cod_barrio <- as.data.frame(datos_barrio_cali$ID_BARRIO)
names(cod_barrio) <- 'Codigo_unico'
data_cali <- left_join(cod_barrio, datos_cali, by = 'Codigo_unico')
data_cali <- data_cali[,-c(2,3)]
names(data_cali)

##----Ingresar los datos al shapefile-------------
##Se reducen los nombres
barrios_Cali$Estrato_moda <- data_cali$Estrato_moda
barrios_Cali$Poblacion_total <- data_cali$Poblacion_total
barrios_Cali$Hogares <- data_cali$Hogares
barrios_Cali$Edad_T_80omas <- data_cali$Edad_T_80omas
barrios_Cali$Proporcion_afiliados <- data_cali$Proporcion_afiliados
summary(barrios_Cali$Hogares)
summary(barrios_Cali$Estrato_moda)
summary(barrios_Cali$Proporcion_afiliados)

### Graficos 
## Grafico de estrato por barrio en Cali
barrios_Cali1 <- sf::st_as_sf(barrios_Cali)

ggplot() + geom_sf(data = barrios_Cali1, aes(fill = as.factor(Estrato_moda)), 
                   color = "grey85") + theme(text = element_text(size = 10))+
  ggtitle ("Estratos en Cali por barrio") + labs(fill = "Estrato")

## Grafico de edad de 80 años o mas en Cali 
barrios_Cali1 <- barrios_Cali
barrios_Cali1$Edad_T_80omas <- as.character(cut(barrios_Cali1$Edad_T_80omas,
                                   breaks = c(0, 175, 350, 525, 695),
                                   labels = c("0-175",
                                              "176-350",
                                              "351-525",
                                              "526-695")))

barrios_Cali1 <- sf::st_as_sf(barrios_Cali1)

ggplot() + geom_sf(data = barrios_Cali1, aes(fill = Edad_T_80omas),
                   color = "grey85") + scale_fill_manual(values = c( "dark blue",
                                                                     "light blue", 
                                                                     "pink", "red"))+
  theme(text = element_text(size = 10))+
  ggtitle ("Población de 80 años o más en Cali") + labs(fill = "Edad de 80 años o más")

## Grafico de poblacion total 
summary(barrios_Cali$Poblacion_total)

tm_shape(barrios_Cali) + tm_fill("Poblacion_total", title=c("Población total 2018"), n=7, 
                                 palette="viridis") + 
  tm_borders(lwd=0.5) + tm_layout(legend.text.size = 0.6, 
                                  legend.position = c(0.7,0.02), 
                                  panel.labels="Población total por barrios en Cali",
                                  bg="grey90")

barrios_Cali1 <- barrios_Cali
barrios_Cali1$Poblacion_total <- as.character(cut(barrios_Cali1$Poblacion_total,
                                               breaks = c(0, 6700, 13400, 20100, 26800, 33500),
                                               labels = c("0-6700",
                                                          "6701-13400",
                                                          "13401-20100",
                                                          "20101-26800",
                                                          "26801-33500")))

barrios_Cali1 <- sf::st_as_sf(barrios_Cali1)
#barrios_Cali1$Poblacion_total <- with(barrios_Cali1, reorder(Poblacion_total, ID_BARRIO))
# Arreglar el orden del grafico
ggplot() + geom_sf(data = barrios_Cali1, aes(fill = Poblacion_total),
                   color = "grey85") + scale_fill_manual(values = c( "dark blue",
                                                                     "light blue",
                                                                     "pink",
                                                                     "red", "blue"))+
  theme(text = element_text(size = 10))+
  ggtitle ("Población por barrio en Cali") + labs(fill = "Población por barrio")

## Graficos afiliacion a salud  
barrios_Cali1 <- barrios_Cali
barrios_Cali1$Proporcion_afiliados <- as.character(cut(barrios_Cali1$Proporcion_afiliados,
                                                breaks = c(0, 0.20, 0.40, 0.60, 0.80, 1),
                                                labels = c("0%-20%",
                                                           "20%-40%",
                                                           "40%-60%",
                                                           "60%-80%",
                                                           "80%-100%")))

barrios_Cali1 <- sf::st_as_sf(barrios_Cali1)


ggplot() + geom_sf(data = barrios_Cali1, aes(fill = Proporcion_afiliados),
                   color = "grey85") + scale_fill_manual(values = c("red1","red2",
                                                                    "red4","rosybrown",
                                                                    "rosybrown1","rosybrown3"))+
  theme(text = element_text(size = 10))+
  ggtitle ("Proporción de afiliados por barrio en Cali") + labs(fill = "Porcentaje de afiliados")

###### graficar el mapa 
barrios_Cali$Covid <- data_cali$Covid

tm_shape(barrios_Cali) + tm_fill("Covid", title=c("Contagios por Covid-19"), n=7, 
                                 palette="viridis") + 
  tm_borders(lwd=0.5) + tm_layout(legend.text.size = 0.6, 
                                  legend.position = c(0.7,0.02), 
                                  panel.labels="Covid-19 por barrios en Cali",
                                  bg="grey90")

library(leaflet)
paleta <- colorRampPalette(c('green', 'orange', 'red'))(100)
q_pal <- colorNumeric(paleta, barrios_Cali$Covid)
#q_pal <- colorNumeric("OrRd", barrios_Cali$Covid)

labels <- sprintf(
  "<strong>%s</strong><br/>%g casos",
  barrios_Cali$NOMBRE, barrios_Cali$Covid
) %>% lapply(htmltools::HTML)

leaflet(barrios_Cali) %>% addTiles() %>% 
  addPolygons(color = '#444444', weight = 1, smoothFactor = 0.5, opacity = 1.0,
              fillOpacity = 0.5, fillColor = q_pal(barrios_Cali$Covid), 
              highlightOptions = highlightOptions(color = "white", weight = 2, 
                                                  bringToFront = TRUE), 
              label = labels, labelOptions = labelOptions(
                style = list("font-weight" = "normal", padding = "3px 8px"),
                textsize = "15px", direction = "auto")) %>% 
  addLegend(pal = q_pal, values = ~Covid)
#https://rstudio.github.io/leaflet/choropleths.html

# Segundo mapa, controlado por la poblacion
barrios_Cali$Covid_poblacion <- data_cali$Covid/data_cali$Poblacion_total

tm_shape(barrios_Cali) + tm_fill("Covid_poblacion", title=c("Contagios por Covid-19"), n=7, 
                                 palette="viridis") + 
  tm_borders(lwd=0.5) + tm_layout(legend.text.size = 0.6, 
                                  legend.position = c(0.7,0.02), 
                                  panel.labels="Covid-19 por barrios en Cali",
                                  bg="grey90")

library(leaflet)
paleta <- colorRampPalette(c('red', 'orange', 'green'))(100)
q_pal <- colorNumeric(paleta, barrios_Cali$Covid_poblacion)
#q_pal <- colorNumeric("OrRd", barrios_Cali$Covid)

labels <- sprintf(
  "<strong>%s</strong><br/>%g casos",
  barrios_Cali$NOMBRE, barrios_Cali$Covid_poblacion
) %>% lapply(htmltools::HTML)

leaflet(barrios_Cali) %>% addTiles() %>% 
  addPolygons(color = '#444444', weight = 1, smoothFactor = 0.5, opacity = 1.0,
              fillOpacity = 0.5, fillColor = q_pal(barrios_Cali$Covid_poblacion), 
              highlightOptions = highlightOptions(color = "white", weight = 2, 
                                                  bringToFront = TRUE), 
              label = labels, labelOptions = labelOptions(
                style = list("font-weight" = "normal", padding = "3px 8px"),
                textsize = "15px", direction = "auto")) %>% 
  addLegend(pal = q_pal, values = ~Covid_poblacion)
#https://rstudio.github.io/leaflet/choropleths.html


# Tercer mapa, controlado por riesgo relativo
barrios_Cali$Riesgo_rel <- data_cali$Covid/mean(data_cali$Covid, na.rm = T)

tm_shape(barrios_Cali) + tm_fill("Riesgo_rel", title=c("Contagios por Covid-19, RR"), n=7, 
                                 palette="viridis") + 
  tm_borders(lwd=0.5) + tm_layout(legend.text.size = 0.6, 
                                  legend.position = c(0.7,0.02), 
                                  panel.labels="RR de Covid-19 por barrios en Cali",
                                  bg="grey90")

library(leaflet)
paleta <- colorRampPalette(c('red', 'orange', 'green'))(100)
q_pal <- colorNumeric(paleta, barrios_Cali$Riesgo_rel)
#q_pal <- colorNumeric("OrRd", barrios_Cali$Covid)

labels <- sprintf(
  "<strong>%s</strong><br/>%g casos",
  barrios_Cali$NOMBRE, barrios_Cali$Riesgo_rel
) %>% lapply(htmltools::HTML)

leaflet(barrios_Cali) %>% addTiles() %>% 
  addPolygons(color = '#444444', weight = 1, smoothFactor = 0.5, opacity = 1.0,
              fillOpacity = 0.5, fillColor = q_pal(barrios_Cali$Riesgo_rel), 
              highlightOptions = highlightOptions(color = "white", weight = 2, 
                                                  bringToFront = TRUE), 
              label = labels, labelOptions = labelOptions(
                style = list("font-weight" = "normal", padding = "3px 8px"),
                textsize = "15px", direction = "auto")) %>% 
  addLegend(pal = q_pal, values = ~Riesgo_rel)
#https://rstudio.github.io/leaflet/choropleths.html


################ Análisis de los datos
#### Planteamiento de modelos con base a la información de Python y R
# Población total <- lineal positiva
# E60 <- lineal positiva
# VIVIENDAS_TOTAL <- lineal positiva
# HOGARES <- lineal positiva
# ADULTO_MAYOR <- lineal positiva


cuadraticas_edad <- datos_cali[,c(230,55:71,89:105,123:139)] # prop edades
colnames(datos_cali[,c(16:29,38:54,72:88,106:122,230)]) # edades
cor_edades <- as.data.frame(cor(datos_cali[,c(15:30,32,37:54,72:88,106:122,230)], 
                                use = 'complete.obs'))
cor_edades <- cor_edades[order(-cor_edades$Covid),]

options(scipen = 999)
p_valores_nolineal <- lapply(datos_cali[,c(230,55:71,89:105,123:139)],
                             function(x) (lmtest::coeftest(glm(formula = datos_cali$Covid~x+I(x^2),
                                                               family = poisson), 
                                                           vcov.=vcovHC(glm(formula = datos_cali$Covid~x+I(x^2),
                                                                            family = poisson), type="const"))))
p_valores_nolineal # De los térnimos cuadráticos el que mejor ajustó fue:
modelo <- glm(Covid ~ Prop_Edad_M_0a4 + I(Prop_Edad_M_0a4^2), family = poisson, 
               data = datos_cali)
summary(modelo) # Mejor?
ss <- lmtest::coeftest(modelo, vcov.=vcovHC(modelo, type="const"))
ss

lmtest::coeftest(glm(Covid ~ Prop_Edad_H_0a4 + I(Prop_Edad_H_0a4^2),
                     family = poisson, data = datos_cali), 
                 vcov.=vcovHC(glm(Covid ~ Prop_Edad_H_0a4 + I(Prop_Edad_H_0a4^2),
                                  family = poisson, data = datos_cali), type="const"))

# Ahora con todas las variables:
glms_covid <- lapply(datos_cali[,c(6:227,230)],
                     function(x) (lmtest::coeftest(glm(formula = datos_cali$Covid~x,
                                                       family = poisson),
                                                   vcov.=vcovHC(glm(formula = datos_cali$Covid~x,
                                                                    family = poisson), type="const"))))

Variables_sig_glm <- subset(x = glms_covid,
                            subset = sapply(datos_cali[,c(6:227,230)],
                                            function(x) (lmtest::coeftest(glm(formula = datos_cali$Covid~x,
                                                                              family = poisson),
                                                                          vcov.=vcovHC(glm(formula = datos_cali$Covid~x,
                                                                                           family = poisson), type="const"))[8]) * 100) <5 &
                              sapply(datos_cali[,c(6:227,230)], 
                                     function(x) (lmtest::coeftest(glm(formula = datos_cali$Covid~x,
                                                                       family = poisson),
                                                                   vcov.=vcovHC(glm(formula = datos_cali$Covid~x,
                                                                                    family = poisson), type="const"))[8]) * 100) != 0)
Variables_sig_glm # Linealmente hay 101 variables que son significativas a menos del 10% 

p_valores <- data.frame()
for(i in 1:length(Variables_sig_glm)){
  p_valores[i,1] <- Variables_sig_glm[[i]][8]
  rownames(p_valores)[i] <- names(Variables_sig_glm[i])
}
ggplot(data = datos_cali, aes(x = Recuperados, y = Covid)) + geom_point()
ggplot(data = datos_cali, aes(x = Asistencia_esc_universi, y = Covid)) + geom_point()


correlaciones <- cor(x = datos_cali[,rownames(p_valores)[1:165]], 
                     use = 'complete.obs')
#Variables:
View(correlaciones[c('Poblacion_total',
                    'Viviendas_total', # Preferible que hogares
                    'Hogares',
                    'Adulto_mayor',
                    'Tam_promedio_hog', # Correlaci?n baja con las de arriba
                    'Personas_con_dificultades', # Correlaci?n moderadamente alta con las de arriba
                    'Asistencia_esc_universi', # Correlaci?n mod baja con las de arriba, pero con Covid es positiva, rara
                    'Tipo_vivienda_cuarto', # Correlaci?n baja con las de arriba, con covid positiva
                    'Cobertura_Internet_fijoomovil', # Correlaci?n baja con las de arriba
                    'Cobertura_Gas_natural', # Tambi?n podria servir
                    'Alimentos_cuarto_usado_cocinar_y_dormir',
                    'Hacinam_dormit', # Importante, correlaci?n baja con las de arriba, con covid positiva
                    'Tasa_Pobre_nbi',
                    'ind_cobert_serv_bajo',
                    'Comorbilidad',
                    'Diabetes',
                    'HTA',
                    'Muertos',
                    'Recuperados',
                    'razon_letalidad'),
                  c('Poblacion_total',
                    'Viviendas_total',
                    'Hogares',
                    'Adulto_mayor',
                    'Tam_promedio_hog', 
                    'Personas_con_dificultades', 
                    'Asistencia_esc_universi',
                    'Tipo_vivienda_cuarto',
                    'Cobertura_Internet_fijoomovil',
                    'Cobertura_Gas_natural',
                    'Alimentos_cuarto_usado_cocinar_y_dormir',
                    'Hacinam_dormit',
                    'Tasa_Pobre_nbi',
                    'ind_cobert_serv_bajo',
                    'Comorbilidad',
                    'Diabetes',
                    'HTA',
                    'Muertos',
                    'Recuperados',
                    'razon_letalidad')])


#'Poblacion_total',
#'Viviendas_total',
#'Hogares',
#'Adulto_mayor',
#'Tam_promedio_hog', #no
#'Personas_con_dificultades', 
#'Asistencia_esc_universi', #si
#'Tipo_vivienda_cuarto',
#'Cobertura_Internet_fijoomovil',
#'Cobertura_Gas_natural',
#'Alimentos_cuarto_usado_cocinar_y_dormir',
#'Tasa_Pobre_nbi',
#'ind_cobert_serv_bajo', #si
#'Comorbilidad',
#'Diabetes',
#'HTA',
#'Muertos',
#'Recuperados',
#'razon_letalidad'
#'Prop_Edad_M_0a4'

datos_cali$'Prop_Edad_M_0a4^2' <- I(datos_cali$Prop_Edad_M_0a4^2)
datos_cali2 <- datos_cali[,c('Poblacion_total','Viviendas_total','Hogares',
                             'Adulto_mayor','Tam_promedio_hog',
                             'Personas_con_dificultades', 
                             'Asistencia_esc_universi',
                             'Tipo_vivienda_cuarto',
                             'Cobertura_Internet_fijoomovil',
                             'Cobertura_Gas_natural',
                             'Alimentos_cuarto_usado_cocinar_y_dormir',
                             'Tasa_Pobre_nbi',
                             'ind_cobert_serv_bajo',
                             'Comorbilidad',
                             'Diabetes',
                             'HTA',
                             'Muertos',
                             'Recuperados',
                             'razon_letalidad',
                             'Prop_Edad_M_0a4','Prop_Edad_M_0a4^2','Covid')]
datos_cali$Estrato_moda <- as.factor(datos_cali$Estrato_moda)
attach(datos_cali)
mod1 <- glm(formula = Covid ~ Adulto_mayor+ind_cobert_serv_bajo+log(Comorbilidad)+
              Cobertura_Internet_fijoomovil%in%Estrato_moda, 
            family = poisson, data = datos_cali, offset = razon_letalidad) #offset = Tasa_Pobre_nbi
summary(mod1)
lmtest::coeftest(mod1, vcov.=vcovHC(mod1, type="const"))
anova(mod1,test="Chisq")

library(ResourceSelection)
#hoslem.test(Covid,fitted(mod1))
pchisq(mod1$deviance, df=mod1$df.residual, lower.tail=FALSE)

library(leaps)
as <- regsubsets(Covid ~ ., data = datos_cali[,c(6:227,230)], nvmax = 20, 
                 method = 'backward')
step(glm(formula = Covid ~ ., 
         family = poisson, data = datos_cali2), direction = 'forward')

################ Regresi?n ridge
library(glmnet)

# Conversi?n a matriz modelo de los datos de train y test
datos.train.mat <- model.matrix(Covid ~ ., data = datos_cali2[1:270,])
datos.test.mat <- model.matrix(Covid ~ ., data = datos_cali2[271:337,])

# Conjunto de valores de lambda 
lambda = seq(from = 100, to = 0.001, length = 2000)
# 10-fold cross validation para obtener el mejor lambda
set.seed(12)
cv.ridge <- cv.glmnet(x = datos.train.mat, y = na.omit(datos_cali2[1:270,])$Covid, alpha = 0, 
                      lambda = lambda, thresh = 1e-12, type.measure="mse")

plot(cv.ridge)

cv.ridge$lambda.min
modelo.ridge.train <- glmnet(x = datos.train.mat, y = na.omit(datos_cali2[1:270,])$Covid, 
                             alpha = 0, lambda = lambda, thresh = 1e-12)

plot(modelo.ridge.train, xvar = "lambda", label = TRUE)

# Predicciones del modelo con los datos de test y el mejor lambda
pred.modelo.ridge <- predict(modelo.ridge.train, s = 93.34674, newx = datos.test.mat)
# Test error (MSE)
test.MSE.ridge <- mean((pred.modelo.ridge - na.omit(datos_cali2[271:337,])$Covid)^2)
test.MSE.ridge

# Se excluye la primera columna
modelo.ridge <- glmnet(x = model.matrix(Covid ~ ., data = na.omit(datos_cali2))[,-1], 
                       y = na.omit(datos_cali2)$Covid, alpha = 93.34674)
# Coeficientes del modelo
predict(modelo.ridge, type = "coefficients", s = 93.34674)
################


################### an?lisis con variable categ?rica
mod1 <- glm(formula = Covid ~ Adulto_mayor+ind_cobert_serv_bajo+log(Comorbilidad)+
              Si_Asistencia_escolar%in%Estrato_moda, 
            family = poisson, data = datos_cali, offset = razon_letalidad) #offset = Tasa_Pobre_nbi
summary(mod1)
lmtest::coeftest(mod1, vcov.=vcovHC(mod1, type="const"))
anova(mod1,test="Chisq")

# Ahora con todas las variables:
glms_covid_cat <- lapply(datos_cali[,c(6:227,230)],
                         function(x) (lmtest::coeftest(glm(formula = datos_cali$Covid~Adulto_mayor+ind_cobert_serv_bajo+log(Comorbilidad)+x%in%Estrato_moda,
                                                           family = poisson),
                                                       vcov.=vcovHC(glm(formula = datos_cali$Covid~Adulto_mayor+ind_cobert_serv_bajo+log(Comorbilidad)+x%in%Estrato_moda,
                                                                        family = poisson), type="const"))))

Variables_sig_glm_cat <- subset(x = glms_covid_cat,
                                subset = sapply(datos_cali[,c(6:227,230)],
                                                function(x) mean((lmtest::coeftest(glm(formula = datos_cali$Covid~Adulto_mayor+ind_cobert_serv_bajo+log(Comorbilidad)+x%in%Estrato_moda,
                                                                                  family = poisson),
                                                                              vcov.=vcovHC(glm(formula = datos_cali$Covid~Adulto_mayor+ind_cobert_serv_bajo+log(Comorbilidad)+x%in%Estrato_moda,
                                                                                               family = poisson), type="const"))[35:40]) * 100)) <10)
Variables_sig_glm_cat 

p_valores_cat <- data.frame()
for(i in 1:length(Variables_sig_glm_cat)){
  p_valores_cat[i,1] <- mean(Variables_sig_glm_cat[[i]][35:40])
  rownames(p_valores_cat)[i] <- names(Variables_sig_glm_cat[i])
}

lmtest::coeftest(glm(formula = datos_cali$Covid~Adulto_mayor+ind_cobert_serv_bajo+log(Comorbilidad)+Proporcion_trabajando%in%Estrato_moda,
                     family = poisson),
                 vcov.=vcovHC(glm(formula = datos_cali$Covid~Adulto_mayor+ind_cobert_serv_bajo+log(Comorbilidad)+Proporcion_trabajando%in%Estrato_moda,
                                  family = poisson), type="const"))#[22:28]

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
variables_continuas <- colnames(data_cali[,c(4:225,228)]) #131

plots <- purrr::map(
  .x = variables_continuas,
  .f = custom_corr_plot,
  variable2 = "Covid",
  df = data_cali[,c(4:225,228)]
)

plots

#ggpubr::ggarrange(plotlist = plots, ncol = 3, nrow = 2) %>%
#  ggpubr::annotate_figure(
#    top = ggpubr::text_grob("Comportamiento variables con Covid", 
#                            face = "bold", size = 16, x = 0.20)
#  )
#

