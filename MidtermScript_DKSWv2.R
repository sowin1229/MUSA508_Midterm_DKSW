#########
#SET UP
#########

#Load map features
library(tidyverse)
library(sf)
library(spdep)
library(caret)
library(ckanr)
library(FNN)
library(grid)
library(gridExtra)
library(ggcorrplot)
library(dplyr)
library(mapview)
library(jtools)     # for regression model plots

#root.dir = "https://raw.githubusercontent.com/urbanSpatial/Public-Policy-Analytics-Landing/master/DATA/"

mapTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle=element_text(face="italic"),
    plot.caption=element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),axis.title = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )
}

plotTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle = element_text(face="italic"),
    plot.caption = element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_line("grey80", size = 0.1),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    strip.background = element_rect(fill = "grey80", color = "white"),
    strip.text = element_text(size=12),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(colour = "black", face = "italic"),
    legend.text = element_text(colour = "black", face = "italic"),
    strip.text.x = element_text(size = 14)
  )
}

palette5 <- c("#25CB10", "#5AB60C", "#8FA108",   "#C48C04", "#FA7800")
palette4blueyellow <- c("#F2E85C", "#F2CA52", "5E7EBF", "4F6573", "#3552F2")

#Set up quintile breaks for map legends
qBr <- function(df, variable, rnd) {
  if (missing(rnd)) {
    as.character(quantile(round(df[[variable]],0),
                          c(.01,.2,.4,.6,.8), na.rm=T))
  } else if (rnd == FALSE | rnd == F) {
    as.character(formatC(quantile(df[[variable]]), digits = 3),
                 c(.01,.2,.4,.6,.8), na.rm=T)
  }
}

q5 <- function(variable) {as.factor(ntile(variable, 5))}

nn_function <- function(measureFrom,measureTo,k) {
  measureFrom_Matrix <- as.matrix(measureFrom)
  measureTo_Matrix <- as.matrix(measureTo)
  nn <-   
    get.knnx(measureTo, measureFrom, k)$nn.dist
  output <-
    as.data.frame(nn) %>%
    rownames_to_column(var = "thisPoint") %>%
    gather(points, point_distance, V1:ncol(.)) %>%
    arrange(as.numeric(thisPoint)) %>%
    group_by(thisPoint) %>%
    summarize(pointDistance = mean(point_distance)) %>%
    arrange(as.numeric(thisPoint)) %>% 
    dplyr::select(-thisPoint) %>%
    pull()
  
  return(output)  
}


#################
#DATA WRANGLING
#################


#read data for zip codes and isolate zipcodes in correct area
zipcodes <- st_read("https://opendata.arcgis.com/datasets/fee863cb3da0417fa8b5aaf6b671f8a7_0.geojson") %>%   
  st_transform('ESRI:102658') %>%  select(OBJECTID, PZIPCODEID, ZIP, ZIPCODE, SHAPE_Length, SHAPE_Area, geometry)%>% 
  filter(ZIPCODE %in% c("33134","33146", "33139", "33140", "33141", "33154", "33138", "33137", "33150","33151", "33153", "33127", "33136", "33128", "33130", "33129", "33133", "33145", "33135", "33125", "33142", "33147", "33126", "33144", "33155", "33143"))
#create shapefile of zipcodes
zipcode_sf <- zipcodes %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
                                      st_transform('ESRI:102658')

#Load housing data 
miamidata <- read_sf("~/GitHub/MUSA508_Midterm_DKSW/studentsData (1).geojson")%>%
  mutate(Price_PerSqFt = Assessed/ActualSqFt)

#create shapefile using Miami data that was just created from the geojson 
miamisf <- miamidata %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')

#plot zipcodes & housing data
ggplot() +
  geom_sf(data = zipcode_sf, fill = "white") +
  geom_sf(data = miamisf, aes(colour = q5(Price_PerSqFt)), 
          show.legend = "point", size = .75) +
  scale_colour_manual(values = palette5,
                      labels=qBr(miamidata,"Price_PerSqFt"),
                      name="Quintile\nBreaks") +
  labs(title="Price/SqFt in Miami & Miami Beach") +
  mapTheme()


#water
water <- read_sf("Water.geojson") %>%
  st_transform('ESRI:102658')
mapview(water, zcol='TYPE')

##Filter out water that is not part of ocean bay
waterbay<- water %>% 
  filter(TYPE == "B")

waterbay <- waterbay %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')

mapview(waterbay)


#Nearest "neighbor" homes to water distance
st_c <- st_coordinates

miamisf <- 
  miamisf %>% 
  mutate(
    bay_nn1 = nn_function(
  st_coordinates(st_centroid(miamisf)), 
  st_coordinates(st_centroid(waterbay)),
  1))


#Trader Joes and Whole Foods
TJWF <- read.csv("TJWF.csv")


TJWF <- TJWF %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658') 

mapview(list(TJWF, miamisf))


#Nearest "neighbor" homes to TJs and WFs distance
## Nearest Neighbor Feature

miamisf <- 
  miamisf %>% 
  mutate(
    tjwf_nn1 = nn_function(
      st_coordinates(st_centroid(miamisf)), 
      st_coordinates((TJWF)),
      1))

summary(miamisf)

#parks
parks <- read_sf("County_Park_Boundary.geojson") %>%
  st_transform('ESRI:102658')
mapview(parks)

#THIS ALSO DID NOT WORK
# Counts of parks per buffer of house sale 
 miamisf$park.buffer =
   st_buffer(miamisf, 660) %>% 
   aggregate(mutate(parks, counter = 1),., sum) %>%
#   pull(counter)



#NEAREST NEIGHBOR DID NOT WORK
# st_c <- st_coordinates
# 
# miamisf <- 
#   miamisf %>% 
#   mutate(
#     parks_nn1 = nn_function(
#       st_coordinates(st_centroid(miamisf)), 
#       st_coordinates(st_centroid(parks)),
#       1))



# beaches
# NOTE Beaches on hold for now.
# beaches <- st_read("https://opendata.arcgis.com/datasets/9e30807e3efd44f3b16ab8d3657249f2_0.geojson")%>%
#   st_transform('ESRI:102658') %>%
#   select(FID, NAME, ADDRESS, CITY, ZIPCODE, PHONE, OWNER, TOTALACRE,TYPE, WEBSITE, LAT, LON, POINT_X, POINT_Y, geometry)%>%
#   filter(CITY %in% c("Miami", "Miami Beach"))
# beachpoints <- beaches %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
#   st_transform('ESRI:ESRI:102658')
# 
# #CreateBuffer
# beachbuffer <-  st_buffer(beachpoints, 2640) %>%
#   mutate(Legend = "Buffer") %>%
#   dplyr::select(Legend) %>% 
#   st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
#   st_transform('ESRI:102658')
# 
# beachbufferunion <-st_union(st_buffer(beachbuffer, 2640)) %>%
#   st_sf() %>%
#   mutate(Legend = "Unioned Buffer")
# 
# ggplot() +
#   geom_sf(data=zipcode_sf)+
#   geom_sf(data=beachbufferunion, fill = "transparent") +
#   geom_sf(data=beachpoints, show.legend = "point") +
#   facet_wrap(~Legend) +
#   mapTheme()
# 
# ggplot() +
#   geom_sf(data = zipcode_sf, fill = "white") +
#   geom_sf(data=miamidata, aes(colour = q5(Price_PerSqFt)))+
#   geom_sf(data = beachpoints) +
#   geom_sf(data = beachbufferunion, fill = "transparent") +
#   geom_sf(data = miamisf, aes(colour = q5(Price_PerSqFt)), 
#           show.legend = "point", size = .75) +
#   scale_colour_manual(values = palette5,
#                       labels=qBr(miamidata,"Price_PerSqFt"),
#                       name="Quintile\nBreaks") +
#   labs(title="Price Per Sq Ft & Beaches") +
#   mapTheme()
# 
# 
# #FROM KENS BOOK: Mapping Variable Denisty ... beaches
# beachpoints2 <-
#   beaches %>%
#   dplyr::select(LAT, LON) %>%
#   na.omit() %>%
#   st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
#   st_transform('ESRI:102658') %>%
#   distinct()
# 
# ggplot() + geom_sf(data = zipcode_sf, fill = "white") +
#   stat_density2d(data = data.frame(st_coordinates(beachpoints2)), 
#                  aes(X, Y, fill = ..level.., alpha = ..level..),
#                  size = 0.01, bins = 40, geom = 'polygon') +
#   scale_fill_gradient(low = "#25CB10", high = "#FA7800", name = "Density") +
#   scale_alpha(range = c(0.00, 0.35), guide = FALSE) +
#   labs(title = "Density of Beaches") +
#   mapTheme()
# 
# #measure distance to beach - house to beach
#   st_c <- st_coordinates
#   miamisf <-
#     miamisf %>% 
#     mutate(
#       beach_nn1 = nn_function(st_c(miamisf), st_c(beachpoints2), 1),
#       beach_nn2 = nn_function(st_c(miamisf), st_c(beachpoints2), 2))
#       
# 
# ##miamisf$beachbuffer =
#   st_buffer(miamisf, 660) %>%
#   aggregate(mutate(beachpoints2, counter = 1),., sum) %>%
#   pull(counter)
# 
# nn_function <- function(measureFrom,measureTo,k) {
#   measureFrom_Matrix <- as.matrix(measureFrom)
#   measureTo_Matrix <- as.matrix(measureTo)
#   nn <-   
#     get.knnx(measureTo, measureFrom, k)$nn.dist
#   output <-
#     as.data.frame(nn) %>%
#     rownames_to_column(var = "thisPoint") %>%
#     gather(points, point_distance, V1:ncol(.)) %>%
#     arrange(as.numeric(thisPoint)) %>%
#     group_by(thisPoint) %>%
#     summarize(pointDistance = mean(point_distance)) %>%
#     arrange(as.numeric(thisPoint)) %>% 
#     dplyr::select(-thisPoint) %>%
#     pull()
#   
#   return(output)  
#}

#Great Schools ratings
allschools <- read_sf("Public_SChool.csv")
 
allschools    <- allschools %>% 
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658') 
mapview(allschools)

ggplot()+
  geom_sf(data = zipcode_sf, fill = "white") +
  geom_sf(data=miamisf, aes(colour = q5(ActualSqFt)),
          show.legend = "point", size = .75)+
  geom_sf(data=allschools, size = 1, shape = 21, fill = "darkred")+
  scale_colour_manual(values = palette5,
                      labels=qBr(miamidata,"ActualSqFt"),
                      name="Quintile\nBreaks") +
  labs(title="Price Per Sq Ft & Beaches") +
  mapTheme()

#make zipcodes into boundary
boundary <- st_union(zipcode_sf)

#Clipped only Miami Schools
mschools <- st_intersection(allschools, boundary)

mapview(list(mschools, miamisf))

ggplot()+
  geom_sf(data = zipcode_sf, fill = "white") +
  geom_sf(data=miamisf, aes(colour = q5(ActualSqFt)),
          show.legend = "point", size = .75)+
  geom_sf(data=mschools, size = 1, shape = 21, fill = "darkred")+
  scale_colour_manual(values = palette5,
                      labels=qBr(miamidata,"ActualSqFt"),
                      name="Quintile\nBreaks") +
  labs(title="Price Per Sq Ft & Beaches") +
  mapTheme()

#JOIN ZIP DATA AND STUDENT DATA
alldatamiami <- st_join(miamisf, zipcode_sf, left = TRUE)

ggplot()+
  geom_sf(data = zipcode_sf, fill = "white") +
  geom_sf(data=alldatamiami, aes(colour = q5(ActualSqFt)),
          show.legend = "point", size = .75)+
  scale_colour_manual(values = palette5,
                      labels=qBr(miamidata,"ActualSqFt"),
                      name="Quintile\nBreaks") +
  labs(title="Price Per Sq Ft & Beaches") +
  mapTheme()


#Join Census Tract and Student Data
tracts <- read_sf("Tract_Pop_2010.geojson") %>%
  st_transform('ESRI:102658')
mapview(tracts)

mtracts <- st_intersection(tracts, boundary)

alldatamiami <- st_join(alldatamiami, mtracts, left = TRUE)

ggplot()+
  geom_sf(data = mtracts, fill = "white") +
  geom_sf(data=alldatamiami, aes(colour = q5(ActualSqFt)),
          show.legend = "point", size = .75)+
  scale_colour_manual(values = palette5,
                      labels=qBr(miamidata,"ActualSqFt"),
                      name="Quintile\nBreaks") +
  labs(title="Price Per Sq Ft & Beaches") +
  mapTheme()

#Opportunity Insights Data 
oppinsights <- read_sf("OppInsights.csv")

#join oppinsights with alldatamiami
alldatamiami <-alldatamiami %>% left_join(oppinsights, by = c("GEOID10" = "GEOID"))

#Can I drop unnecessary variables?

