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

#Multiringbuffer
multipleRingBuffer <- function(inputPolygon, maxDistance, interval) 
{
  #create a list of distances that we'll iterate through to create each ring
  distances <- seq(0, maxDistance, interval)
  #we'll start with the second value in that list - the first is '0'
  distancesCounter <- 2
  #total number of rings we're going to create
  numberOfRings <- floor(maxDistance / interval)
  #a counter of number of rings
  numberOfRingsCounter <- 1
  #initialize an otuput data frame (that is not an sf)
  allRings <- data.frame()
  
  #while number of rings  counteris less than the specified nubmer of rings
  while (numberOfRingsCounter <= numberOfRings) 
  {
    #if we're interested in a negative buffer and this is the first buffer
    #(ie. not distance = '0' in the distances list)
    if(distances[distancesCounter] < 0 & distancesCounter == 2)
    {
      #buffer the input by the first distance
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #different that buffer from the input polygon to get the first ring
      buffer1_ <- st_difference(inputPolygon, buffer1)
      #cast this sf as a polygon geometry type
      thisRing <- st_cast(buffer1_, "POLYGON")
      #take the last column which is 'geometry'
      thisRing <- as.data.frame(thisRing[,ncol(thisRing)])
      #add a new field, 'distance' so we know how far the distance is for a give ring
      thisRing$distance <- distances[distancesCounter]
    }
    
    
    #otherwise, if this is the second or more ring (and a negative buffer)
    else if(distances[distancesCounter] < 0 & distancesCounter > 2) 
    {
      #buffer by a specific distance
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #create the next smallest buffer
      buffer2 <- st_buffer(inputPolygon, distances[distancesCounter-1])
      #This can then be used to difference out a buffer running from 660 to 1320
      #This works because differencing 1320ft by 660ft = a buffer between 660 & 1320.
      #bc the area after 660ft in buffer2 = NA.
      thisRing <- st_difference(buffer2,buffer1)
      #cast as apolygon
      thisRing <- st_cast(thisRing, "POLYGON")
      #get the last field
      thisRing <- as.data.frame(thisRing$geometry)
      #create the distance field
      thisRing$distance <- distances[distancesCounter]
    }
    
    #Otherwise, if its a positive buffer
    else 
    {
      #Create a positive buffer
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #create a positive buffer that is one distance smaller. So if its the first buffer
      #distance, buffer1_ will = 0. 
      buffer1_ <- st_buffer(inputPolygon, distances[distancesCounter-1])
      #difference the two buffers
      thisRing <- st_difference(buffer1,buffer1_)
      #cast as a polygon
      thisRing <- st_cast(thisRing, "POLYGON")
      #geometry column as a data frame
      thisRing <- as.data.frame(thisRing[,ncol(thisRing)])
      #add teh distance
      thisRing$distance <- distances[distancesCounter]
    }  
    
    #rbind this ring to the rest of the rings
    allRings <- rbind(allRings, thisRing)
    #iterate the distance counter
    distancesCounter <- distancesCounter + 1
    #iterate the number of rings counter
    numberOfRingsCounter <- numberOfRingsCounter + 1
  }
  
  #convert the allRings data frame to an sf data frame
  allRings <- st_as_sf(allRings)
}


#read data for zip codes and isolate zipcodes in correct area
zipcodes <- st_read("https://opendata.arcgis.com/datasets/fee863cb3da0417fa8b5aaf6b671f8a7_0.geojson") %>%   
  st_transform('ESRI:102286') %>%  select(OBJECTID, PZIPCODEID, ZIP, ZIPCODE, SHAPE_Length, SHAPE_Area, geometry)%>% 
  filter(ZIPCODE %in% c("33134","33146", "33139", "33140", "33141", "33154", "33138", "33137", "33150","33151", "33153", "33127", "33136", "33128", "33130", "33129", "33133", "33145", "33135", "33125", "33142", "33147", "33126", "33144", "33155", "33143"))
#create shapefile of zipcodes
zipcode_sf <- zipcodes %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
                                      st_transform('ESRI:102286')

#Load housing data 
miamidata <- read_sf("~/GitHub/MUSA508_Midterm_DKSW/studentsData (1).geojson")%>%
  mutate(Price_PerSqFt = Assessed/ActualSqFt)
#Pricepersquarefoot



#create shapefile using Miami data that was just created from the geojson 
miamisf <- miamidata %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102286')

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


#other data: beaches
beaches <- st_read("https://opendata.arcgis.com/datasets/9e30807e3efd44f3b16ab8d3657249f2_0.geojson")%>%
  st_transform('ESRI:102286') %>%
  select(FID, NAME, ADDRESS, CITY, ZIPCODE, PHONE, OWNER, TOTALACRE,TYPE, WEBSITE, LAT, LON, POINT_X, POINT_Y, geometry)%>%
  filter(CITY %in% c("Miami", "Miami Beach"))
beachpoints <- beaches %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102286')

beachbuffer <-  st_buffer(beachpoints, 5280) %>%
  mutate(Legend = "Buffer") %>%
  dplyr::select(Legend)

ggplot() +
  geom_sf(data=beachbuffer) +
  geom_sf(data=beachpoints, show.legend = "point") +
  facet_wrap(~Legend) +
  mapTheme()

ggplot() +
  geom_sf(data = zipcode_sf, fill = "white") +
  geom_sf(data = beachpoints, fill = "blue") +
  geom_sf(data = miamisf, aes(colour = q5(Price_PerSqFt)), 
          show.legend = "point", size = .75) +
  scale_colour_manual(values = palette5,
                      labels=qBr(miamidata,"Price_PerSqFt"),
                      name="Quintile\nBreaks") +
  labs(title="Price Per Sq Ft & Beaches") +
  mapTheme()


#JOIN ZIP DATA AND STUDENT DATA

