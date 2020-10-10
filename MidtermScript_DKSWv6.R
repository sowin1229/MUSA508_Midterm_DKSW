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
library(ggmap)

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

#Convert given miamisf data into points
miamipts <-   st_centroid(miamisf)
mapview(miamipts)

#Overwrite miamipts to our miamisf file
miamisf <- miamipts
rm(miamipts)

###############
#DATA CLEANING
###############


#Removal of unncessary miamisf variables
colnames(miamisf)
summary(miamisf)

miamisf$HEX <- NULL
miamisf$GPAR<- NULL
miamisf$Legal1 <- NULL
miamisf$Legal2 <- NULL
miamisf$Legal3 <- NULL
miamisf$Legal4 <- NULL
miamisf$Legal5 <- NULL
miamisf$Legal6 <- NULL
miamisf$County.Senior <- NULL
miamisf$County.Other.Exempt <- NULL
miamisf$County.2nd.HEX <- NULL
miamisf$County.LongTermSenior <- NULL
miamisf$City.2nd.HEX <- NULL
miamisf$City.LongTermSenior <- NULL
miamisf$City.Senior <- NULL
miamisf$City.Other.Exempt <- NULL
miamisf$Owner2 <- NULL


###############
#NEW FEATURES
##############

#Water
water <- read_sf("Water.geojson") %>%
  st_transform('ESRI:102658')
mapview(water, zcol='TYPE')

  ##Filter out water that is not part of ocean bay
  waterbay<- water %>% 
    filter(TYPE == "B")
  
  waterbay <- waterbay %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
    st_transform('ESRI:102658')
  
  mapview(waterbay)
  
  
  #Remove original water file from Global Env
  rm(water)
  
  #Nearest "neighbor" homes to water distance
  
  
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
  miamisf <- 
    miamisf %>% 
    mutate(
      tjwf_nn1 = nn_function(
        st_coordinates(st_centroid(miamisf)), 
        st_coordinates((TJWF)),
        1))



#parks
parks <- read_sf("County_Park_Boundary.geojson") %>%
  st_transform('ESRI:102658')
mapview(parks)
  
  #Nearest Neighbor for parks
   miamisf <- 
     miamisf %>% 
     mutate(
       parks_nn1 = nn_function(
         st_coordinates(st_centroid(miamisf)), 
         st_coordinates(st_centroid(parks)),
         1),
      parks_nn2 = nn_function(
       st_coordinates(st_centroid(miamisf)), 
      st_coordinates(st_centroid(parks)),
      2))
   
   
 # beaches
 beaches <- st_read("https://opendata.arcgis.com/datasets/9e30807e3efd44f3b16ab8d3657249f2_0.geojson") %>% 
   st_transform('ESRI:102658') %>%
   select(FID, NAME, ADDRESS, CITY, ZIPCODE, PHONE, OWNER, TOTALACRE,TYPE, WEBSITE, LAT, LON, POINT_X, POINT_Y, geometry)%>%
  filter(CITY %in% c("Miami", "Miami Beach"))
 
   
   miamisf <- 
     miamisf %>% 
     mutate(
       beach_nn1 = nn_function(
         st_coordinates(st_centroid(miamisf)), 
         st_coordinates((beaches)),
         1))
 
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
    scale_colour_manual(values = palette4blueyellow,
                        labels=qBr(miamidata,"ActualSqFt"),
                        name="Quintile\nBreaks") +
    labs(title="Price Per Sq Ft & Beaches") +
    mapTheme()

  #Making GSRATINGS Numeric
  mschools <- mschools %>% mutate(GSRATING = as.numeric(GSRATING))


  #Developing High and Low  Rating School Set
  mschoolsHI <- filter (mschools, mschools $ GSRATING > 6)
  mschoolsLO <- filter (mschools, mschools $ GSRATING < 5)


  #Nearest Neighbor Schools HIGH
  miamisf <- 
    miamisf %>% 
    mutate(
      mschoolsHI_nn1 = nn_function(
        st_coordinates(st_centroid(miamisf)), 
        st_coordinates((mschoolsHI)),
        1))

  miamisf <- 
    miamisf %>% 
    mutate(
      mschoolsHI_nn3 = nn_function(
        st_coordinates(st_centroid(miamisf)), 
        st_coordinates((mschoolsHI)),
        3))

  #Nearest Neighbor Schools LOW
  
  miamisf <- 
    miamisf %>% 
    mutate(
      mschoolsLO_nn1 = nn_function(
        st_coordinates(st_centroid(miamisf)), 
        st_coordinates((mschoolsLO)),
        1))
  
  miamisf <- 
    miamisf %>% 
    mutate(
      mschoolsLO_nn3 = nn_function(
        st_coordinates(st_centroid(miamisf)), 
        st_coordinates((mschoolsLO)),
        3))


#Transit stops
Transit <- st_read("https://opendata.arcgis.com/datasets/ee3e2c45427e4c85b751d8ad57dd7b16_0.geojson")

Transitstops <- Transit %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')

transitmiami <- st_intersection(Transitstops, boundary)%>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')

ggplot()+
  geom_sf(data = zipcode_sf, fill = "white") +
  geom_sf(data = waterbay, fill = "lightblue", xlim = c(-80, -82), ylim = c(30, 35), boundary = TRUE,)
  geom_sf(data=miamisf, aes(colour = q5(ActualSqFt)),
          show.legend = "point", size = .75)+
  geom_sf(data=transitmiami, size = 2, shape = 21, fill = "blue")+
  scale_colour_manual(values = palette5,
                      labels=qBr(miamidata,"ActualSqFt"),
                      name="Quintile\nBreaks") +
  labs(title="Price Per Sq Ft & Transit") +
  mapTheme()


mapview(list(transitmiami, miamisf))

  # Nearest neighbor transitstops
  miamisf <- 
    miamisf %>% 
    mutate(
      transitstops_nn1 = nn_function(
        st_coordinates(st_centroid(miamisf)), 
        st_coordinates((transitmiami)),
        1))


################################
#JOIN ZIP DATA AND STUDENT DATA
###############################


#Joining miami data and zipcodes
alldatamiami <- st_join(miamisf, zipcode_sf, left = TRUE)


#Join Census Tract and Student Data
tracts <- st_read("Tract_Pop_2010.geojson")%>%
  st_transform('ESRI:102658') %>% 
  mutate(tracts, pct_white = (WHITENH/POP2010*100)) %>% 
  mutate(tracts, pct_black = (BLACKNH/POP2010*100)) %>% 
  mutate(tracts, pct_hisp = (HISPAN/POP2010*100))

mapview(tracts)


mtracts <- st_intersection(tracts, boundary)

#Adding in income and renaming fields
income <- st_read("Income.csv") %>% 
rename(HHless10000=8, 
       HHmore200000=9, 
       MedHHinc=10, 
       AvgHHinc=11, 
       MedFinc=12, 
       AvgFinc=13,
       MedNFinc=14,
       AVgNFinc=15,
       MedHHincOwnOcc=16,
       MedHHincRentOcc=17,
       HHPubAs=18,
       HHNoPubAs=19,
       PerCapInc=20)

#fixing character columns to numeric
income <- income %>% mutate(HHless10000 = as.numeric(HHless10000))
income <- income %>% mutate(HHmore200000 = as.numeric(HHmore200000))
income <- income %>% mutate(MedHHinc = as.numeric(MedHHinc))
income <- income %>% mutate(AvgHHinc = as.numeric(AvgHHinc))
income <- income %>% mutate(MedFinc = as.numeric(MedFinc))
income <- income %>% mutate(AvgFinc = as.numeric(AvgFinc))
income <- income %>% mutate(MedHHincOwnOcc = as.numeric(MedHHincOwnOcc))
income <- income %>% mutate(MedHHincRentOcc = as.numeric(MedHHincRentOcc))
income <- income %>% mutate(HHPubAs = as.numeric(HHPubAs))
income <- income %>% mutate(HHNoPubAs = as.numeric(HHNoPubAs))
income <- income %>% mutate(PerCapInc = as.numeric(PerCapInc))

glimpse(income)


#joining Income to tracts
mtracts <-mtracts %>% left_join(income, by = c("GEOID10" = "FIPS"))


#Opportunity Insights Data 
oppinsights <- st_read("OppInsights.csv")
oppinsights <- oppinsights %>% mutate(HHIncomeOpp = as.numeric(HHIncomeOpp))

#join oppinsights with alldatamiami
mtracts <-mtracts %>% left_join(oppinsights, by = c("GEOID10" = "GEOID"))

#joining miami tracts to all
alldatamiami <- st_join(alldatamiami, mtracts, left = TRUE)

ggplot()+
  geom_sf(data = mtracts, fill = "white")+
  geom_sf(data=alldatamiami, aes(colour = q5(WHITENH)),
          show.legend = "point", size = .55)+
  scale_colour_manual(values = palette5,
                      labels=qBr(alldatamiami,"ActualSqFt"),
                      name="Quintile\nBreaks") +
  labs(title="Price Per Sq Ft & Beaches") +
  mapTheme()




#Creating exploratory data frame

miamiexplor <- select (alldatamiami, Folio, SalePrice, Assessed, AdjustedSqFt, 
                       LotSize, Bed, Bath, Stories, Units, YearBuilt, LivingSqFt,
                       ActualSqFt, toPredict, Price_PerSqFt, bay_nn1, tjwf_nn1, 
                       parks_nn1, parks_nn2, beach_nn1,mschoolsHI_nn1, mschoolsHI_nn3,
                       mschoolsLO_nn1, transitstops_nn1, mschoolsLO_nn3, ZIPCODE, 
                       GEOID10, POP2010, HU2010, HISPAN, WHITENH, BLACKNH, HHless10000,
                       HHmore200000, HHIncomeOpp, 88, 89, 90, 91, 92, 93, 94, 95,
                       96, 97, 98)

#############
#Correlations
#############

colnames(miamiexplor)
## Home Features cor
#ROUND 1 Features
st_drop_geometry(miamiexplor) %>% 
  dplyr::select(SalePrice, LivingSqFt, Price_PerSqFt, bay_nn1) %>%
  filter(SalePrice <= 1000000) %>%
  gather(Variable, Value, -SalePrice) %>% 
  ggplot(aes(Value, SalePrice)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se=F, colour = "#FA7800") +
  facet_wrap(~Variable, ncol = 3, scales = "free") +
  labs(title = "Price as a function of continuous variables") +
  plotTheme()
#Round 2 Features
st_drop_geometry(miamiexplor) %>% 
  dplyr::select(SalePrice, tjwf_nn1, parks_nn1, parks_nn2) %>%
  filter(SalePrice <= 1000000) %>%
  gather(Variable, Value, -SalePrice) %>% 
  ggplot(aes(Value, SalePrice)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se=F, colour = "#FA7800") +
  facet_wrap(~Variable, ncol = 3, scales = "free") +
  labs(title = "Price as a function of continuous variables") +
  plotTheme()

#Round 3 Features
st_drop_geometry(miamiexplor) %>% 
  dplyr::select(SalePrice, beach_nn1, transitstops_nn1, HISPAN) %>%
  filter(SalePrice <= 1000000) %>%
  gather(Variable, Value, -SalePrice) %>% 
  ggplot(aes(Value, SalePrice)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se=F, colour = "#FA7800") +
  facet_wrap(~Variable, ncol = 3, scales = "free") +
  labs(title = "Price as a function of continuous variables") +
  plotTheme()

#Round 4 Features
st_drop_geometry(miamiexplor) %>% 
  dplyr::select(SalePrice, WHITENH, BLACKNH, HHless10000) %>%
  filter(SalePrice <= 1000000) %>%
  gather(Variable, Value, -SalePrice) %>% 
  ggplot(aes(Value, SalePrice)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se=F, colour = "#FA7800") +
  facet_wrap(~Variable, ncol = 3, scales = "free") +
  labs(title = "Price as a function of continuous variables") +
  plotTheme

#Round 5 Features
st_drop_geometry(miamiexplor) %>% 
  dplyr::select(SalePrice, HHmore200000, HHIncomeOpp, MedHHinc) %>%
  filter(SalePrice <= 1000000) %>%
  gather(Variable, Value, -SalePrice) %>% 
  ggplot(aes(Value, SalePrice)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se=F, colour = "#FA7800") +
  facet_wrap(~Variable, ncol = 3, scales = "free") +
  labs(title = "Price as a function of continuous variables") +
  plotTheme()

#Round 6 Features
st_drop_geometry(miamiexplor) %>% 
  dplyr::select(SalePrice, AvgHHinc, MedHHincOwnOcc, MedHHincRentOcc) %>%
  filter(SalePrice <= 1000000) %>%
  gather(Variable, Value, -SalePrice) %>% 
  ggplot(aes(Value, SalePrice)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se=F, colour = "#FA7800") +
  facet_wrap(~Variable, ncol = 3, scales = "free") +
  labs(title = "Price as a function of continuous variables") +
  plotTheme()

#Round 7 Features
st_drop_geometry(miamiexplor) %>% 
  dplyr::select(SalePrice, mschoolsHI_nn1, mschoolsLO_nn3, mschoolsLO_nn1) %>%
  filter(SalePrice <= 1000000) %>%
  gather(Variable, Value, -SalePrice) %>% 
  ggplot(aes(Value, SalePrice)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se=F, colour = "#FA7800") +
  facet_wrap(~Variable, ncol = 3, scales = "free") +
  labs(title = "Price as a function of continuous variables") +
  plotTheme()

## Waterbay cor
st_drop_geometry(miamiexplor) %>% 
  dplyr::select(SalePrice, starts_with("park")) %>%
  filter(SalePrice <= 1000000) %>%
  gather(Variable, Value, -SalePrice) %>% 
  ggplot(aes(Value, SalePrice)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se=F, colour = "#FA7800") +
  facet_wrap(~Variable, nrow = 1, scales = "free") +
  labs(title = "Price as a function of continuous variables") +
  plotTheme()

### Corr matrix
numericVars <- 
  select_if(st_drop_geometry(miamiexplor), is.numeric) %>% na.omit()

ggcorrplot(
  round(cor(numericVars), 1), 
  p.mat = cor_pmat(numericVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across numeric variables") 


## Univarite correlation vs. regression
cor.test(miamiexplor$bay_nn1, alldatamiami$SalePrice, method = "pearson")

ggplot(filter(miamiexplor, SalePrice <= 2000000), aes(y=SalePrice, x = bay_nn1)) +
  geom_point() +
  geom_smooth(method = "lm")

## Univarite Regression
parksReg <- lm(SalePrice ~ MedHHinc, data = miamiexplor)

summary(MedHHincReg)
summ(MedHHincReg)



## Multivariate Regression
reg1 <- lm(SalePrice ~ ., data = st_drop_geometry(miamiexplor) %>% 
             dplyr::select(SalePrice, tjwf_nn1, 
                           Bed, parks_nn2, HHIncomeOpp,
                          mschoolsLO_nn3, ActualSqFt, MedHHinc,
                          Land, Bath, Stories, LivingSqFt, bay_nn1,
                          mschoolsHI_nn3, transitstops_nn1,  
                          HHPubAs, HHNoPubAs)) 



# R^2 = 0.59
summ(reg1)

## Plot of marginal response
effect_plot(reg1, pred = HHIncomeOpp, interval = TRUE, plot.points = TRUE)

## Plot coefficients
plot_summs(reg1)

## plot multiple model coeffs
plot_summs(reg1, parksReg)