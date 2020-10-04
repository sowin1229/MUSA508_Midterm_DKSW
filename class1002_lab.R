library(tidyverse)
library(sf)
library(spdep)
library(caret)
library(ckanr)
library(FNN)
library(grid)
library(gridExtra)
library(ggcorrplot)
library(jtools)     # for regression model plots

# functions
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



data_src <- "https://raw.githubusercontent.com/urbanSpatial/Public-Policy-Analytics-Landing/master/DATA/Chapter3_4/"

nhoods <-
  st_read("http://bostonopendata-boston.opendata.arcgis.com/datasets/3525b0ee6e6b427f9aab5d0a1d0a1a28_0.geojson") %>%
  st_transform('ESRI:102286')

boston.sf   <- st_read(file.path(data_src,"boston_sf_Ch1_wrangled.geojson"))

glimpse(boston.sf)


reg1 <- lm(SalePrice ~ ., data = st_drop_geometry(boston.sf) %>% 
             dplyr::select(SalePrice, LivingArea, 
                           GROSS_AREA, R_TOTAL_RM, NUM_FLOORS,
                           R_BDRMS, R_FULL_BTH, R_HALF_BTH, 
                           R_KITCH, R_FPLACE))
summary(reg1)

paste0("Adjust R-Squared: ", round(summary(reg1)$adj.r.squared,2))
broom::glance(reg1)

# set random seed
set.seed(31357)

# get index for training sample
inTrain <- caret::createDataPartition(
  y = boston.sf$SalePrice, 
  p = .60, list = FALSE)
# split data into training and test
boston.training <- boston.sf[inTrain,] 
boston.test     <- boston.sf[-inTrain,]  

# Regression  
reg2 <- lm(SalePrice ~ ., data = st_drop_geometry(boston.training) %>% 
             dplyr::select(SalePrice, LivingArea, 
                           GROSS_AREA, R_TOTAL_RM, NUM_FLOORS,
                           R_BDRMS, R_FULL_BTH, R_HALF_BTH, 
                           R_KITCH, R_FPLACE))

# Run this a number of times to see Adjusted R2
summary(reg2)


reg2_predict <- predict(reg2, newdata = boston.test)
reg2_predict

## Mean Square Error train and test
rmse.train <- caret::MAE(predict(reg2), boston.training$SalePrice)
rmse.test  <- caret::MAE(reg2_predict, boston.test$SalePrice)

cat("Train MAE: ", as.integer(rmse.train), " \n","Test MAE: ", as.integer(rmse.test))

#plotting accuracy
preds.train <- data.frame(pred   = predict(reg2),
                          actual = boston.training$SalePrice,
                          source = "training data")
preds.test  <- data.frame(pred   = reg2_predict,
                          actual = boston.test$SalePrice,
                          source = "testing data")
preds <- rbind(preds.train, preds.test)

ggplot(preds, aes(x = pred, y = actual, color = source)) +
  geom_point() +
  geom_smooth(method = "lm", color = "green") +
  geom_abline(color = "orange") +
  coord_equal() +
  theme_bw() +
  facet_wrap(~source, ncol = 2) +
  labs(title = "Comparing predictions to actual values",
       x = "Predicted Value",
       y = "Actual Value") +
  theme(
    legend.position = "none"
  )

# use caret package cross-validation method
fitControl <- trainControl(method = "cv", 
                           number = 10,
                           # savePredictions differs from book
                           savePredictions = TRUE)
set.seed(717)
# crimes.buffer feature added
# for k-folds CV
reg.cv <- 
  train(SalePrice ~ ., data = st_drop_geometry(boston.sf) %>% 
          dplyr::select(SalePrice, LivingArea,  
                        GROSS_AREA, R_TOTAL_RM, NUM_FLOORS,
                        R_BDRMS, R_FULL_BTH, R_HALF_BTH, 
                        R_KITCH, R_FPLACE, crimes.Buffer), 
        method = "lm", 
        trControl = fitControl, 
        na.action = na.pass)

reg.cv

reg.cv$resample
RMSE  Rsquared      MAE Resample
1  675941.7 0.2945504 276416.7   Fold01
2  369841.9 0.3410855 250108.0   Fold02
3  368167.9 0.5866561 246425.4   Fold03
4  380621.9 0.3786402 240403.3   Fold04
5  305656.3 0.6074219 229799.8   Fold05
6  335957.9 0.3715664 241791.9   Fold06
7  627726.9 0.7405234 250344.2   Fold07
8  354013.4 0.6453627 239917.3   Fold08
9  404586.6 0.5956115 247275.3   Fold09
10 329898.6 0.5909977 230215.4   Fold10


reg.cv$resample %>% 
  pivot_longer(-Resample) %>% 
  mutate(name = as.factor(name)) %>% 
  ggplot(., aes(x = name, y = value, color = name)) +
  geom_jitter(width = 0.1) +
  facet_wrap(~name, ncol = 3, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "none"
  )

reg.cv$results
str(reg.cv$results)
> reg.cv$results
intercept     RMSE  Rsquared      MAE   RMSESD RsquaredSD    MAESD
1      TRUE 415241.3 0.5152416 245269.7 128249.2  0.1534117 13144.64
> str(reg.cv$results)
'data.frame':	1 obs. of  7 variables:
  $ intercept : logi TRUE
$ RMSE      : num 415241
$ Rsquared  : num 0.515
$ MAE       : num 245270
$ RMSESD    : num 128249
$ RsquaredSD: num 0.153
$ MAESD     : num 13145

names(reg.cv)


# extract predictions from CV object
cv_preds <- reg.cv$pred
# compare number of observations between data sets
nrow(boston.sf)


#spatial correlation
inTrain <- createDataPartition(
  y = paste(boston.sf$Name, boston.sf$NUM_FLOORS.cat, 
            boston.sf$Style, boston.sf$R_AC), 
  p = .60, list = FALSE)
boston.training <- boston.sf[inTrain,] 
boston.test <- boston.sf[-inTrain,]  

reg.training <- 
  lm(SalePrice ~ ., data = as.data.frame(boston.training) %>% 
       dplyr::select(SalePrice, LivingArea, Style, 
                     GROSS_AREA, NUM_FLOORS.cat,
                     R_BDRMS, R_FULL_BTH, R_HALF_BTH, 
                     R_KITCH, R_AC, R_FPLACE, crimes.Buffer))

boston.test <-
  boston.test %>%
  mutate(Regression = "Baseline Regression",
         SalePrice.Predict = predict(reg.training, boston.test),
         SalePrice.Error = SalePrice.Predict - SalePrice,
         SalePrice.AbsError = abs(SalePrice.Predict - SalePrice),
         SalePrice.APE = (abs(SalePrice.Predict - SalePrice)) / SalePrice.Predict)%>%
  filter(SalePrice < 5000000) 

k_nearest_neighbors = 5
#prices
coords <- st_coordinates(boston.sf) 
# k nearest neighbors
neighborList <- knn2nb(knearneigh(coords, k_nearest_neighbors))
spatialWeights <- nb2listw(neighborList, style="W")
boston.sf$lagPrice <- lag.listw(spatialWeights, boston.sf$SalePrice)

#errors
coords.test <-  st_coordinates(boston.test) 
neighborList.test <- knn2nb(knearneigh(coords.test, k_nearest_neighbors))
spatialWeights.test <- nb2listw(neighborList.test, style="W")
boston.test$lagPriceError <- lag.listw(spatialWeights.test, boston.test$SalePrice.AbsError)

ggplot(boston.sf, aes(x=lagPrice, y=SalePrice)) +
  geom_point(colour = "#FA7800") +
  geom_smooth(method = "lm", se = FALSE, colour = "#25CB10") +
  labs(title = "Price as a function of the spatial lag of price",
       caption = "Public Policy Analytics, Figure 6.6",
       x = "Spatial lag of price (Mean price of 5 nearest neighbors)",
       y = "Sale Price") +
  plotTheme()
#is homeprice spatially correlated? based on the plot, yes!

#morans I
moranTest <- moran.mc(boston.test$SalePrice.AbsError, 
                      spatialWeights.test, nsim = 999)

ggplot(as.data.frame(moranTest$res[c(1:999)]), aes(moranTest$res[c(1:999)])) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(aes(xintercept = moranTest$statistic), colour = "#FA7800",size=1) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title="Observed and permuted Moran's I",
       subtitle= "Observed Moran's I in orange",
       x="Moran's I",
       y="Count",
       caption="Public Policy Analytics, Figure 6.8") +
  plotTheme()