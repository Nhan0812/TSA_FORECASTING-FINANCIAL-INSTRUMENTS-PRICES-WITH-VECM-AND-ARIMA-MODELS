#1. Load Librabry

library(xts)
library(tidyverse)
library(quantmod)
library(lmtest)
library(dygraphs)


getwd()
getwd()
source("testdf.R")
Data <- read.csv("TSA_2023_project_data_1.csv",header = TRUE, dec = ".")

#2. Summary Data


Data %>% glimpse()
structure(Data)
glimpse(Data)

#Notice that the first column is the Date column, which is under "Character" type
#We need to transform it into "Date" type.

Data$X <- as.Date(Data$X, 
                  format = "%Y-%m-%d") 
glimpse(Data)

#Until now the class is "Data.Frame" object
class(Data) 

#Create xts objects
Data <- xts(Data[,-1], Data$X, order.by = Data$X)

#After creating xts objects, the class is "xts", "zoo"
class(Data) 

#Checking for missing data
summary(Data)
#There is no missing data for the table.

#3. Checking for Cointegration
#a/ Visualize the chart for all 10 time series.

#Plot 10 Time Series together
#First Method
plot(Data,  
     main = "Line Chart for 10 Time Series",  
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     legend.loc = "bottomleft",
     type="l")

#Second Method
Data %>%
  dygraph() %>%
  dyRangeSelector(height = 40)


#Plot 10 Time Series separately

plot(Data,
     main = "Line Chart for 10 Time Series",  
     major.ticks = "years", 
     grid.ticks.on = "years",
     legend.loc = "bottomleft",
     multi.panel = TRUE,
     yaxis.same = FALSE,
     type="l")

plot(Data,
     main = list("Line Chart for 10 Time Series", 
                 cex = 1.5, 
                 line = 2,
                 col = "darkblue",
                 font = 2),
     major.ticks = "years", 
     grid.ticks.on = "years",
     legend.loc = "bottomleft",
     multi.panel = TRUE,
     yaxis.same = FALSE,
     type = "l")

plot(Data$x2,type="l")
plot(Data$x3,type="l")
plot(Data$x4,type="l")
plot(Data$x5,type="l")
plot(Data$x6,type="l")
plot(Data$x7,type="l")
plot(Data$x8,type="l")
plot(Data$x9,type="l")
plot(Data$x10,type="l")


#Time series 1 & 2:
#4. Create first difference
Data$dx1 <- diff.xts(Data$x1)
Data$dx2 <- diff.xts(Data$x1)
Data$dx3 <- diff.xts(Data$x3)
Data$dx4 <- diff.xts(Data$x4)
Data$dx5 <- diff.xts(Data$x5)
Data$dx6 <- diff.xts(Data$x6)
Data$dx7 <- diff.xts(Data$x7)
Data$dx8 <- diff.xts(Data$x8) 
Data$dx9 <- diff.xts(Data$x9)
Data$dx10 <- diff.xts(Data$x10) 

#Plot both variables on the graph:
plot(Data[, 1:2],
     col = c("red", "blue"),
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 2,
     main = "Time Series X1 & X2",
     legend.loc = "topright")

#2 Testing cointegration
#We perform the tests of integration order.

#Testing the order of the time series 1 and its difference:

testdf(variable = Data$x1,
       max.augmentations = 3)

#Interprete result: The 
#Beusch-Godfrey: Ho: There is no auto-correlation in Residuals

testdf(variable = Data$dx1, 
       max.augmentations = 3)

#Since we can reject the null in the case of the first differences, 
#we can conclude that the Time Series X1 is integrated of order 1.

#Testing the order of the time series 2 and its difference:

testdf(variable = Data$x2,
       max.augmentations = 3)

testdf(variable = Data$dx2, 
       max.augmentations = 3)


#Since we can reject the null in the case of the first differences, 
#we can conclude that the Time Series X2 is integrated of order 1.


#As a result, both variables has the same order 1: ∼I(1)
#so in the next step we can check whether they are cointegrated or not.

#3. Estimating the Cointegrated Vector:

model.coint <- lm(x1 ~ x2, data = Data)

#Examine the model summary:
summary(model.coint)

#The model is significantly explained by x2.
#We further  test the stationary of the residuals:
testdf(variable = residuals(model.coint), max.augmentations = 3)

#The ADF test with no augmentations can be used. 
#The result is that non-stationarity of residuals is STRONGLY REJECTED, 
#so residuals are stationary, which means that x1 and x2 are cointegrated.

#The cointegrating vector is [1, -29.890 , -0.713]

#which defines the cointegrating relationship as: 1 * x1 - 29.890 - 0.713 * x2.