#1. Load Librabry

library(xts)
library(tidyverse)
library(quantmod)
library(lmtest)
library(forecast)
library(dygraphs)


getwd()
setwd("D://UW//2. Summer 22-23//1. Time Series Analysis//Project//TSA")


source("testdf.R")
source("function_plot_ACF_PACF_resids.R")

Data <- read.csv("TSA_2023_project_data_1.csv",header = TRUE, dec = ".")

#2. Summary Data


Data %>% glimpse()
str(Data)
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

# 4. Applying Box-Jenkins procedure for Time Series X1:

# Step 1: Model Parameters:  
par(mfrow = c(2, 1)) 
acf(Data$dx1,
    lag.max = 36, # max lag for ACF
    ylim = c(-0.1, 0.1),   # limits for the y axis - we give c(min, max)
    lwd = 5,               
    col = "dark green",
    na.action = na.pass)   # do not stop if there are missing values in the data
pacf(Data$dx1, 
     lag.max = 36, 
     lwd = 5, col = "dark green",
     na.action = na.pass)

par(mfrow = c(1, 1)) # restore the original single panel

#The PACF shown is suggestive of an AR(5) model. 
#So an initial candidate model is an ARIMA(5,1,0). 
#We also have some variations of this model: ARIMA(5,1,1), ARIMA(4,1,0),ARIMA(3,1,0).


#Step 2 - model estimation

arima510 <- Arima(Data$x1,  # variable	
                  order = c(5, 1, 0)  # (p,d,q) parameters	
)	
coeftest(arima510)
summary(arima510)

#The above model is not included constant. Let's include constant into the model: 
arima510_2 <- Arima(Data$x1,  # variable	
                  order = c(5, 1, 0),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant
)
coeftest(arima510_2)
#Adding the constant did not change the result much, so we keep the model without constant

#Step 3 - model diagnostics
#Method 1:
plot(resid(arima510))


#Method 2:
tibble(
  date = index(Data),
  resid = arima510 %>% resid() %>% as.numeric()
) %>%
  ggplot(aes(date, resid)) +
  geom_line(col = "royalblue3") +
  theme_bw()

#Lets check the ACF and the PACF of the Residual values:

#Method 1:
par(mfrow = c(2, 1)) 
acf(resid(arima510), 
    lag.max = 36,
    ylim = c(-0.1, 0.1), 
    lwd = 5, col = "dark green",
    na.action = na.pass)
pacf(resid(arima510), 
     lag.max = 36, 
     lwd = 5, col = "dark green",
     na.action = na.pass)
par(mfrow = c(1, 1))

#Method 2:
plot_ACF_PACF_resids(arima510)

# The Ljung-Box test (for a maximum of 10 lags):

Box.test(resid(arima510), type = "Ljung-Box", lag = 10)
Box.test(resid(arima510), type = "Ljung-Box", lag = 15)
Box.test(resid(arima510), type = "Ljung-Box", lag = 20)
Box.test(resid(arima510), type = "Ljung-Box", lag = 25)

#-> Very large p-values: The residual is white-noise  

#Plot graph for Ljung-Box test:
bj_pvalues = c()

for(i in c(1:100)){
  bj = Box.test(resid(arima510), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}

plot(bj_pvalues, type='l')

abline(h=0.05, col='red', type='l')


#Model ARIMA(5,1,1)

arima511 <- Arima(Data$x1,  # variable	
                    order = c(5, 1, 1),  # (p,d,q) parameters
                    include.constant = TRUE)  # including a constant
)
coeftest(arima511)

#ARIMA(4,1,0)
arima410 <- Arima(Data$x1,  # variable	
                  order = c(4, 1, 0),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant
)
coeftest(arima410)

#ARIMA(3,1,0).
arima310 <- Arima(Data$x1,  # variable	
                  order = c(3, 1, 0),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant
)
coeftest(arima310)

#Step 4. Evaluate Model:

AIC(arima510,arima510_2, arima511, arima410, arima310)
BIC(arima510, arima510_2, arima511, arima410, arima310)

#-> Suggestion ARIMA(5,1,0)

#Cross check with Auto Correlation:

arima.best.AIC <- 
  auto.arima(Data$x1,
             d = 1,             # parameter d of ARIMA model
             max.p = 6,         # Maximum value of p
             max.q = 6,         # Maximum value of q
             max.order = 12,    # maximum p+q
             start.p = 1,       # Starting value of p in stepwise procedure
             start.q = 1,       # Starting value of q in stepwise procedure
             ic = "aic",        # Information criterion to be used in model selection.
             stepwise = FALSE,  # if FALSE considers all models
             allowdrift = TRUE, # include a constant
             trace = TRUE)      # show summary of all models considered

#Return: Best model: ARIMA(3,1,4)
arima314 <- Arima(Data$x1,  # variable	
                  order = c(3, 1, 4),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant
)
coeftest(arima314)

AIC(arima510,arima314)

BIC(arima510,arima314)

arima.best.BIC <- 
  auto.arima(Data$x1,
             d = 1,             # parameter d of ARIMA model
             max.p = 6,         # Maximum value of p
             max.q = 6,         # Maximum value of q
             max.order = 12,    # maximum p+q
             start.p = 1,       # Starting value of p in stepwise procedure
             start.q = 1,       # Starting value of q in stepwise procedure
             ic = "bic",        # Information criterion to be used in model selection.
             stepwise = FALSE,  # if FALSE considers all models
             allowdrift = TRUE, # include a constant
             trace = TRUE)      # show summary of all models considered
#Return: Best model: ARIMA(3,1,2)
arima312 <- Arima(Data$x1,  # variable	
                  order = c(3, 1, 2),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant
)
coeftest(arima312)

AIC(arima510,arima314,arima312)

BIC(arima510,arima314,arima312)

#5. Applying Box-Jenkins procedure for Time Series X2:

# Step 1: Model Parameters:  
par(mfrow = c(2, 1)) 
acf(Data$dx2,
    lag.max = 36, # max lag for ACF
    ylim = c(-0.1, 0.1),   # limits for the y axis - we give c(min, max)
    lwd = 5,               
    col = "dark green",
    na.action = na.pass)   # do not stop if there are missing values in the data
pacf(Data$dx2, 
     lag.max = 36, 
     lwd = 5, col = "dark green",
     na.action = na.pass)

par(mfrow = c(1, 1)) # restore the original single panel

#The PACF shown is suggestive of an AR(5) model. 
#So an initial candidate model is an ARIMA(5,1,0). 
#We also have some variations of this model: ARIMA(5,1,1), ARIMA(4,1,0),ARIMA(3,1,0).

#Step 2 - model estimation

arima510_x2 <- Arima(Data$x2,  # variable	
                  order = c(5, 1, 0)  # (p,d,q) parameters
)	
coeftest(arima510_x2)
summary(arima510_x2)

#The above model is not included constant. Let's include constant into the model: 
arima510_2_x2 <- Arima(Data$x2,  # variable	
                    order = c(5, 1, 0),  # (p,d,q) parameters
                    include.constant = TRUE)  # including a constant
)
coeftest(arima510_2_x2)
#Adding the constant did not change the result much, so we keep the model without constant

#Step 3 - model diagnostics
#Method 1:
plot(resid(arima510_x2))

# The Ljung-Box test (for a maximum of 10 lags):

Box.test(resid(arima510_x2), type = "Ljung-Box", lag = 10)
Box.test(resid(arima510_x2), type = "Ljung-Box", lag = 15)
Box.test(resid(arima510_x2), type = "Ljung-Box", lag = 20)
Box.test(resid(arima510_x2), type = "Ljung-Box", lag = 25)

#-> Very large p-values: The residual is white-noise  

#Plot graph for Ljung-Box test:
bj_pvalues = c()

for(i in c(1:100)){
  bj = Box.test(resid(arima510_x2), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}

plot(bj_pvalues, type='l')

abline(h=0.05, col='red', type='l')


#Model ARIMA(5,1,1)

arima511_x2 <- Arima(Data$x2,  # variable	
                  order = c(5, 1, 1),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant

coeftest(arima511_x2)

#ARIMA(4,1,0)
arima410_x2 <- Arima(Data$x2,  # variable	
                  order = c(4, 1, 0),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant

coeftest(arima410_x2)

#ARIMA(3,1,0).
arima310_x2 <- Arima(Data$x2,  # variable	
                  order = c(3, 1, 0),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant

coeftest(arima310_x2)

#Step 4. Evaluate Model:

AIC(arima510_x2,arima510_2_x2, arima511_x2, arima410_x2, arima310_x2)
BIC(arima510_x2, arima510_2_x2, arima511_x2, arima410_x2, arima310_x2)

#-> Suggestion model: arima510_x2 - ARIMA(5,1,0)


#Cross check with Auto Correlation:

arima.best.AIC <- 
  auto.arima(Data$x2,
             d = 1,             # parameter d of ARIMA model
             max.p = 6,         # Maximum value of p
             max.q = 6,         # Maximum value of q
             max.order = 12,    # maximum p+q
             start.p = 1,       # Starting value of p in stepwise procedure
             start.q = 1,       # Starting value of q in stepwise procedure
             ic = "aic",        # Information criterion to be used in model selection.
             stepwise = FALSE,  # if FALSE considers all models
             allowdrift = TRUE, # include a constant
             trace = TRUE)      # show summary of all models considered

#Return: Best model: ARIMA(3,1,3)
arima313_x2 <- Arima(Data$x2,  # variable	
                  order = c(3, 1, 3),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant

AIC(arima510_x2,arima313_x2)

BIC(arima510_x2,arima313_x2)

arima.best.BIC <- 
  auto.arima(Data$x2,
             d = 1,             # parameter d of ARIMA model
             max.p = 6,         # Maximum value of p
             max.q = 6,         # Maximum value of q
             max.order = 12,    # maximum p+q
             start.p = 1,       # Starting value of p in stepwise procedure
             start.q = 1,       # Starting value of q in stepwise procedure
             ic = "bic",        # Information criterion to be used in model selection.
             stepwise = FALSE,  # if FALSE considers all models
             allowdrift = TRUE, # include a constant
             trace = TRUE)      # show summary of all models considered

#Return: Best model: ARIMA(3,1,2)
arima312_x2 <- Arima(Data$x2,  # variable	
                  order = c(3, 1, 2),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant
AIC(arima510_x2,arima313_x2, arima312_x2)

BIC(arima510_x2,arima313_x2,arima312_x2 )
