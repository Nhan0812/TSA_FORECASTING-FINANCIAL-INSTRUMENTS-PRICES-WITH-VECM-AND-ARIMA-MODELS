#1. Load Librabry

library(xts)
library(tidyverse)
library(quantmod)
library(lmtest)
library(forecast)
library(dygraphs)
library(vars)
library(knitr)
library(kableExtra)
library(ggplot2)


getwd()
setwd("D://UW//2. Summer 22-23//1. Time Series Analysis//Project//TSA")


source("testdf.R")
source("function_plot_ACF_PACF_resids.R")

Data <- read.csv("TSA_2023_project_data_1.csv",header = TRUE, dec = ".")

#2. Summary Data

str(Data)


#Notice that the first column is the Date column, which is under "Character" type
#We need to transform it into "Date" type.

Data$X <- as.Date(Data$X, 
                  format = "%d-%m-%y") 

#Until now the class is "Data.Frame" object
class(Data) 

#Create xts objects
Data <- xts(Data[,-1], Data$X)

#After creating xts objects, the class is "xts", "zoo"
class(Data) 

#Checking for missing data
summary(Data)
head(Data,6)
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

plot(Data$x1,type="l")
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

#The cointegrating vector is [1, -29.853 , -0.713]

#which defines the cointegrating relationship as: 1 * x1 - 29.853 - 0.713 * x2.

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

abline(h = 0.05, col='red')


#Model ARIMA(5,1,1)

arima511 <- Arima(Data$x1,  # variable	
                    order = c(5, 1, 1),  # (p,d,q) parameters
                    include.constant = TRUE)  # including a constant

coeftest(arima511)

#ARIMA(4,1,0)
arima410 <- Arima(Data$x1,  # variable	
                  order = c(4, 1, 0),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant

coeftest(arima410)

#ARIMA(3,1,0).
arima310 <- Arima(Data$x1,  # variable	
                  order = c(3, 1, 0),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant

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

#Return: Best model: ARIMA(4,1,2)
arima412 <- Arima(Data$x1,  # variable	
                  order = c(4, 1, 2),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant

coeftest(arima412)

AIC(arima510,arima412)

BIC(arima510,arima412)

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

coeftest(arima312)

AIC(arima510,arima412,arima312)

BIC(arima510,arima412,arima312)

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

abline(h=0.05, col='red')


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


#FORECASTING
#Import Data Forecast:

out_of_sample <-read.csv("Out_of_sample.csv",header = TRUE, dec = ".")

class(out_of_sample)

#Change Date
out_of_sample$X <- as.Date(out_of_sample$X, 
                           format = "%d-%m-%y") 

#Create xts objects
out_of_sample <- xts(out_of_sample[,-1], out_of_sample$X)

#Assign forecast X1 to object oos_x1
oos_x1 <- out_of_sample$x1

oos_x1

#FORECAST X1
tail(Data, 12)

forecasts_x1 <- forecast(arima510, # model for prediction
                      h = 30) # how many periods outside the sample

forecasts_x1

#Extract the forecast points
forecasts_x1$mean
class(forecasts_x1$mean)
#It is a ts object, not xts!
#However, the xts objects are more convenient and modern.
#In terms of plotting the real data and forecast data to compare.

forecasts_x1$lower
forecasts_x1$upper

#Create xts object: (We use the second column (95% confidence interval))
forecasts_x1_data <- data.frame(f_mean = as.numeric(forecasts_x1$mean),
                             f_lower = as.numeric(forecasts_x1$lower[, 2]),
                             f_upper = as.numeric(forecasts_x1$upper[, 2]))

head(forecasts_x1_data,30)


#Adding real value X1 below the current Dataset:
Data_x1 <- rbind(Data[, "x1"], oos_x1)
tail(Data_x1, n = 40)


#Turn Forecast with Index into Forecast with Date
forecasts_x1_xts <- xts(forecasts_x1_data,
                     order.by = index(oos_x1))
forecasts_x1_xts

#Merge it together with the original data
Data_x1_combined <- merge(Data_x1, forecasts_x1_xts)
head(Data_x1_combined)
tail(Data_x1_combined,40)


#Plot Chart

plot(Data_x1_combined [, c("x1", "f_mean", "f_lower", "f_upper")], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 day forecast of x1",
     col = c("black", "blue", "red", "red"))


plot(Data_x1_combined ["2020-11/", c("x1", "f_mean", "f_lower", "f_upper")], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 day forecast of x1",
     col = c("black", "blue", "red", "red"))

#Extract Data for Evaluating the Forecast:
Data_x1_Eva <- tail(Data_x1_combined, 30)  
Data_x1_Eva

Data_x1_Eva$mae   <-  abs(Data_x1_Eva$x1 - Data_x1_Eva$f_mean)
Data_x1_Eva$mse   <-  (Data_x1_Eva$x1 - Data_x1_Eva$f_mean) ^ 2
Data_x1_Eva$mape  <-  abs((Data_x1_Eva$x1 - Data_x1_Eva$f_mean)/Data_x1_Eva$x1)
Data_x1_Eva$amape <-  abs((Data_x1_Eva$x1 - Data_x1_Eva$f_mean)/(Data_x1_Eva$x1 + Data_x1_Eva$f_mean))
Data_x1_Eva

colMeans(Data_x1_Eva[, c("mae", "mse", "mape", "amape")])
apply(Data_x1_Eva[, c("mae", "mse", "mape", "amape")], 2, FUN = median)



##########################################################################################
#FORECAST X2



#Assign forecast X2 to object oos_x2
oos_x2 <- out_of_sample$x2

oos_x2

#FORECAST X2
tail(Data, 12)

forecasts_x2 <- forecast(arima510_x2, # model for prediction
                         h = 30) # how many periods outside the sample

forecasts_x2

#Extract the forecast points
forecasts_x2$mean
class(forecasts_x2$mean)
#It is a ts object, not xts!
#However, the xts objects are more convenient and modern.
#In terms of plotting the real data and forecast data to compare.

forecasts_x2$lower
forecasts_x2$upper

#Create xts object: (We use the second column (95% confidence interval))
forecasts_x2_data <- data.frame(f_mean = as.numeric(forecasts_x2$mean),
                                f_lower = as.numeric(forecasts_x2$lower[, 2]),
                                f_upper = as.numeric(forecasts_x2$upper[, 2]))

head(forecasts_x2_data,30)


#Adding real value X1 below the current Dataset:
Data_x2 <- rbind(Data[, "x2"], oos_x2)
tail(Data_x2, n = 40)


#Turn Forecast with Index into Forecast with Date
forecasts_x2_xts <- xts(forecasts_x2_data,
                        order.by = index(oos_x2))
forecasts_x2_xts

#Merge it together with the original data
Data_x2_combined <- merge(Data_x2, forecasts_x2_xts)
head(Data_x2_combined)
tail(Data_x2_combined,40)


#Plot Chart

plot(Data_x2_combined ["2020-11/", c("x2", "f_mean", "f_lower", "f_upper")], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 day forecast of x2",
     col = c("black", "blue", "red", "red"))

#Extract Data for Evaluating the Forecast:
Data_x2_Eva <- tail(Data_x2_combined, 30)  
Data_x2_Eva

Data_x2_Eva$mae   <-  abs(Data_x2_Eva$x2 - Data_x2_Eva$f_mean)
Data_x2_Eva$mse   <-  (Data_x2_Eva$x2 - Data_x2_Eva$f_mean) ^ 2
Data_x2_Eva$mape  <-  abs((Data_x2_Eva$x2 - Data_x2_Eva$f_mean)/Data_x2_Eva$x2)
Data_x2_Eva$amape <-  abs((Data_x2_Eva$x2 - Data_x2_Eva$f_mean)/(Data_x2_Eva$x2 + Data_x2_Eva$f_mean))
Data_x2_Eva

colMeans(Data_x2_Eva[, c("mae", "mse", "mape", "amape")])
apply(Data_x2_Eva[, c("mae", "mse", "mape", "amape")], 2, FUN = median)

#7. Johansen cointegration test

#We will assume the K=5 lag structure:
  
johan.test.trace <- 
ca.jo(Data[,1:2], # data 
        ecdet = "const", # "none" for no intercept in cointegrating equation, 
        # "const" for constant term in cointegrating equation and 
        # "trend" for trend variable in cointegrating equation
        type = "trace",  # type of the test: trace or eigen
        K = 5           # lag order of the series (levels) in the VAR
)


summary(johan.test.trace) 

cbind(summary(johan.test.trace)@teststat, summary(johan.test.trace)@cval)
  
#Test Statistic < Critical Value: CANNOT reject the null
#Test Statistic > Critical Value: reject the null
#First we start with r=0, (no cointegrating vector): t-test statistic > Critical value: we reject the null hypothesis about NO cointerating vector.
#Next, we test the hypothesis of r=1: Test Statistic < Critical Value: we CANNOT reject the null about 1 cointegrating vector.
#The model has ONLY 1 cointegrating vector.

summary(johan.test.trace)@V

#Weights W:
  
summary(johan.test.trace)@W

##Check for another type of test: for Eigen:

johan.test.eigen <- 
  ca.jo(Data[,1:2], # data 
        ecdet = "const", # "none" for no intercept in cointegrating equation, 
        # "const" for constant term in cointegrating equation and 
        # "trend" for trend variable in cointegrating equation
        type = "eigen",  # type of the test: trace or eigen
        K = 5           # lag order of the series (levels) in the VAR
) 
summary(johan.test.eigen) 

#The conclusion stays the same: There is only one cointegrating vector.

#8. The VECM model

#Define the specification of the VECM model, with cointegrating vector from either trace or eigen test from Johansen test.

Data.vec5 <- cajorls(johan.test.eigen, # defined specification
                        r = 1) # number of cointegrating vectors

summary(Data.vec5)

#Summary results:
summary(Data.vec5$rlm)


#extract the cointegrating vector:
Data.vec5$beta

#Reparametrize the VEC model into VAR:
Data.vec5.asVAR <- vec2var(johan.test.eigen, r = 1)

#Check result:
Data.vec5.asVAR

#Calculate and plot Impulse Response Functions:
plot(irf(Data.vec5.asVAR, n.ahead = 36), ask = FALSE)

#Perform variance decomposition:
#ERROR HEREEEE!
#plot(fevd(Data.vec5.asVAR, n.ahead = 36), ask = FALSE)

#Check if model residuals are autocorrelated or not:
#Residuals can be extracted only from the VAR reparametrized model.

head(residuals(Data.vec5.asVAR))
serial.test(Data.vec5.asVAR)

#p-value = 0.0831 > p-critical = 5%
#The null about no-autocorrelation is Not Rejected 
#-> There is auto-correlated.????

#Plot ACF and PACF for the model:
#plot(serial.test(Data.vec5.asVAR))

#9. Forecasting based on the VECM

Data.vec5.fore <- 
  predict(
    vec2var(
      johan.test.eigen, 
      r = 1),     # no of cointegrating vectors 
    n.ahead = 30, # forecast horizon
    ci = 0.95)    # confidence level for intervals

summary(Data.vec5.fore)

#VEC forecasts for x1
Data.vec5.fore$fcst$x1

#VEC forecasts for x2
Data.vec5.fore$fcst$x2

#Lets store it as an xts object. The correct set of dates (index) can be extracted from the out_of_sample xts data object.
tail(index(out_of_sample), 30)

x1_forecast <- xts(Data.vec5.fore$fcst$x1[,-4], 
                    # we exclude the last column with CI
                    tail(index(out_of_sample), 30))

#Correction of the variable names:
names(x1_forecast) <- c("x1_fore", "x1_lower", "x1_upper")

#Apply similarly for x2:
x2_forecast <- xts(Data.vec5.fore$fcst$x2[, -4],
                    # we exclude the last column with CI
                    tail(index(out_of_sample), 30))

names(x2_forecast) <- c("x2_fore", "x2_lower", "x2_upper")

#Merge forecast into orignial data:
Data.fore <- merge(out_of_sample[,1:2], 
                 x1_forecast,
                 x2_forecast)

tail(Data.fore,40) 

plot(Data.fore ["2020-11/", c("x1", "x1_fore", "x1_lower", "x1_upper")], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 day forecast of x1",
     col = c("black", "blue", "red", "red"))


plot(Data.fore ["2020-11/", c("x2", "x2_fore", "x2_lower", "x2_upper")], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 day forecast of x2",
     col = c("black", "blue", "red", "red"))


#9. Calculate forecast accuracy measures

#Extract the out-of-sample data to evaluate:
Data.fore2 <- Data.fore[,-30]

Data.fore2$mae.x1   <-  abs(Data.fore2$x1 - Data.fore2$x1_fore)
Data.fore2$mse.x1   <-  (Data.fore2$x1 - Data.fore2$x1_fore)^2
Data.fore2$mape.x1  <-  abs(Data.fore2$x1 - Data.fore2$x1_fore)/Data.fore2$x1
Data.fore2$amape.x1 <-  abs(Data.fore2$x1 - Data.fore2$x1_fore)/(Data.fore2$x1 + Data.fore2$x1_fore)

Data.fore2$mae.x2   <-  abs(Data.fore2$x2 - Data.fore2$x2_fore)
Data.fore2$mse.x2   <-  (Data.fore2$x2 - Data.fore2$x2_fore)^2
Data.fore2$mape.x2  <-  abs(Data.fore2$x2 - Data.fore2$x2_fore)/Data.fore2$x2
Data.fore2$amape.x2 <-  abs(Data.fore2$x2 - Data.fore2$x2_fore)/(Data.fore2$x2 + Data.fore2$x2_fore)

# and calculate its averages

colMeans(Data.fore2[, c("mae.x1", 
                      "mse.x1",
                      "mape.x1",
                      "amape.x1")], na.rm = TRUE)

colMeans(Data.fore2[, c("mae.x2", 
                      "mse.x2",
                      "mape.x2",
                      "amape.x2")], na.rm = TRUE)  
