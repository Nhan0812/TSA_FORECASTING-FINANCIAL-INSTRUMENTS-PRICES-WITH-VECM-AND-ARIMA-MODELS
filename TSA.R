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
plot(Data,  
     main = "Line Chart for 10 Time Series",  
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     legend.loc = "bottomleft",
     type="l")

#Plot 10 Time Series separately
plot(Data,
     main = "Graph for 10 Time Series",
     major.ticks = "years", 
     grid.ticks.on = "years",
     legend.loc = "bottomleft",
     multi.panel = TRUE,
     yaxis.same = FALSE,
     type="l")

#So our group decide to choose Time Series 1 and 2 for analysis.

#4. Create first difference
Data$dx1 <- diff.xts(Data$x1, na.pad = FALSE)
Data$dx2 <- diff.xts(Data$x2, na.pad = FALSE)
head(Data)
#Plot both variables on the graph:
plot(Data[, 1:2],
     col = c("black", "blue"),
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

#Interpret result: The p-bg is very small -> there are auto-correlation in residuals, even when we add the augmentations, we still cant get rid of auto-correlation.
#Hence, we test for its 1st difference. 

testdf(variable = Data$dx1, 
       max.augmentations = 3)

#The p-bg is greater than 5% -> there are no auto-correlations in residual
#Next, we test the stationary of its 1st lag, by checking p-adf. 
#p-adf is smaller than 5%. We can reject the null in the case of the first differences about non-stationary.
#Its 1st lag is stationary.
#We can conclude that the Time Series X1 is integrated of order 1.

#Testing the order of the time series 2 and its difference:

testdf(variable = Data$x2,
       max.augmentations = 3)

#Interpret result: Similar to X1, the p-bg here is smaller than 5% -> there are auto-correlation in residuals, even when we add the augmentations, we still cant get rid of auto-correlation.
#Hence, we test for its 1st difference. 

testdf(variable = Data$dx2, 
       max.augmentations = 3)

#By adding 1 augmentation, there will be no auto correlation in Residuals.
#p-adf is greater than 5%, we can reject the null hypothesis in the case the non-stationary of first differences. 
#Its 1st lag is stationary.
#we can conclude that the Time Series X2 is integrated of order 1.

#Conclusion:
#As a result, both variables has the same order 1: ∼I(1)
#so in the next step we can check whether they are cointegrated or not.

#3. Estimating the Cointegrated Vector:

model.coint <- lm(x1 ~ x2, data = Data)

#Examine the model summary:
summary(model.coint)

#Both the intercept and Coefficient for x2 are statistically significant.
#The model is significantly explained by x2.

#We further test the stationary of the residuals:
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

#The PACF shown is suggestive of an AR(5) model or AR(7). 
#So an initial candidate model is an ARIMA(5,1,0). 
#We also have some variations of this model: ARIMA(5,1,1), ARIMA(4,1,0), ARIMA(3,1,0), ARIMA(7,1,0)

#Step 2 - model estimation

arima510 <- Arima(Data$x1,  # variable	
                  order = c(5, 1, 0)  # (p,d,q) parameters	
)

arima511 <- Arima(Data$x1,  # variable	
                  order = c(5, 1, 1)  # (p,d,q) parameters	
)

arima410 <- Arima(Data$x1,  # variable	
                  order = c(4, 1, 0)  # (p,d,q) parameters	
)

arima310 <- Arima(Data$x1,  # variable	
                  order = c(3, 1, 0)  # (p,d,q) parameters	
)

arima313 <- Arima(Data$x1,  # variable	
                  order = c(3, 1, 3),  # (p,d,q) parameters	
)

arima710 <- Arima(Data$x1,  # variable	
                  order = c(7, 1, 0),  # (p,d,q) parameters	
)	

#The above model is not included constant. Let's include constant into the model: 
arima510_2 <- Arima(Data$x1,  # variable	
                  order = c(5, 1, 0),  # (p,d,q) parameters
                  include.constant = TRUE)  # including a constant

coeftest(arima510_2)

#Adding the constant did not change the result much for ALL the above models, 
#so we keep the model without constant


#Summary All The Models:

coeftest(arima510)
# All the terms/coeffs are significant.

coeftest(arima511)
# All the terms/coeffs are significant, except MA1.

coeftest(arima410)
# All the terms/coeffs are significant.

coeftest(arima310)
# All the terms/coeffs are significant.

coeftest(arima313)
# All the terms/coeffs are significant, except MA1, MA3.

coeftest(arima710)
# All the terms/coeffs are significant, except AR5, AR7.

#Step 3 - model diagnostics
plot(resid(arima510),col = "royalblue")
plot(resid(arima511),col = "royalblue")
plot(resid(arima410),col = "royalblue")
plot(resid(arima310),col = "royalblue")
plot(resid(arima313),col = "royalblue")
plot(resid(arima710),col = "royalblue")

#Lets check the ACF and the PACF of the Residual values:
plot_ACF_PACF_resids(arima510)
plot_ACF_PACF_resids(arima511)
plot_ACF_PACF_resids(arima410)
plot_ACF_PACF_resids(arima310)
plot_ACF_PACF_resids(arima313)
plot_ACF_PACF_resids(arima710)

#The Ljung-Box test:

Box.test(resid(arima510), type = "Ljung-Box", lag = 10)
Box.test(resid(arima510), type = "Ljung-Box", lag = 15)
Box.test(resid(arima510), type = "Ljung-Box", lag = 20)
Box.test(resid(arima510), type = "Ljung-Box", lag = 25)

Box.test(resid(arima511), type = "Ljung-Box", lag = 10)
Box.test(resid(arima511), type = "Ljung-Box", lag = 15)
Box.test(resid(arima511), type = "Ljung-Box", lag = 20)
Box.test(resid(arima511), type = "Ljung-Box", lag = 25)

Box.test(resid(arima410), type = "Ljung-Box", lag = 10)
Box.test(resid(arima410), type = "Ljung-Box", lag = 15)
Box.test(resid(arima410), type = "Ljung-Box", lag = 20)
Box.test(resid(arima410), type = "Ljung-Box", lag = 25)

Box.test(resid(arima310), type = "Ljung-Box", lag = 10)
Box.test(resid(arima310), type = "Ljung-Box", lag = 15)
Box.test(resid(arima310), type = "Ljung-Box", lag = 20)
Box.test(resid(arima310), type = "Ljung-Box", lag = 25)
#We have very low p-values, greater than 5% 
#-> Reject Ho about no autocorrelation 
#Hence, the residual is auto-correlated.

Box.test(resid(arima313), type = "Ljung-Box", lag = 10)
Box.test(resid(arima313), type = "Ljung-Box", lag = 15)
Box.test(resid(arima313), type = "Ljung-Box", lag = 20)
Box.test(resid(arima313), type = "Ljung-Box", lag = 25)

Box.test(resid(arima710), type = "Ljung-Box", lag = 10)
Box.test(resid(arima710), type = "Ljung-Box", lag = 15)
Box.test(resid(arima710), type = "Ljung-Box", lag = 20)
Box.test(resid(arima710), type = "Ljung-Box", lag = 25)

#We have very large p-values, greater than 5% for all models (except ARIMA(3,1,0))
#-> fail to reject Ho about no autocorrelation 
#Hence, the residual is white-noise.


#Plot graph for Ljung-Box test:
#ARIMA(5,1,0)
bj_pvalues = c()

for(i in c(1:100)){
  bj = Box.test(resid(arima510), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}

plot(bj_pvalues, type='l')

abline(h = 0.05, col='red')


#ARIMA(3,1,3)
bj_pvalues = c()

for(i in c(1:100)){
  bj = Box.test(resid(arima313), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}

plot(bj_pvalues, type='l')

abline(h = 0.05, col='red')


#Step 4. Evaluate Model:

AIC(arima510, arima511,
    arima410, 
    arima310, arima313,
    arima710)
#Model ARIMA(3,1,3) returns the lowest AIC test.

BIC(arima510, arima511,
    arima410,
    arima310, arima313,
    arima710)
#Model ARIMA(5,1,0) returns the lowest BIC test; Model ARIMA(3,1,3) comes second.

#However, we should prefer AIC over BIC.

#From the perspective of sensibility, the ARIMA(3,1,3) seems to be the most attractive one: 
#all terms are significant, residuals are white noise and we observe low values of information criteria (IC).


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

arima511_x2 <- Arima(Data$x2,  # variable	
                     order = c(5, 1, 1)  # (p,d,q) parameters
)

arima410_x2 <- Arima(Data$x2,  # variable	
                     order = c(4, 1, 0)  # (p,d,q) parameters
)  

arima310_x2 <- Arima(Data$x2,  # variable	
                     order = c(3, 1, 0)  # (p,d,q) parameters
)

arima313_x2 <- Arima(Data$x2,  # variable	
                     order = c(3, 1, 3)  # (p,d,q) parameters
)  

#The above model is not included constant. Let's include constant into the model: 
arima510_2_x2 <- Arima(Data$x2,  # variable	
                    order = c(5, 1, 0),  # (p,d,q) parameters
                    include.constant = TRUE)  # including a constant

coeftest(arima510_2_x2)
#Adding the constant did not change the result much for ALL Models, 
#so we keep the model without constant.


#Summary Results:

coeftest(arima510_x2)
# All parameters are significant.

coeftest(arima511_x2)
# All parameters are significant, except MA1.

coeftest(arima410_x2)
# All parameters are significant.

coeftest(arima310_x2)
# All parameters are significant.

coeftest(arima313_x2)
# All parameters are significant, except MA3.

#Step 3 - model diagnostics
plot(resid(arima510_x2))
plot(resid(arima511_x2))
plot(resid(arima410_x2))
plot(resid(arima310_x2))
plot(resid(arima313_x2))

# The Ljung-Box test:
Box.test(resid(arima510_x2), type = "Ljung-Box", lag = 10)
Box.test(resid(arima510_x2), type = "Ljung-Box", lag = 15)
Box.test(resid(arima510_x2), type = "Ljung-Box", lag = 20)
Box.test(resid(arima510_x2), type = "Ljung-Box", lag = 25)
#-> Very large p-values: The residual is white-noise  

Box.test(resid(arima511_x2), type = "Ljung-Box", lag = 10)
Box.test(resid(arima511_x2), type = "Ljung-Box", lag = 15)
Box.test(resid(arima511_x2), type = "Ljung-Box", lag = 20)
Box.test(resid(arima511_x2), type = "Ljung-Box", lag = 25)
#-> Very large p-values: The residual is white-noise  

Box.test(resid(arima410_x2), type = "Ljung-Box", lag = 10)
Box.test(resid(arima410_x2), type = "Ljung-Box", lag = 15)
Box.test(resid(arima410_x2), type = "Ljung-Box", lag = 20)
Box.test(resid(arima410_x2), type = "Ljung-Box", lag = 25)
#-> Not so large p-values: The residual is not really a white-noise for lag20. 

Box.test(resid(arima310_x2), type = "Ljung-Box", lag = 10)
Box.test(resid(arima310_x2), type = "Ljung-Box", lag = 15)
Box.test(resid(arima310_x2), type = "Ljung-Box", lag = 20)
Box.test(resid(arima310_x2), type = "Ljung-Box", lag = 25)
#-> Very small p-values: The residual is auto-correlated

Box.test(resid(arima313_x2), type = "Ljung-Box", lag = 10)
Box.test(resid(arima313_x2), type = "Ljung-Box", lag = 15)
Box.test(resid(arima313_x2), type = "Ljung-Box", lag = 20)
Box.test(resid(arima313_x2), type = "Ljung-Box", lag = 25)
#-> Very large p-values: The residual is white-noise  

#Plot graph for Ljung-Box test:
bj_pvalues = c()

for(i in c(1:100)){
  bj = Box.test(resid(arima510_x2), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}

plot(bj_pvalues, type='l')

abline(h=0.05, col='red')


#Step 4. Evaluate Model:

AIC(arima510_x2,arima511_x2,
    arima410_x2,  
    arima310_x2, arima313_x2)
#ARIMA(3,1,3)

BIC(arima510_x2,arima511_x2,
    arima410_x2,  
    arima310_x2, arima313_x2)
#ARIMA(5,1,0)


#From the perspective of sensibility, the ARIMA(5,1,0) seems to be the most attractive one: 
#all terms are significant, residuals are white noise and we observe low values of information criteria (IC).

#-> Suggestion model: arima510_x2 - ARIMA(5,1,0)


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

forecasts_x1 <- forecast(arima313, # model for prediction
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

ARIMA_x1 <- colMeans(Data_x1_Eva[, c("mae", "mse", "mape", "amape")])
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

ARIMA_x2 <- colMeans(Data_x2_Eva[, c("mae", "mse", "mape", "amape")])
apply(Data_x2_Eva[, c("mae", "mse", "mape", "amape")], 2, FUN = median)

#7. Johansen cointegration test

#Determine the lag for Johansen test:
VARselect(Data[,1:2], 
          lag.max = 7)
#Three out of four test suggest K = 6.

#We will choose the K=6 lag structure:
  
johan.test.trace <- 
ca.jo(Data[,1:2], # data 
        ecdet = "const", # "none" for no intercept in cointegrating equation, 
        # "const" for constant term in cointegrating equation and 
        # "trend" for trend variable in cointegrating equation
        type = "trace",  # type of the test: trace or eigen
        K = 6           # lag order of the series (levels) in the VAR
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
        K = 6           # lag order of the series (levels) in the VAR
) 
summary(johan.test.eigen) 

#The conclusion stays the same: There is only one cointegrating vector.

#8. The VECM model

#Define the specification of the VECM model, with cointegrating vector from either trace or eigen test from Johansen test.

Data.vec6 <- cajorls(johan.test.eigen, # defined specification
                        r = 1) # number of cointegrating vectors

summary(Data.vec6)

#Summary results:
summary(Data.vec6$rlm)

#extract the cointegrating vector:
Data.vec6$beta

#extract the adjustment coefficients (check for sign to determine whether ECM works or not):
johan.test.eigen@W
#-> The adjustment Coeff has different sign here -> ECM works.

#Reparametrize the VEC model into VAR:
Data.vec6.asVAR <- vec2var(johan.test.eigen, r = 1)

#Check result:
Data.vec6.asVAR

#Calculate and plot Impulse Response Functions:
plot(irf(Data.vec6.asVAR, n.ahead = 36))
#The residuals seem to be stable: first increases then decreases.


#Perform variance decomposition:
plot(fevd(Data.vec6.asVAR, n.ahead = 36))


#Check if model residuals are autocorrelated or not:
#Residuals can be extracted only from the VAR reparametrized model.

head(residuals(Data.vec6.asVAR))
serial.test(Data.vec6.asVAR)

#p-value = 0.6616 > p-critical = 5%
#The null about no-autocorrelation is fail to Reject.
#=> There is no auto-correlation in Residuals. 

#Plot ACF and PACF for the model:
plot(serial.test(Data.vec6.asVAR))

#Checking the Nomarlity for x1 and x2 by creating Histogram:

#For Timeseries x1:
Data.vec6.asVAR %>%
  residuals() %>%
  as_tibble() %>%
  ggplot(aes(`resids of x1`)) +
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "pink") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(residuals(Data.vec6.asVAR)[, 1]), 
                            sd = sd(residuals(Data.vec6.asVAR)[, 1]))) +
  theme_bw() + 
  labs(
    title = "Density of x1 residuals", 
    y = "", x = "",
    caption = "source: own calculations"
  )

#For Timeseries x2:
Data.vec6.asVAR %>%
  residuals() %>%
  as_tibble() %>%
  ggplot(aes(`resids of x2`)) +
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "pink") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(residuals(Data.vec6.asVAR)[, 2]), 
                            sd = sd(residuals(Data.vec6.asVAR)[, 2]))) +
  theme_bw() + 
  labs(
    title = "Density of x2 residuals", 
    y = "", x = "",
    caption = "source: own calculations"
  )

#We can also check it formally by using the Jarque-Bera (JB) test.

normality.test(Data.vec6.asVAR)

#p-value > 0.05 => Fail to reject Ho about the normality
#Conclusion: the residual has a normal distribution.

#9. Forecasting based on the VECM

Data.vec6.fore <- 
  predict(
    vec2var(
      johan.test.eigen, 
      r = 1),     # no of cointegrating vectors 
    n.ahead = 30, # forecast horizon
    ci = 0.95)    # confidence level for intervals

summary(Data.vec6.fore)

#VEC forecasts for x1
Data.vec6.fore$fcst$x1

#VEC forecasts for x2
Data.vec6.fore$fcst$x2

#Lets store it as an xts object. The correct set of dates (index) can be extracted from the out_of_sample xts data object.
tail(index(out_of_sample), 30)

x1_forecast <- xts(Data.vec6.fore$fcst$x1[,-4], 
                    # we exclude the last column with CI
                    tail(index(out_of_sample), 30))

#Correction of the variable names:
names(x1_forecast) <- c("x1_fore", "x1_lower", "x1_upper")

#Apply similarly for x2:
x2_forecast <- xts(Data.vec6.fore$fcst$x2[, -4],
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


#10. Calculate forecast accuracy measures

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

VECM_x1 <- colMeans(Data.fore2[, c("mae.x1", 
                      "mse.x1",
                      "mape.x1",
                      "amape.x1")], na.rm = TRUE)

VECM_x2 <- colMeans(Data.fore2[, c("mae.x2", 
                      "mse.x2",
                      "mape.x2",
                      "amape.x2")], na.rm = TRUE)  

#6. Comparing VECM model’s forecasts with ARIMAs:

result <- rbind(ARIMA_x1,VECM_x1, ARIMA_x2,VECM_x2)

result %>%
  knitr::kable(digits = 4) %>%
  kableExtra::kable_styling(full_width = F,
                            bootstrap_options = c("striped", 
                                                  "hover", 
                                                  "condensed"))

#Conclusion: 

#For Timeseries x1:
#The forecasts from VECM model outperforms that of ARIMA.

#For Timeseries x2:
#The forecasts from ARIMA model outperforms that of VECM.

#####################################################################################
## Additional Task

str(Data)
head(Data)
Data[,1:2]

VARselect(Data[,1:2], 
          lag.max = 7)

library(formattable)
VARselect(Data[,1:2],      
          lag.max = 7,     
          season = 12) 

Data.var6 <- VAR(Data[,1:2],
                    p = 6)
summary(Data.var6)
plot(Data.var6)
serial.test(Data.var6)
serial.test(Data.var6, type = "BG")



########
Data.var5 <- VAR(Data[,1:2],
                 p = 5)
summary(Data.var5)
plot(Data.var5)
serial.test(Data.var5)
serial.test(Data.var6, type = "BG")

AIC(Data.var6,Data.var5)
BIC(Data.var6,Data.var5)
-> VAR6 is still better. Chosse to be our model.

plot(irf(Data.var6, n.ahead = 36))
