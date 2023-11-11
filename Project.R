### Import Library
install.packages("knitr")
# Load useful libraries we will use for time series.
library(AER)
library(readxl)
library(xts)
library(zoo)
library(dynlm)
library(stargazer)
library(urca)
library(orcutt)
library(vars)
library(fGarch)
library(quantmod)

# Set Working Directory where the data is stored on your computer
#(Steps: Session --> Set Working Directory --> Choose Directory)

### Load and Format Data
# load UK macroeconomic data
UK_Manufacturing_AWE_Quaterly <- read_excel(path= "AWE_SeasonallyAdjusted_Manufacturing_CleanedQuarterly.xls")
UK_Manufacturing_AWE_Quaterly$Date <- as.yearqtr(UK_Manufacturing_AWE_Quaterly$Date, format = "%Y Q%q") #Format the "Date" column as year-quarter data
UK_Manufacturing_AWE_Quaterly_xts <- xts(UK_Manufacturing_AWE_Quaterly$AWE,UK_Manufacturing_AWE_Quaterly$Date)

### TASK 1.1.2 - Time Series Plots
plot(as.zoo(UK_Manufacturing_AWE_Quaterly_xts),col = "steelblue",
     lwd = 2,
     ylab = "Quarterly Manufacturing AWE",
     xlab = "Date",
     main = "U.K. Average Weekly Earnings(Quaterly)")
# It seems that there are outliers in 2020. We are going to look at data up to 2019 instead.
UK_M_AWE_Q_xts <- UK_Manufacturing_AWE_Quaterly_xts["2000::2019"]  #Take subsets of our data.
plot(as.zoo(UK_M_AWE_Q_xts),col = "steelblue",
     lwd = 2,
     ylab = "Quarterly Manufacturing AWE",
     xlab = "Date",
     main = "U.K. Average Weekly Earnings (Quaterly)")

### 1.2.1 - Estimate an Autoregression Model
UK_AWE_lm <- lm(UK_M_AWE_Q_xts ~ lag(UK_M_AWE_Q_xts,1))  #run an AR(1) regression
summary(UK_AWE_lm)
#Compute BIC for AR models using a convenient function.
IC <- function(model) {
  ssr <- sum(model$residuals^2) 
  t <- length(model$residuals)  
  p <- length(model$coef) - 1 
  return(
    round(c("p" = p, # 'p' here is just the length of the vector of coefficients
            "BIC" = log(ssr/t) + ((p + 1) * log(t)/t), 
            "R2" = summary(model)$r.squared), 
          4) #round the output to 4 decimal points.
  )
}
for (p in 1:4) { 
  print(IC(lm(UK_M_AWE_Q_xts ~ lag(UK_M_AWE_Q_xts,1:p))))   
}
#Based on the above we select p = 1

# test unit root around a time trend
UK_trend_AWE_lm <- lm(diff(UK_M_AWE_Q_xts) ~ time(UK_M_AWE_Q_xts) + lag(UK_M_AWE_Q_xts,1))
coeftest(UK_trend_AWE_lm)
AR1_UK_trend_AWE_manual = coeftest(UK_trend_AWE_lm)[3,3]
#automatic test unit root around a time trend
unit_root_test_UK_M_AWE_Q_xts = ur.df(UK_M_AWE_Q_xts, lags=0, type='trend')
unit_root_test_UK_M_AWE_Q_xts@teststat[1]
unit_root_test_UK_M_AWE_Q_xts@cval[1,]

# Taking the first difference
transformedUK_AWE_lm <- lm(diff(UK_M_AWE_Q_xts) ~ lag(diff(UK_M_AWE_Q_xts),1))
summary(transformedUK_AWE_lm)

# plot transformed data
plot(as.zoo(diff(UK_M_AWE_Q_xts)),col = "steelblue",
          lwd = 2,
          ylab = "Quarterly difference change in AWE",
          xlab = "Date",
          main = "U.K. Difference change in AWE (Quaterly)")

#Compute BIC for AR transformed models
for (p in 1:4) { 
  print(IC(lm(diff(UK_M_AWE_Q_xts) ~ lag(diff(UK_M_AWE_Q_xts),1:p))))   
}

### QLR Test 
#  set up the interval 
num_periods <- length(UK_M_AWE_Q_xts)
# Then, we use these lines to approximate the first fifteen percent and last 
# fifteen percent of the time series:
tau_zero =round(0.15*num_periods,digits=0)
tau_one =round((1-0.15)*num_periods,digits=0)
# Using these tau's to approximate the start and end of the testing interval,
# we can figure out how many tests we need to run, which tells us how many loops
# we need to iterate over:
num_tests <- tau_one - tau_zero + 1
tau <- seq(tau_zero, tau_one)

chow_test_statistics <- array(num_tests)
for (i in 1:num_tests) {
  # for each of the tests we run, we set up D with the appropriate tau
  D <- 1*(time(UK_M_AWE_Q_xts) > time(UK_M_AWE_Q_xts)[tau[i]])  
  # Use the appropriate regression model, including D interacted with all the regressors
  Chow_test_model = lm(UK_M_AWE_Q_xts ~ lag(UK_M_AWE_Q_xts,1) + D + (D*lag(UK_M_AWE_Q_xts,1)))
  # Use the output in your coeftest table above to set up the F-test hypothesis below:
  chow_test = linearHypothesis(Chow_test_model, c("D=0", "lag(UK_M_AWE_Q_xts, 1):D"), test="F", white.adjust = FALSE)
  # Store the Chow test statistic:
  chow_test_statistics[i] = chow_test$F[2]
}
data.frame("Level" = tau, "F-stat" = chow_test_statistics)

# Now, we are going to get R to pick the maximum Chow test stat.
# This is the test statistic you need to compare against the appropriate critical values.
QLR_test_stat = max(chow_test_statistics)
cat("QLR Test Statistic: ", QLR_test_stat, "\n")
# This line tells us at WHICH tau we get the highest Chow test stat
tau_hat <- tau[which.max(chow_test_statistics)]
# Now, let's see which observation the likely break is at:
UK_M_AWE_Q_xts[tau_hat]

# Now, we can see what happens AT the period where our QLR test thinks the break 
# was most likely by running the Chow test model and interpreting specific 
# coefficients.
D <- 1*(time(UK_M_AWE_Q_xts) > time(UK_M_AWE_Q_xts)[tau_hat])  
Chow_test_model = lm(UK_M_AWE_Q_xts ~ lag(UK_M_AWE_Q_xts,1) + D + (D*lag(UK_M_AWE_Q_xts,1)))
coeftest(Chow_test_model)

# Finally, we can see what this looks like by plotting our data with a line at
# break period to see if the results look like what we see on the plot
plot(as.zoo(UK_M_AWE_Q_xts),col = "steelblue",
     lwd = 2,
     ylab = "Quarterly Manufacturing AWE",
     xlab = "Date",
     main = "U.K. Average Weekly Earnings (Quaterly)")
abline(v=as.yearqtr("2010 Q1"),col="red")

stargazer(UK_AWE_lm, transformedUK_AWE_lm, digit =3, header = F, type ='html',out = 'Assignment_Table3.doc')

### 1.2.2 Estimate an Autoregressive Distributed Lag Model

#Load Xt data
UK_weekly_hour <- read_excel(path= "WeeklyHoursWorked_UK_All_clean.xlsx")
UK_weekly_hour$Date <- as.yearqtr(UK_weekly_hour$Date, format = "%Y Q%q") #Format the "Date" column as year-quarter data
UK_weekly_hour_xts <- xts(UK_weekly_hour$Hours,UK_weekly_hour$Date)
UK_weekly_hour_2000_2019 <- UK_weekly_hour_xts["2000::2019"]


# plot Xt
plot(as.zoo(UK_weekly_hour_2000_2019),col = "steelblue",
     lwd = 2,
     ylab = "Weekly hours (million)",
     xlab = "Date",
     main = "U.K. Weekly hours (Quaterly)")

#Test Unit root of Xt
UnitRoot_UK_hour_lm <- lm(diff(UK_weekly_hour_2000_2019) ~ lag(UK_weekly_hour_2000_2019,1))
UK_hour_teststat_manual = coeftest(UnitRoot_UK_hour_lm)[2,3]

unit_root_test_hour = ur.df(UK_weekly_hour_2000_2019, lags=0, type='drift')

# These two commands retrieve the parameters of interest we care about:
DickeyFuller_teststat_automatic_hour = unit_root_test_hour@teststat[1]
DickeyFuller_CriticalValues_hour = unit_root_test_hour@cval[1,]

# ADL(1,1)
ADL11 <- lm(diff(UK_M_AWE_Q_xts) ~ lag(diff(UK_M_AWE_Q_xts),1) + lag(diff(UK_weekly_hour_2000_2019),1))
summary(ADL11)  

## BIC computation for ADL(1,1)
for (p in 1:4) {
  print(IC(lm(diff(UK_M_AWE_Q_xts) ~ lag(diff(UK_M_AWE_Q_xts),1:p) + lag(diff(UK_weekly_hour_2000_2019),1:p)))) 
}
  
## Granger Causality
#ADL(2,2)
GC_ADL22 <- lm(diff(UK_M_AWE_Q_xts) ~ lag(diff(UK_M_AWE_Q_xts),1) + lag(diff(UK_M_AWE_Q_xts),2) 
               + lag(diff(UK_weekly_hour_2000_2019),1) + lag(diff(UK_weekly_hour_2000_2019),2))
summary(GC_ADL22)
linearHypothesis(GC_ADL22, c("lag(diff(UK_weekly_hour_2000_2019), 1)=0", 
                             "lag(diff(UK_weekly_hour_2000_2019), 2)=0"), test ="F", white.adjust=FALSE)  

#ADL(4,4)
GC_ADL44 <- lm(diff(UK_M_AWE_Q_xts) ~ lag(diff(UK_M_AWE_Q_xts),1) + 
                 lag(diff(UK_M_AWE_Q_xts),2) + lag(diff(UK_M_AWE_Q_xts),3) + 
                 lag(diff(UK_M_AWE_Q_xts),4) + lag(diff(UK_weekly_hour_2000_2019),1) + 
                 lag(diff(UK_weekly_hour_2000_2019),2) +lag(diff(UK_weekly_hour_2000_2019),3) + 
                 lag(diff(UK_weekly_hour_2000_2019),4))
summary(GC_ADL44)
linearHypothesis(GC_ADL44, c("lag(diff(UK_weekly_hour_2000_2019), 1)=0", 
                             "lag(diff(UK_weekly_hour_2000_2019), 2)=0", 
                             "lag(diff(UK_weekly_hour_2000_2019), 3)=0", 
                             "lag(diff(UK_weekly_hour_2000_2019), 4)=0"), test ="F", white.adjust=FALSE)

stargazer(ADL11, GC_ADL44, digit =3, header = F, type ='html',out = 'Assignment_Table5.doc')

### 1.2.3 Check Out-Of-Sample Forecast Performance  
## Test for excluded observations
end_dates <- seq(2015.00, 2019.75, 0.25)     
P <- length(end_dates)  

forecasts <- array(c(0),dim = c(P))
true_outcomes <- array(c(0),dim = c(P))
PsOOSF_errors <- array(c(0),dim = c(P))
SER <- array(c(0),dim = c(P))

for (i in 1:P) {
  EndDate_YearQuarter = as.yearqtr(end_dates[i])  # First, we convert our data back into the Year-Quarter format we are used to.
  # Now, run our regression on the limited sample.
  # First, we limit our sample to observations whose index (i.e the DATE of that observation) comes before the designated EndDate_YearQuarter
  AWE_diff_sample <- diff(UK_M_AWE_Q_xts)[index(diff(UK_M_AWE_Q_xts)) <  EndDate_YearQuarter] # Now we use this to limit our estimating sample
  Hour_diff_sample <- diff(UK_weekly_hour_2000_2019)[index(diff(UK_weekly_hour_2000_2019)) < EndDate_YearQuarter]
  ADL11_lm <- lm(AWE_diff_sample ~ lag(AWE_diff_sample,1) + lag(Hour_diff_sample,1))
  SER[i] <- summary(ADL11_lm)$sigma
  # Save Coefficients from the ADL(1,1) above:
  beta_0_hat = ADL11_lm$coefficients[1]
  beta_1_hat = ADL11_lm$coefficients[2]
  beta_2_hat = ADL11_lm$coefficients[3]
  
  true_outcome <- diff(UK_M_AWE_Q_xts)[EndDate_YearQuarter]    # Find the true value of the outcome:
  pseudo_forecast <- (beta_0_hat + (beta_1_hat %*% diff(UK_M_AWE_Q_xts)[EndDate_YearQuarter - 0.25]) 
                      + (beta_2_hat %*% diff(UK_weekly_hour_2000_2019)[EndDate_YearQuarter - 0.25]))
  
  pseudo_error <- true_outcome  - pseudo_forecast  
  true_outcomes[i] <- true_outcome    
  forecasts[i] <- pseudo_forecast
  PsOOSF_errors[i] <- pseudo_error
}

SER_within_sample <- SER[1] # Here we take the first value from our array of SER's. This is the original SER of our first within-sample regression.
Estimated_RMSFE_OutOfSample <- sd(PsOOSF_errors) # Here we estimate the RMSFE by using the standard error of our pseudo out of sample forecast errors

cat("Within-Sample Errors: ", SER_within_sample, "\n")
cat("Estimated RMSFE: ", Estimated_RMSFE_OutOfSample, "\n")

#### Use the mean(), sd() and sqrt() functions to fill in the formula for the t-statistic for our t-test:
t_stat_PsOOSF_errors = mean(PsOOSF_errors)/(sd(PsOOSF_errors)*sqrt(P))
print(t_stat_PsOOSF_errors)

t.test(PsOOSF_errors)

true_outcomes_xts <- xts(true_outcomes, as.yearqtr(end_dates))
forecasts_xts <- xts(forecasts, as.yearqtr(end_dates))

# plot the AWE difference series
plot(as.zoo(true_outcomes_xts),
     col = "purple",
     lwd = 4,
     ylab = "Pounds",
     xlab = "Date",
     main = "AR2: Pseudo-Out-Of-Sample Forecasts of AWE difference")
# add the series of pseudo-out-of-sample forecasts
lines(as.zoo(forecasts_xts), 
      lwd = 4,
      lty = 2)
# shade area between curves (the pseudo forecast error)
polygon(x= c(time(true_outcomes_xts), rev(time(forecasts_xts))),
        y= c(true_outcomes, rev(forecasts)),
        col = "grey85",
        border = NA)
# add a legend
legend("bottomright",
       lty = c(1, 2, 1),
       lwd = c(2, 2, 10),
       cex = 0.5,
       col = c("purple", "black", "grey85"),
       legend = c("Actual AWE difference",
                  "Forecasted AWE difference",
                  "Pseudo AWE difference"))

## test for the non excluded observations

UK_Manufacturing_AWE_Quaterly_xts
UK_weekly_hour_2000_2020 <- UK_weekly_hour_xts["2000::2020"]

end_dates <- seq(2015.75, 2020.50, 0.25)     
P <- length(end_dates)  

forecasts <- array(c(0),dim = c(P))
true_outcomes <- array(c(0),dim = c(P))
PsOOSF_errors <- array(c(0),dim = c(P))
SER <- array(c(0),dim = c(P))

for (i in 1:P) {
  EndDate_YearQuarter = as.yearqtr(end_dates[i])  # First, we convert our data back into the Year-Quarter format we are used to.
  # Now, run our regression on the limited sample.
  # First, we limit our sample to observations whose index (i.e the DATE of that observation) comes before the designated EndDate_YearQuarter
  AWE_diff_sample <- diff(UK_Manufacturing_AWE_Quaterly_xts)[index(diff(UK_Manufacturing_AWE_Quaterly_xts)) <  EndDate_YearQuarter] # Now we use this to limit our estimating sample
  Hour_diff_sample <- diff(UK_weekly_hour_2000_2020)[index(diff(UK_weekly_hour_2000_2020)) < EndDate_YearQuarter]
  ADL11_lm <- lm(AWE_diff_sample ~ lag(AWE_diff_sample,1) + lag(Hour_diff_sample,1))
  SER[i] <- summary(ADL11_lm)$sigma
  # Save Coefficients from the ADL(1,1) above:
  beta_0_hat = ADL11_lm$coefficients[1] 
  beta_1_hat = ADL11_lm$coefficients[2]
  beta_2_hat = ADL11_lm$coefficients[3]
  
  true_outcome <- diff(UK_Manufacturing_AWE_Quaterly_xts)[EndDate_YearQuarter]    # Find the true value of the outcome:
  pseudo_forecast <- (beta_0_hat + (beta_1_hat %*% diff(UK_Manufacturing_AWE_Quaterly_xts)[EndDate_YearQuarter - 0.25]) 
                      + (beta_2_hat %*% diff(UK_weekly_hour_2000_2020)[EndDate_YearQuarter - 0.25]))
  
  pseudo_error <- true_outcome  - pseudo_forecast    
  true_outcomes[i] <- true_outcome    
  forecasts[i] <- pseudo_forecast
  PsOOSF_errors[i] <- pseudo_error
}

SER_within_sample <- SER[1] # Here we take the first value from our array of SER's. This is the original SER of our first within-sample regression.
Estimated_RMSFE_OutOfSample <- sd(PsOOSF_errors) # Here we estimate the RMSFE by using the standard error of our pseudo out of sample forecast errors

cat("Within-Sample Errors: ", SER_within_sample, "\n")
cat("Estimated RMSFE: ", Estimated_RMSFE_OutOfSample, "\n")

#### Use the mean(), sd() and sqrt() functions to fill in the formula for the t-statistic for our t-test:
t_stat_PsOOSF_errors = mean(PsOOSF_errors)/(sd(PsOOSF_errors)*sqrt(P))
print(t_stat_PsOOSF_errors)

t.test(PsOOSF_errors)

true_outcomes_xts <- xts(true_outcomes, as.yearqtr(end_dates))
forecasts_xts <- xts(forecasts, as.yearqtr(end_dates))

# plot the AWE differene series
plot(as.zoo(true_outcomes_xts),
     col = "purple",
     lwd = 4,
     ylab = "Pounds",
     xlab = "Date",
     main = "ADL(1,1): Pseudo-Out-Of-Sample Forecasts of AWE difference")
# add the series of pseudo-out-of-sample forecasts
lines(as.zoo(forecasts_xts),
      lwd = 4,
      lty = 2)
# shade area between curves (the pseudo forecast error)
polygon(x= c(time(true_outcomes_xts), rev(time(forecasts_xts))),
        y= c(true_outcomes, rev(forecasts)),
        col = "grey85",
        border = NA)
# add a legend
legend("bottomleft",
       lty = c(1, 2, 1),
       lwd = c(2, 2, 10),
       cex = 0.8,
       col = c("purple", "black", "grey85"),
       legend = c("Actual AWE difference",
                  "Forecasted AWE difference",
                  "Pseudo AWE difference"))

### 1.3 Dynamic Causal Effects

DL3_lm <- lm(diff(UK_M_AWE_Q_xts) ~ diff(UK_weekly_hour_2000_2019) + lag(diff(UK_weekly_hour_2000_2019),1) +lag(diff(UK_weekly_hour_2000_2019),2))
# Store residuals:
DL3_residuals <- DL3_lm$residuals
# Regress residuals on lag of residuals:
residuals_lm <- lm(DL3_residuals ~ lag(DL3_residuals))
# Store the coefficient on the lag of residuals -- this is our estimated phi_1
coeftest(residuals_lm)
phi1_hat <- residuals_lm$coefficients[2]
cat("Phi_1_hat from manual Cochrane-Orcutt: ", phi1_hat, "\n")

# Now, we generate our quasi-diff variables using the estimated phi1:
Y_tilde_xts <- diff(UK_M_AWE_Q_xts) - (phi1_hat* lag(diff(UK_M_AWE_Q_xts)))
X_tilde_xts <- diff(UK_weekly_hour_2000_2019) - (phi1_hat* lag(diff(UK_weekly_hour_2000_2019),1))
Xlag1_tilde_xts <- lag(diff(UK_weekly_hour_2000_2019),1) - (phi1_hat* lag(diff(UK_weekly_hour_2000_2019),2))
Xlag2_tilde_xts <- lag(diff(UK_weekly_hour_2000_2019),2) - (phi1_hat* lag(diff(UK_weekly_hour_2000_2019),3))
# Using these quasi-differenced variables, we run the GLS model:
fGLS_DL_lm <- lm(Y_tilde_xts ~ X_tilde_xts + Xlag1_tilde_xts+Xlag2_tilde_xts)
# we can recover the estimates from this, which are our dynamic multipliers:
coeftest(fGLS_DL_lm)
acf(fGLS_DL_lm$residuals, main = "GLS Residuals: Autocorrelations")

stargazer(DL3_lm, fGLS_DL_lm, digit =3, header = F, type ='html',out = 'Assignment_Table6.doc')

### 1.4 Multiperiod Forecasts

prediction_start_date = 2020.00
prediction_end_date = 2022.25
prediction_dates <- seq(prediction_start_date, prediction_end_date, 0.25)
P <- length(prediction_dates)

# Iterated Multiperiod Forecasts

iterated_forecasts <- array(c(0),dim = c(P))
AWE_diff_sample <- diff(UK_M_AWE_Q_xts)

# We will see how an AR(2) model performs, so we define this model here:
AR2_lm <- lm(AWE_diff_sample ~ lag(AWE_diff_sample) + lag(AWE_diff_sample,2))
print(coeftest(AR2_lm))

beta_0_hat = AR2_lm$coefficients[1]    # The first coefficient is the intercept, etc.
beta_1_hat = AR2_lm$coefficients[2]
beta_2_hat = AR2_lm$coefficients[3]

# Here we identify the first time period where we are forecasting.
FirstForecast_YearQuarter = as.yearqtr(prediction_start_date)  # First, we convert our data back into the Year-Quarter format we are used to.

# For the first two periods, our forecast involves true values,
#     so I have set these up for you.
# Here, FirstForecast_YearQuarter is the first period we predict, so it is
#     the same as T+1 in our usual notation. This means that to get the observation for
#     period T, we use [FirstForecast_YearQuarter - (0.25)].
iterated_forecast <- (beta_0_hat + (beta_1_hat %*% diff(UK_M_AWE_Q_xts)[FirstForecast_YearQuarter- 0.25]) 
                      + (beta_2_hat %*% diff(UK_M_AWE_Q_xts)[FirstForecast_YearQuarter - 0.50]))  
iterated_forecasts[1] <- iterated_forecast

iterated_forecast <- (beta_0_hat + (beta_1_hat %*% iterated_forecasts[1]) 
                      + (beta_2_hat %*% diff(UK_M_AWE_Q_xts)[FirstForecast_YearQuarter - 0.25]))  
iterated_forecasts[2] <- iterated_forecast

# For each observation i that we forecast, we will use the previous two
#     forecasts as Y_{t-1} and Y_{t-1}, so we use the  iterated_forecasts[i-1]
#     and iterated_forecasts[i-2] terms to capture this.
for (i in 3:P) {
  iterated_forecast <- (beta_0_hat + (beta_1_hat %*% iterated_forecasts[i-1]) 
                        + (beta_2_hat %*% iterated_forecasts[i-2]))  
  iterated_forecasts[i] <- iterated_forecast
}

# Now, let's save this as an xts so we can process it as a time series.
iterated_forecasts_xts <- xts(iterated_forecasts, as.yearqtr(prediction_dates))

### Direct Multiperiod Forecasts
direct_forecasts <- array(c(0),dim = c(P))

for (i in 1:P) {
  # Here, for every period we forecast in the future, we need to run regressions
  #     with the appropriate lag. So to predict 10 periods in the future, the 
  #     appropriate AR(2) includes regressors Y_(t-10) and Y_(t-11).
  AR2_DirectForecast_lm <- lm(AWE_diff_sample ~ lag(AWE_diff_sample,(i)) + lag(AWE_diff_sample,(1+i)))   
  
  # Every forecasting period is a separate regression, so we need to use
  #     new parameter estimates each time.
  beta_0_hat = AR2_DirectForecast_lm$coefficients[1]    # The first coefficient is the intercept, etc.
  beta_1_hat = AR2_DirectForecast_lm$coefficients[2]
  beta_2_hat = AR2_DirectForecast_lm$coefficients[3]
  
  # Here, FirstForecast_YearQuarter is the first period we predict, so it is
  #     the same as T+1 in our usual notation. This means that to get the observation for
  #     period T, we use [FirstForecast_YearQuarter - (0.25)].
  direct_forecasts[i] <- (beta_0_hat + (beta_1_hat %*% diff(UK_M_AWE_Q_xts)[FirstForecast_YearQuarter - (0.25)]) 
                          + (beta_2_hat %*% diff(UK_M_AWE_Q_xts)[FirstForecast_YearQuarter - (0.50)]))  
}
direct_forecasts_xts <- xts(direct_forecasts, as.yearqtr(prediction_dates))

### Plotting Forecasts

# Here we are just selecting the part of the series we will plot.
forecasting_region <- xts(array(c(NA),dim = c(P)), as.yearqtr(prediction_dates))
AWE_2022 <- rbind(AWE_diff_sample["2015::2019"],forecasting_region)
# plot
plot(as.zoo(AWE_2022), 
     col = "purple",
     lwd = 4,
     ylab = "Pounds",
     main = "Multiperiod AR(2) Forecasts of difference AWE")
# add the series of pseudo-out-of-sample forecasts
lines(as.zoo(direct_forecasts_xts),
      lwd = 4,
      col= "green",
      lty = 2)
lines(as.zoo(iterated_forecasts_xts),
      lwd = 4,
      col= "blue",
      lty = 2)
legend("bottomright",
       lty = c(1, 2, 2),
       lwd = c(2, 2, 10),
       cex = 0.8,
       col = c("purple", "green", "blue"),
       legend = c("Actual", "Direct Forecast", "Iterated Forecast"))

### 1.5 Cointergration 

# plot 
plot(as.zoo(UK_M_AWE_Q_xts),
    plot.type = "single",
    lty = 2,
    lwd = 2,
    col = "blue",
    xlab = "Date",
    ylab = "Pounds",
    ylim = c(300, 1100),
    main = "Earnings and Hours")
lines(as.zoo(UK_weekly_hour_2000_2019),
      col = "orange",
      lwd = 2,
      xlab = "Date",
      ylab = "Hours",
      main = "Term Spread")
legend("bottomright",
       lty = c(2, 1, 1),
       lwd = c(2, 2, 2),
       cex = 0.6,
       col = c("blue", "orange"),
       legend = c("Weekly Average Earnings", "Total hours"))

# Automatic GLS Dickey-Fuller Test
GLS_AWE <- ur.ers(UK_M_AWE_Q_xts,
                  lag.max = 1,
                  type = "DF-GLS",
                  model = "trend")
summary(GLS_AWE)

GLS_hours <- ur.ers(UK_weekly_hour_2000_2019,
                    lag.max = 1,
                    type = "DF-GLS", 
                    model="trend")
summary(GLS_hours)

## test cointegration
CointegrationFirstStep_lm <- lm(UK_M_AWE_Q_xts ~ UK_weekly_hour_2000_2019)
coeftest(CointegrationFirstStep_lm)
theta_hat = CointegrationFirstStep_lm$coefficients[2]

# Now use this estimated theta_hat to generate z_hat:
z_hat = UK_M_AWE_Q_xts - (theta_hat * UK_weekly_hour_2000_2019)

# Now, run the ADF test on z_hat, to test for stationarity:
ADFtest_z_hat <- ur.df(z_hat, lags = 1, selectlags = "BIC", type = "none")
summary(ADFtest_z_hat)

## 1.6 Volatility Clustering Analysis

GARCH_AWE <- garchFit( ~ arma(1,1) + garch(1,1), data = diff(UK_M_AWE_Q_xts)[-1], trace = F)
summary(GARCH_AWE)

stargazer(GARCH_AWE, digit =3, header = F, type ='html',out = 'Assignment_Table8.doc')

N <- length(UK_M_AWE_Q_xts)
dates <- index(diff(UK_M_AWE_Q_xts))[-N]

AWE_fitted_xts = xts(GARCH_AWE@fitted, dates)
AWE_fitted_xts

plot(as.zoo(diff(UK_M_AWE_Q_xts)), 
     type = "l", 
     col = "steelblue",
     ylab = "pounds", 
     xlab = "Date",
     main = "ARMA(1,1)-GARCH(1,1)",
     lwd = 0.25)
# add the model-fitted values of the percentage change to the plot:
lines(as.zoo(AWE_fitted_xts), 
      col = "forestgreen", 
      lwd = 1.5)
legend("topright",
       lty = c(1, 1),
       lwd = c(0.2, 1),
       cex = 0.8,
       col = c("steelblue", "forestgreen"),
       legend = c("Actual difference Change"))

demean_AWE_diff <- diff(UK_M_AWE_Q_xts)[-1] - GARCH_AWE@fit$coef[1]
# The demean_AWE_diff_Pcsd time series indicates the mean estimated by our model,
#     plus the conditional standard error: this is our model fit of the
#     conditional standard error at a given time, based on other values
demean_AWE_diff_Pcsd <- GARCH_AWE@sigma.t
# The demean_AWE_diff_Mcsd time series indicates the mean estimated by our model,
#     minus the conditional standard error: this is our model fit of the
#     conditional standard error at a given time, based on other values
demean_AWE_diff_Mcsd <- - GARCH_AWE@sigma.t
# Together, the plus/minus one conditional standard error gives a confidence
#     interval for the predicted volatility of the time series.

# Now, let' convert these into xts objects so we can use our usual tools:
demean_AWE_diff_xts <- xts(demean_AWE_diff, dates)
demean_AWE_diff_Pcsd_xts <- xts(demean_AWE_diff_Pcsd, dates)
demean_AWE_diff_Mcsd_xts <- xts(demean_AWE_diff_Mcsd, dates)

# plot deviation of percentage changes from mean
plot(as.zoo(demean_AWE_diff_xts), 
     type = "l", 
     col = "steelblue",
     ylab = "Pounds", 
     xlab = "Date",
     main = "ARMA(1,1)-GARCH(1,1) \n Estimated Bands of +- One Conditional Standard Deviation",
     lwd = 0.2)
# add horizontal line at y = 0
abline(0, 0)
# add GARCH confidence bands (one standard deviation) to the plot
lines(as.zoo(demean_AWE_diff_Pcsd_xts), 
      col = "darkred", 
      lwd = 0.5)
lines(as.zoo(demean_AWE_diff_Mcsd_xts), 
      col = "darkred", 
      lwd = 0.5)
legend("topright",
       lty = c(1, 1),
       lwd = c(0.5, 0.5),
       cex = 0.7,
       col = c("steelblue", "darkred"),
       legend = c("Actual Pounds Change", "+/- Conditional Standard Deviation"))

knitr::stitch_rhtml('Empirical Assignment.r')




