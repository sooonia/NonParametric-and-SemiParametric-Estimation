library(np)
library(sm)
library(lubridate)
source('Helper2.R')
rd <- read.csv('df.csv')
colnames(rd)[1] <- 'idx'

seasons <- c('winter', 'spring', 'summer', 'fall', 'winter', 'spring', 'summer', 'fall', 'winter')
rd$round_date <- as.factor(round_date(mdy(rd$date), "season"))
rd$season <- "none"
for(i in seq(1:length(levels(rd$round_date)))){
  s <- levels(rd$round_date)[i]
  rd[rd$round_date == s,'season'] <- seasons[i]
}
rd$season <- as.factor(rd$season)

# Question 1
bw <- Silverman(rd$SO2_level)$bw

plot(density(rd$SO2_level, 1.5^-2 * bw), col = 'purple', main="S02 Density Estimation")
lines(density(rd$SO2_level, 1.5^-1 * bw), col = 'blue')
lines(density(rd$SO2_level, bw), col = 'black')
lines(density(rd$SO2_level, 1.5 * bw), col = 'orange')
lines(density(rd$SO2_level, 1.5^2 * bw), col = 'red')
bws <- round(c(1.5^-2 * bw, 1.5^-1 * bw, bw, 1.5^2 * bw, 1.5 * bw), digits = 2)
legend(75, .05, legend = bws, col = c('purple', 'blue', 'black', 'red', 'orange'), 
       lty = 1, title = "BWs")


#SEASONS
bw <- Silverman(rd[rd$season == 'winter', 'SO2_level'])$bw
plot(density(rd[rd$season == 'winter', 'SO2_level'], bw, kernel = 'ep'), 
      col = 'blue', main = 'Seasonal Subsets', xlab = 'S02 Level')
bw <- Silverman(rd[rd$season == 'spring', 'SO2_level'])$bw
lines(density(rd[rd$season == 'spring', 'SO2_level'], bw, kernel = 'ep'), 
      col = 'green4', main = 'Seasonal Subsets')
bw <- Silverman(rd[rd$season == 'fall', 'SO2_level'])$bw
lines(density(rd[rd$season == 'fall', 'SO2_level'], bw, kernel = 'ep'), 
     col = 'coral', main = 'Seasonal Subsets')
bw <- Silverman(rd[rd$season == 'summer', 'SO2_level'])$bw
lines(density(rd[rd$season == 'summer', 'SO2_level'], bw, kernel = 'ep'), 
     col = 'gold', main = 'Seasonal Subsets')
legend(50, .06, legend=c("Winter", "Spring", 'Summer', 'Fall'),
       col=c("blue", "green4", 'gold', 'coral'), lty=1)

ks.test(rd[rd$season == 'winter', 'SO2_level'], rd[rd$season == 'spring', 'SO2_level'])
ks.test(rd[rd$season == 'winter', 'SO2_level'], rd[rd$season == 'fall', 'SO2_level'])
ks.test(rd[rd$season == 'winter', 'SO2_level'], rd[rd$season == 'summer', 'SO2_level'])
ks.test(rd[rd$season == 'spring', 'SO2_level'], rd[rd$season == 'fall', 'SO2_level'])
ks.test(rd[rd$season == 'spring', 'SO2_level'], rd[rd$season == 'summer', 'SO2_level'])
ks.test(rd[rd$season == 'fall', 'SO2_level'], rd[rd$season == 'summer', 'SO2_level'])

# Question 2
loc.lin <- loess(log_hosp ~ SO2_level, data=rd)
lin <- glm(log_hosp ~ SO2_level, data=rd)

plot(rd$SO2_level, rd$log_hosp, xlab = "SO2", ylab = "log(hospital admits)", cex = .1,
     main = 'Linear vs NP Model')
lines(seq(2,100,.5), predict(loc.lin, seq(2,100,.5)), lty=1, col='blue')
lines(rd$SO2_level, fitted(lin), lty = 1, col = "red")
legend(60, 3, legend=c("GLM", "Local Lin."),
       col=c("red", "blue"), lty=1)

plot(rd$SO2_level, rd$log_hosp, xlab = "SO2", ylab = "log(hospital admits)", cex = .1,
     main = 'Zoomed: Linear vs NP Model', ylim=c(4.5, 6))
lines(seq(2,100,.5), predict(loc.lin, seq(2,100,.5)), lty=1, col='blue')
lines(rd$SO2_level, fitted(lin), lty = 1, cex=.1, col = "red")
legend(70, 5.1, legend=c("GLM", "Local Lin."),
       col=c("red", "blue"), lty=1)
plot(rd$SO2_level, loc.lin$residuals, main = 'NP Model Residuals')
plot(rd$SO2_level, lin$residuals, main = "Linear Model Residuals")

preds <- predict(loc.lin, newdata= rd)
print(SSE(preds,rd$log_hosp))


# Question 3
s <- sm.poisson(x= rd$SO2_level, y= rd$hospital_admissions, h= 30, col='red')
lin <- glm(hospital_admissions ~ SO2_level, family = poisson, data=rd)
lines(rd$SO2_level, fitted(lin), lty = 1, cex=.1, col = "red")

preds <- predict(s, newdata= rd)
print(SSE(preds,rd$log_hosp))

# Question 4
rd$high <- 0
rd[rd$hospital_admissions >= 300, 'high'] <- 1
f= formula('high ~ SO2_level')
loc.lin <- loess(formula = f, data = rd)
plot(rd$SO2_level, rd$high, xlab = "SO2", ylab = "log(hospital admits)", cex = .5, 
     main="Predicting High (>=300) Hospital Admits")
plot(loc.lin$residuals)


# Question 5
np <- sm.binomial(x=rd$SO2_level, y=rd$high, h=10)
lin <- glm(high ~ SO2_level, family=binomial(link = "logit"), rd)
lines(seq(2,100,.5), predict(loc.lin, seq(2,100,.5)), lty=1, col='blue')
lines(rd$SO2_level, fitted(lin), lty = 1, cex=.1, col = "red")

legend(65, .9, legend=c("Logistic NP", "Local Lin.", "Linear Logistic"),
       col=c("black", "blue", "red"), lty=1)


