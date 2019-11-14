#loading the package
require(PerformanceAnalytics)
require(ISLR)
# Dataset
data <- Auto
# Display the dataset
cat("\n Auto dataset: \n")
print(head(data))
# 1---------------------------------------------------------------------------------
print("-------------------------Problem 4.1-------------------------")

# Create scatter plot of all variables
svg(filename = "figure1.svg")
pairs(data[1:8], main="Auto Dataset Scatterplots", col = c("red"))
dev.off()
# 2---------------------------------------------------------------------------------
print("-------------------------Problem 4.2-------------------------")

# Correlation between variables
cat("\n Correlation between the variables: \n")
print(cor(data[1:8]))
svg(filename = "figure2.svg")
chart.Correlation(data[1:8], histogram = FALSE)
dev.off()

# 3---------------------------------------------------------------------------------
print("-------------------------Problem 4.3-------------------------")

# Simple linear regression
#mpg as the response using variable cylinder
linearMod_cylinder <- lm(mpg ~ cylinders, data=Auto)
print(linearMod_cylinder)
print(summary(linearMod_cylinder))

#mpg as the response using variable displacement
linearMod_displacement <- lm(mpg ~ displacement, data=Auto)
print(linearMod_displacement)
print(summary(linearMod_displacement))

#mpg as the response using variable horsepower
linearMod_horespower <- lm(mpg ~ horsepower, data=Auto)
print(linearMod_horespower)
print(summary(linearMod_horespower))
      
#mpg as the response using variable year
linearMod_year <- lm(mpg ~ year, data=Auto)
print(linearMod_year)
print(summary(linearMod_year))

# 4---------------------------------------------------------------------------------
print("-------------------------Problem 4.4-------------------------")

#Mulitple regression model with mpg as response to all variables except names
mulregModel <- lm(mpg~.-name, data=Auto)
print(mulregModel)
print(summary(mulregModel))

# 5---------------------------------------------------------------------------------
print("-------------------------Problem 4.5-------------------------")

par(mfrow=c(2,2))
svg(filename = "figure3.svg")
plot(mulregModel)
dev.off()

# 6---------------------------------------------------------------------------------
print("-------------------------Problem 4.6-------------------------")

##Pairwise interaction terms (X1X2) for cylinders, weight, and year
lm_cw <- lm(mpg~cylinders+weight+cylinders*weight, data=Auto)
print(lm_cw)
print(summary(lm_cw))

lm_wy <- lm(mpg~weight+year+weight*year, data=Auto)
print(lm_wy)
print(summary(lm_wy))

lm_yc <- lm(mpg~year+cylinders+year*cylinders, data=Auto)
print(lm_yc)
print(summary(lm_yc))

# Non-linear transformations log (X), sqrt(X), X^2 for the displacement variable
lm4 <- lm(mpg~sqrt(displacement), data=Auto)
print(lm4)
print(summary(lm4))

lm5 <- lm(mpg~(log(displacement)), data=Auto)
print(lm5)
print(summary(lm5))

lm6 <- lm(mpg~(displacement^2), data=Auto)
print(lm6)
print(summary(lm6))