pdf("P1_Programming_Solutions.pdf") 
load("ozone.RData")
cat("\n The objects in the dataset are:", ls(), "\n")

print("Structure of the ozone object is:")
str(ozone)
print("Structure of the trainset is:")
str(trainset)
print("Structure of the testset is:")
str(testset)
print("Ozone Summary")
print(summary(ozone))

cat("\n", dim(ozone), dim(trainset), dim(testset), "\n")
#NOTE dim of values is NULL

cat("\n Length: ",length(ozone), length(trainset), length(testset), "\n")
cat("\n Range: ",range(ozone), range(trainset), range(testset), "\n")
print(colnames(ozone))

# 2---------------------------------------------------------------------------------
#How many observations do you have (in total, in the training set, in the testset)
#Answer: ozone: 111; trainset: 80; testset: 31; total: 111

# 3---------------------------------------------------------------------------------
print("Problem 4.3")
#Range
print(apply(ozone, 2, range))
#Mean
print(apply(ozone, 2, mean))
#Standard Deviation
print(apply(ozone, 2, sd))

# 4---------------------------------------------------------------------------------
print("Problem 4.4")
# Create scatterplots for every pair of features in the dataset
pairs(ozone[1:4], main="Ozone Data Scatterplots", col = c("red"))

# Calculate the Pearson correlation
print(cor(ozone[1:4], method = "pearson" ))

# Coefficients for each pair of datatypes. In general, what is the range of the Pearson correlation coefficient?
#Answer: Range for pearson coefficient is +1 to -1

# What does a correlation coefficient of 0 tell you about the relationship between two variables? 
# Answer: A correlation of zero indicates that the variables do not have a linear association between them.

# What trends do you observe in the data according to the correlation coefficient? 
# There are postive as well as negative trends among the different variables. For instance,
# a) ozone has a strong positive (0.7) correlation with temperature, moderate correlation with radiation and 
# negatively correlated with wind.
# b) The second variable, radiation, postive correlation with temperature and a small negative correlation with wind.
# c) The third variable, temperature is negatively correlated with wind.

# Can you see them directly from the plot (visually)?
library("PerformanceAnalytics")
chart.Correlation(ozone, histogram = FALSE)

# 5---------------------------------------------------------------------------------
print("Problem 4.5")
rss <- function(true_value, predicted_value){
    print("Calculating RSS ...")
    return(sum((true_value - predicted_value)^2))
}
print(rss(testset,testset+1))

# 6---------------------------------------------------------------------------------
print("Problem 4.6")
train_data <- ozone[trainset,]
test_data <- ozone[testset,]

ozone_model <- lm(ozone ~ radiation+temperature+wind, data = train_data, graph = TRUE)
print(summary(ozone_model))
pred_data <- predict.lm(ozone_model, test_data)

dev.off()