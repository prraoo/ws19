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
require("PerformanceAnalytics")
chart.Correlation(ozone, histogram = FALSE)

# 5---------------------------------------------------------------------------------
print("Problem 4.5")
rss <- function(true_value, predicted_value){
    return(sum((true_value - predicted_value)^2))
}
print(rss(testset,testset+1))

# 6---------------------------------------------------------------------------------
print("Problem 4.6")
train_data <- ozone[trainset,]
test_data <- ozone[testset,]

ozone_model <- lm(ozone ~ radiation+temperature+wind, data = train_data)
print(summary(ozone_model))
pred_data <- predict.lm(ozone_model, test_data)

test_lm_rss <- rss(test_data["ozone"], pred_data)
test_lm_cor <- cor(test_data["ozone"], pred_data, method = "pearson")
plot(test_data$ozone, pred_data, main = "Predicted and True Values Scatterplot",col = c("Blue"), pch=16)

# 7---------------------------------------------------------------------------------
print("Problem 4.7")
require("FNN")
knn_model <- function(train, test, n){
   knn_pred <- knn.reg(train = train[2:4],test =  test[2:4],y = train[1],k = n)
   knn_rss = rss(test$ozone, knn_pred$pred)
   return(list(pred=knn_pred$pred, rss=knn_rss))
}

# Train Data
i = 1
train_knn_rss <- c()
knn_predictions <- knn_model(train_data, train_data, i)
train_knn_rss <- c(train_knn_rss, knn_predictions$rss)

while(i<30){
    knn_predictions <- knn_model(train_data, train_data, i+2)
    train_knn_rss <- c(train_knn_rss, knn_predictions$rss)
    i <- i+1
}
grid(1,1)
plot(train_knn_rss, col="green", ylim = c(0,70000), main = "Train vs Pred RSS", pch=16, panel.first = grid())

# Test Data
i = 1
test_knn_rss <- c()
knn_predictions <- knn_model(train_data, test_data, i)
test_knn_rss <- c(test_knn_rss, knn_predictions$rss)

while(i<30){
    knn_predictions <- knn_model(train_data, test_data, i+2)
    test_knn_rss <- c(test_knn_rss, knn_predictions$rss) 
    i <- i+1
}
points(test_knn_rss, col="blue", pch=16)

# Assumptions in KNN:
# https://saravananthirumuruganathan.wordpress.com/2010/05/17/a-detailed-introduction-to-k-nearest-neighbor-knn-algorithm/
# 8---------------------------------------------------------------------------------
print("Problem 4.8")
# Chosen nearest neighbour = 6
knn_predictions_best <- knn_model(train_data, test_data, 6)
test_knn_cor_best <- cor(knn_predictions_best$pred, test_data$ozone, method = "pearson")
cat("\n The RSS value of best kNN model: ", knn_predictions_best$rss,"\n Pearson Coeff for best kNN model is: ", test_knn_cor_best, "\n")
cat("\n The RSS value of lm  model: ", test_lm_rss,"\n Pearson Coeff for lm model is: ", test_lm_cor[1], "\n")

dev.off()