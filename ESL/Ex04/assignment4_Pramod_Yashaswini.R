require(glmnet )
require(Metrics)
require(boot)
set.seed(4)

# load data
raw_data <- read.table("prostate.txt")
#data_mean <- apply(raw_data[,1:8],2, mean)
#data_sd <- apply(raw_data[,1:8],2, sd)
#
#data_zero <- raw_data[,1:8] - data_mean
#raw_data[,1:8] <- data_zero/data_sd
raw_data[,1:8] <- apply(raw_data[,1:8],2, scale)

# 1---------------------------------------------------------------------------------
print("-------------------------Problem 4.1-------------------------")
# train-test split
train.data <- (raw_data[,1:9])[raw_data$train == "TRUE",]
test.data <- (raw_data[,1:9])[raw_data$train == "FALSE",]

# 2---------------------------------------------------------------------------------
print("-------------------------Problem 4.2-------------------------")
# Linear Regression model 

# LOOCV
LOOCV.fit <- glm(lpsa ~., data = train.data)
LOOCV.err <- cv.glm(train.data, LOOCV.fit)
print("LOOCV Avg MSE: "); print(LOOCV.err$delta[1])

# 5-Fold
CV_5fold.fit <- glm(lpsa ~., data = train.data)
CV_5fold.err <- cv.glm(train.data, CV_5fold.fit, K=5)
print("5-Fold CV Avg MSE: "); print(CV_5fold.err$delta[1])

# 10-Fold
CV_10fold.fit <- glm(lpsa ~., data = train.data)
CV_10fold.err <- cv.glm(train.data, CV_10fold.fit, K=10)
print("10-Fold CV Avg MSE: "); print(CV_10fold.err$delta[1])

#Complete train data
no_CV.fit <- lm(lpsa ~ ., data = train.data)
no_CV.err <- predict(no_CV.fit, test.data)
print("Full data Avg MSE: "); print(mean((no_CV.err-test.data$lpsa )^2))

# 3---------------------------------------------------------------------------------
print("-------------------------Problem 4.3-------------------------")
matrix_data <- model.matrix(lpsa ~ ., raw_data[,1:9])

train.feat <- matrix_data[raw_data$train,]
train.response <- train.data$lpsa

test.feat <- matrix_data[!raw_data$train,]
test.response <- test.data$lpsa
#RIDGE Regression
grid =10^seq(10 , -2, length = 100)

ridge.mod <- glmnet(train.feat, train.response, alpha = 0, lambda = grid)
ridge.pred <-predict(ridge.mod, s = 10, newx = test.feat)

plot(ridge.mod, xvar = "lambda", label="TRUE")
legend("bottomright", lwd = 1,col = 1:9, legend = colnames(train.feat), cex = .7)

# 4---------------------------------------------------------------------------------
print("-------------------------Problem 4.4-------------------------")

Ridge_CV_10fold.fit <- cv.glmnet(train.feat, train.response, alpha = 0, K=10)
plot(Ridge_CV_10fold.fit)

ridge_bestlam <- Ridge_CV_10fold.fit$lambda.min
Ridge_CV_10fold.train_pred <- predict(Ridge_CV_10fold.fit, s = ridge_bestlam, newx = train.feat)
Ridge_CV_10fold.test_pred <- predict(Ridge_CV_10fold.fit, s = ridge_bestlam, newx = test.feat)

print("Ridge Coeffs: "); print(coef(Ridge_CV_10fold.fit, Ridge_CV_10fold.fit$lambda.min))
print("Best Ridge lambda: "); print(ridge_bestlam)
print("Ridge 10fold-CV Train MSE: "); print(mean((Ridge_CV_10fold.train_pred -train.response)^2))
print("Ridge 10fold-CV Test MSE: "); print(mean((Ridge_CV_10fold.test_pred -test.response)^2))

# 5---------------------------------------------------------------------------------
print("-------------------------Problem 4.5-------------------------")
#LASSO Regression
grid =10^seq(10 , -2, length = 100)

lasso.mod <- glmnet(train.feat, train.response, alpha = 1, lambda = grid)
lasso.pred <-predict(lasso.mod, s = 10, newx = test.feat)

plot(lasso.mod, xvar = "lambda", label="TRUE")
legend("bottomright", lwd = 1,col = 1:9, legend = colnames(train.feat), cex = .7)

# 6---------------------------------------------------------------------------------
print("-------------------------Problem 4.6-------------------------")

Lasso_CV_10fold.fit <- cv.glmnet(train.feat, train.response, alpha = 1, K=10)
plot(Lasso_CV_10fold.fit)

lasso_bestlam <- Lasso_CV_10fold.fit$lambda.min
Lasso_CV_10fold.train_pred <- predict(Lasso_CV_10fold.fit, s = lasso_bestlam, newx = train.feat)
Lasso_CV_10fold.test_pred <- predict(Lasso_CV_10fold.fit, s = lasso_bestlam, newx = test.feat)

print("Lasso Coeffs: "); print(coef(Lasso_CV_10fold.fit, Lasso_CV_10fold.fit$lambda.min))
print("Best Lasso Lambda: ");print(lasso_bestlam)
print("Lasso 10fold-CV Train MSE: "); print(mean((Lasso_CV_10fold.train_pred -train.response)^2))
print("Lasso 10fold-CV Test MSE: "); print(mean((Lasso_CV_10fold.test_pred -test.response)^2))