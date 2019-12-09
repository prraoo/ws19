require(caret)
#require()

# load data
raw_data <- read.table("prostate.txt")
# 1---------------------------------------------------------------------------------
print("-------------------------Problem 4.1-------------------------")
# train-test split
train <- raw_data[raw_data$train == "TRUE",]
train <- (dplyr::select(train, -c("train")))
train.data <- data.frame(apply(train, 2, scale))

test <- raw_data[raw_data$train == "FALSE",]
test <- dplyr::select(test, -c("train"))
test.data <- data.frame(apply(test, 2, scale))

# 2---------------------------------------------------------------------------------
print("-------------------------Problem 4.2-------------------------")
# Linear Regression model
#train.LOOCV <- lm(lpsa ~ .,data = train.data)
# LOOCV
train.control = trainControl(method = "LOOCV")
train.LOOCV <- caret::train(lpsa ~ ., data= train.data, method = "lm", trControl=train.control)
print("LOOCV Results: "); print(train.LOOCV$results)

# 5-Fold
train.control = trainControl(method = "cv", number = 5)
train.5fold_CV <- caret::train(lpsa ~ ., data= train.data, method = "lm", trControl=train.control)
print("5-Fold CV Results: "); print(train.5fold_CV$results)
# 10-Fold
train.control = trainControl(method = "cv", number = 10)
train.10fold_CV <- caret::train(lpsa ~ ., data= train.data, method = "lm", trControl=train.control)
print("10-Fold CV Results: "); print(train.10fold_CV$results)
#Full train data
train.no_CV <- caret::train(lpsa ~ ., data= train.data, method = "lm")
print("No CV Results: "); print(train.no_CV$results)
# Complete Train