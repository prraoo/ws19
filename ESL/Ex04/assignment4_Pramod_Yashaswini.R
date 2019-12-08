# load data
raw_data <- read.table("prostate.txt")
# 1---------------------------------------------------------------------------------
print("-------------------------Problem 4.1-------------------------")
# train-test split
train <- raw_data[raw_data$train == "TRUE",]
train.data <- (dplyr::select(train, -c("train", "lpsa")))
train.data <- (train.data - colMeans(train.data))/(apply(train.data, 2, sd))
train.labels <- (train$lpsa - mean(train$lpsa))/(sd(train$lpsa))

test <- raw_data[raw_data$train == "FALSE",]
test.data <- dplyr::select(test, -c("train", "lpsa"))
test.data <- (test.data - colMeans(test))/apply(test.data, 2, sd)
test.labels <- (test$lpsa - mean(test$lpsa))/sd(test$lpsa)

# 2---------------------------------------------------------------------------------
print("-------------------------Problem 4.2-------------------------")
