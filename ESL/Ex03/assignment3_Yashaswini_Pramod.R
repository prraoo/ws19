#-------Authors----
# 1.  Pramod Ramesh Rao - 2581261
# 2. Yashaswini Mysuru Udaya Kumar - 2581353
# ------------------

#loading the package
require(dplyr)
require(MASS)
require(ggplot2)
require(caret)
require(RColorBrewer)

# Dataset
raw_data <- read.csv("phoneme.csv", header = TRUE,sep = ",")

# 1---------------------------------------------------------------------------------
print("-------------------------Problem 4.1-------------------------")
raw_data <- raw_data %>% mutate(index=gsub("\\..*","",speaker))
data <- dplyr::select(raw_data, -c("speaker", "row.names"))

train.data <- data[data$index == "train",]
train.responses <- data[data$index == "train",]$g
train.data <- dplyr::select(train.data,-c("index", "g"))

test.data <- data[data$index == "test",]
test.responses <- data[data$index == "test",]$g
test.data <- dplyr::select(test.data,-c("index", "g"))

# 2---------------------------------------------------------------------------------
print("-------------------------Problem 4.2-------------------------")
#LDA MODEL
train.lda_model <- lda(train.data, grouping = train.responses)

train.predictions <- predict(train.lda_model, train.data)
test.predictions <- predict(train.lda_model, test.data)

#accuracy
cat("The train error is: \t", 1 - mean(train.predictions$class==train.responses))
cat("\nThe test error is: \t", 1 - mean(test.predictions$class==test.responses))
cat("\n")

# 3---------------------------------------------------------------------------------
print("-------------------------Problem 4.3-------------------------")
# projection of the training data onto the first two canonical coordinates of the LDA
plot.df <- data.frame(train.predictions$x , "Outcome" = train.responses)
ggplot(plot.df, aes(x = LD1, y = LD2, color = Outcome)) + geom_point()

#TODO Plot various color plotS
colors <- brewer.pal(8, "Accent")
my.cols <- colors[match(train.predictions$class, levels(data$g))]
plot(train.lda_model, dimen = 4, col = my.cols)
# 4---------------------------------------------------------------------------------
print("-------------------------Problem 4.4-------------------------")
#Use only "aa" and "ao" samples for LDA

train.data <- data[data$index == "train",]
train.data_aa_ao <- train.data[train.data$g %in% c("aa", "ao"),]

train.responses_aa_ao <- train.data_aa_ao[train.data_aa_ao$index == "train",]$g
train.responses_aa_ao <- droplevels(train.responses_aa_ao)

train.data_aa_ao <- dplyr::select(train.data_aa_ao,-c("index", "g"))

test.data <- data[data$index == "test",]
test.data_aa_ao <- test.data[test.data$g %in% c("aa", "ao"),]

test.responses_aa_ao <- test.data_aa_ao[test.data_aa_ao$index == "test",]$g
test.responses_aa_ao <- droplevels(test.responses_aa_ao)

test.data_aa_ao <- dplyr::select(test.data_aa_ao,-c("index", "g"))

#fit the LDA MODEL with new data
train.lda_model_aa_ao <- lda(train.data_aa_ao, grouping = train.responses_aa_ao)

train.lda_predictions_aa_ao <- predict(train.lda_model_aa_ao, train.data_aa_ao)
test.lda_predictions_aa_ao <- predict(train.lda_model_aa_ao, test.data_aa_ao)

#accuracy
cat("The train error is: \t", 1 - mean(train.lda_predictions_aa_ao$class==train.responses_aa_ao))
cat("\nThe test error is: \t", 1 - mean(test.lda_predictions_aa_ao$class==test.responses_aa_ao))
cat("\n")


# 5---------------------------------------------------------------------------------
print("-------------------------Problem 4.5-------------------------")
#QDA MODEL
train.data <- data[data$index == "train",]
train.responses <- data[data$index == "train",]$g
train.data <- dplyr::select(train.data,-c("index", "g"))

test.data <- data[data$index == "test",]
test.responses <- data[data$index == "test",]$g
test.data <- dplyr::select(test.data,-c("index", "g"))
train.qda_model <- qda(train.data, grouping = train.responses)

train.predictions <- predict(train.qda_model, train.data)
test.predictions <- predict(train.qda_model, test.data)

#accuracy
cat("The train error is: \t", 1 - mean(train.predictions$class==train.responses))
cat("\nThe test error is: \t", 1 - mean(test.predictions$class==test.responses))
cat("\n")

#Use only "aa" and "ao" samples for QDA
train.data <- data[data$index == "train",]
train.data_aa_ao <- train.data[train.data$g %in% c("aa", "ao"),]

train.responses_aa_ao <- train.data_aa_ao[train.data_aa_ao$index == "train",]$g
train.responses_aa_ao <- droplevels(train.responses_aa_ao)

train.data_aa_ao <- dplyr::select(train.data_aa_ao,-c("index", "g"))

test.data <- data[data$index == "test",]
test.data_aa_ao <- test.data[test.data$g %in% c("aa", "ao"),]

test.responses_aa_ao <- test.data_aa_ao[test.data_aa_ao$index == "test",]$g
test.responses_aa_ao <- droplevels(test.responses_aa_ao)

test.data_aa_ao <- dplyr::select(test.data_aa_ao,-c("index", "g"))

#fit the QDA MODEL with new data
train.qda_model_aa_ao <- qda(train.data_aa_ao, grouping = train.responses_aa_ao)

train.qda_predictions_aa_ao <- predict(train.qda_model_aa_ao, train.data_aa_ao)
test.qda_predictions_aa_ao <- predict(train.qda_model_aa_ao, test.data_aa_ao)

#accuracy
cat("The train error is: \t", 1 - mean(train.qda_predictions_aa_ao$class==train.responses_aa_ao))
cat("\nThe test error is: \t", 1 - mean(test.qda_predictions_aa_ao$class==test.responses_aa_ao))
cat("\n")

# 6---------------------------------------------------------------------------------
print("-------------------------Problem 4.6-------------------------")
#Confusion Matrices
lda_CM <- confusionMatrix(test.lda_predictions_aa_ao$class, test.responses_aa_ao)
qda_CM <- confusionMatrix(test.qda_predictions_aa_ao$class, test.responses_aa_ao)
cat("\n The LDA Confusion matrix is:\n")
print(lda_CM$table)
cat("\n The QDA Confusion Matrix is:\n")
print(qda_CM$table)
