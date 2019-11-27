#-------Authors----
# 1.  Pramod Ramesh Rao - 2581261
# 2. Yashaswini Mysuru Udaya Kumar - 2581353
# ------------------
#loading the package
require(dplyr)
require(MASS)
require(ggplot2)
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
train.model <- lda(train.data, grouping = train.responses)

train.predictions <- predict(train.model, train.data)
test.predictions <- predict(train.model, test.data)

#accuracy
print(mean(train.predictions$class==train.responses))
print(mean(test.predictions$class==test.responses))

# 3---------------------------------------------------------------------------------
print("-------------------------Problem 4.3-------------------------")
# projection of the training data onto the first two canonical coordinates of the LDA
plot.df <- data.frame(train.predictions$x , "Outcome" = train.responses)
ggplot(plot.df, aes(x = LD1, y = LD2, color = Outcome)) + geom_point()

# 4---------------------------------------------------------------------------------
print("-------------------------Problem 4.4-------------------------")

train.data <- data[data$index == "train",]
train.data_aa_ao <- train.data[data$g %in% c("aa", "ao"),]
#train.responses <- data[data$index == "train",]$g
#train.data <- dplyr::select(train.data,-c("index", "g"))

#test.data <- data[data$index == "test",]
#test.responses <- data[data$index == "test",]$g
#test.data <- dplyr::select(test.data,-c("index", "g"))