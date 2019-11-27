#-------Authors----
# 1.  Pramod Ramesh Rao - 2581261
# 2. Yashaswini Mysuru Udaya Kumar - 2581353
# ------------------
#loading the package
require(dplyr)
require(MASS)
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

