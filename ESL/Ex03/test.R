df <- read.csv("phoneme.csv")
#https://www.datascienceblog.net/post/machine-learning/linear-discriminant-analysis/
    
train <- grepl("^train", df$speaker)
# remove non-feature columns
to.exclude <- c("row.names", "speaker", "g")
feature.df <- df[, !colnames(df) %in% to.exclude]
test.set <- subset(feature.df, !train)
train.set <- subset(feature.df, train)
train.responses <- subset(df, train)$g
test.responses <- subset(df, !train)$g

lda.model <- lda(train.set, grouping = train.responses)

means <- colSums(lda.model$prior * lda.model$means)
train.mod <- scale(train.set, center = means, scale = FALSE) %*% 
    lda.model$scaling

lda.prediction.train <- predict(lda.model, train.set)

# visualize the features in the two LDA dimensions
plot.df <- data.frame(train.mod, "Outcome" = train.responses)
library(ggplot2)
ggplot(plot.df, aes(x = LD1, y = LD2, color = Outcome)) + geom_point()
