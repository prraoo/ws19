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