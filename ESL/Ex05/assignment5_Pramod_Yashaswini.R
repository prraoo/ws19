require(leaps)
require(Metrics)
require(pls)

# load data
load("prostate.Rdata")
# 1---------------------------------------------------------------------------------
print("-------------------------Problem 4.1-------------------------")
regfit.train <- regsubsets(lpsa~.,prostate.train)
regfit.summary <- summary(regfit.train)

train.mat = model.matrix(lpsa~., data = prostate.train) 
test.mat = model.matrix(lpsa~., data = prostate.test) 

par(mfrow=c(2,2))
plot(regfit.summary$rsq , xlab =" Number of Variables ", ylab =" RSq ",type ="l")
plot(regfit.summary$adjr2 , xlab =" Number of Variables ", ylab =" Adj RSq ",type ="l")
plot(regfit.summary$cp , xlab =" Number of Variables ", ylab =" CP ",type ="l")
plot(regfit.summary$bic , xlab =" Number of Variables ", ylab =" BIC ",type ="l")

print(regfit.summary)
coef.best <- coef(regfit.train,3)

pred.regfit_train <- train.mat[,names(coef.best)]%*%coef.best
cat("Regfit Train MSE loss: ", mse(prostate.train$lpsa, pred.regfit_train), "\n")
pred.regfit_test <- test.mat[,names(coef.best)]%*%coef.best
cat("Regfit Test MSE loss: ", mse(prostate.test$lpsa, pred.regfit_test), "\n")

# 2---------------------------------------------------------------------------------
print("-------------------------Problem 4.2-------------------------")
pcr.fit <- pcr(lpsa~., data = prostate.train, scale = TRUE)

for (comp in 1:8){
  cat("ncomp = ", comp, "\n")
  pred.pcr_train <- predict(pcr.fit, prostate.train, ncomp = comp)
  cat("PCR Train MSE loss is :", mse(prostate.train$lpsa, pred.pcr_train), "\n")
  pred.pcr_test <- predict(pcr.fit, prostate.test, ncomp = comp)
  cat("PCR Test MSE loss is :", mse(prostate.test$lpsa, pred.pcr_test), "\n")
}

# 3---------------------------------------------------------------------------------
print("-------------------------Problem 4.3-------------------------")
par(mfrow=c(1,1))
pls.fit <- plsr(lpsa~., data = prostate.train, scale = TRUE)

pls_tr_error_list <- c()
pls_te_error_list <- c()

for (comp in 1:8){
  cat("ncomp = ", comp, "\n")
  pred.pls_train <- predict(pls.fit, prostate.train, ncomp = comp)
  pls_tr_error_list <- c(pls_tr_error_list, mse(prostate.train$lpsa, pred.pls_train))
  cat("PLS Train MSE loss is :", mse(prostate.train$lpsa, pred.pls_train), "\n")
  pred.pls_test <- predict(pls.fit, prostate.test, ncomp = comp)
  pls_te_error_list <- c(pls_te_error_list, mse(prostate.test$lpsa, pred.pls_test))
  cat("PLS Test MSE loss is :", mse(prostate.test$lpsa, pred.pls_test), "\n")
}

# 3---------------------------------------------------------------------------------
print("-------------------------Problem 4.3-------------------------")
par(mfrow=c(1,1))

pca.train <- prcomp(prostate.train, scale = TRUE)
Cols= function (vec ){
  cols= rainbow(length(unique(vec)))
  return (cols[as.numeric(as.factor(vec))])
}

#plot(pca.train$x[,2:3])

biplot(pca.train)


