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
plot(regfit.summary$rsq , xlab =" Number of Variables ", ylab =" R-Sq ",type ="l")
points(which.max(regfit.summary$rsq), regfit.summary$rsq[which.max(regfit.summary$rsq)], col="red", pch=15)
title("R-Sq vs n-comp")

plot(regfit.summary$adjr2 , xlab =" Number of Variables ", ylab =" Adj R-Sq ",type ="l")
points(which.max(regfit.summary$adjr2), regfit.summary$adjr2[which.max(regfit.summary$adjr2)], col="red", pch=15)
title("Adj R-Sq vs n-comp")

plot(regfit.summary$cp , xlab =" Number of Variables ", ylab =" Cp ",type ="l")
points(which.min(regfit.summary$cp), regfit.summary$cp[which.min(regfit.summary$cp)], col="red", pch=15)
title("Cp vs n-comp")

plot(regfit.summary$bic , xlab =" Number of Variables ", ylab =" BIC ",type ="l")
points(which.min(regfit.summary$bic), regfit.summary$bic[which.min(regfit.summary$bic)], col="red", pch=15)
title("BIC vs n-comp")

print(regfit.summary)
coef.best <- coef(regfit.train,3)

pred.regfit_train <- train.mat[,names(coef.best)]%*%coef.best
cat("Regfit Train MSE loss: ", mse(prostate.train$lpsa, pred.regfit_train), "\n")
pred.regfit_test <- test.mat[,names(coef.best)]%*%coef.best
cat("Regfit Test MSE loss: ", mse(prostate.test$lpsa, pred.regfit_test), "\n")

# 2---------------------------------------------------------------------------------
print("-------------------------Problem 4.2-------------------------")

pcr.fit <- pcr(lpsa~., data = prostate.train, validation="CV")

pcr_tr_error_list <- c()
pcr_te_error_list <- c()

for (comp in 1:8){
  cat("ncomp = ", comp, "\n")
  
  pred.pcr_train <- predict(pcr.fit, prostate.train, ncomp = comp)
  pcr_tr_error_list <- c(pcr_tr_error_list, mse(prostate.train$lpsa, pred.pcr_train))
  cat("PCR Train MSE loss is :", mse(prostate.train$lpsa, pred.pcr_train), "\n")
  
  pred.pcr_test <- predict(pcr.fit, prostate.test, ncomp = comp)
  pcr_te_error_list <- c(pcr_te_error_list, mse(prostate.test$lpsa, pred.pcr_test))
  cat("PCR Test MSE loss is :", mse(prostate.test$lpsa, pred.pcr_test), "\n")
}

par(mfrow=c(1,1))

#plot PCR
plot(pcr_tr_error_list, col="green", ylim = c(0.2,1), pch = 16,
     panel.first = grid(), 
     xlab = "n-Components", ylab = "MSE Error", type = "b")
points(pcr_te_error_list, col = "blue", pch=16, type = "b")
title("PCR Train vs Test Error")
legend("topright", legend = c("Train", "Test"), col = c("green", "blue"), pch = 16)

# 3---------------------------------------------------------------------------------
print("-------------------------Problem 4.3-------------------------")

pls.fit <- plsr(lpsa~., data = prostate.train, validation="CV")

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

#plot PLS
plot(pls_tr_error_list, col="green", ylim = c(0.3,.7), pch = 16,
     panel.first = grid(), 
     xlab = "n-Components", ylab = "MSE Error", type = "b")
points(pls_te_error_list, col = "blue", pch=16, type = "b")
title("PLS Train vs Test Error")
legend("bottomright", legend = c("Train", "Test"), col = c("green", "blue"), pch = 16)

# 4---------------------------------------------------------------------------------
print("-------------------------Problem 4.4-------------------------")

par(mfrow=c(1,1))
# Train Data
pca.train <- prcomp(prostate.train)
pairs(pca.train$x[,1:4], main="PCA Plot - Train Data", 
      col = ifelse(prostate.train$lpsa >= 2.5, "red",  "green"), pch = 16)

# Complete data <- Train + Test 
prostate.data <- rbind(prostate.train, prostate.test)
pca.data <- prcomp(prostate.data)
pairs(pca.data$x[,1:4], main="PCA Plot - Whole Data", 
      col = ifelse(prostate.data$lpsa >= 2.5, "red",  "green"), pch = 16)

# 5---------------------------------------------------------------------------------
print("-------------------------Problem 4.5-------------------------")

par(mfrow=c(1,1))
# Train Data
pls.train <- plsr(lpsa~., data = prostate.train)
pairs(pls.train$scores[,1:4], main="PLS Plot - Train Data", 
      col = ifelse(prostate.train$lpsa >= 2.5, "blue",  "orange"), pch = 10)

# Complete data <- Train + Test 
prostate.data <- rbind(prostate.train, prostate.test)
pls.data <- plsr(lpsa~., data = prostate.data)
pairs(pls.data$scores[,1:4], main="PLS Plot - Whole data", 
      col = ifelse(prostate.data$lpsa >= 2.5, "blue",  "orange"), pch = 10)
