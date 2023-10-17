##Load Data
data <- read.table("https://hastie.su.domains/ElemStatLearn/datasets/prostate.data")
data <- as.matrix(data)

##Split to Train and Test Sets
train = data[data[,ncol(data)]==1,]
test = data[data[,ncol(data)]==0,]

##Split into X , Y Train and test
X_train <- train[, c(1:8)]
Y_train <- train[,ncol(train)-1]

X_test <- test[, c(1:8)]
Y_test = test[, (ncol(data)-1)]

###### Estimator's Functions ######

##OLS Estimator
OLS <- function(X , Y){
  return(solve(t(X)%*%X)%*%t(X)%*%Y);
}

##Ridge Estimator
Ridge <- function(X,Y,lam){
  return (solve((t(X)%*%X)+(lam*diag(ncol(X))))%*%t(X)%*%Y);
}

##Principal Components Estimator
BCP <- function(X, Y, K) {
  X <- svd(X)
  return(X$v[, 1:K] %*% solve(diag(X$d)[1:K, 1:K]) %*% t(X$u[, 1:K]) %*% Y)
}

##Lasso Estimator
install.packages("quadprog")
library(quadprog)
LASSO<- function(X, Y, t){
  p <- ncol(X)
  D <- 2*t(X)%*%X
  d <- 2*t(X)%*%Y
  R <- as.matrix(expand.grid(rep(list(c(-1,1)),p)))
  b <- rep(-t,2^p)
  Beta_LASSO <- solve.QP(Dmat=D,dvec=d,Amat=t(R),bvec=b)
  return(Beta_LASSO$solution)
}

##Elastic Net Estimator 
EN<-function(X,Y,lambda1,lambda2){
  beta1<-rep(0,dim(X)[2])
  beta2<-rep(1,dim(X)[2])
  while(sqrt(t(beta1-beta2)%*%(beta1-beta2))>0.0000001){
    
    if(dim(X)[2]==2)
    { beta2<-beta1
    for(j in 1:2)
    {Rj<-t(X[,j])%*%(Y-X[,-j]*beta1[-j])
    betaj<-(1/(1+(lambda2/(t(X[,j])%*%X[,j]))))*Rj*max(1/(t(X[,j])%*%X[,j])-lambda1/(2*abs(Rj)*(t(X[,j])%*%X[,j])),0)
    beta1[j]<-betaj}
    }
    else
    {  beta2<-beta1
    for(j in 1:dim(X)[2])
    {Rj<-t(X[,j])%*%(Y-X[,-j]%*%(beta1[-j]))
    betaj<-(1/(1+(lambda2/(t(X[,j])%*%X[,j]))))*Rj*max(1/(t(X[,j])%*%X[,j])-lambda1/(2*abs(Rj)*(t(X[,j])%*%X[,j])),0)
    beta1[j]<-betaj}}
    
  }
  return(beta1)
}

###### MSE Function ######
MSE<- function(X, Y, Beta){
  Y_hat <- (X%*%Beta)
  return(mean((Y - Y_hat)^2))
}

###### Tuning Parameters selection : Leave One Out Cross Validation Algorithm (LOOCV) ######

##For Ridge
LOOCV1 <- function(X, Y, lambdas , method) {
  n <- nrow(X)
  erreur_lambda <- Inf
  lambda_opti <- NULL
  
  for (lambda in lambdas) {
    erreur_lambda_new <- 0
    
    for (i in 1:n) {
      training_set_X <- X[-i, ]
      test_X <- X[i, ]
      training_set_Y <- Y[-i]
      test_Y <- Y[i]
      
      Beta <- method(training_set_X, training_set_Y, lambda)
      erreur_i <- abs(test_Y - test_X %*% Beta)
      erreur_lambda_new <- erreur_lambda_new + erreur_i
    }
    
    erreur_lambda_new <- erreur_lambda_new / n
    
    if (erreur_lambda_new < erreur_lambda) {
      erreur_lambda <- erreur_lambda_new
      lambda_opti <- lambda
    }
  }
  
  return(lambda_opti)
}

##For 2 Parameter (ELASTIC NET)

LOOCV2 <- function(X, Y, lambdas, alphas) {
  n <- nrow(X)
  erreur_alpha_lambda <- Inf
  alpha_opti <- NULL
  lambda_opti <- NULL
  
  for (alpha in alphas) {
    for (lambda in lambdas) {
      erreur_alpha_lambda_new <- 0
      
      for (i in 1:n) {
        training_set_X <- X[-i, ]
        test_X <- X[i, ]
        training_set_Y <- Y[-i]
        test_Y <- Y[i]
        
        Beta_EN <- EN(training_set_X, training_set_Y, alpha, lambda)
        erreur_i <- abs(test_Y - test_X %*% Beta_EN)
        erreur_alpha_lambda_new <- erreur_alpha_lambda_new + erreur_i
      }
      
      erreur_alpha_lambda_new <- erreur_alpha_lambda_new / n
      
      if (erreur_alpha_lambda_new < erreur_alpha_lambda) {
        erreur_alpha_lambda <- erreur_alpha_lambda_new
        alpha_opti <- alpha
        lambda_opti <- lambda
      }
    }
  }
  
  return(list(alpha_opti = alpha_opti, lambda_opti = lambda_opti))
}

############# COPARAISON ###############
#I. Estimator Calculators
#Ridge
lambdas <- seq(0.1 , 10 , 0.1)
lam_Ridge <- LOOCV1(X_train , Y_train , lambdas , Ridge)
Beta_Ridge <- Ridge(X_train , Y_train , lam_Ridge) 

#OLS
Beta_OLS <- OLS(X_train , Y_train) 

#Principal Components 
K <- seq(1 , 8 , 1)
K_opt <- LOOCV1(X_train , Y_train , K , BCP)
Beta_CP <- BCP(X_train , Y_train , K_opt) 

#LASSO
lambdas <- seq(0.1 , 10 , 0.1)
lam_LASSO <- LOOCV1(X_train , Y_train , lambdas , LASSO)
Beta_LASSO <- LASSO(X_train , Y_train , lam_LASSO) 

#EN
lambdas1 <- seq(0.1 , 10 , 0.1)
lambdas2 <- seq(0.1 , 10 , 0.1)
#lam_EN <- LOOCV2(X_train , Y_train , lambdas1 ,lambdas2 )
Beta_EN <- EN(X_train , Y_train , 0.00001 , 0.00001)
# II. MSE Calculators 

#Ridge 
MSE_Ridge <- MSE(X_test , Y_test , Beta_Ridge)

#OLS 
MSE_OLS <- MSE(X_test , Y_test , Beta_OLS)

#ACP 
MSE_ACP <- MSE(X_test , Y_test , Beta_CP)

#LASSO 
MSE_LASSO <- MSE(X_test , Y_test , Beta_LASSO)

#ELastic Net 
MSE_EN <- MSE(X_test , Y_test , Beta_EN)



