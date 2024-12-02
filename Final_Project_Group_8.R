#create the function

bank <- read.csv('bank.csv', sep = ';')
bank

#must identify x columns and y variable before using function
x <- bank[,1:16]
bank$y <- as.factor(bank[,17])
y <- as.factor(bank[,17])

logistic_regression <- function(x, y, B = 20, conf_level = 0.95){
  x_matrix <- model.matrix(~., data = data.frame(x))
  y_matrix <- ifelse(y == unique(y)[1], 0, 1)
  
  betaval <- solve(t(x_matrix)%*%x_matrix)%*%t(x_matrix)%*%y_matrix
  
  beta_function <- function(x_matrix, y_matrix, betaval){
    p <- 1/(1+exp((-x_matrix)%*%betaval))
    sum(-y_matrix %*% log(p)-(1-y_matrix)%*%log(1-p))
  }
  
  coef_ests <- optim(par = betaval, fn = beta_function, x_matrix = x_matrix, y_matrix = y_matrix)
  
  boot_samp_mat <- matrix(NA, nrow = ncol(x_matrix), ncol = B)
  
  rownames(boot_samp_mat) <- colnames(x_matrix)
  
  for(i in 1:B){
    bootrows <- sample(nrow(x_matrix), nrow(x_matrix), replace = TRUE)
    bootsamp_x <- x_matrix[bootrows,]
    bootsamp_y <- y_matrix[bootrows]
    
    betaval_boot <- solve(t(bootsamp_x)%*%bootsamp_x)%*%t(bootsamp_x)%*%bootsamp_y
    
    boot_samp_mat[,i] <- optim(par = betaval_boot, fn = beta_function, x_matrix = bootsamp_x, y_matrix = bootsamp_y)$par
  }
  
  bootstrap_intervals <- apply(boot_samp_mat, 1, quantile, probs = c((1-conf_level)/2, 1-(1-conf_level)/2))
  
  log_odds <- x_matrix %*% coef_ests$par
  
  prob_ests <- exp(log_odds)/(1+exp(log_odds))
  
  predictions <- ifelse(prob_ests > 0.5, 1, 0)
  
  confusion <- table(y_matrix, predictions)
  
  confusion_metrics <- data.frame(
    prevalence = sum(confusion[2,])/sum(confusion),
    accuracy = sum(diag(confusion))/sum(confusion),
    sensitivity = confusion[2,2]/sum(confusion[2,]),
    specificity = confusion[1,1]/sum(confusion[1,]),
    false_discovery_rate = confusion[1,2]/sum(confusion[,2]),
    diagnostic_odds_ratio = (confusion[2,2]/confusion[2,1])/(confusion[1,2]/confusion[1,1]))
  
  list(coefficient_estimates = coef_ests$par,
       bootstrap_intervals = bootstrap_intervals,
       confusion_matrix = confusion,
       confusion_metrics = confusion_metrics)
  
}



logistic_regression(x, y)

run_function <- logistic_regression(x,y)

run_function$



model.matrix(~., data = data.frame(x))





#test
x <- matrix(data = rnorm(100, 2, 3), nrow = 20, ncol = 5)
y <- sample(c(0,1), size = 20, replace = TRUE, prob = c(0.5, 0.5))

x_matrix <- cbind(rep(1, 20), x)
y_matrix <- y

betaval <- solve(t(x_matrix)%*%x_matrix)%*%t(x_matrix)%*%y_matrix

optim(par = betaval, fn = beta_function, x_matrix = x_matrix, y_matrix = y_matrix)

summary(glm(y~x[,1]+x[,2]+x[,3]+x[,4]+x[,5], family = binomial))