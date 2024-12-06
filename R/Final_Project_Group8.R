#' @title Logistic Regression for Binary Classification
#'
#' @description Calculate the estimated logistic regression coefficients for a binary
#' classification problem. The function also produces bootstrap confidence intervals
#' for each estimated coefficient. Additionally, the function will use the estimated
#' coefficients to produce predictions for the class for each data entry.
#' These predictions will be used to create a confusion matrix that compares the
#' estimated classes to the true classes.
#' @param x A \code{matrix} or \code{data frame} of the x-variables (independent variables)
#' that will be used in your logistic regression.
#' @param y A \code{vector} containing your y (dependent) variable. This should be a
#' binary variable (two categories), and it must be a a \code{factor} to use in this
#' function.
#' @param B A \code{numeric} (integer) used to denote the number of bootstrap
#' simulations to perform in order to compute bootstrap confidence intervals for
#' the estimated coefficients.
#' @param conf_level A \code{numeric} used to specify the confidence level for
#' the bootstrap confidence intervals.
#' @return A \code{list} containing the following attributes:
#' \describe{
#'      \item{coefficient_estimates}{Estimated values of the logistic regression coefficients}
#'      \item{bootstrap_intervals}{Boostrap confidence intervals for each of the estimated coefficients}
#'      \item{confusion_matrix}{A confusion matrix comparing the estimated classes based on your logistic regression model and the true classes}
#'      \item{confusion_metrics}{Metrics from your confusion matrix, including prevalence, accuracy, sensitivity, specificity, false discovery rate, and diagnostic odds ratio}
#' }
#' @author John Skarmeas, Kwabena Bayity, Alireza Bidar
#' @importFrom stats optim
#' @examples
#' set.seed(123)
#' x <- data.frame(X1 = rnorm(20), X2 = rgamma(20, 5, 1), X3 = rbeta(20, 1, 1), X4 = sample(c('yes', 'no'), 20, TRUE, prob = c(0.7, 0.3)))
#' y <- as.factor(sample(c(0,1), 20, replace = TRUE, prob = c(0.3, 0.7)))
#' logistic_regression(x, y, B = 50, conf_level = 0.95)
#'@export
logistic_regression <- function(x, y, B = 20, conf_level = 0.95){

  #create a matrix for the x-variables that turns the
  #categorical variables into dummy variables
  x_matrix <- model.matrix(~., data = data.frame(x))

  #if the y-variables are not already in terms of 1 and 0,
  #rewrite the y-variable so that it is only 0s and 1s
  y_matrix <- ifelse(y == unique(y)[1], 0, 1)

  #this is our starting point when optimizing
  #we use the least squares estimates as a starting point
  betaval <- solve(t(x_matrix)%*%x_matrix)%*%t(x_matrix)%*%y_matrix

  #this is the function we need to minimize
  beta_function <- function(x_matrix, y_matrix, betaval){
    p <- 1/(1+exp((-x_matrix)%*%betaval))
    sum(-y_matrix %*% log(p)-(1-y_matrix)%*%log(1-p))
  }

  #using the optim function and the least squares estimates as a starting point
  #we minimize the beta_function above
  coef_ests <- optim(par = betaval, fn = beta_function, x_matrix = x_matrix, y_matrix = y_matrix)

  #this empty matrix is used to store the B number of coefficient estimates we get
  #from our bootstrap sampling procedure
  boot_samp_mat <- matrix(NA, nrow = ncol(x_matrix), ncol = B)

  #this just names the rows in the above matrix for convenience
  rownames(boot_samp_mat) <- colnames(x_matrix)

  #this is the bootstrap sampling procedure which we do B times
  for(i in 1:B){
    #randomly sample row numbers (based on the number of rows in the original data set)
    #with replacement
    bootrows <- sample(nrow(x_matrix), nrow(x_matrix), replace = TRUE)

    #extract the "x"-data corresponding to the above random sample of rows
    bootsamp_x <- x_matrix[bootrows,]

    #extract the "x"-data corresponding to the above random sample of rows
    bootsamp_y <- y_matrix[bootrows]

    #the starting point for our bootstrap estimated coefficients
    betaval_boot <- solve(t(bootsamp_x)%*%bootsamp_x)%*%t(bootsamp_x)%*%bootsamp_y

    #find our bootstrap estimates for the coefficients and store them in the empty matrix
    boot_samp_mat[,i] <- optim(par = betaval_boot, fn = beta_function, x_matrix = bootsamp_x, y_matrix = bootsamp_y)$par
  }

  #find our bootstrap intervals by getting the lower and upper quantiles from each
  #row in our bootstrap matrix
  bootstrap_intervals <- apply(boot_samp_mat, 1, quantile, probs = c((1-conf_level)/2, 1-(1-conf_level)/2))

  #make the format of the intervals look nicer
  bootstrap_intervals <- t(bootstrap_intervals)

  #calculate the log-odds based on our coefficient estimates for the full data
  log_odds <- x_matrix %*% coef_ests$par

  #transform the log-odds into probabilities
  prob_ests <- exp(log_odds)/(1+exp(log_odds))

  #predict either 1 or 0 based on the probability estimates above.
  #this is done for every data row in the original data set.
  predictions <- ifelse(prob_ests > 0.5, 1, 0)

  #using our predictions and the true data, create a confusion matrix.
  confusion <- table(y_matrix, predictions)
  rownames(confusion) <- unique(y)
  colnames(confusion) <- unique(y)

  #this calculates all the metrics from the confusion matrix
  confusion_metrics <- data.frame(
    prevalence = sum(confusion[2,])/sum(confusion),
    accuracy = sum(diag(confusion))/sum(confusion),
    sensitivity = confusion[2,2]/sum(confusion[2,]),
    specificity = confusion[1,1]/sum(confusion[1,]),
    false_discovery_rate = confusion[1,2]/sum(confusion[,2]),
    diagnostic_odds_ratio = (confusion[2,2]/confusion[2,1])/(confusion[1,2]/confusion[1,1]))

  #these results are returned in a list when the function is run
  list(coefficient_estimates = coef_ests$par,
       bootstrap_intervals = bootstrap_intervals,
       confusion_matrix = confusion,
       confusion_metrics = confusion_metrics)

}
