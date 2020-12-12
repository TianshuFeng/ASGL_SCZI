library(MASS)
library(lsgl)


#-----------------------------------------#
#### Data generation ####
#-----------------------------------------#

#' Generate simulation data
#'
#' @param n_obs Int, number of observations
#' @param n_var Int, number of variables
#' @param n_fullZero Int, number of variables that have no influence on the
#'   entire model
#' @param n_partialZero Int, number of variables that have partial influence on
#'   either zero or nonzero part
#' @param mean_beta Numeric, mean of beta for random generation, for ZERO part
#' @param mean_tau Numeric, mean of tau for random generation, for NONZERO part
#' @param sd_beta Numeric, standard deviation of beta for random generation
#' @param sd_tau Numeric, standard deviation of tau for random generation
#' @param nu Numeric, set nu (shape) for gamma distribution and sd for lognormal
#'   distribution
#' @param seed Int for set.seed
#'
#' @return A list. design matrix `X`, response for ZIG `y_gamma`, response for
#'   ZILN `y_LN`, coefficients `beta`, `tau` , vetor of probability of y > 0
#'   `p_1`, binary vector indicating if y = 0 `y_1`, expectation of ZIG and ZILN
#'   `mu`, scalar controling the shape of gamma and lognormal `nu`=.
#' @export
#'
#' @examples
dataSim <- function(n_obs, n_var, 
                    n_fullZero = 0, n_partialZero = 0, 
                    mean_beta = 1, mean_tau = 0.4,
                    sd_beta = 0.7, sd_tau = 0.2,
                    nu = 3,
                    seed = 1) {
  set.seed(seed)
  
  ## Generate the coefficients
  beta <- c(rnorm(n = ceiling(n_var/2), 
                  mean = mean_beta, 
                  sd = sd_beta), 
            rnorm(n = n_var - ceiling(n_var/2), 
                  mean = -mean_beta, 
                  sd = sd_beta))
  tau <- rnorm(n = n_var,
               mean = beta/2,
               sd = sd_tau)
  true_zero_beta <- rep(1, n_var)
  true_zero_tau <- rep(1, n_var)

  
  if(n_fullZero > n_var) {
    stop("n_fullZero is larger than the number of variables")
  } else {
    # Get the index of coefficients that are set to zero for both beta and tau
    fullZeroIdx <- sample(1:n_var, 
                          size = n_fullZero)  
    beta[fullZeroIdx] <- rnorm(n = n_fullZero, 0, 0.005)
    true_zero_beta[fullZeroIdx] <- 0
    tau[fullZeroIdx] <- rnorm(n = n_fullZero, 0, 0.005)
    true_zero_tau[fullZeroIdx] <- 0
  }
  
  if(n_partialZero > n_var - n_fullZero) {
    stop("n_partialZero is larger than n_var - n_fullZero")
  } else {
    # Get the index of coefficients that are set to zero for both beta and tau
    partialZeroIdx <- sample((1:n_var)[-fullZeroIdx], 
                             size = n_partialZero)  
    midIdx <- ceiling(n_partialZero/2)
    
    # First half indices are for beta, second half are for tau
    beta[partialZeroIdx[1:midIdx]] <- rnorm(n = midIdx, 0, 0.005)
    true_zero_beta[partialZeroIdx[1:midIdx]] <- 0
    tau[partialZeroIdx[(midIdx+1):n_partialZero]] <- rnorm(n = n_partialZero-midIdx, 0, 0.005)
    true_zero_tau[partialZeroIdx[(midIdx+1):n_partialZero]] <- 0
  }
  
  
  ## Generate the data X
  X <- matrix(rnorm(n_obs*n_var, 0, 1), 
              nrow = n_obs)
  
  p_1 <- 1/(1 + exp(- X %*% beta))  # pobability of y>0
  y_1 <- rbinom(n_obs, 1, p_1)  # y_1 = 0 if y = 0, y_1 = 1 if y > 0
  
  y_gamma <- y_1  # For ZIG gamma

  pos_obs <- which(y_1 == 1)
  mu <- exp(X[pos_obs, ] %*% tau)
  
  # Shape is nu and scale is mu
  nu <- 3
  y_gamma[pos_obs] <- rgamma(length(pos_obs), 
                             shape = nu, 
                             scale = mu)/nu
  
  return(list(X = X, y_gamma = y_gamma, 
              beta = beta, tau = tau, 
              p_1 = p_1, y_1 = y_1,
              mu = mu, nu = nu,
              true_zero_beta = true_zero_beta, true_zero_tau = true_zero_tau))
}


#-----------------------------------------#
#### Model ####
#-----------------------------------------#

#' Adaptive sparse group lasso for zero inflated model
#'
#' @param y response matrix, matrix of size n.
#' @param X design matrix, matrix of size n \times p.
#' @param omage scalar, tuning parameter for adaptive lasso
#' @param alpha scalar, balancing the l1 and l2 penalty.
#' @param lambda lambda.min relative to lambda.max or the lambda sequence for
#'   the regularization path.
#' @param family a description of the error distribution and link function to be
#'   used in the model. See \code{\link[stats]{family}} for details of family
#'   functions.
#' @param grouping a numeric vector specifying the group of variables. Each
#'   element of the factor/vector specifying the group of the feature.
#' @param if_weight if a vector of ‘prior weights’ to be used in the fitting
#'   process. Decided based on the percentage of zero reponses.
#'
#' @return A list. \describe { \item{model_fit}{a `lsgl` object with the fitted
#'   model.} \item{cv_fit}{a `lsgl` object with the cross validation results.}
#'   \item{cv_index}{the index of the best model identified by cross
#'   validation.} \item{cv_beta}{a vector of coefficients from the model
#'   selected by cv.} \item{beta_mle}{a vector of mle estimate of beta.}
#'   \item{alpha}{the sparse group lasso penalty}
#'   \item{sigmaInverseSqrt}{\hat{\Sigma}^{1/2}} \item{family}{an object of
#'   class "family"} \item{weight_logistic}{weight for logistic regression} }
#'
#' @export
#'
#' @examples
ASGL_ZI <- function(y, X, omega = 1,
                    alpha = 0.5, lambda = 0.01,
                    family = Gamma(link = "log"),
                    grouping = NULL,
                    if_weight = FALSE) { # For ZIG model, omage: tuning parameter
  require(MASS)
  require(lsgl)
  
  n_obs <- nrow(X)
  n_var <- ncol(X)
  y_binary <- as.numeric(y > 0)
  
  #### ___ Fit the glm model to get the mle ####
  # binary part
  if(if_weight) {
    num_zero <- sum(y_binary == 0)
    num_nonzero <- length(y_binary) - num_zero
    weight_logistic <- rep(1, length(y_binary))
    weight_logistic[y_binary != 0] <- ceiling(num_zero/(num_nonzero + 1))
  } else {
    weight_logistic <- NULL
  }
  
  model_logit <- glm(y_binary ~ X, family = binomial(link="logit"), weights = weight_logistic, maxit = 100)
  sigma_logit <- vcov(model_logit)
  
  # gamma part
  model_nonzero <- glm(y[y>0] ~ X[y>0,], family = family)
  # MASS::gamma.shape(model_nonzero)  # To get the nu and its sd
  sigma_nonzero <- vcov(model_nonzero)
  
  # Final sigma matrix
  sigma <- matrix(0, nrow = n_var*2, ncol = n_var*2)
  sigma[1:n_var, 1:n_var] <- sigma_logit[-1,-1]
  sigma[(n_var+1):(2*n_var), (n_var+1):(2*n_var)] <- sigma_nonzero[-1,-1]
  
  #### ___ Construct the linear regression for lasso/group lasso ####
  sigmaInverseSqrt <- inverseSqrtMatrix(sigma)
  beta_mle <- c(model_logit$coefficients[-1], model_nonzero$coefficients[-1])
  names(beta_mle) <- NULL
  y_lm <- sigmaInverseSqrt %*% beta_mle
  
  
  #### ___ Construct the adaptive lasso weights for parameters ####
  adaPara <- matrix(1/abs(beta_mle)^omega, nrow = 1)
  
  #### ___ Construct the groups and associated weights ####
  
  if(is.null(grouping)){
    # right now, just for the variables shared by the two parts
    grouping <- c(1:n_var, 1:n_var)
  } 
  groupingWeight <- rep(1, length.out = length(unique(grouping)))
  for(i in unique(grouping)) {
    idx <- which(grouping == i)
    groupingWeight[i] <- 1/norm(matrix(beta_mle[idx]), type = "F")^omega
  }
  
  
  #### ___ Fit the model with lsgl ####
  fit <- lsgl::fit(x = sigmaInverseSqrt, 
                   y = y_lm, 
                   alpha = alpha, 
                   lambda = lambda,
                   intercept=FALSE,
                   grouping = grouping,
                   groupWeights = groupingWeight,
                   parameterWeights = adaPara
  )
  
  return(list(model_fit = fit, 
              # cv_fit = fit.cv, 
              # cv_index = idx_cv, 
              # cv_beta = fit$beta[[idx_cv]],
              beta_mle = beta_mle, 
              alpha = alpha,
              omega = omega,
              sigmaInverseSqrt = sigmaInverseSqrt,
              grouping = grouping,
              n_obs = n_obs,
              family = family,
              weight_logistic = weight_logistic))
}


#-----------------------------------------#
#### Variable selection ####
#-----------------------------------------#

#' Calculate the Bayesian information criterion for a fitted model
#'
#' @param fit The fitted model from ASGL_ZI
#'
#' @return
#' @export
#'
#' @examples
bic_ASGL_ZI_poor_approx <- function(fit) {
  # Claculate BIC for each lambda
  
  bic <- c()
  df_list <- c()
  
  i <- 0
  for (beta_lambda in fit$model_fit$beta) {
    i <- i + 1
    
    # Calculate the weighted beta
    weighted_beta <- fit$sigmaInverseSqrt%*%(beta_lambda - fit$beta_mle)
    
    # Calculate the degree of freedom
    df_lasso <- sum(beta_lambda != 0) # df for the lasso part
    df_group <- 0  # df for the group lasso part
    for(j in unique(fit$grouping)) {
      idx <- which(fit$grouping == j)
      beta_lambda_group <- beta_lambda[idx]
      beta_mle_group <-  fit$beta_mle[idx]
      
      if(norm(matrix(beta_lambda_group), 'f') > 0){
        df_group <- df_group + 1 + norm(matrix(beta_lambda_group), 'f')/norm(matrix(beta_mle_group), 'f') * (length(idx) - 1)
      }
    }
    df <- fit$alpha * df_group + (1-fit$alpha) * df_lasso
    
    df_list <- c(df_list,
                 df)
    bic <- c(bic, 
             sum(weighted_beta^2) + log(fit$n_obs) * df)
  }
  return(list(bic = bic, df_list = df_list))
}



bic_ASGL_ZI <- function(fit) {
  
  bic <- c()
  df_list <- c()
  i <- 0
  alpha <- fit$alpha
  
  for(lambda in fit$model_fit$lambda) {
    i <- i + 1
    
    # Calculate the weighted beta
    beta_lambda <- fit$model_fit$beta[[i]]
    weighted_beta <- fit$sigmaInverseSqrt%*%(beta_lambda - fit$beta_mle)
    
    # Calculate the degree of freedom
    df_hat <- 0
    for(j in unique(fit$grouping)) {
      
      idx <- which(fit$grouping == j)
      beta_lambda_group <- beta_lambda[idx]
      beta_mle_group <-  fit$beta_mle[idx]
      
      tmp_soft_thres <- soft_thres(beta_mle_group, alpha * lambda)
      
      if(norm(matrix(tmp_soft_thres), 'f') > (1 - alpha) * lambda){
        for(idx_coef in 1:length(tmp_soft_thres)) {
          if(abs(beta_mle_group[idx_coef]) > alpha * lambda) {
            df_hat = df_hat + 
              1 - (1 - alpha) * lambda * (sum(tmp_soft_thres^2) - 
                                            tmp_soft_thres[idx_coef]^2)/sum(tmp_soft_thres^2)^(3/2)
          }
        }
      }
    }
    df_list <- c(df_list,
                 df_hat)
    bic <- c(bic, 
             sum(weighted_beta^2) + log(fit$n_obs)/2 * df_hat)
    
  }
  return(list(bic = bic, df_list = df_list))
}


soft_thres <- function(z, alpha_lambda) {
  last_part <- pmax(abs(z) - alpha_lambda, 0)
  return(sign(z) * last_part)
}


#-----------------------------------------#
#### Utility functions ####
#-----------------------------------------#

#' Calculate the square root of the inverse of a matrix
#'
#' Use eigenvalue decomposition
#'
#' @param x The matrix
#'
#' @return The square root of the inverse of the given matrix.
#'
#' @examples
inverseSqrtMatrix <- function(x) {
  eigen_decomp <- eigen(x)
  
  V <- eigen_decomp$vectors
  
  B <- V %*% diag(1/sqrt(eigen_decomp$values)) %*% t(V)
  
  return(B)
}


#' Calculate the the percentage of correctly identified coefficients for given betas
#'
#' @param beta_est Estimated coefficients
#' @param beta_true True coefficients
#'
#' @return A scalalr. The percentage of correctly identified coefficients
#'
#' @examples
percent_true_coef <- function(beta_est, beta_true) {
  
  # num_discover <- sum(beta_est != 0 & beta_true != 0)
  # 
  # num_true <- sum(beta_true != 0)
  # num_est <- sum(beta_est != 0)
  
  # return(min(num_discover/num_est, num_discover/num_true))
  
  return(1 - sum( abs(as.numeric(beta_est != 0) - as.numeric(beta_true != 0)) )/ length(beta_true))
}

