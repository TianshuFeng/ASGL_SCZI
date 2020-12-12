# ASGL_SCZI
Adaptive sparse group lasso regularized semi-continuous zero-inflated model

This repo includes R codes to demonstrate the adaptive sparse group lasso regularized semi-continuous zero-inflated model. 

```
source("SCZI_ASGL.R")
set.seed(1)
```

## Data simulation

First, we use `dataSim()` to simulate data from zero-inflated gamma distribution.
```
dat <- dataSim(n_obs = 300, n_var = 15,n_fullZero = 7,n_partialZero = 3)
```

## Model fitting -----------

Simulated dataset is then fitted with `ASGL_ZI()`. Here, the tuning parameter `alpha` determines the convex combination of the lasso and group lasso penalties, and `lambda` is lambda.min relative to lambda.max or the lambda sequence for the regularization path. Details of the tuning paramters can be found in R package `lsgl`. Details of the return values of `ASGL_ZI` can be found in the R code document.
```
model_gamma <- ASGL_ZI(y = dat$y_gamma, X = dat$X, 
                       family = Gamma(link = "log"), 
                       alpha = 0.1, lambda = 0.001)
```

## Model selection

After fitting the model, `bic_ASGL_ZI()` can be used to calculate BIC for different lambdas.

```
bic_res <- bic_ASGL_ZI(model_gamma)
bic_selected <- which.min(bic_res$bic)
```

# For alpha

We can also select alpha using the following codes.

```
model_gamma_list <- list()
bic_list <- list()
smallest_bic <- c()

for (alpha in 1:9/10) {
  model_gamma_list[[alpha*10]] <- ASGL_ZI(y = dat$y_gamma, X = dat$X, 
                         family = Gamma(link = "log"), 
                         alpha = alpha, lambda = 0.001)
  bic_res <- bic_ASGL_ZI(model_gamma_list[[alpha*10+1]])
  bic_list[[alpha*10]] <- bic_res
  smallest_bic <- c(smallest_bic, min(bic_res$bic))
  smallest_bic_idx <- c(smallest_bic, which.min(bic_res$bic))
}
```

## Result evaluation
Finally, we calculate the mean squared error and percentage of coefficients correctly selected with the following codes.
beta_true <- c(dat$true_zero_beta, dat$true_zero_tau)
B <- matrix(beta_true, ncol = 1)

mse <- sapply(model_gamma_list[[which.min(smallest_bic)]]$model_fit$beta, 
              function(beta) sum((B - beta)^2))


percent_correctly_selected <- sapply(model_gamma_list[[which.min(smallest_bic)]]$model_fit$beta, 
                                     function(beta) percent_true_coef(beta_est = beta,  B))
