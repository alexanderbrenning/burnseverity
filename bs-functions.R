############################################################################
## Code related to the paper:
## "Benchmarking burn severity predictors using Sentinel-2 data and
## penalized regression: the 2016-17 Chilean fire-storm"
## (c) M.A. Pena and A. Brenning
############################################################################
## bs-functions.R - functions for fitting tuned ridge regression models
##                  and estimating model performance
############################################################################


#' Customized error function for burn severity paper
#'
#' @param obs numeric vector: observed values
#' @param pred numeric vector: predicted values
#'
#' @return a named list with several error measures, including the RMSE
#' @export
#'
#' @examples
error <- function(obs, pred) {
  list(
    bias = mean(obs - pred), 
    stddev = sd(obs - pred),
    rmse = sqrt(mean((obs - pred)^2)),
    dev05 = mean(abs(obs - pred) > 0.5),
    dev1 = mean(abs(obs - pred) > 1.0),
    outside = mean((pred < 0) | (pred > 3)),
    count = length(obs)
  )
}


#' Prediction function with outputs restricted to an interval
#'
#' @param object a fitted regression model whose `predict()` method will be called
#' @param newdata data frame with new data to be passed to `object`'s `predict()` method
#' @param ylim interval to which the predictions are to be limited by bottom- 
#'   and top-coding the raw predictions
#' @param ... additional arguments to be passed to the `predict()` method
#'
#' @return a numeric vector of predictions that are strictly limited to the 
#'   interval given by the `ylim` argument
#' @export
#'
#' @examples
mod_predict <- function(object, newdata, ylim = c(0,3), ...) {
  pred <- predict(object, newdata, ...)
  pred[pred < ylim[1]] <- ylim[1]
  pred[pred > ylim[2]] <- ylim[2]
  pred
}


#' Fit a `glmnet` regression model
#'
#' @param formula model formula
#' @param data a data frame with training data
#' @param lambda the `lambda` penalty, passed to [glmnet::glmnet]
#' @param prefit 
#' @param ... 
#' 
#' @details This function provides a convenient wrapper for the [glmnet::glmnet]
#'   function for penalized regression models. It comes with a formula 
#'   interface, which is needed for model assessment with the 
#'   [sperrorest::sperrorest] function.
#'   This function is used by [tuned_glmnet()] as a work horse.
#'
#' @return a fitted `fglmnet` model; this class inherits from class `glmnet`,
#'   which is extended by storing additional metadata and the training sample
#' @export
#'
#' @examples
fglmnet <- function(formula, data, lambda = 1, prefit = NULL, ...) {
  xvars <- all.vars(formula)[-1]
  yvar <- all.vars(formula)[1]
  if (is.null(prefit)) {
    prefit <- glmnet::glmnet(x = as.matrix(data[,xvars]), 
                     y = data[,yvar], lambda = lambda, ...)
  }
  
  prefit$meta <- list(
    formula = formula,
    xvars = xvars,
    yvar = yvar,
    lambda = lambda,
    training_data = data
  )
  class(prefit) <- c("fglmnet", class(prefit))
  return(prefit)
}


#' `predict` method for `fglmnet` models
#'
#' @param object Fitted `fglmnet` model.
#' @param newdata Data frame with predictor variables
#' @param lambda Chosen \eqn{\lambda} value.
#' @param type `type` argument of [glmnet::predict.glmnet()], e.g. `"response"` (default)
#'
#' @return Model predictions, depends on `"type"`; see [glmnet::predict.glmnet()] for details.
#' @export
predict.fglmnet <- function(object, newdata = NULL, lambda = NULL, type = "response") {
  if (is.null(lambda)) lambda <- object$meta$lambda
  if (is.null(newdata)) newdata <- object$meta$training_data
  xvars <- object$meta$xvars
  object$meta <- NULL
  
  # couldn't get NextMethod call to work...  
  class(object) <- class(object)[ class(object) != "fglmnet" ]
  as.vector(predict(object, newx = as.matrix(newdata[,xvars]), 
                    s = lambda, type = type))
}



#' Self-tuning penalized regression using `glmnet`
#'
#' @inheritParams fglmnet
#' @param alpha Elastic net parameter; defaults to 0 for ridge regression.
#' @param lambdas Candidate regularization parameter values
#' @param smp_fun,smp_args,... Arguments for [sperrorest::sperrorest()]
#' @param verbose Numeric; silent if 0.
#'
#' @return Fitted `fglmnet` object.
#' @export
tuned_glmnet <- function(formula, data, alpha = 0, 
                         lambdas = 10^seq(-3, 5, length = 30),
                         smp_fun = partition_cv,
                         smp_args = list(repetition = 1:1, 
                                         nfold = 10, seed1 = 321),
                         verbose = 0, ...) {
  # Gather RMSEs for all tested lambda values in this vector:
  rmses <- rep(NA, length(lambdas))

  # Estimate RMSEs for all candidate lambda values:  
  # I guess there was a reason why I didn't parallelize this for loop and
  # the sperrorest call inside the for loop; it likely has to do with the
  # handling of global objects in R, with namespaces and the like - but
  # performance was not a critical issue in this analysis...
  for (i in 1:length(lambdas)) {
    res <- sperrorest::sperrorest(formula = formula, data = data,
                      model_fun = fglmnet, 
                      model_args = list(alpha = alpha, 
                                        lambda = lambdas[i]),
                      smp_fun = smp_fun,
                      smp_args = smp_args,
                      mode_rep = "sequential",
                      mode_fold = "sequential",
                      progress = verbose > 1,
                      ...)
    # extract cross-validation RMSE:
    rmses[i] <- summary(res$error_rep)["test_rmse",1]
  }
  
  # Identify optimal lambda value:
  # Note that which.min() returns the first minimum in the unlikely case that
  # two RMSEs are identical; this is an unintended side effect in this case
  # since we'd rather prefer a larger lambda penalty, but tied values are
  # SO unlikely in this context...
  opt.lambda <- lambdas[ which.min(rmses) ]
  
  if (verbose) {
    cat("Optimal lambda: ", round(opt.lambda,4),
        " = 10^", round(log10(opt.lambda),2), 
        "   RMSE: ", round(min(rmses), 3), "\n", sep = "")
  }

  # Fit tuned model with optimal lambda parameter, and return it:  
  fglmnet(formula, data = data, 
          alpha = alpha, lambda = opt.lambda)
}

