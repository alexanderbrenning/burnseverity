############################################################################
## Code related to the paper:
## "Benchmarking burn severity predictors using Sentinel-2 data and
## penalized regression: the 2016-17 Chilean fire-storm"
## (c) M.A. Pena and A. Brenning
############################################################################
## Data analysis:
## - Exploratory analysis
## - Model assessment of linear regression and ridge regression models
## - Model interpretation of selected ridge regression model
############################################################################

## Load required packages

library("purrr")        # for data handling
library("stringr")      # for operations on character strings
library("MASS")         # what for??

library("glmnet")       # for ridge regression
library("sperrorest")   # for spatial cross-validation

library("plotly")       # for PCA visualization (optional)

# Optional, not required for the paper:
library("iml")          # for model-agnostic ML model interpretation
library("wiml")         # for model interpretation from transformed perspective


## Load burn severity functions, including self-tuning glmnet model

source("bs-functions.R", encoding = "UTF-8")


############################################################################
## Initialization
############################################################################


## Load preprocessed data

d <- readRDS("burnseverity.rds")


## Feature sets:

xvars_bands <- stringr::str_subset(colnames(d), "^B[0-9]*")
xvars_ndvis <- stringr::str_subset(colnames(d), "^NDVI*")
xvars_ndwis <- stringr::str_subset(colnames(d), "^NDWI*")
xvars_nbrs <- stringr::str_subset(colnames(d), "^NBR*")
xvars_nbr <- c("NBR2017") 
xvars_indxs <- c(xvars_nbrs, xvars_ndvis, xvars_ndwis)
xvars_indxs2017 <- stringr::str_subset(xvars_indxs, "2017")
xvars_indxs2018 <- stringr::str_subset(xvars_indxs, "2018")
xvars_indxs2019 <- stringr::str_subset(xvars_indxs, "2019")
xvars_indxs2020 <- stringr::str_subset(xvars_indxs, "2020")

xvars <- c(xvars_bands, xvars_ndvis, xvars_ndwis, xvars_nbrs)

yvar <- "CBI.Dosel"


## Prepare model formulas for the various settings:

fo_null <- CBI.Dosel ~ 1
fo <- as.formula(paste(yvar, "~", paste(xvars, collapse = " + ")))
fo_nbr <- as.formula(paste(yvar, "~", paste(xvars_nbr, collapse = " + ")))
fo_nbrs <- as.formula(paste(yvar, "~", paste(xvars_nbrs, collapse = " + ")))
fo_ndvis <- as.formula(paste(yvar, "~", paste(xvars_ndvis, collapse = " + ")))
fo_ndwis <- as.formula(paste(yvar, "~", paste(xvars_ndwis, collapse = " + ")))
fo_indxs <- as.formula(paste(yvar, "~", paste(xvars_indxs, collapse = " + ")))
fo_indxs2017 <- as.formula(paste(yvar, "~", paste(xvars_indxs2017, collapse = " + ")))
fo_indxs2018 <- as.formula(paste(yvar, "~", paste(xvars_indxs2018, collapse = " + ")))
fo_indxs2019 <- as.formula(paste(yvar, "~", paste(xvars_indxs2019, collapse = " + ")))
fo_indxs2020 <- as.formula(paste(yvar, "~", paste(xvars_indxs2020, collapse = " + ")))
fo_bands <- as.formula(paste(yvar, "~", paste(xvars_bands, collapse = " + ")))


############################################################################
## Exploratory data analysis
############################################################################


## Principal component analysis:

pca_fo <- as.formula(paste("~", paste(xvars, collapse = " + ")))
pca <- prcomp(pca_fo, d[,xvars], scale. = TRUE, center = TRUE)
plot(pca)
biplot(pca)
summary(pca) # first 6 PCs explain 90% of variance

# Visualize PCA results:
qplot(x = 1:length(xvars), y = cumsum(pca$sdev)/sum(pca$sdev), geom = "line")
plot_ly(as.data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color = ~d[,yvar]) %>% add_markers()

# For curiosity, how strongly are the first PCs correlated with the response:
round(cor(cbind(CBI.Dosel = d[,yvar], pca$x[,1:30]))*100)[-1,1]
# mainly first PC seems to be related to CBI.Dosel

# For curiosity, what are the weights of the predictors in the first PCs:
round(100*pca$rotation[,1:4])


## Correlation analysis:

# Print pairs of variables with strongest correlations:
for (ixv in 1:(length(xvars) - 1)) {
  xv <- xvars[ixv]
  for (iyv in (ixv+1):length(xvars)) {
    yv <- xvars[iyv]
    co <- cor(d[,xv], d[,yv], method = "spearman")
    if (co > 0.9) cat(xv, yv, round(co, 2), "\n")
  }
}



############################################################################
## Fit linear regression and ridge models on entire data set
############################################################################

# LM and ridge models fitted on entire sample
# using only NBR features:

# Correlations among NBR predictor variables (for table in the paper):
round(cor(d[,xvars_nbrs]),2)

# Fit multiple linear regression model:
fit.lm <- lm(fo_nbrs, data = d)
summary(fit.lm)
coef.lm <- coef(fit.lm)

fit.ridge <- tuned_glmnet(formula = fo_nbrs, data = d, verbose = 1)
fit.ridge$lambda  # optimal lambda
fit.ridge$a0      # intercept
fit.ridge$beta    # coefficients
coef.ridge <- c(unname(fit.ridge$a0), fit.ridge$beta@x)
names(coef.ridge) <- names(coef.lm)

# Display model coefficients for table 8 in the paper:
coef.lm
coef.ridge

# Just for curiosity, models with all variables:
fit.lm.all <- lm(fo, data = d)
coef.lm.all <- coef(fit.lm.all)
fit.ridge.all <- tuned_glmnet(formula = fo, data = d, verbose = 1)
coef.ridge.all <- c(unname(fit.ridge.all$a0), fit.ridge.all$beta@x)
names(coef.ridge.all) <- names(coef.lm.all)
coef.ridge.all

# Correlation of predictions:
cor(mod_predict(fit.ridge, newdata = d), mod_predict(fit.lm, newdata = d))
# This is just to show how biased a performance assessment on the training set
# would be! -> we must cross-validate tuned_glmnet!

# Just for fun, some plots that are not meant for publication:
par(mfrow = c(1,3))
# Response vs. linear regression prediction on the training set:
plot(d$CBI.Dosel ~ mod_predict(fit.lm, newdata = d), pch = "+", xlab = "LM prediction", ylab = "Overstory CBI")
# Response vs. ridge regression prediction on the training set:
plot(d$CBI.Dosel ~ mod_predict(fit.ridge, newdata = d), pch = "+", xlab = "Ridge prediction", ylab = "Overstory CBI")
# Linear versus ridge regression prediction on the training set:
plot(mod_predict(fit.ridge, newdata = d) ~ mod_predict(fit.lm, newdata = d), pch = "+", xlab = "Ridge prediction", ylab = "LM prediction")
abline(c(0,1))



############################################################################
## Model assessment using spatial cross-validation
############################################################################

# For reproducibility of sperrorest:
smp_args <- list(nfold = 10, repetition = 1:10, seed1 = 123)

### Mean prediction distance in spatial CV resampling:
parti <- partition_kmeans(data = d, nfold = smp_args$nfold,
                          repetition = smp_args$repetition,
                          seed1 = smp_args$seed1)
parti <- add.distance(parti, data = d)
get_distances <- function(x) {
  x %>% map("distance") %>% flatten() %>% unlist()
}
parti %>% map(get_distances) %>% flatten() %>% unlist() %>% mean()
# 16.7 km


## Model assessment using all indices and bands:

res <- sperrorest::sperrorest(formula = fo, data = d, 
                  model_fun = tuned_glmnet,
                  pred_fun = mod_predict,
                  smp_fun = partition_kmeans, 
                  smp_args = smp_args,
                  mode_rep = "loop",
                  mode_fold = "loop",
                  err_fun = error, verbose = 1)
smry <- summary(res$error_rep)[c("train_bias","train_stddev","train_rmse", "test_bias","test_stddev","test_rmse"), ]
round(smry[,c(1,3)], 3)
# for comparison:
sd(d$CBI.Dosel)


## Model assessments for all relevant model / predictor combinations:

analyses <- list(
  null = list(
    model = lm,
    formula = fo_null
  ),
  lmNBR1 = list(
    model = lm,
    formula = fo_nbr
  ),
  lmNBRs = list(
    model = lm,
    formula = fo_nbrs
  ),
  NBRs = list(
    model = tuned_glmnet,
    formula = fo_nbrs
  ),
  NDVIs = list(
    model = tuned_glmnet,
    formula = fo_ndvis
  ),
  NDWIs = list(
    model = tuned_glmnet,
    formula = fo_ndwis
  ),
  Indexes2017 = list(
    model = tuned_glmnet,
    formula = fo_indxs2017
  ),
  Indexes2018 = list(
    model = tuned_glmnet,
    formula = fo_indxs2018
  ),
  Indexes2019 = list(
    model = tuned_glmnet,
    formula = fo_indxs2019
  ),
  Indexes2020 = list(
    model = tuned_glmnet,
    formula = fo_indxs2020
  ),
  Indexes = list(
    model = tuned_glmnet,
    formula = fo_indxs
  ),
  Bands = list(
    model = tuned_glmnet,
    formula = fo_bands
  ),
  All = list(
    model = tuned_glmnet,
    formula = fo
  )
)


for (i in 1:length(analyses)) {
  cat("\n\n", names(analyses)[i], "\n\n")
  vnms <- unique(c(all.vars(analyses[[i]]$fo), c("x","y")))
  res <- sperrorest::sperrorest(formula = analyses[[i]]$fo, 
                    data = d[,vnms], 
                    model_fun = analyses[[i]]$model,
                    pred_fun = mod_predict,
                    smp_fun = partition_kmeans, 
                    smp_args = smp_args,
                    mode_rep = "loop",
                    mode_fold = "loop",
                    err_fun = error, 
                    verbose = 1, progress = FALSE)
  smry <- summary(res$error_rep)[c("train_bias","train_stddev","train_rmse", "test_bias","test_stddev","test_rmse"), ]
  print(round(smry[,c(1,3)], 3))
  analyses[[i]]$result <- res
  analyses[[i]]$summary <- smry
}

## Save results of time-consuming model assessment to file:

# saveRDS(analyses, file = "cvresults.rds")



## Summarize cross-validation results

analyses <- readRDS("cvresults.rds")

train_rmse <- sapply(analyses, function(x) x$summary["train_rmse", "mean"])
train_sd <- sapply(analyses, function(x) x$summary["train_rmse", "sd"])
test_rmse <- sapply(analyses, function(x) x$summary["test_rmse", "mean"])
test_sd <- sapply(analyses, function(x) x$summary["test_rmse", "sd"])

options(digits = 3)
data.frame(train_rmse = train_rmse, 
           train_sd = train_sd,
           test_rmse = test_rmse,
           test_sd = test_sd)



############################################################################
## Model interpretation from a multivariate angle using 'wiml'
## (This is not included in the paper.)
############################################################################

### variable importance of model with 12 indices

library("iml")
library("wiml")

# Fit the model to be analyzed:
fit <- tuned_glmnet(formula = fo, data = d, verbose = 1)

# What kind of plot, at what resolution:
plot_type <- "ale"
grid_size <- 20

# Prepare PCA transformation object:
wrp <- pca_warper(d[,all.vars(fo)], xvars = xvars, yvar = yvar, uvars = NULL)
wfit <- warp_fitted_model(fit, warper = wrp)
wd <- warp(d[,all.vars(fo)], warper = wrp)
wvars <- wrp$wvars[1:14]

# Feature effects and importance, PCA-transformed:
wpredictor <- Predictor$new(wfit, data = wd, y = yvar)

# ALE plot:
weffs <- FeatureEffects$new(wpredictor, features = wvars, method = plot_type, grid.size = grid_size)
plot(weffs)

# Permutation variable importance on training sample, PC features:
wimp <- FeatureImp$new(wpredictor, loss = "rmse", n.repetitions = 30)
plot(wimp)

# Feature effects and importance, untransformed:
predictor <- Predictor$new(fit, data = d[,all.vars(fo)], y = yvar)
effs <- FeatureEffects$new(predictor, features = xvars, method = plot_type, grid.size = grid_size)
plot(effs)

# Permutation variable importance on training sample, original features:
imp <- FeatureImp$new(predictor, loss = "rmse", n.repetitions = 30)
plot(imp)
