############################################################################
## Code related to the paper:
## "Benchmarking burn severity predictors using Sentinel-2 data and
## penalized regression: the 2016-17 Chilean fire-storm"
## (c) M.A. Pena and A. Brenning
############################################################################
## Import and prepare data for modeling
############################################################################

d <- read.csv("tablaCBI_2019_20_indices y difbanda_filled data.csv", 
              header = TRUE, sep = ";", na.strings = c("","na"))

# Select only used columns:
d <- d[, c(6,13,14,17:76)]

# rename x/y columns:
colnames(d)[2:3] <- c("x", "y")

# 10 out of 76 observations have missing data:
sum(!complete.cases(d)) 
# one observation has no spectral information:
d[1,] 
# the other 9 have no CBI.Dosel:
d$CBI.Dosel[!complete.cases(d)] 

# Remove cases with missing data:
d <- d[complete.cases(d), ] # 75 observations

# Save cleaned data:
saveRDS(d, file = "burnseverity.rds")

# Here's what we've got:
head(d)
nrow(d)
