## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 3.5,
  fig.align = "center",
  warning = FALSE,
  message = FALSE
)

## ----library------------------------------------------------------------------
library(BarmaRidgeBJPS2025)
library(forecast)   # time series
library(ggplot2)    # plotting
library(zoo)        # yearmon handling
library(dplyr)      # data manipulation
library(gridExtra)  # grid.arrange for plots

# Bootstrap libraries
library(doMC)
library(doRNG)

## ----config_data_visualization------------------------------------------------
# Series size and end time
sample_size <- length(brasilia_ts)
end_time_series <- time(brasilia_ts)[sample_size] + 1

# Size font of the plot
ggplot_size_font <- 9

# ggplot theme
ggplot_theme <- theme(
  legend.position = "bottom",
  title = element_text(size = ggplot_size_font),
  axis.text = element_text(size = ggplot_size_font),
  axis.title = element_text(size = ggplot_size_font),
  legend.text = element_text(size = ggplot_size_font),
  legend.title = element_text(size = ggplot_size_font)
  )

# y-axis scale
ggplot_scale_y <-
  scale_y_continuous(
    breaks = seq(0.0, 1.00, 0.10)
  )

# x-axis scale
ggplot_scale_x <- scale_x_continuous(
  breaks = seq(1999, end_time_series, 2),
  limits = c(1999, end_time_series)
)

## ----ggplot_brasilia----------------------------------------------------------
# -------------------------------------------------------------------------- #
# Time series plot
# -------------------------------------------------------------------------- #
ggplot_brasilia <- ggplot(brasilia_df, aes(time, y)) +
  geom_point(size = 0.5) +
  geom_line(aes(group = 1)) +
  ggplot_scale_x +
  ggplot_scale_y +
  ggplot_theme

# -------------------------------------------------------------------------- #
# Autocorrelation function (ACF) plot
# -------------------------------------------------------------------------- #
ggacf_plot <- ggAcf(brasilia_ts) +
  ggtitle(" ") +
  xlab("Lag") +
  ylab("ACF")

# -------------------------------------------------------------------------- #
# Partial autocorrelation function (PACF) plot
# -------------------------------------------------------------------------- #
ggpacf_plot <- ggPacf(brasilia_ts) +
  ggtitle(" ") +
  xlab("Lag") +
  ylab("PACF")

## ----print_ggplot_brasilia, message=FALSE, echo=FALSE, fig.height=5, fig.cap="Relative humidity in Brasília (top), with its corresponding autocorrelation function (bottom left) and partial autocorrelation function (bottom right)."----
# -------------------------------------------------------------------------- #
# Grid all plots
# -------------------------------------------------------------------------- #
grid.arrange(
  ggplot_brasilia,
  ggacf_plot,
  ggpacf_plot,
  layout_matrix = rbind(c(1, 1), c(2, 3))
)

## ----descriptive--------------------------------------------------------------
# -------------------------------------------------------------------------- #
# Descriptive statistics: Useful volume of the brasilia
# -------------------------------------------------------------------------- #
descriptive_df <- data.frame(
  min             = min(brasilia_ts),
  max             = max(brasilia_ts),
  median          = median(brasilia_ts),
  mean            = mean(brasilia_ts),
  sd              = sd(brasilia_ts),
  skewness        = moments::skewness(brasilia_ts),
  excess_kurtosis = moments::kurtosis(brasilia_ts)
)

# Round to two decimals
descriptive_df <- round(descriptive_df, 2)

## ----print_descriptive, echo=FALSE--------------------------------------------
knitr::kable(
  descriptive_df,
  caption = 
    "Descriptive statistics of the brasilia relative humidity in Brasília."
)

## ----sample-------------------------------------------------------------------
# The length of the training sample
sample_size <- 150

# The lag in months (how far back to go) 
# In this example specifies the data from April 2006 to September 2018.
lag_months <- 69   

# Define n_obs before using it.
n_obs <- length(brasilia_ts)

# The end of the primary training period (our anchor point)
end_train <- time(brasilia_ts)[n_obs]

# The end of our lagged sample window
end_train_lag <- end_train - (lag_months / 12)

# The start of our lagged sample window
start_train_lag <- end_train_lag - deltat(brasilia_ts) * (sample_size - 1)

# Sample
# ----------------------------------------------------------------------- #
y_sample_ts <- window(brasilia_ts,
                      start = start_train_lag,
                      end = end_train_lag)

## ----ggplot_brasilia_subsample------------------------------------------------
len_y_sample <- length(y_sample_ts)

start_subsample_time <- time(y_sample_ts)[1]
end_subsample_time <- time(y_sample_ts)[len_y_sample]

# -------------------------------------------------------------------------- #
# Time series plot
# -------------------------------------------------------------------------- #
ggplot_brasilia <- ggplot(brasilia_df, aes(time, y)) +
  geom_point(size = 0.5) +
  geom_line(aes(group = 1)) +
  ggplot_scale_x +
  ggplot_scale_y +
  ggplot_theme +
  guides(colour = guide_legend(title = "Relative humidity in Brasília")) +

  # Vertical line for the start of the sub-sample
  geom_vline(
    xintercept = start_subsample_time,
    linetype = "dashed",
    linewidth = 0.5,
    color = "red"
  ) +

  # Vertical line for the end of the sub-sample
  geom_vline(
    xintercept = end_subsample_time,
    linetype = "dashed",
    linewidth = 0.5,
    color = "red"
  )

# -------------------------------------------------------------------------- #
# Autocorrelation function (ACF) plot
# -------------------------------------------------------------------------- #
ggacf_plot <- ggAcf(brasilia_ts) +
  ggtitle(" ") +
  xlab("Lag") +
  ylab("ACF")

# -------------------------------------------------------------------------- #
# Partial autocorrelation function (PACF) plot
# -------------------------------------------------------------------------- #
ggpacf_plot <- ggPacf(brasilia_ts) +
  ggtitle(" ") +
  xlab("Lag") +
  ylab("PACF")

## ----print_ggplot_brasilia_subsample, message=FALSE, echo=FALSE, fig.height=5, fig.cap=" Relative humidity in Brasília (top), with its corresponding autocorrelation function (bottom left) and partial autocorrelation function (bottom right)."----
# -------------------------------------------------------------------------- #
# Grid all plots
# -------------------------------------------------------------------------- #
grid.arrange(
  ggplot_brasilia,
  ggacf_plot,
  ggpacf_plot,
  layout_matrix = rbind(c(1, 1), c(2, 3))
)

## ----regressor----------------------------------------------------------------
# Get Time Series Frequency
freq <- frequency(y_sample_ts)

# The number of observations for forecast
n_test <- 6        

# Create Trend Regressors
trend_index <- seq_along(y_sample_ts)
trend_index_hat <- (max(trend_index) + 1):(max(trend_index) + n_test)

# Create Harmonic Regressors 
hs <- sin(2 * pi * trend_index / freq)
hc <- cos(2 * pi * trend_index / freq)
hs_hat <- sin(2 * pi * trend_index_hat / freq)
hc_hat <- cos(2 * pi * trend_index_hat / freq)

# Final
regressors_list <- list(
    trend_index = trend_index,
    trend_index_hat = trend_index_hat,
    hs = hs,
    hc = hc,
    hs_hat = hs_hat,
    hc_hat = hc_hat
)


## ----estimation_setup---------------------------------------------------------
# Model Specification
ar_vec <- 1
ma_vec <- 1:4

## ----estimates_mle------------------------------------------------------------
# The 'penalty' parameter controls the ridge penalization. Setting it to 0 
# effectively disables the penalty, resulting in a standard unpenalized model.
penalty <- FALSE

# Fit the unpenalized BARMA model (standard CMLE)
fit_default <- BarmaRidgeBJPS2025::barma(
  y = y_sample_ts,
  ar = ar_vec,
  ma = ma_vec,
  penalty = penalty,
  X = cbind(
    hs = regressors_list$hs,
    hc = regressors_list$hc
  ), 
  X_hat = cbind(
    hs = regressors_list$hs_hat,
    hc = regressors_list$hc_hat
  )
)

## ----estimates_mle_df---------------------------------------------------------
# Data frame with CMLE estimates
estimates <- cbind(
  CMLE = as.numeric(fit_default$coef)
)

estimates_df <- data.frame(estimates)

# Row names of the LaTeX table
ar_names <- paste0("$\\varphi_", ar_vec, "$")
ma_names <- paste0("$\\theta_", ma_vec, "$")
regressors_names <- c("$\\beta_1$ (sin)", "$\\beta_2$ (cos)")

barma_names <- c("$\\alpha$", 
                 ar_names, 
                 ma_names, 
                 regressors_names, 
                 "$\\phi$")


rownames(estimates_df) <- barma_names
estimates_df <- round(estimates_df, 4)

## ----print_estimates_mle_df, echo=FALSE---------------------------------------
knitr::kable(
  t(estimates_df),
  caption = "**CMLE**; Parameter estimates for the $\\beta\\text{ARMA(1,4)}$ model; relative humidity in Brasília, data from April 2006 to September 2018."
)

## ----penalty------------------------------------------------------------------
# Next, define the small penalty value (lambda) for the ridge regression.
a_max <- max(ar_vec, ma_vec)
penalty <- 1 / (len_y_sample - a_max)^(0.9)

## ----estimates_pmle-----------------------------------------------------------
# Fit the penalized BARMA model (PCMLE)
fit_ridge <- BarmaRidgeBJPS2025::barma(
  y = y_sample_ts,
  ar = ar_vec,
  ma = ma_vec,
  penalty = penalty,
  X = cbind(
    hs = regressors_list$hs,
    hc = regressors_list$hc
  ), 
  X_hat = cbind(
    hs = regressors_list$hs_hat,
    hc = regressors_list$hc_hat
  )
)

## ----estimates_df-------------------------------------------------------------
# Data frame with CMLE and PCMLE estimates
estimates <- cbind(
  CMLE = as.numeric(fit_default$coef),
  PCMLE = as.numeric(fit_ridge$coef)
)

estimates_df <- data.frame(estimates)

# Row names of the LaTeX table
rownames(estimates_df) <- barma_names
estimates_df <- round(estimates_df, 4)

## ----print_estimates_df, echo=FALSE-------------------------------------------
knitr::kable(
  t(estimates_df),
  caption = "**CMLE and PCMLE**; Parameter estimates for the $\\beta\\text{ARMA(1,4)}$ model; relative humidity in Brasília, data from April 2006 to September 2018."
)

## ----bootstrap_conf-----------------------------------------------------------
# -------------------------------------------------------------------------- #
# Bootstrap Configuration
# -------------------------------------------------------------------------- #
# Set a seed for the random number generator to ensure reproducibility.
seed <- 14
# Define the number of CPU cores to use for parallel execution.
n_cores <- 2
# Specify the total number of bootstrap replications to be performed.
n_boot_rep <- 200

## ----bootstrap12--------------------------------------------------------------
# -------------------------------------------------------------------------- #
# Run Bootstrap Procedures
# -------------------------------------------------------------------------- #

# Perform the seasonal block bootstrap using a block length of 12.
block_length <- 12

# This value corresponds to the annual seasonal cycle (12 months).
ridge_boot_seasonal_block12 <-
  BarmaRidgeBJPS2025::BARMAX_ridge_bootstrap_seasonal_parallel(
  n_cores = n_cores,
  y = y_sample_ts,
  seed = seed,
  ar_vec = ar_vec,
  ma_vec = ma_vec,
  fit_ridge = fit_ridge,
  penalty = penalty,
  regressors_list = regressors_list,
  block_length = block_length,
  n_boot_rep = n_boot_rep
)

## ----bootstrap24--------------------------------------------------------------
# Repeat the bootstrap procedure with a different block length of 24.
block_length <- 24

ridge_boot_seasonal_block24 <- 
  BarmaRidgeBJPS2025::BARMAX_ridge_bootstrap_seasonal_parallel(
  n_cores = n_cores,
  y = y_sample_ts,
  ar_vec = ar_vec,
  seed = seed,
  ma_vec = ma_vec,
  fit_ridge = fit_ridge,
  penalty = penalty,
  regressors_list = regressors_list,
  block_length = block_length,
  n_boot_rep = n_boot_rep
)

## ----estimates_final_df-------------------------------------------------------
# -------------------------------------------------------------------------- #
# Format and Print Bootstrap Estimates
# -------------------------------------------------------------------------- #

# Combine original model estimates with the mean bootstrap results.
# The new columns are from the runs with block lengths 12 and 24.
estimates_final_df <- cbind(
  estimates_df,
  PCMLE12 = ridge_boot_seasonal_block12$mean_bootstrap_estimates,
  PCMLE24 = ridge_boot_seasonal_block24$mean_bootstrap_estimates
)

# Rename the new columns using LaTeX for better table formatting in the report.
# PCMLE* is the Penalized Conditional Maximum Likelihood Estimator.
colnames(estimates_final_df)[3] <- "$\\text{PCMLE}^*_{12}$"
colnames(estimates_final_df)[4] <- "$\\text{PCMLE}^*_{24}$"

# Round all numerical estimates to four decimal places for presentation.
estimates_final_df <- round(estimates_final_df, 4)

print_estimates_final_df <- cbind(
  Parameter = rownames(estimates_final_df), estimates_final_df
  )
rownames(print_estimates_final_df) <- NULL

## ----print_estimates_final_df, echo=FALSE-------------------------------------
# Generate a formatted table for the report using the kable() function.
knitr::kable(
  print_estimates_final_df,
  caption = "Parameter estimates for the $\\beta\\text{ARMA(1,4)}$ model;
  relative humidity in Brasília, data from April 2006 to September 2018."
)

## ----session_info, echo=FALSE-------------------------------------------------
cat("=================================================================", "\n")
cat("Session Information", "\n")
cat("=================================================================", "\n")

cat("This report was generated at:", 
    format(Sys.time(), "%B %d, %Y at %I:%M %p"), "\n")
cat("\n")

print(sessionInfo())

