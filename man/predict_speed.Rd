% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_changes.R
\name{predict_speed}
\alias{predict_speed}
\title{Predict Speed of Change Based on Initial Value}
\usage{
predict_speed(
  data_model,
  min = NULL,
  max = NULL,
  lambdas = NULL,
  granularity = 0.1,
  floor = 0,
  ceiling = 100,
  verbose = TRUE
)
}
\arguments{
\item{data_model}{A data frame with the training data. Obtained through `prep_data()`. Must contain:
- `initialvalue`: numeric, the initial level of the indicator.
- `change`: numeric, the change in the indicator.
- `fold_id`: numeric or factor, used for cross-validation.}

\item{min}{Optional. Minimum value of `initialvalue` to predict. If `NULL`,
it is set to `floor` (if provided) or `min(data_model$initialvalue)`.}

\item{max}{Optional. Maximum value of `initialvalue` to predict. If `NULL`,
it is set to `ceiling` (if provided) or `max(data_model$initialvalue)`.}

\item{lambdas}{Optional. A vector of lambda values, i.e., levels of flexibility that should be tried and evaluated through cross validation
Defaults to `0.1 * 1.148^(0:50)`.}

\item{granularity}{Numeric. Granularity in outcome variable. Default is `0.1`.}

\item{floor}{Numeric or `NULL`.Minimum value of indicator.
If `NULL`, predictions are unrestricted on the lower end.}

\item{ceiling}{Maximum value of indicator (NA if none).
If `NULL`, predictions are unrestricted on the upper end.}

\item{verbose.}{Logical. If TRUE, display messages in console. Default is TRUE}
}
\value{
A data frame with:
  - `initialvalue`: sequence of values used for prediction.
  - `predictions_speed`: predicted change.
}
\description{
Fits a generalized constrained regression quantile model (`gcrq`) to estimate
the change in a variable as a function of its initial value
The function can restrict predictions to a range between a specified `floor` and `ceiling`
(e.g., targets), or default to the observed data range if these are not specified.
}
