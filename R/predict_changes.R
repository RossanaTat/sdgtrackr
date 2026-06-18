#' Predict Speed of Change Based on Initial Value
#'
#' Fits a generalized constrained regression quantile model (`gcrq`) to estimate
#' the change in a variable as a function of its initial value
#' The function can restrict predictions to a range between a specified `min` and `max`
#' (e.g., targets), or default to the observed data range if these are not specified.
#'
#' @param data_model A data frame with the training data. Obtained through `prep_data()`. Must contain:
#'   - `initialvalue`: numeric, the initial level of the indicator.
#'   - `change`: numeric, the change in the indicator.
#'   - `fold_id`: numeric or factor, used for cross-validation.
#' @param min Optional. Minimum value of `initialvalue` to predict.
#' @param max Optional. Maximum value of `initialvalue` to predict.
#' @param granularity Numeric. Granularity in outcome variable. Default is `0.1`.
#'   If `NULL`, predictions are unrestricted on the upper end.
#' @param verbose. Logical. If TRUE, display messages in console. Default is TRUE
#'
#' @return A data frame with:
#'   - `initialvalue`: sequence of values used for prediction.
#'   - `changes_speed`: predicted change.
#'
predict_speed <- function(data_model,
                          min         = NULL,
                          max         = NULL,
                          granularity = 0.1,
                          verbose     = TRUE) {

  # Early exit if no data
  if (nrow(data_model) == 0L) {

    cli::cli_alert_warning("Input data_model is empty. Returning NA data.table.")

    return(data.table(y      = NA_real_,
                      change = NA_real_))
  }

  # Fit model
  fit_speed <- gcrq(change ~ ps(initialvalue,
                                lambda = 0.1 * 1.148^seq(0, 50, by=2)),
                    foldid = data_model$fold_id,
                    tau    = 0.5,
                    data   = data_model)

  # Generate prediction sequence
  x_seq <- seq(min,
               max,
               by = granularity)

  # Predict (charts() returns a named vector or data.frame)
  # Predict changes over the grid
  # changes_speed <- charts(fit_speed,
  #                         k = x_seq)
  #
  # # Updating to change with granularity
  # predictions_dt <- data.table(
  #   y      = round(x_seq / granularity) * granularity,
  #   change = round(as.numeric(changes_speed) / granularity) * granularity
  # )

  # changes_speed <- charts(fit_speed,
  #                         k = x_seq,
  #                         digits=-log10(granularity))
  
  changes_speed <- charts(fit_speed, k = x_seq)

  #fittedvalues <- charts(fit, k=seq(min,max,granularity),digits=-log10(granularity))

  # predictions_dt <- data.table(
  #   y      = round(x_seq / granularity) * granularity,
  #   change = round(as.numeric(changes_speed) / granularity) * granularity
  # )

  predictions_dt <- data.table(
    y      = round(x_seq / granularity) * granularity,
    change = as.numeric(changes_speed)
  )


  if (verbose) cli::cli_alert_success("Changes speed successfully calculated")

  return(predictions_dt)
}


#' Generate Speed Path from Predicted Changes
#'
#' Computes the expected trajectory (or "speed path") of a variable over time,
#' based on predicted annualized changes as a function of the variable's initial value.
#' This is typically used to simulate progress from the worst to the best value
#' of an indicator
#'
#'
#' @param changes_speed A data.table containing two columns:
#'   \code{y} (the starting level of the indicator) and
#'   \code{change} (the expected annualized change from that level).
#' @param granularity A numeric value indicating the step size used when
#'   computing the trajectory over time. Smaller values produce smoother curves.
#'   Defaults to \code{0.1}.
#' @param best A character string indicating whether higher or lower values of thei indicator are better. Use "high" if higher values of the indicator reflect better outcomes, or "low" if lower values are better.
#' @param verbose Logical. If \code{TRUE} (default), a success message is displayed after computation.
#'
#' @return A \code{data.table} with two columns:
#'   \item{time}{Cumulative time (in years) needed to reach each step in the path.}
#'   \item{y}{The corresponding value of the indicator at each time step.}
#'
#'
#'
#' @export
get_speed_path <- function(changes_speed,
                           granularity = 0.1,
                           best        = "high",
                           verbose     = TRUE) {

  # Ensure changes_speed is a data.table
  path_speed <- copy(changes_speed)  # don't overwrite original

  if (best == "low") {

    setorder(path_speed,
             -y)       # reorder rows by descending y

    path_speed[, change := -change] # mutate/change 'change' column by negating its values


  } else {
      setorder(path_speed,
        #changes_speed,
               y)       # reorder rows by descending y
  }

  # Step 1: Rename 'change' to 'time'
  setnames(path_speed, "change", "time")

# filter out zero and negative time increments (wrong-direction changes)
 path_speed <- path_speed[time > 0]

  # Step 2: Compute lagged time and cumulative transformed time

  path_speed[, time := {
    ltime        <- shift(time,
                          type = "lag")                      # lag of 'time'
    time_new     <- 1 / ltime * granularity                     # transformed time
    time_new[1L] <- 0                                       # first value is 0
    cumsum(time_new)                                        # cumulative sum
  }]


  # Step 3: Filter out non-finite values
  path_speed <- path_speed[is.finite(time)]

  # Step 4: Keep only 'time' and 'y' columns
  path_speed <- path_speed[,
                           .(y, time)]

  # time is cumulative time it takes to reach the value of y from the starting point
  # example:  If you start at y = 0.30, and follow the speed of progress described in changes_speed,
  #           it will take 2 years to reach y = 0.40, 4.5 years to reach y = 0.50, etc.

  if (verbose) cli::cli_alert_success("Path speed successfully calculated")

  return(path_speed)

}

#' Predict percentiles
#'
#' This function calculates the changes as a function of initial values by percentile
#'
#' @inheritParams predict_speed
#' @param verbose Logical. If TRUE, prints info messages in console
#' @param sequence_pctl Numeric vector of percentile paths to calculate.
#'                      The percentile score is determined by the granularity of the percentiles chosen here. For example, if seq(20,80,20) is chosen, then a country’s progress will fall in the 0th-20th percentile, 20th-40th percentile etc.
#'                      This argument is only relevant if `percentiles=TRUE`.
#' @return A data frame with calculated changes by percentiles.
#' @export
predict_pctls <- function(data_model,
                          min           = NULL,
                          max           = NULL,
                          granularity   = 0.1,
                          sequence_pctl = NULL,
                          verbose       = TRUE) {

  # __________________________________ #
  # Validate input ~~~~~ ####
  # __________________________________ #

  # Prediction grid
  x_seq <- seq(min,
               max,
               by = granularity)

  # Uses cross-validation to find the optimal smoothing of percentile-curves.
  # gcrq automatically ensures that the percentile-curves do not cross

  fit_pctl <- gcrq(change ~ ps(initialvalue,
                               lambda = 0.1 * 1.148^seq(0, 50, by=2)),
                   foldid = data_model$fold_id,
                   tau    = sequence_pctl / 100,
                   data   = data_model)


  # Predict changes at each initial value for each percentile
  changes_mat <- charts(fit_pctl,
                        k = x_seq)  # matrix with cols = percentiles

  # Convert to data.table and reshape long
  changes_dt <- qDT(changes_mat)

  changes_dt[, y := round(x_seq / granularity) * granularity]

  changes_dt <- melt(changes_dt,
                     id.vars       = "y",
                     variable.name = "pctl",
                     value.name    = "change")

  # Convert factor or character percentile labels to numeric percentiles
  changes_dt[,
             pctl := as.numeric(gsub("[^0-9.]", "", pctl)) * 100]

  # Optional: sort for nicer output
  setorder(changes_dt, y, pctl)

  if (verbose) {
    cli::cli_alert_success("Changes as function of initial values by percentiles successfully calculated")
  }

  return(changes_dt)

}



#' Predict Changes in an indicator Over Time
#'
#' This function calculates predicted changes in an indicator based on its initial value. Two methods are available: `speed` and `percentiles`.
#'
#' @param data A data frame or list representing the input data model to be used for the prediction.
#' This should be the output of the [`prep_data()`] function, which prepares the data in the required format.
#' @param speed Logical. If TRUE, calculate speed of progress scores. Default is FALSE.
#' @param percentiles Logical. If TRUE, calculate percentile scores. Default is TRUE
#' @inheritParams predict_speed
#' @inheritParams predict_pctls
#' @inheritParams get_speed_path
#' @return A list containing prediction results. The content depends on the selected method:
#' \describe{
#'   \item{`changes_speed`}{Predicted changes using the `"speed"` method (if selected).}
#'   \item{`path_speed`}{The path or trajectory of predictions based on the `"speed"` method (if selected).}
#'   \item{`changes_pctls`}{Predicted changes using the `"percentiles"` method (if selected).}
#' }
#'
#' @details
#' - On The `"speed"` method: for indicators without overall progress, the speed of progress cannot be computed.
#' - The `"percentiles"` method can be computed for all indicators but ideally require rich training data.
#'
#' @seealso [predict_speed()], [get_speed_path()], [predict_pctls()]
#'
#' @export
predict_changes <- function(data,
                            min           = NULL,
                            max           = NULL,
                            granularity   = 0.1,
                            best          = "high",
                            speed         = FALSE,
                            percentiles   = TRUE,
                            sequence_pctl = NULL,
                            support       = 1,
                            verbose       = TRUE) {

  # Add validations on inputs

  # At least one between speed and percentiles should be true

  if(speed == FALSE && percentiles == FALSE) {

    cli::cli_abort("At least one between speed and percentiles should be TRUE")
  }

  if (percentiles == TRUE & is.null(sequence_pctl)) {

    cli::cli_abort("If percentiles method is selected, need to provide pctls sequence")
  }

  # Validate 'best'
  if (!best %in% c("high", "low")) {

    cli::cli_abort("`best` must be either 'high' or 'low'")

  }

  if (verbose) {cli::cli_alert_info("Speed method for computations: {.strong {speed}}")}
  if (verbose) {cli::cli_alert_info("Percentiles method for computations: {.strong {percentiles}}")}


  res_list <- list()

  if (speed) {

    # Predict speed
    changes_speed <- predict_speed(data_model  = data,
                                   min         = min,
                                   max         = max,
                                   granularity = granularity,
                                   verbose     = verbose)

    # Get speed paths
    path_speed <- get_speed_path(changes_speed     = changes_speed,
                                 granularity       = granularity,
                                 best              = best,
                                 verbose           = verbose)

    res_list$changes_speed <- changes_speed
    res_list$path_speed    <- path_speed

  }

  if (percentiles) {

    changes_pctl <- predict_pctls(data_model    = data,
                                  min           = min,
                                  max           = max,
                                  granularity   = granularity,
                                  sequence_pctl = sequence_pctl,
                                  verbose       = verbose)

    res_list$changes_pctl <- changes_pctl
  }

  return(
    invisible(res_list)
  )

}

