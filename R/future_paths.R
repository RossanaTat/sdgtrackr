#' Prepare data for future target projections
#'
#'
#' @inheritParams prep_data
#' @return A `data.table` with the latest observation per country, including:
#' \describe{
#'   \item{code}{Standardized country code (from `code_col`)}
#'   \item{year}{Year of the latest available observation}
#'   \item{y}{Value of the indicator}
#'   \item{y_fut}{Rounded value of `y`, based on the specified granularity}
#' }
#' @export
prep_data_fut <- function(data               = NULL,
                          indicator          = "EG.ELC.ACCS.ZS",
                          granularity        = 0.1,
                          code_col           = "iso3c",
                          year_col           = "date",
                          support            = 1,
                          extreme_percentile = getOption("trackr.extreme_pctl"),
                          startyear_data     = 2000,
                          endyear_data       = 2023,
                          verbose            = TRUE) {

  # Convert to data.table
  dt <- qDT(data)

  # Standardize column names
  data.table::setnames(dt,
                       old         = c(code_col, year_col, indicator),
                       new         = c("code", "year", "y"),
                       skip_absent = FALSE)

  # Keep relevant columns and latest observation per country
  # Filter to intended future horizon

  dt <- dt[!is.na(y), .(code, year, y)]
  dt <- dt[year >= startyear_data & year <= endyear_data]

  setorder(dt, code, year)

  dt <- dt[, .SD[.N], by = code]

  if (nrow(dt) == 0) {
    cli::cli_warn("No data remains within startyear_data and endyear_data window.")
    return(NULL)
  }


  # ____________________________
  # Tail-aware support filtering on raw y
  # ____________________________

  # Compute percentiles of y
  # --------------------------------------
  # Tail cutoff logic: support asymmetric percentiles ####
  # --------------------------------------
  if (length(extreme_percentile) == 1) {
    lower_pct <- extreme_percentile
    upper_pct <- extreme_percentile
  } else if (length(extreme_percentile) == 2) {
    lower_pct <- extreme_percentile[1]
    upper_pct <- extreme_percentile[2]
  } else {
    cli::cli_abort("extreme_percentile must be a numeric vector of length 1 or 2.")
  }

  # Compute percentile-based cutoffs
  lower_cutoff <- quantile(dt$y, probs = lower_pct, na.rm = TRUE)
  upper_cutoff <- quantile(dt$y, probs = 1 - upper_pct, na.rm = TRUE)



  # Compute bin on y for support checking
  dt[, y_bin := round(y / granularity) * granularity]

  # Support counts by raw y bin
  support_table <- dt[, .(n_countries = uniqueN(code)), by = y_bin]
  support_table[, region := fifelse(y_bin < lower_cutoff, "low",
                                    fifelse(y_bin > upper_cutoff, "high", "middle"))]
  support_table[, keep := region == "middle" | n_countries >= support]

  supported_bins <- support_table[keep == TRUE, y_bin]

  # Filter based on support in y bins
  dt <- dt[y_bin %in% supported_bins]

  # Now compute y_fut from filtered y
  dt[, y_fut := round(y / granularity) * granularity]
  dt[, y_bin := NULL]  # clean up

  if (verbose && support >= 1) {
    cli::cli_alert_info("{length(supported_bins)} starting levels retained after applying support >= {support} to tails.")
  }

  return(dt)
}


# -------------------- #
# Percentiles ~~ ####
# -------------------- #

#' Generate Future Percentile Paths
#'
#' Computes projected future paths for given percentiles based on initial levels and predicted changes.
#'
#' @param data_fut A data frame containing historical data with columns `code`, `year`, and `y`. This serves as the baseline for projections.
#' @param target_year An integer specifying the year to project to. Default is '2030'. Only relevant if `future = TRUE`.
#' @inheritParams path_historical
#' @param verbose Logical; if `TRUE`, messages will be printed for each processing year. Default is `TRUE`.
#'
#' @return A data.table with the projected values (`y_fut`) by `code`, `year`, and `pctl`, including observed values where available.
#'
#' @examples
#' \dontrun{
#' data_fut <- data.frame(code = c("A", "A"), year = c(2020, 2021), y = c(50, 52))
#' predictions_pctl <- data.frame(initialvalue = 52, pctl = 20, change = 1)
#' future_path_pctls(data_fut = data_fut, predictions_pctl = predictions_pctl, target_year = 2025)
#' }
#'
#' @export
future_path_pctls <- function(data_fut,
                              sequence_pctl     = seq(20,80,20),
                              changes_pctl,
                              target_year = 2030,
                              granularity = 0.1,
                              min         = NULL,
                              max         = NULL,
                              verbose     = TRUE) {

  # Create a new dataset which will contain the predicted path from the last observation to targetyear at all selected percentiles

  path_fut_pctl <- expand.grid(code = unique(data_fut$code),
                               year = seq(min(data_fut$year),
                                        target_year,
                                        1),
                               pctl = sequence_pctl) |>
    # Merge in the actual data from WDI
    joyn::joyn(data_fut,
               by = c("code","year"),
               match_type ="m:1",
               keep = "left",
               reportvar = FALSE,
               verbose = FALSE)

  # Year-by-year, calculate the percentile paths for the last value observed. For future target creation.

  cli::cli_alert_info("Calculating future percentile paths")

  startyear_target = min(path_fut_pctl$year)

  n = startyear_target - min(path_fut_pctl$year) + 1

  # Continue this iteratively until the target year
  while (n + min(path_fut_pctl$year) - 1 <= target_year) {
    # Year processed concurrently
    path_fut_pctl <- path_fut_pctl |>
      # Merge in data with predicted changes based on initial levels
      joyn::joyn(changes_pctl,
                 match_type = "m:1",
                 keep = "left",
                 by = c("y_fut = y",
                        "pctl"),
                 reportvar = FALSE,
                 verbose = FALSE) |>
      group_by(code,
               pctl) |>
      arrange(year) |>
      # Calculate new level based on the predicted changes.
      mutate(y_pctl = if_else(row_number() == n & !is.na(lag(y_fut)),
                        round((lag(y_fut) + lag(change)) / granularity) * granularity,
                        y_fut),
       y_fut = dplyr::coalesce(y_pctl, y_fut)) |>
      ungroup() |>
      select(-change)

    n=n+1
  }

  path_fut_pctl <- qDT(path_fut_pctl)[
    !is.na(y_pctl) & (is.na(y) | (y >= min & y <= max))
  ][
    , y_fut := NULL
  ][
    , setorder(.SD, code, year)
  ]

  # Return ####
  return(path_fut_pctl)
}

# -------------------- #
# Speed ~~ ####
# -------------------- #

#' Project future indicator paths at different speeds
#'
#' This function simulates future indicator trajectories under different speeds of progress.
#' It generates yearly projections from the latest available data point up to a given target year.
#'
#' @param data_fut A data frame or data table containing the historical and/or
#'   baseline data. Must include columns: `code`, `year`, and `y` (the indicator),
#'   and optionally `y_fut` and `time` for projections.
#' @param sequence_speed Vector of numeric speeds to simulate.
#' @param path_speed A data frame or table containing the modeled speed paths.
#' @param best Character: "high" or "low", the direction of progress.
#' @param target_year Numeric, final year for projections.
#' @param min Optional numeric lower bound for indicator values.
#' @param max Optional numeric upper bound for indicator values.
#'
#' @return A data table with projected indicator values (`y_fut`) for each
#'   country-year-speed combination from the last available year through
#'   the `target_year`.
#' @export
future_path_speed <- function(data_fut,
                              sequence_speed = c(0.25, 0.5, 1, 2, 4),
                              path_speed,
                              best = "high",
                              target_year = 2030,
                              min = NULL,
                              max = NULL) {

  # Step 1: Create country-year-speed grid
  data_fut <- expand.grid(
    code = unique(data_fut$code),
    year = seq(min(data_fut$year), target_year, 1),
    speed = sequence_speed
  ) |>
    rename("yeartemp" = "year") |>
    joyn::joyn(
      data_fut,
      match_type = "m:1",
      by = "code",
      reportvar = FALSE,
      verbose = FALSE,
      y_vars_to_keep = "year"
    ) |>
    filter(yeartemp >= year) |>
    select(-year) |>
    rename("year" = "yeartemp") |>
    joyn::joyn(
      data_fut,
      match_type = "m:1",
      by = c("code", "year"),
      reportvar = FALSE,
      verbose = FALSE
    )

  # Step 2: Prepare dataset for projected paths
  path_fut_speed <- data_fut |>
    select(-y) |>
    filter(!is.na(y_fut)) |>
    dplyr::cross_join(path_speed) |>
    mutate(best = best) |>
    filter(if_else(best == "high",
                   y_fut <= y,
                   y_fut >= y))

  # Step 3: Debug / safety check
  #message("nrow(path_fut_speed) before group_by: ", nrow(path_fut_speed))
  if (nrow(path_fut_speed) == 0) {
    warning("No valid future speed paths available after filtering. Returning empty table.")
    return(data.table(
      code = character(),
      year = numeric(),
      speed = numeric(),
      y_fut = numeric()
    ))
  }

  # Step 4: Build projected paths
  path_fut_speed <- path_fut_speed |>
    group_by(code, speed) |>
    arrange(time) |>
    mutate(year = year + (time - time[1]) / speed) |>
    ungroup() |>
    select(-c(y_fut, time, best)) |>
    rename("y_speed" = "y") |>
    joyn::joyn(
      data_fut,
      match_type = "1:1",
      by = c("code", "year", "speed"),
      reportvar = FALSE,
      verbose = FALSE,
      y_vars_to_keep = "y"
    ) |>
    group_by(code, speed) |>
    arrange(year) |>
    mutate(y_fut = zoo::na.approx(y_speed,
                                  year,
                                  na.rm = FALSE,
                                  rule = 2)) |>
    filter(year %in% seq(min(year), target_year, 1)) |>
    ungroup() |>
    select(-y_speed, -y)

  # Step 5: Apply bounds and ordering
  path_fut_speed <- as.data.table(path_fut_speed)[
    is.na(y_fut) | (y_fut >= min & y_fut <= max)
  ][order(code, year)]

  return(path_fut_speed)
}


#' Future paths. Wrapper to Compute future paths based on either percentiles or speeds method
#'
#' Computes future trajectories for each country either based on percentile projections
#' (`future_path_pctls()`) or speed of progress (`future_path_speed()`), depending on user input.
#'
#' @inheritParams future_path_pctls
#' @inheritParams future_path_speed
#' @param speed Logical. If `TRUE`, calls `future_path_speed()` and always augments
#'   projections with historically observed speeds of progress.
#' @param percentiles Logical. If `TRUE`, calls `future_path_pctls()`.
#'
#' @return A list containing projected future paths based on selected method(s).
#'
#' @seealso [future_path_pctls()], [future_path_speed()], [path_future_his_speed()]
#' @export
future_path <- function(data_fut,
                        target_year      = 2030,
                        min              = 0,
                        max              = 100,
                        granularity      = 0.1,
                        sequence_pctl    = seq(20, 80, 20),
                        changes_pctl     = NULL,
                        sequence_speed   = c(0.25, 0.5, 1, 2, 4),
                        path_speed       = NULL,
                        best             = "high",
                        verbose          = TRUE,
                        speed            = FALSE,
                        percentiles      = FALSE,
                        scores           = NULL) {

  result <- list(pctl = NULL, speed = NULL)

  # --------------------------------------------------
  # Percentile-based future paths
  # --------------------------------------------------

  if (percentiles) {
    if (is.null(changes_pctl)) {
      cli::cli_abort("`changes_pctl` must be provided when `percentiles = TRUE`.")
    }

    result$pctl <- future_path_pctls(
      data_fut         = data_fut,
      sequence_pctl    = sequence_pctl,
      changes_pctl     = changes_pctl,
      target_year      = target_year,
      granularity      = granularity,
      min              = min,
      max              = max,
      verbose          = verbose
    )
  }

  # --------------------------------------------------
  # Speed-based future paths (model + historical)
  # --------------------------------------------------

  if (speed) {
    if (is.null(path_speed)) {
      cli::cli_abort("`path_speed` must be provided when `speed = TRUE`.")
    }
    if (is.null(scores)) {
      cli::cli_abort("`scores` must be provided when `speed = TRUE` to compute historical-speed projections.")
    }

    # --- Model-based speed projections
    path_speed_model <- future_path_speed(
      data_fut       = data_fut,
      sequence_speed = sequence_speed,
      path_speed     = path_speed,
      best           = best,
      target_year    = target_year,
      min            = min,
      max            = max
    )

    if (!is.null(path_speed_model) && nrow(path_speed_model) > 0) {
      path_speed_model[, speed_source := "model"]
    }

    # --- Historical-speed projections
    path_speed_his <- path_future_his_speed(
      data_fut    = data_fut,
      scores      = scores,
      path_speed  = path_speed,
      best        = best,
      target_year = target_year,
      min         = min,
      max         = max
    )

    # --- Combine
    result$speed <- rbind(
      path_speed_model,
      path_speed_his,
      fill = TRUE
    )

    # --- Preserve diagnostics
    attr(result$speed, "diagnostics") <- attr(path_speed_his, "diagnostics")
  }

  if (is.null(result$pctl) && is.null(result$speed)) {
    cli::cli_abort(
      "At least one of `percentiles = TRUE` or `speed = TRUE` must be specified."
    )
  }

  return(result)
}


# --------------------------- #
# Hist. Speed Scores ~~ ####
# --------------------------- #


#' Project future indicator paths using historical country-specific speeds
#'
#' This function simulates future indicator trajectories by extending each
#' country along the global speed path using its own historically observed
#' speed score.
#'
#' Missing values (`NA`) in projected outcomes indicate that a historically
#' anchored future path could not be computed, typically due to missing or
#' non-positive historical speed estimates or infeasible extrapolation beyond
#' the modeled indicator range.
#'
#' @param data_fut Output of `prep_data_fut()`
#' @param scores A data.table containing historical speed scores.
#'   Must include columns: `code`, `year`, and `score`
#' @param path_speed Output of `get_speed_path()`
#' @param best Indicator direction ("high" or "low")
#' @param target_year Projection horizon
#' @param min Optional lower bound on projected values
#' @param max Optional upper bound on projected values
#'
#' @return A data.table with projected indicator values (`y_fut`) by country-year.
#'   A diagnostics table explaining missing projections is attached as an
#'   attribute: `attr(x, "diagnostics")`.
#'
#' @export
path_future_his_speed <- function(data_fut,
                                  scores,
                                  path_speed,
                                  best        = "high",
                                  target_year = 2030,
                                  min         = NULL,
                                  max         = NULL) {

  # --------------------------------------------------
  # 1. Keep latest valid speed score per country
  # --------------------------------------------------

  setorder(scores, code, -year)
  scores_latest <- scores[, .SD[1], by = code]

  scores_latest[, has_valid_speed :=
                  is.finite(score) & score > 0]

  # --------------------------------------------------
  # 2. Create country-year grid for projections
  # --------------------------------------------------

  data_fut <- expand.grid(
    code = unique(data_fut$code),
    year = seq(min(data_fut$year), target_year, 1)
  ) |>
    rename("yeartemp" = "year") |>
    joyn::joyn(
      data_fut,
      match_type     = "m:1",
      by             = "code",
      reportvar      = FALSE,
      verbose        = FALSE,
      y_vars_to_keep = "year"
    ) |>
    filter(yeartemp >= year) |>
    select(-year) |>
    rename("year" = "yeartemp") |>
    joyn::joyn(
      data_fut,
      match_type = "m:1",
      by = c("code", "year"),
      reportvar = FALSE,
      verbose   = FALSE
    )

  # --------------------------------------------------
  # 3. Attach historical speed
  # --------------------------------------------------

  data_fut <- joyn::joyn(
    data_fut,
    scores_latest[, .(code, speed = score, has_valid_speed)],
    match_type = "m:1",
    by = "code",
    reportvar = FALSE,
    verbose   = FALSE
  )

  # --------------------------------------------------
  # 4. Build future path using historical speed
  # --------------------------------------------------

  path_fut_speed <- data_fut |>
    select(-y) |>
    filter(!is.na(y_fut), has_valid_speed) |>
    cross_join(path_speed) |>
    mutate(best = best) |>
    filter(if_else(best == "high",
                   y >= y_fut,
                   y <= y_fut)) |>
    group_by(code) |>
    arrange(time) |>
    mutate(year = year + (time - time[1]) / speed) |>
    ungroup() |>
    select(-c(y_fut, time, best, speed, has_valid_speed)) |>
    rename("y_speed" = "y") |>
    joyn::joyn(
      data_fut,
      match_type = "1:1",
      by = c("code", "year"),
      reportvar = FALSE,
      verbose   = FALSE,
      y_vars_to_keep = "y"
    ) |>
    group_by(code) |>
    arrange(year) |>
    mutate(
      y_fut = zoo::na.approx(
        y_speed,
        year,
        na.rm = FALSE,
        rule = 2
      )
    ) |>
    filter(year %in% seq(min(year), target_year, 1)) |>
    ungroup() |>
    select(-y_speed, -y)

  # --------------------------------------------------
  # 5. Apply bounds and ordering
  # --------------------------------------------------

  path_fut_speed <- as.data.table(path_fut_speed)[
    is.na(y_fut) | (y_fut >= min & y_fut <= max)
  ][
    order(code, year)
  ]

  # --------------------------------------------------
  # 6. Add schema alignment for merging with model-based paths
  # --------------------------------------------------

  path_fut_speed[, `:=`(
    speed        = NA_real_,
    speed_source = "historical"
  )]

  # --------------------------------------------------
  # 7. Build diagnostics for missing projections
  # --------------------------------------------------

  diagnostics <- path_fut_speed[, .(
    has_any_projection = any(!is.na(y_fut)),
    all_na             = all(is.na(y_fut))
  ), by = code]

  diagnostics <- joyn::joyn(
    diagnostics,
    scores_latest[, .(code, has_valid_speed)],
    match_type = "m:1",
    by = "code",
    reportvar = FALSE,
    verbose   = FALSE
  )

  diagnostics[, reason := fifelse(
    !has_valid_speed,
    "no_valid_historical_speed",
    fifelse(all_na,
            "infeasible_or_outside_speed_path",
            NA_character_)
  )]

  # Attach diagnostics as attribute
  attr(path_fut_speed, "diagnostics") <- diagnostics

  return(path_fut_speed)
}



## Helper ####
#' Diagnose speed-path domain violations
#'
#' Identifies countries whose starting values lie outside the
#' speed-path domain and therefore cannot be projected.
#'
#' @param data_fut Output of `prep_data_fut()`
#' @param path_speed Output of `get_speed_path()`
#'
#' @return A data.table with one row per country and a diagnostic flag
#' @keywords internal
diagnose_speed_domain <- function(data_fut, path_speed) {

  data_fut  <- as.data.table(data_fut)
  path_min  <- min(path_speed$y, na.rm = TRUE)
  path_max  <- max(path_speed$y, na.rm = TRUE)

  out <- data_fut[, .(
    code,
    year,
    y_fut,
    status = fifelse(
      y_fut < path_min, "below_speed_domain",
      fifelse(y_fut > path_max, "above_speed_domain", "inside_speed_domain")
    )
  )]

  return(out)
}
## Check which countries would be undefined BEFORE projecting
#diagnose_speed_domain(data_fut, path_speed)[status != "inside_speed_domain"]

# Summary table of violations
# diagnose_speed_domain(data_fut, path_speed)[,
#                                            .N, by = status
# ]


