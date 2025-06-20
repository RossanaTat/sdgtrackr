% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_data.R
\name{prep_data}
\alias{prep_data}
\title{Prepares data for analysis}
\usage{
prep_data(
  indicator = "EG.ELC.ACCS.ZS",
  data = wbstats::wb_data(indicator = indicator, lang = "en", country = "countries_only"),
  startyear_data = 2000,
  floor = 0,
  ceiling = 100,
  granularity = 0.1,
  code_col = "iso3c",
  year_col = "date",
  verbose = TRUE
)
}
\arguments{
\item{indicator}{Character. Indicator code (e.g., `"EG.ELC.ACCS.ZS"`). Defaults to access to electricity.}

\item{data}{Optional. A data frame with indicator data. If NULL, data is downloaded via `wbstats::wb_data()`.}

\item{startyear_data}{Integer. Minimum year to include in the data. Defaults to 2000.}

\item{code_col}{Character. Name of the column with country codes. Defaults to `"iso3c"`.}

\item{year_col}{Character. Name of the column with years. Defaults to `"date"`.}

\item{verbose}{Logical. If TRUE print messages in console. Default is TRUE}
}
\value{
A `list` with 3 elements: 1. data prepared for estimation, 2. min and 3. max. Min and Max are range limits for expected changes, based on floor/ceiling if provided, otherwise on observed values. Rounded to nearest granularity.
}
\description{
Prepares indicator data for estimation by computing annualized changes over 5–10 year periods, selecting the shortest available spell, balancing countries by number of rows, and assigning fold IDs for cross-validation.
}
