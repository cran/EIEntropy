#' Randomly Generated Data
#'
#' This dataset contains 100 observations of 6 variables. The data was
#' generated randomly for the purpose of exemplifying a database that could
#' potentially be used with this function.
#'
#' @docType data
#' @name financial
#' @format A data frame with 100 observations of 6 variables called financial
#' @keywords dataset example
#' @return A data frame containing the loaded "financial" data from the .rds file.
#' \itemize{
#'   \item \code{Dcollege}: Dummy variable indicating college education.
#'   \item \code{Dunemp}: Dummy variable indicating unemployment.
#'   \item \code{Totalincome}: Total income of each observation.
#'   \item \code{poor_liq}: Dummy variable indicating liquid poverty.
#'   \item \code{w}: Weights for each observation.
#'   \item \code{n}: Identifier for observations.
#' }
#'@examples
#' data(financial)
#' head(financial)
#'
#'
#' @export
financial <- function() {
  file_path <- system.file("extdata", "financial.rds", package = "EIEntropy")
  financial_data <- readRDS(file_path)
  return(financial_data)
}

