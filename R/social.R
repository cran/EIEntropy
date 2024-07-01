#' Randomly Generated Data
#'
#' This dataset contains 200 observations of 6 variables. The data was
#' generated randomly for the purpose of exemplifying a database that could
#' potentially be used with this function.
#'
#'
#' @format A data frame with 200 observations of 6 variables called social
#' @keywords dataset example
#' @return A data frame containing the loaded "social" data from the .rds file.
#' \itemize{
#'   \item \code{Dcollege}: Dummy variable indicating college education.
#'   \item \code{Dunemp}: Dummy variable indicating unemployment.
#'   \item \code{Totalincome}: Total income of each observation.
#'   \item \code{reg}:Variable indicating the region of the observation.
#'   \item \code{w}: Weights for each observation.
#'   \item \code{n}: Identifier for observations.
#' }
#' #' @examples
#' data(social)
#' head(social)
#'

#' @export
social<- function() {
  file_path <- system.file("extdata", "social.rds", package = "EIEntropy")
  social_data <- readRDS(file_path)
  return(social_data)
}







