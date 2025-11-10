#' Split a string
#'
#' @param x A character vector with one element
#' @param split What to split on
#'
#' @return A character vector
#' @export
#'
#' @examples
#' x <- "alpha,beta,charlie"
#' strsplit1(x, ",")
strsplit1 <- function(x, split){
  strsplit(x, split = split)[[1]]
}
