#' @title Mean.in.log2space
#'
#' @description Applies the log2 to the mean of ((2^x - pseudo count) + pseudo
#' count).
#'
#' @param x Data
#' @param pseudo.count A pseudocount value
#'
#' @return Values
#'
#' @examples
#'
#' Mean.in.log2space(dataBulk, 0.1)
#'
#'@export Mean.in.log2space
#'
#'
#' @importFrom dplyr "%>%"
#'
#'
Mean.in.log2space=function(x,pseudo.count) {
  return(log2(mean(2^(x)-pseudo.count)+pseudo.count))
}

