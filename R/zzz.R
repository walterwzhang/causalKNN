# zzz.R
# Contains the package documentation



#' causalKNN
#'
#' The causalKNN package provides causal KNN regression and treatment effect projection
#'     algorithms following Hitsch, Misra, and Zhang (2023)
#'
#' @docType package
#' @author Walter Zhang <walterwzhang@chicagobooth.edu>
#' @import Rcpp data.table glmnet
#' @importFrom Rcpp evalCpp
#' @importFrom KernelKnn knn.index.dist
#' @useDynLib causalKNN, .registration = TRUE
#' @name causalKNN
#' @noRd
NULL
