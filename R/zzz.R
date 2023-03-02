# zzz.R
# Contains the package documentation



#' causalKNN
#'
#' The causalKNN package provides causal KNN regression and treatment effect projection
#'     algorithms following Hitsch, Misra, and Zhang (2023),
#'     Heterogeneous Treatment Effects and Optimal Targeting Policy Evaluation,
#'     <https://dx.doi.org/10.2139/ssrn.3111957>
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
