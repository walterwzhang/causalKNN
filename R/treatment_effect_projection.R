# treatment_effect_projection.R
# Projects the Treatment Effect to the Validation Set
# Currently uses the Elastic-Net to project the treatment effect
# Key Functions:
#   - tep_enet



# Functions =======================================================================================

# tep_enet ----------------------------------------------------------------------------
#' Performs the TEP step for the Elastic-Net family of models
#'
#' We search over a grid of values for the elastic-net to find the optimal
#' alpha parameter. Set alpha to NULL to find the optimal alpha value.
#' The Elastic-Net family of regressions is estimated with `glmnet` package.
#' Parallelization is done with the `parallel` package.
#'
#' @param training_X A data frame of the training set features (data.frame)
#' @param training_TE A vector of the individually estimated treatment effects (numeric)
#' @param test_X A data frame of the test set features (data.frame)
#' @param alpha The alpha value for the elastic net (numeric)
#' @param alpha_seq The a vector alpha value to consider (numeric)
#' @param threads The value of number of threads to use (integer)
#' @param cv_folds The number of folds for the cross-validation step (integer)
#' @param est_seed A value for CV estimation seed (integer)
#' @param ... Additional parameters to be passed to \code{glmnet::cv.glmnet}
#' @return A vector of the TEP projection values for the test set (numeric)
#' @importFrom stats predict
#' @export

tep_enet   <-   function(training_X, training_TE, test_X,
                                     alpha       =   1,
                                     alpha_seq   =   seq(0, 1, by = 0.01),
                                     threads     =   1L,
                                     cv_folds    =   10L,
                                     est_seed    =   1234L,
                                     ...)
{
    # Checks
    if (threads < 0 | threads %% 1 != 0)
        stop("Incorrect threads specification")
    if (!requireNamespace("glmnet", quietly = TRUE))
        stop("Requires `glmnet` to be installed for the Elastic-Net.")
    if (!is.null(alpha))
    {
        if (alpha < 0 & alpha > 1)
            stop("Specified alpha value is out of range for the Elastic-Net.")
    }


    training_X_mat   <-   as.matrix(training_X)

    # Parallelization Pre-processing
    if (threads > 1L)
    {
        if ("parallel" %in% (.packages()))
        {
            parallel_flag   <-   TRUE
            cl              <-   parallel::makeCluster(threads)
            parallel::setDefaultCluster(cl)
            parallel::clusterExport(cl, list("training_X_mat", "training_TE", "alpha"), envir = environment())
            invisible(parallel::clusterEvalQ(cl, c(library("glmnet"), library("data.table"))))
        } else
        {
            warning("`parallel` package unavailable. Parallelization will not occur.", immediate. = TRUE)
            parallel_flag   <-   FALSE
        }

    } else
    {
        parallel_flag   <-   FALSE
    }


    # Training Step
    set.seed(est_seed)
    folds_index   <-   sample(1:cv_folds, nrow(training_X), replace = TRUE)

    if (is.null(alpha))
    {
        # Loop over alpha values to find the best value
        cv_alpha_flag    <-   TRUE

        if (!parallel_flag)
        {
            MSE_DT   <-   rbindlist(lapply(alpha_seq, function(alpha_value) # mclapply
            {
                ENet_fit_train   <-   cv.glmnet(x             =   training_X_mat,
                                                y             =   training_TE,
                                                parallel      =   FALSE,
                                                alpha         =   alpha_value,
                                                foldid        =   folds_index,
                                                ...)
                data.table(alpha = alpha_value, CV_Error = min(ENet_fit_train$cvm))
            }))
        } else
        {
            # Parallel search for alpha
            parallel::clusterExport(cl, list("alpha_seq", "folds_index"), envir = environment())
            MSE_DT   <-   rbindlist(parallel::parLapply(cl, alpha_seq, function(alpha_value)
            {
                ENet_fit_train   <-   cv.glmnet(x             =   training_X_mat,
                                                y             =   training_TE,
                                                parallel      =   FALSE,
                                                alpha         =   alpha_value,
                                                foldid        =   folds_index,
                                                ...)
                data.table(alpha = alpha_value, CV_Error = min(ENet_fit_train$cvm))
            }))
        }

        # Save optimal alpha value
        alpha   <-   MSE_DT[[1]][which.min(MSE_DT[[2]])]
    } else
    {
        cv_alpha_flag    <-   FALSE
    }

    # Estimate the Elastic-Net with pre-specified alpha values
    training_fit   <-   cv.glmnet(x             =   training_X_mat,
                                  y             =   training_TE,
                                  parallel      =   parallel_flag,
                                  nfolds        =   cv_folds,
                                  alpha         =   alpha,
                                  foldid        =   folds_index,
                                  ...)

    test_prediction   <-   predict(training_fit,
                                   newx   =   as.matrix(test_X),
                                   s      =   training_fit$lambda.min,
                                   type   =   'response'
                                   )[,1]


    # Parallelization Post-processing
    if (parallel_flag) parallel::stopCluster(cl)

    if (cv_alpha_flag)
    {
        return(list(test_TE        =   test_prediction,
                    alpha          =   alpha,
                    CV_alpha_MSE   =   MSE_DT))
                    # training_fit   =   training_fit))
    }
    return(list(test_TE        =   test_prediction,
                alpha          =   alpha))
                # training_fit   =   training_fit))
}

# -------------------------------------------------------------------------------------------------

# =================================================================================================
