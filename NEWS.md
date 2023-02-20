# causalKNN 0.3

* Added direct causal KNN estimation and a wrapper function `causalKNN_direct` for off-the-shelf estimation. The back end `bootstrapIndexMatrixProcessorTest` function added for estimating direct causal KNN treatment effects for the test set. The package documentation and vignette has been updated to reflect new functionality

* Eliminated unneeded parameters in previous functions and renamed `tep_projection_enet` to `tep_enet` for documentation consistency

* Added bootstrap functionality for the Optimal K and causal KNN treatment effect estimation procedures. A bootstrap vignette has been added

# causalKNN 0.2

* Added Separate K among untreated NN and treated NN functionality with `knn_optimal_k_separate`

* Added a wrapper function `causalKNN_TEP` that wraps the keys functions from the last commit and provides an off-the-shelf estimation

* Added Travis-CI build support

* Added Unit Test functionality with the `testthat` package
    - Added a test to ensure `causalKNN_TEP` yields the same results as running the sequential functions
    - Added tests to ensure the separate KNN approach works and yields the same results when K is set to be equal for treated and untreated



# causalKNN 0.1

* Initialized package using `devtools` and documentation is built with `roxygen2`

* Added key functions `causalknn_treatment_effect`, `knn_index_mat`, `knn_optimal_k`, `tep_projection_enet` for a four step sequential estimation of the causal KNN regression and then a treatment effect projection with the Elastic Net

* Added helper functions `bootstrapIndexMatrixProcessor`, `standardizeVector`, `bootstrapKNNStepProcessor`, `cumWeightedSum`, and `cumWeightedSumStep`. These functions are documented with `@keyword internal` in accordance to `roxygen2` so they are still directly accessible to the user but hidden from manual

* Added a the `causalKNN-Intro` vignette

