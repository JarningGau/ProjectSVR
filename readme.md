## Log

- `2022-07-12`

  `TODO` `Enhancement` Try to speed up `MergeCells()` function.

- `2022-07-11` 

  `Done` import `mlr3verse` properly

  `Done` Add normalization method for `FitEnsemblSVM()` `FitEnsemblMultiClassif()` `PredictNewdata()` and `ProjectNewdata()`

- `2022-06-11`

  `Enhancement` Plot for the VD results.

- `2022-06-09`

  `Done` `Enhancement` Normalize the feature.mat by cells (sum to 1) before model construction and prediction.

- `2022-05-30`

  `Done` `BUG` FitEnsembleSVM() | feature.mat should be a data.frame, when input a matrix. error occurs. Want a clear error information.

- `2022-05-28`

  `Enhancement` Two prediction function is too confused for new users. Try to merge them via generic method.

  `Enhancement` Add a plot function for visualizing the coords before and after refinement.

  `Done` `BUG` AddProjQual() function will generated duplicated results in @cellmeta slot when execute twice or more times. Want codes for colnames check.

- `2022-05-24`

  `Done` `LabelTransfer()` SVM model.

  ​	`FitEnsemblMultiClassif()`

  ​	`PredictNewdata()`

  ​	`MajorityVote()`

  `Done` Test these functions.

- `2022-05-23`

  `Done` Differential expression analysis `VarDecompose()` and `GetExpressionShift()`

  `TODO` `ClusterGenes()`

  `Done` Test `VarDecompose()` and `GetExpressionShift()`

- `2022-05-19`

  `Done` Add `CellProject`  object to store results generated during cell projection.

  `Done` New function: `RefineProjection`

- `2022-04-30`

  `Done` `ProjectNewdata()` add warning info and fill up 0 for not enough feature case. (case2 Ythdc2-KO data)

  ~~`TODO` Add a function to refine low confidence mapping coords.~~

  ~~`TODO` `LabelTransfer()` Add SVM method.~~

- `2022-03-27`

  `Done` ggplot2 in scProject package. `import ggplot2` is not proper for a package.

- `2022-03-26`

  `TODO` Add new signature scoring method -- [UCell](https://github.com/carmonalab/UCell). 

  `Done` RNA velocity on `scenario 2` task

- `2022-03-25`

  ~~`TODO` Differential expressed genes. `cacoa` package solved the problem.~~

- `2022-03-24`

  `Done` Optimize the mapQ function:

  1. Markov affinity matrix replaces the distance matrix to calculate the PCCs. `Done` `Not good`
  2. NN-based mean distance performs better. `Done`

- `2022-03-23`

  `Done` mTCA, levels' cellular annotation.

- `2022-03-22`

  `Done` Add two mapQ functions

  1. Global (random N) rather than local (KNN) topological similarity
  2. Jaccard index

- `2022-03-21`

  scenario2: time-course gonocyte-spermatogonia-transition

  `Done` optimized model training step. 

  1) Drop outliers cells before model training; `TODO`
  2) redo cNMF on cluster-based rather than grid-based pseudocells; `Done`
  3) The number of cells from different cell type should be balanced (balanced sampling replaces random sampling). `TODO`
  4) Only use protein coding genes. `Done` `Boost the results`

- `2022-03-20`

  `Done` Test `MakeMetaCells()`

  `TODO` Drop the `ml3verse` dependence. Use `e1071` for SVR directly.

- `2022-03-18`

  Speed up `MergeCells()`

  Add `LabelTransfer()` for cell type annotation after projection.

  Add `MakeMetaCell()` function for group transcriptome similar cells into meta cells to speed up cNMF.

  `TODO` Add `AddProjQaul()`, cluster-based map quality.

  `TODO` Implementation of DE analysis model.

- `2022-03-17`

  Add `NNHelper()` function for calculating nearest neightbors.

  Modifiy `EstimateKnnDensity()`. Test it again.  Drop @param cores.

  Implementation of `AddProjQual()`, evaluating projection/mapping quality.

- `2022-03-15`

  test cNMF.R worked

  implementation of `FitEnsembleSVM()`, `ComputeModuleScore()`

- `2022-03-14`

  `add function` merge by bins, or merge by mesh points. 

- `2022-03-12`

  `EstimateKnnDensity()`注意数字在R中默认为`float`，这在利用reticulate和python交互时容易出现类型错误

  函数中需要一个类型检查。
