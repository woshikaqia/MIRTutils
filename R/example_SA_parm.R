#' Example standalone item parameter tables
#'
#' @description Example standalone item parameter dataframe for the \code{SA_parm}
#'
#' @docType data
#'
#' @usage data(example_SA_parm)
#'
#' @format An object of class data.frame with columns below\cr\cr
#'  a: slope \cr\cr
#'  b1, b2, ..., b_k: difficulty or step difficulty. \cr
#'  For GPCM items, b1, b2, ..., b_k are step difficulty parameters where k is the maximum possible score
#'                            of the item that has the largest number of score categories in the test (Score starts from 0).
#'                            For example, for two items with 3 and 4 categories (0,1,2 & 0,1,2,3), there should be column b1, b2 and b3, where
#'                            for the first item, b3 should be set to NA \cr
#'  For 3PL items, b1 is the difficulty parameter. All other b should be set to NA \cr\cr
#'  g: guessing parameters for 3PL items. Set to NA for GPCM items \cr\cr
#'  ItemID: Item IDs for SA items \cr\cr
#'  AssertionID: Assertion IDs for SA items \cr\cr
#'  ***\code{SA_parm} columns must follow the above order
"example_SA_parm"

