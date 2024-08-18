#' Item score residual
#' @description
#' Compute item score residual (i.e., item score \code{-} expected item score. Item score for a testlet model item is
#' the summed raw score of all assertions in the testlet, and expected item score is
#' expected raw score calculated with Lord-Wingersky algorithm (See \code{?lord_wing} for details about the algorithm).
#'
#' @param theta a scalar or a vector of student ability
#' @param SA_dat For one student, a vector of response to standalone items.
#' For more than one student, a matrix or dataframe of response to standalone items. One assertion per column.
#' Column order must match row order in \code{SA_parm}. Use NA for missing responses
#' @param Cluster_dat For one student, a vector of response to cluster items.
#' For more than one student, a matrix or dataframe of response to cluster items. One assertion per column.
#' Column order must match row order in \code{Cluster_parm}. Use NA for missing responses.
#' @param SA_parm a matrix or dataframe of item parameters for standalone items, where columns are
#'  a (slope), b1, b2, ..., b_k (difficulty or step difficulty), g (guessing), ItemID, and AssertionID.
#'  Columns must follow the above order.
#'  See \code{example_SA_parm} for an example. Use \code{?example_SA_parm} for detailed column descriptions
#' @param Cluster_parm a matrix or dataframe of item parameters for cluster items, where columns are
#'  a (slope), b (difficulty), cluster variance, cluster position, ItemID, and AssertionID.
#'  Columns must follow the above order.
#'  See \code{example_Cluster_parm} for an example. Use \code{?example_Cluster_parm} for detailed column descriptions
#' @param Dv scaling factor for IRT model (usually 1 or 1.7)
#' @param n.nodes number of nodes used when integrating out the specific dimension
#' @param missing_as_incorrect by default, missings (NAs) are treated as missing; if TRUE, missings are treated as incorrect
#'
#' @note If the test does not have SA items or Cluster items, use default (NULL) for the corresponding data and parameter arguments \cr\cr
#'
#' @author Zhongtian Lin <lzt713@gmail.com>
#' @examples
#' data(example_SA_parm)
#' data(example_Cluster_parm)
#' sigma <- diag(c(1, sqrt(unique(example_Cluster_parm$cluster_var))))
#' mu <- rep(0, nrow(sigma))
#' thetas <- MASS::mvrnorm(7,mu,sigma)
#' thetas[,1] <- seq(-3,3,1) #overall dimension theta values
#' itmDat <- sim_data(thetas = thetas, SA_parm = example_SA_parm, Cluster_parm = example_Cluster_parm)
#' SA_dat <- itmDat[,1:20]
#' Cluster_dat <- itmDat[,-1:-20]
#' rst <- item_residual(thetas[,1], SA_dat, Cluster_dat, example_SA_parm, example_Cluster_parm, n.nodes = 11)
#' @export
item_residual = function(theta, SA_dat=NULL, Cluster_dat=NULL, SA_parm=NULL, Cluster_parm=NULL, Dv=1, n.nodes = 21, missing_as_incorrect = F) {
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")}
  if(is.null(SA_dat) & is.null(Cluster_dat)) {stop("No data found!!!")}

  if(!is.null(SA_dat)) SA_dat = matrix(as.matrix(SA_dat), nrow = length(theta))
  if(!is.null(Cluster_dat)) Cluster_dat = matrix(as.matrix(Cluster_dat), nrow = length(theta))
  combined_dat = cbind(SA_dat, Cluster_dat)
  if(any(rowSums(is.na(combined_dat)) == ncol(combined_dat))) {stop("one or more students did not respond to any item!!!")}
  if(missing_as_incorrect == T) combined_dat[is.na(combined_dat)] = 0  # recode if missing is treated as incorrect

  rst = utility(theta=theta, SA_parm=SA_parm, Cluster_parm=Cluster_parm, Dv=Dv, n.nodes = n.nodes,
                what = c("escore"))
  if (!is.null(Cluster_dat)) {
    Cluster_dat2 = as.data.frame(Cluster_dat)
    names(Cluster_dat2) = Cluster_parm$ItemID
    Cluster_dat2 = split.default(Cluster_dat2, names(Cluster_dat2))
    Cluster_dat2 = Cluster_dat2[as.character(unique(Cluster_parm$ItemID))]
    Cluster_dat2 = do.call(cbind, lapply(Cluster_dat2, rowSums))
    combined_dat = cbind(SA_dat, Cluster_dat2)
  }
  cbind(theta = rst$escore[,1], combined_dat - rst$escore[,-1])
}
