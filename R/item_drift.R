#' Item drift analysis
#' @description Conduct item difficulty drift analysis and flag items with potential drift
#' Model currently supported are 1-3 PL, GPCM, and Rasch Testlet, and a mix of these models.
#'
#' @param theta a vector of examinee ability. The length of the vector should be equal to the number of examinees
#' @param SA_dat For one examinee, a vector of response to standalone items.
#' For more than one examinee, a matrix or dataframe of response to standalone items. One assertion per column.
#' Column order must match row order in \code{SA_parm}. Use NA for missing responses
#' @param Cluster_dat For one examinee, a vector of response to cluster items.
#' For more than one examinee, a matrix or dataframe of response to cluster items. One assertion per column.
#' Column order must match row order in \code{Cluster_parm}. Use NA for missing responses.
#' @param SA_parm a matrix or dataframe of item parameters for standalone items, where columns are
#'  a (slope), b1, b2, ..., b_k (difficulty or step difficulty), g (guessing), ItemID, and AssertionID.
#'  Columns must follow the above order.
#'  See \code{example_SA_parm} for an example. Use \code{?example_SA_parm} for detailed column descriptions
#' @param Cluster_parm a matrix or dataframe of item parameters for cluster items, where columns are
#'  a (slope), b (difficulty), cluster variance, cluster position, ItemID, and AssertionID.
#'  Columns must follow the above order.
#'  See \code{example_Cluster_parm} for an example. Use \code{?example_Cluster_parm} for detailed column descriptions
#' @param drift a numeric scalar for the amount of item difficulty parameter drift tested
#' @param Alpha a numeric vector for one or more nominal type I error rates for flagging aberrant responses.
#' @param Dv scaling factor for IRT model (usually 1 or 1.7)
#' @param n.nodes number of nodes used when integrating out the nuisance dimension
#' @param missing_as_incorrect by default, missings (NAs) are treated as missing; if TRUE, missings are treated as incorrect
#'
#' @note If the test does not have SA items or Cluster items, use default (NULL) for the corresponding data and parameter arguments. \cr\cr
#'
#' @return a list of
#' * \code{RR} = N (number of examinees) by n (number of assertions) dataframe of assertion residuals \cr
#' * \code{LRR} = N by n dataframe of the lower bound assertion residuals\cr
#' * \code{URR} = N by n dataframe of the upper bound assertion residuals\cr
#' If run for more than one examinees, returns an additional \code{drift_output} table in the list \cr
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
#' out_scoring <- list()
#' for (i in 1:nrow(itmDat)) {
#'   out_scoring[[i]] <- scoring(SA_dat[i,], Cluster_dat[i,], example_SA_parm, example_Cluster_parm, n.nodes = 11, SE=TRUE)
#' }
#' est_theta <- sapply(out_scoring, function(x) x$par)
#' rst <- item_drift(est_theta, SA_dat=SA_dat, Cluster_dat=Cluster_dat, SA_parm=example_SA_parm, Cluster_parm=example_Cluster_parm,
#'                   drift = 0.3, Alpha=0.05, Dv=1, n.nodes = 21, missing_as_incorrect = FALSE)
#' @export
item_drift <- function(theta, SA_dat=NULL, Cluster_dat=NULL, SA_parm=NULL, Cluster_parm=NULL, drift = 0.3, Alpha=0.05, Dv=1, n.nodes = 21, missing_as_incorrect = FALSE) {
  # -------------------------------------------------
  # ------------ Data Cleansing ---------------------
  # -------------------------------------------------
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")}
  if(is.null(SA_dat) & is.null(Cluster_dat)) {stop("No data found!!!")}
  if((is.null(SA_parm) | is.null(SA_dat)) && (is.null(Cluster_parm) | is.null(Cluster_dat)))  {stop("Data or Parameter file missing!!!")}

  if(!is.null(SA_dat)) {
    if(is.vector(SA_dat)) SA_dat = matrix(SA_dat, nrow = 1) else SA_dat = as.matrix(SA_dat)
    if (ncol(SA_dat) != nrow(SA_parm)) {stop("Number of items does not match between data and parameter file for standalone items")}
  }
  if(!is.null(Cluster_dat)) {
    if(is.vector(Cluster_dat)) Cluster_dat = matrix(Cluster_dat, nrow = 1) else Cluster_dat = as.matrix(Cluster_dat)
    if (ncol(Cluster_dat) != nrow(Cluster_parm)) {stop("Number of items does not match between data and parameter file for cluster items")}
  }
  if(!is.null(SA_dat) & !is.null(Cluster_dat) & nrow(Cluster_dat)!=nrow(SA_dat)) {stop("Cluster and SA data should have the same number of rows!")}

  combined_dat = cbind(SA_dat, Cluster_dat)
  if (any(rowSums(is.na(combined_dat)) == ncol(combined_dat))) {stop("one or more examinees did not respond to any item!!!")}
  if (nrow(combined_dat)!=1 & nrow(combined_dat)<30) warning("Number of examinees less than 30. Results may not be accurate!")

  if(!is.null(SA_parm)) {
    SA_parm = as.data.frame(SA_parm)
    names(SA_parm)[1] = c("a")
    names(SA_parm)[(ncol(SA_parm)-2):ncol(SA_parm)] = c("g", "ItemID", "AssertionID")
    names(SA_parm)[2:(ncol(SA_parm)-3)] = paste0("b",2:(ncol(SA_parm)-3)-1)
    colnames(SA_dat) = paste0(SA_parm$ItemID, "_", SA_parm$AssertionID)

  }
  if(!is.null(Cluster_parm)) {
    Cluster_parm = as.data.frame(Cluster_parm)
    colnames(Cluster_parm) = c("a","b","cluster_var","position", "ItemID", "AssertionID")
    Cluster_parm$position = match(Cluster_parm$position, sort(unique(Cluster_parm$position))) # reorder the position so it starts from 1
    colnames(Cluster_dat) = paste0(Cluster_parm$ItemID, "_", Cluster_parm$AssertionID)
  }
  # -------------------------------------------------
  # ------------ Calculate Residuals ----------------
  # -------------------------------------------------
  SA_parm_L <- SA_parm_U <- SA_parm
  SA_parm_L[,grepl("b",names(SA_parm))] <- SA_parm_L[,grepl("b",names(SA_parm))] - drift
  SA_parm_U[,grepl("b",names(SA_parm))] <- SA_parm_U[,grepl("b",names(SA_parm))] + drift

  Cluster_parm_L <- Cluster_parm_U <- Cluster_parm
  Cluster_parm_L[,grepl("b",names(Cluster_parm))] <- Cluster_parm_L[,grepl("b",names(Cluster_parm))] - drift
  Cluster_parm_U[,grepl("b",names(Cluster_parm))] <- Cluster_parm_U[,grepl("b",names(Cluster_parm))] + drift

  RR <- item_residual(theta, SA_dat, Cluster_dat, SA_parm, Cluster_parm,
                      Dv = Dv, n.nodes = n.nodes, missing_as_incorrect = missing_as_incorrect)
  LRR <- item_residual(theta, SA_dat, Cluster_dat, SA_parm_L, Cluster_parm_L,
                       Dv = Dv, n.nodes = n.nodes, missing_as_incorrect = missing_as_incorrect)
  URR <- item_residual(theta, SA_dat, Cluster_dat, SA_parm_U, Cluster_parm_U,
                       Dv = Dv, n.nodes = n.nodes, missing_as_incorrect = missing_as_incorrect)
  # Another way is to use utility() with what = "mprob" and compute assertion level residual and then sum over a cluster

  # if adaptive, need to run the function one person at a time and do further analysis
  if (nrow(combined_dat)==1) {
    return(list(RR=RR, LRR=LRR, URR=URR))
  } else {
    lower_bound_mean <- colMeans(LRR[,-1])
    lower_bound_se <- sqrt(apply(LRR[,-1], 2, var) / nrow(LRR))
    lower_bound_index <- lower_bound_mean/lower_bound_se
    lower_bound_flag <- ifelse(lower_bound_index < qnorm(Alpha), 1, 0)
    upper_bound_mean <- colMeans(URR[,-1])
    upper_bound_se <- sqrt(apply(URR[,-1], 2, var) / nrow(URR))
    upper_bound_index <- upper_bound_mean/upper_bound_se
    upper_bound_flag <- ifelse(upper_bound_index > qnorm(1 - Alpha), 1, 0)
    drift_output <- cbind(lower_bound_mean,lower_bound_se,lower_bound_index,lower_bound_flag,
                          upper_bound_mean,upper_bound_se,upper_bound_index,upper_bound_flag)
    return(list(RR=RR, LRR=LRR, URR=URR, drift_output=drift_output))
  }
}
