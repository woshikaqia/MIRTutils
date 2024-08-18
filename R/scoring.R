#' IRT person scoring
#' @description
#' Compute IRT latent scores (theta estimates). Scoring method currently supported are: \cr
#' * MMLE: marginal maximum likelihood estimation. When no Rasch testlet model item is involved,
#'  it reduces to regular MLE.
#' * More methods to be added in the near future...
#'
#'
#' Use '?MIRTutils-package' for more details, such as the context of the current package and models supported.
#'
#' @param SA_dat For one examinee, a vector or row matrix/dataframe of response to standalone items.
#' Responses must follow the row order in \code{SA_parm}. Use NA for missing responses.
#' @param Cluster_dat For one examinee, a vector or row matrix/dataframe of response to cluster items.
#' Responses must follow the row order in \code{Cluster_parm}. Use NA for missing responses.
#' @param SA_parm a matrix or dataframe of item parameters for standalone items, where columns are
#'  a (slope), b1, b2, ..., b_k (difficulty or step difficulty), g (guessing), ItemID, and AssertionID.
#'  Columns must follow the above order.
#'  See \code{example_SA_parm} for an example. Use \code{?example_SA_parm} for detailed column descriptions
#' @param Cluster_parm a matrix or dataframe of item parameters for cluster items, where columns are
#'  a (slope), b (difficulty), cluster variance, cluster position, ItemID, and AssertionID.
#'  Columns must follow the above order.
#'  See \code{example_Cluster_parm} for an example. Use \code{?example_Cluster_parm} for detailed column descriptions
#' @param Dv scaling factor for IRT model (usually 1 or 1.7)
#' @param n.nodes number of nodes used when integrating out the nuisance dimension
#' @param censor when there's perfect or all wrong score, this is a 2-element vector of lower and upper limit of
#'      estimated theta for such score patterns.\cr
#'      NOTE: This argument only censor the estimated theta for the perfect or all wrong scores. It dose not censor
#'      possible extreme theta values estimated from other score patterns. Perform censoring after running this code
#'      if a global censorship is desired.
#' @param correction_val a value to add or subtract when there's perfect or all wrong score to avoid extremely large theta estimate.
#' @param SE if TRUE, returns standard error
#'
#' @return a list of scoring results, where the first element is the estimated theta.
#'
#' @note
#' This function should be run for one examinee at a time. Use parallel processing for higher processing speed.\cr\cr
#' If the test does not have SA items or Cluster items, use default (NULL) for the corresponding data and parameter arguments \cr\cr
#'
#' @author Zhongtian Lin lzt713@gmail.com
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
#' # Scoring for the first examinee
#' rst <- scoring(SA_dat[1,], Cluster_dat[1,], example_SA_parm, example_Cluster_parm, n.nodes = 11, SE=TRUE)
#' rst$par  # estimated theta
#' rst$SE  # estimated standard error
#' @export

##########################################################################################################
scoring = function(SA_dat=NULL, Cluster_dat=NULL, SA_parm=NULL, Cluster_parm=NULL, Dv=1, n.nodes = 21, censor=c(-6,6), correction_val = 0.5, SE=FALSE) {
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")}
  if(is.null(SA_dat) & is.null(Cluster_dat)) {stop("No data found!!!")}

  if(!is.null(SA_dat)) {
    if(is.vector(SA_dat)) SA_dat = matrix(SA_dat, nrow = 1) else SA_dat = as.matrix(SA_dat)
  }
  if(!is.null(Cluster_dat)) {
    if(is.vector(Cluster_dat)) Cluster_dat = matrix(Cluster_dat, nrow = 1) else Cluster_dat = as.matrix(Cluster_dat)
  }
  combined_dat = cbind(SA_dat, Cluster_dat)

  if(!is.null(SA_dat) && !is.null(SA_parm)) {
    temp_SA_parm = SA_parm
    temp_SA_parm$cluster_var = 0 # this is for all wrong/perfect score correction with only SA items (see highest_cluster_var_pos)
  } else {
    temp_SA_parm = NULL
  }

  if(correction_val == 0) {warning("No correction value used for perfect or all wrong scores. Theta estimates could be very extreme at these scores")}

  # -----starting values and extreme handling------
  # maxmium possible score  = number of item + extra score points from polytomous item
  if (!is.null(SA_parm)) {
    possible = ncol(combined_dat) + sum(!is.na(SA_parm[,c(-1,-2,-(ncol(SA_parm)-2):-ncol(SA_parm))]))
  } else {
    possible = ncol(combined_dat)
  }
  rawscore = rowSums(combined_dat, na.rm = TRUE)
  allperfect = which(rawscore==possible)
  allwrong = which(rawscore==0)
  # if all wrong or perfect, add or subtract 0.5 from the cluster with largest variance
  highest_cluster_var_pos = which.max(c(temp_SA_parm$cluster_var, Cluster_parm$cluster_var))
  if (length(allperfect)>0) {
    combined_dat[allperfect,highest_cluster_var_pos] = combined_dat[allperfect,highest_cluster_var_pos] - correction_val
  }
  if (length(allwrong)>0) {
    combined_dat[allwrong,highest_cluster_var_pos] = combined_dat[allwrong,highest_cluster_var_pos] + correction_val
  }
  # starting value
  start_val = ((log((rawscore + .00001)/(possible-rawscore + .00001))))

  #objective function
  obj_fn = function(theta) {
    ll = data_loglik(theta, SA_dat, Cluster_dat, SA_parm, Cluster_parm, Dv, n.nodes, return_additional=FALSE)
    return(-sum(ll$loglik))
  }

  if(prod(rowSums(combined_dat, na.rm = TRUE) == 0)) {
    out = nlminb(censor[1], obj_fn, control = list(maxit = 0))
  } else if (prod(rowSums(combined_dat, na.rm = TRUE) == possible)) {
    out = nlminb(censor[2], obj_fn, control = list(maxit = 0))
  } else {
    out = nlminb(start_val, obj_fn, control = list(rel.tol = 1e-10))
  }

  ### Return standard errors
  if(out$par < censor[1]){
    se = 1/sqrt(optimHess(censor[1], obj_fn))
  } else if (out$par > censor[2]){
    se = 1/sqrt(optimHess(censor[2], obj_fn))
  } else {
    se = 1/sqrt(optimHess(out$par, obj_fn))
  }

  if(SE==TRUE) out$SE = se
  out
}
