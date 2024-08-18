#'Marginal loglikelihood of the observed data
#' @description Compute observed data marginal loglikihood.
#' * When there are only standalone items, this computes the regular loglik of the observed data.
#' * When there are cluster items, this computes the marginal loglik where the word marginal means
#' integrating out the nuisance dimension from the conditional likelihood of the cluster items. \cr\cr
#' Use '?MIRTutils-package' for more details, such as the context of the current package and models supported.
#'
#' @param theta a scalar or a vector of examinee ability
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
#' @param Dv scaling factor for IRT model (usually 1 or 1.7)
#' @param n.nodes number of nodes used when integrating out the nuisance dimension
#' @param return_additional if TRUE, returns data loglikelihood with some additional by-product of the function in a list. See \strong{Value} section for details\cr
#' @param missing_as_incorrect by default, missings (NAs) are treated as missing; if TRUE, missings are treated as incorrect
#'
#' @return
#' If \code{return_additional} is `FALSE`, returns a dataframe with two columns: theta and marginalized data loglikelihood. \cr\cr
#' If \code{return_additional} is `TRUE`, returns the dataframe of loglikelihood with following additional tables in a list
#' * \code{probs.SA}: probability of correct response for standalone items \cr
#' * \code{probs.cluster}: (conditional) probability table of correct response for clusters at each given nodes \cr
#' * \code{parms}: parameter tables in a list
#'
#' @note If the test does not have SA items or Cluster items, use default (NULL)
#' for the corresponding data and parameter arguments \cr\cr
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
#' rst <- data_loglik(thetas[,1], SA_dat, Cluster_dat, example_SA_parm, example_Cluster_parm, n.nodes = 11, return_additional = TRUE)
#' rst$loglik
#' rst$prob.SA
#' length(rst$probs.cluster) # a list conditional probabilities. The length of the list = number of clusters
#' @export
data_loglik = function(theta, SA_dat=NULL, Cluster_dat=NULL, SA_parm=NULL, Cluster_parm=NULL, Dv=1, n.nodes = 21, return_additional=TRUE, missing_as_incorrect = F) {
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
  combined_dat = cbind(SA_dat, Cluster_dat)
  if(any(rowSums(is.na(combined_dat)) == ncol(combined_dat))) {stop("one or more examinees did not respond to any item!!!")}

  # -------------------------------------------------
  # ------------ Standalone Part --------------------
  # -------------------------------------------------
  if(is.null(SA_parm) | is.null(SA_dat)) {
    loglik_SA = 0
    probs.SA = NA
    SA_parm_3pl = SA_parm_gpc = NA
  } else {
    # +++ Add SA_parm validity checks later
    SA_dat = matrix(as.matrix(SA_dat), nrow = length(theta))
    if(missing_as_incorrect == T) SA_dat[is.na(SA_dat)] = 0  # recode if missing is treated as incorrect

    SA_parm = as.data.frame(SA_parm)
    names(SA_parm)[1] = c("a")
    names(SA_parm)[(ncol(SA_parm)-2):ncol(SA_parm)] = c("g", "ItemID", "AssertionID")
    names(SA_parm)[2:(ncol(SA_parm)-3)] = paste0("b",2:(ncol(SA_parm)-3)-1)

    # --- Separate out 3PL and GPC item par---
    dich.pos = which(apply(as.matrix(is.na(SA_parm[grepl("^b",names(SA_parm))][,-1])), 1, prod) == 1)
    poly.pos = which(apply(as.matrix(is.na(SA_parm[grepl("^b",names(SA_parm))][,-1])), 1, prod) == 0)
    SA_parm_3pl = SA_parm[dich.pos,]
    SA_parm_gpc = SA_parm[poly.pos,]
    # --- Separate out 3PL and GPC data ---
    SA_dat_3pl = SA_dat[,dich.pos]
    SA_dat_gpc = SA_dat[,poly.pos]

    if (nrow(SA_parm_3pl) != 0) {
      a = SA_parm_3pl$a
      b = SA_parm_3pl$b1
      g = rep(1,length(theta)) %o% SA_parm_3pl$g
      lin_pred_SA = Dv * t(apply(outer(theta, b, "-"), 1, function(x) a*x))
      probs.SA.3pl = g + (1 - g) * plogis(lin_pred_SA)
      colnames(probs.SA.3pl) = SA_parm_3pl$AssertionID
      lik_table_observed.3pl = matrix(as.vector(probs.SA.3pl)^as.vector(SA_dat_3pl) * as.vector(1 - probs.SA.3pl)^as.vector(1-SA_dat_3pl), nrow = length(theta))
      # lik_table_observed = matrix(dbinom(SA_dat,1,probs.SA.3pl), nrow = length(theta))
      colnames(lik_table_observed.3pl) = SA_parm_3pl$AssertionID
      loglik_SA_3pl = as.matrix(rowSums(log(lik_table_observed.3pl), na.rm = T))
    } else {
      probs.SA.3pl = NA
      loglik_SA_3pl = 0
    }

    if (nrow(SA_parm_gpc) != 0) {
      a.gpc = SA_parm_gpc[,1]
      b.gpc = SA_parm_gpc[,grepl("^b",names(SA_parm_gpc))]
      maxscr = rowSums(!is.na(b.gpc))
      b.gpc.list = split(b.gpc, seq(nrow(b.gpc)))
      b.gpc.list = lapply(b.gpc.list,function(x) x[!is.na(x)])
      # theta is a vector of examinee ability
      # a is a scalar of the a parameter for one single item
      # bs is a vector of all the step difficulty parameters of one single item
      # maxscr is the maximum possible score (where possible score starts from 0)
      # a dataframe or matrix of step difficulty, where row is item, col is steps
      gpcm = function (theta, a, bs, maxscr, Dv=1) {
        p = matrix(nrow = length(theta), ncol = length(bs)+1)
        denom = sapply(theta, function(x) 1 + sum(exp(cumsum(unlist(Dv * a * (x - bs))))))
        for (k in 1:maxscr) {
          numer = exp(Dv * a * (k * theta - sum(bs[1:floor(k)])))
          p[,k+1] = numer/denom
        }
        p[,1] = 1 - rowSums(p,na.rm = TRUE)
        p
      }
      probs.SA.gpc = mapply(gpcm, rep(list(theta),nrow(SA_parm_gpc)), a.gpc, b.gpc.list, maxscr, Dv, SIMPLIFY = FALSE)
      names(probs.SA.gpc) = SA_parm_gpc$AssertionID
      SA_dat_gpc = matrix(SA_dat_gpc, nrow = length(theta))
      lik_table_observed.gpc = sapply(1:length(probs.SA.gpc), function(x) probs.SA.gpc[[x]][cbind(1:nrow(SA_dat_gpc),SA_dat_gpc[,x]+1)])
      lik_table_observed.gpc = matrix(lik_table_observed.gpc, nrow=length(theta))
      colnames(lik_table_observed.gpc) = SA_parm_gpc$AssertionID
      loglik_SA_gpc = as.matrix(rowSums(log(lik_table_observed.gpc), na.rm = T))
    } else {
      probs.SA.gpc = NA
      loglik_SA_gpc = 0
    }
    probs.SA = list(probs.SA.3pl = probs.SA.3pl, probs.SA.gpc = probs.SA.gpc)
    loglik_SA = loglik_SA_3pl + loglik_SA_gpc
  }

  # -------------------------------------------------
  # --------------- Cluster Part --------------------
  # -------------------------------------------------
  if(is.null(Cluster_parm) | is.null(Cluster_dat)) {
    loglik_cluster = 0
    probs = NA
  } else {
    gq = statmod::gauss.quad.prob(n.nodes, dist = 'normal', sigma = 1)
    nodes = gq$nodes
    whts = gq$weights

    Cluster_dat = matrix(as.matrix(Cluster_dat), nrow = length(theta))
    if(missing_as_incorrect == T) Cluster_dat[is.na(Cluster_dat)] = 0  # recode if missing is treated as incorrect
    Cluster_parm = as.data.frame(Cluster_parm)
    colnames(Cluster_parm) = c("a","b","cluster_var","position", "ItemID", "AssertionID")
    Cluster_parm$position = match(Cluster_parm$position, sort(unique(Cluster_parm$position))) # reorder the position so it starts from 1
    colnames(Cluster_dat) = paste0(Cluster_parm$ItemID, "_", Cluster_parm$AssertionID)

    loglik_cluster = probs = list()
    for (k in 1:length(unique(Cluster_parm$position))) {
      one_cluster_parm = Cluster_parm[Cluster_parm$position == k,]
      one_cluster_dat = Cluster_dat[,grep(paste0("^",unique(one_cluster_parm$ItemID)), colnames(Cluster_dat))]
      mvars = one_cluster_parm$cluster_var
      rescaled.nodes = nodes %o% sqrt(mvars) # rescaled nodes for these assertions
      All_thetas = outer(theta, rescaled.nodes, "+")
      ma = one_cluster_parm$a # a parameters for these assertions
      mb = one_cluster_parm$b # b parameters for these assertions
      a_long = rep(ma, each = dim(All_thetas)[1] * dim(All_thetas)[2])  # each value of b parameter vector repeated [person by nodes] times
      b_long = rep(mb, each = dim(All_thetas)[1] * dim(All_thetas)[2])
      lin_pred_cluster = Dv * a_long * (All_thetas - b_long) # a*(theta + u - b) for all theta, all nodes of u, and all assertions [person by nodes by assertion]
      probs[[k]] = plogis(lin_pred_cluster)
      x = matrix(one_cluster_dat, nrow = length(theta))
      x = aperm(replicate(n.nodes, x), c(1,3,2))
      p.x = probs[[k]]^x
      q.x = (1 - probs[[k]])^(1-x)
      pqx =  p.x * q.x
      pqx_prod = apply(pqx, c(1,2), prod, na.rm=T)
      loglik_cluster[[k]] = log(pqx_prod %*% whts)
    }
    loglik_cluster = as.matrix(rowSums(matrix(simplify2array(loglik_cluster), nrow = length(theta))))
  }

  # -------------------------------------------------
  # -------------------- Output ---------------------
  # -------------------------------------------------
  loglik = loglik_SA + loglik_cluster
  colnames(loglik) = "loglik"
  if (return_additional == F) {
    output = data.frame(theta=theta, loglik=loglik)
  } else {
    output = list(loglik = data.frame(theta=theta, loglik=loglik),
                  probs.SA = probs.SA, probs.cluster = probs,
                  parms = list(SA_parm_3pl=SA_parm_3pl,SA_parm_gpc=SA_parm_gpc,Cluster_parm=Cluster_parm))
  }
  return(output)
}
