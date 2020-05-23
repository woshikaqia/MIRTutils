#' Marginal probability of the observed data
#'
#' @description Computes the marginal probability of the observed data
#'
#' @param theta a scalar or a vector of student ability
#' @param SA_dat For one student: a vector of responses to standalone items \cr\cr
#' For two or more students: a matrix or dataframe of responses to standalone items. One assertion per column. Column order must match row order in SA_parm \cr\cr
#' Use NA for missing reponses
#' @param Cluster_dat For one student: a vector of responses to cluster items \cr\cr
#' For two or more students: a matrix or dataframe of responses to cluster items. One assertion per column. Column order must match row order in Cluster_parm \cr\cr
#' Use NA for missing reponses
#' @param SA_parm a matrix or dataframe of a,b,g parameters (column must be in this order), ItemID, and Assertion_ID for SA items
#' @param Cluster_parm a matrix or dataframe with columns of a, b and variance parameters for each assertion, as well as a column of cluster position, a column of cluster ItemID, and a column of Assertion_ID for Cluster items
#' @param Dv scaling factor for IRT model [use 1 or 1.7]
#' @param n.nodes number of nodes used when integrating out the specific dimension
#' @param return_additional if TRUE, return a list of the loglik plus some additional by-product of the function such as: \cr\cr
#' probs.SA: \tabular{ll}{\tab probablity table of correct response for standalones}
#' probs.cluster \tabular{ll}{\tab (conditional) probablity table of correct response for clusters at each given node}
#' @param missing_as_incorrect by default, missings (NAs) are treated as missing; if TRUE, missings are treated as incorrect
#' @return If \emph{return_additional} is FALSE, returns a dataframe with two columns: theta and marginalized data loglikelihood; \cr
#' If \emph{return_additional} is TRUE, returns the dataframe of loglikelihood plus additional tables \emph{probs.SA} and \emph{probs.cluster} in a list
#' @export
#'
#' @note If the test does not have SA items or Cluster items, use default (NULL) for the corresponding data and parameter arguments \cr\cr
#' SA items can be treated as clusters. To do so, store SA item parameters in the "Cluster_parm" argument with 0 variances, and store all student responses in "Cluster_dat"
#' @examples
#'
#'
obs.data.loglik = function(theta, SA_dat=NULL, Cluster_dat=NULL, SA_parm=NULL, Cluster_parm=NULL, Dv=1, n.nodes = 21, return_additional=F, missing_as_incorrect = F) {
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")}
  if(is.null(SA_dat) & is.null(Cluster_dat)) {stop("No data found!!!")}
  if((is.null(SA_parm) | is.null(SA_dat)) && (is.null(Cluster_parm) | is.null(Cluster_dat)))  {stop("Data or Parameter file missing!!!")}

  if(!is.null(SA_dat)) {
    if(is.vector(SA_dat)) SA_dat = matrix(SA_dat, nrow = 1) else SA_dat = as.matrix(SA_dat)
  }
  if(!is.null(Cluster_dat)) {
    if(is.vector(Cluster_dat)) Cluster_dat = matrix(Cluster_dat, nrow = 1) else Cluster_dat = as.matrix(Cluster_dat)
  }
  combined_dat = cbind(SA_dat, Cluster_dat)
  if(any(rowSums(is.na(combined_dat)) == ncol(combined_dat))) {stop("one or more students did not respond to any item!!!")}

  if(is.null(SA_parm) | is.null(SA_dat)) {
    loglik_SA = 0
    probs.SA = NA
  } else {
    SA_dat = matrix(as.matrix(SA_dat), nrow = length(theta))
    if(missing_as_incorrect == T) SA_dat[is.na(SA_dat)] = 0  # recode if missing is treated as incorrect
    SA_parm = as.data.frame(SA_parm)
    names(SA_parm) = c("a","b","g", "ItemID", "Assertion_ID")
    a = SA_parm[,1]
    b = SA_parm[,2]
    g = SA_parm[,3]
    a.parm = rep(1,length(theta)) %o% a
    b.parm = rep(1,length(theta)) %o% b
    g.parm = rep(1,length(theta)) %o% g
    theta.parm = theta %o% rep(1,length(b))
    lin_pred_SA = Dv * (theta.parm * a.parm - b.parm)
    probs.SA = g.parm + (1 - g.parm) * plogis(lin_pred_SA)
    lik_table_observed = as.data.frame(matrix(dbinom(SA_dat,1,probs.SA), nrow = length(theta)))
    colnames(lik_table_observed) = SA_parm$Assertion_ID
    loglik_SA = as.matrix(rowSums(log(lik_table_observed), na.rm = T))
  }

  if(is.null(Cluster_parm) | is.null(Cluster_dat)) {
    loglik_cluster = 0
    probs = NA
  } else {
    gq = statmod::gauss.quad.prob(n.nodes, dist = 'normal', sigma = 1)
    nodes = gq$nodes
    whts = gq$weights

    Cluster_dat = matrix(as.matrix(Cluster_dat), nrow = length(theta))
    if(missing_as_incorrect == T) Cluster_dat[is.na(Cluster_dat)] = 0  # recode if missing is treated as incorrect
    colnames(Cluster_parm) = c("a","b","cluster_var","position", "ItemID", "Assertion_ID")
    Cluster_parm$position = match(Cluster_parm$position, sort(unique(Cluster_parm$position))) # reorder the position so it starts from 1
    colnames(Cluster_dat) = paste0(Cluster_parm$ItemID, "_", Cluster_parm$Assertion_ID)

    loglik_cluster = probs = list()
    for (k in 1:length(unique(Cluster_parm$position))) {
      one_cluster_parm = filter(Cluster_parm, position == k)
      one_cluster_dat = Cluster_dat[,grep(paste0("^",unique(one_cluster_parm$ItemID)), colnames(Cluster_dat))]
      mvars = one_cluster_parm$cluster_var
      rescaled.nodes = nodes %o% sqrt(mvars) # rescaled nodes for these assertions
      All_thetas = outer(theta, rescaled.nodes, "+")
      ma = one_cluster_parm$a # a parameters for these assertions
      mb = one_cluster_parm$b # b parameters for these assertions
      mtheta = theta  # students' theta in a vector
      a_long = rep(ma, each = dim(All_thetas)[1] * dim(All_thetas)[2])  # each value of 5-element b parameter vector repeated [person by nodes] times
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
                  probs.SA = probs.SA, probs.cluster = probs)
  }
  return(output)
}


# now with git
