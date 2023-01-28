#' Lord-Wingersky algorithm for computing the marginal probability of the raw scores
#' @description
#' For N persons (theta values), this function makes use of the recursive algorithm by Lord & Wingersky (1984) to compute
#' the marginal probability of the raw scores for
#' a cluster item (with J dichotomously scored assertions) modeled by the Rasch testlet model. The word "marginal" means
#' integrating out the nuisance dimension from the conditional likelihood of the cluster items.
#' @param cluster_var a vector of  of length J with repeated values of the cluster
#' variance (e.g., rep(0.8, 9) for a 9-assertion cluster item).
#' Alternativaly, a scalar value of the cluster variance for the item.
#' @param a a vector of length J for the a (slope) parameters
#' @param b a vector of of length J for the b (difficulty) parameters
#' @param theta a vector of length N for the thetas
#' @param n.nodes number of nodes used when integrating out the specific dimension
#' @param return_additional if TRUE, return a list containing the marginal probability as well
#' as some additional by-product of the function such as the conditional probability tables. See \strong{Value} section for details.
#' @param Dv scaling factor for IRT model (usually 1 or 1.7)
#'
#' @return When \code{return_additional = FALSE}, returns the marginal probability of raw scores, which is
#' a J+1 by N matrix, where J+1 is the number of possible raw scores \cr
#' When \code{return_additional = TRUE}, returns a list containing the marginal probability (\code{prk.marginal})
#' as well as \code{prob} and \code{qrob}, where \code{qrob} is a N by \code{n.nodes} by J array containing the conditional probability
#' of correct response for each theta at each node of the nuisance dimension for each assertion, and \code{qrob = 1 - prob}
#' @references Lord, F. M., & Wingersky, M. S. (1984). Comparison of IRT true-score and equipercentile observed-score "equatings".
#' \emph{Applied Psychological Measurement}, 8(4), 453-461.
#' @author Zhongtian Lin lzt713@gmail.com
#' @examples
#' data(example_Cluster_parm)
#' # Compute on the first cluster, for 5 students
#' one_cluster_parm <- example_Cluster_parm[example_Cluster_parm$position == 1,]
#' rst <- lord_wing(one_cluster_parm$cluster_var , one_cluster_parm$a, one_cluster_parm$b,
#'  theta <- seq(-2,2,1), return_additional = TRUE)
#' @export

lord_wing = function(cluster_var, a, b, theta, n.nodes=21, return_additional=FALSE, Dv=1) {
  gq = statmod::gauss.quad.prob(n.nodes, dist = 'normal', sigma = 1)
  nodes = gq$nodes
  whts = gq$weights
  if(length(cluster_var)==1) cluster_var = rep(cluster_var, length(b))
  rescaled.nodes = nodes %o% sqrt(cluster_var)
  All_thetas = outer(theta, rescaled.nodes, "+") #  theta+u; all theta with all nodes of u [person by nodes]
  a_long = rep(a, each = dim(All_thetas)[1] * dim(All_thetas)[2])  # each value of b parameter vector repeated [person by nodes] times
  b_long = rep(b, each = dim(All_thetas)[1] * dim(All_thetas)[2])  # each value of b parameter vector repeated [person by nodes] times
  lin_pred = Dv * a_long * (All_thetas - b_long) # a*(theta + u - b) for all theta, all nodes of u, and all assertions [person by nodes by assertion]
  probs = plogis(lin_pred)
  qrobs = 1 - probs
  qrobs[qrobs==0] = .Machine$double.xmin
  n.ass = length(b)

  # run for all persons
  Ps = aperm(probs, c(3,2,1)) # permute the conditional probabilities array so it's [assertion by nodes by person]
  Qs = aperm(qrobs, c(3,2,1)) # 1-p

  # Three lines below prepare a 5-dimensional array for LW calculation [assertion+2 by assertion+2 by assertion by nodes by person]
  PPs = rbind(cbind(rep(0,n.ass+1),diag(n.ass+1)),rep(0,n.ass+2)) %o% Ps
  QQs = rbind(rep(0,n.ass+2),cbind(rep(0,n.ass+1),diag(n.ass+1))) %o% Qs
  PQs = PPs + QQs

  # prk.marginal: marginal probability of raw scores for the item "k" (notation k stems from upper level function notation)
  prk.marginal = sapply(1:dim(PQs)[5], function(z) { # sapply over the fifth dimension (the person dimension)
    prk.cond.u = matrix(c(c(0,1),rep(0,n.ass)), dim(PQs)[1], dim(PQs)[4])
    for (i in 1:n.ass) {
      prk.cond.u = matrix(c(c(rbind(prk.cond.u[,1:(dim(PQs)[4]-1)],matrix(0, dim(PQs)[1]*dim(PQs)[4], dim(PQs)[4]-1))), prk.cond.u[,dim(PQs)[4]]), dim(PQs)[1]*dim(PQs)[4], dim(PQs)[4])
      PQs_u = PQs[,,i,1:dim(PQs)[4],z]
      PQs_u = matrix(aperm(PQs_u, c(1,3,2)), dim(PQs)[1]*dim(PQs)[4], dim(PQs)[1])
      prk.cond.u = t(PQs_u) %*% prk.cond.u
    }
    as.vector(prk.cond.u[-1,] %*% whts)
  })

  # return additional by-product or not?
  if(return_additional == TRUE) {
    return(list(prk.marginal=prk.marginal, probs=probs, qrobs=qrobs))
  } else {
    return(prk.marginal)
  }
}
