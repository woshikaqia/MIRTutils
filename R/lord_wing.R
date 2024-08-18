#' Lord-Wingersky algorithm for computing the marginal probability of the raw scores
#' @description
#' For N persons (theta values), this function makes use of the recursive algorithm
#' by Lord & Wingersky (1984) and its extensions (Lin et.al, 2024) to compute
#' the marginal probability of the raw scores for a cluster
#' item (with J dichotomously scored assertions) modeled by the Rasch testlet model.
#' The word "marginal" here means integrating out the nuisance dimension from the
#' conditional likelihood of the cluster items.
#' @param cluster_var a vector of  of length J with repeated values of the cluster
#' variance (e.g., rep(0.8, 9) for a 9-assertion cluster item).
#' Alternatively, a scalar value of the cluster variance for the item.
#' @param a a vector of length J for the a (slope) parameters
#' @param b a vector of of length J for the b (difficulty) parameters
#' @param theta a vector of length N for the thetas
#' @param n.nodes number of nodes used when integrating out the specific dimension
#' @param return_additional if TRUE, return a list containing the marginal probability as well
#' as some additional by-product of the function such as the conditional probability tables. See \strong{Value} section for details.
#' @param Dv scaling factor for IRT model (usually 1 or 1.7)
#'
#' @return
#' When \code{return_additional = FALSE}, returns the marginal probability of raw scores, which is
#' a J+1 by N matrix, where J+1 is the number of possible raw scores \cr\cr
#' When \code{return_additional = TRUE}, returns a list containing
#' * \code{prk}: the marginal probability
#' * \code{prob}: N by \code{n.nodes} by J array containing the conditional probability
#' of correct response for each theta at each node of the nuisance dimension for each assertion
#' * \code{qrob}: \code{1 - prob}
#' * \code{nodes} and \code{whts}: nodes and weights used in the calculation.
#'
#' @references
#' Lord, F. M., & Wingersky, M. S. (1984). Comparison of IRT true-score and equipercentile observed-score "equatings".
#' \emph{Applied Psychological Measurement}, 8(4), 453-461. \cr\cr
#' Lin, Z., Jiang, T., Rijmen, F. et al. (2024). Asymptotically Correct Person Fit z-Statistics For the Rasch Testlet Model.
#' \emph{Psychometrika}, https://doi.org/10.1007/s11336-024-09997-y.
#' @author Zhongtian Lin <lzt713@gmail.com>
#' @examples
#' data(example_Cluster_parm)
#' # Compute on the first cluster, for 5 students
#' one_cluster_parm <- example_Cluster_parm[example_Cluster_parm$position == 1,]
#' rst <- lord_wing(one_cluster_parm$cluster_var , one_cluster_parm$a, one_cluster_parm$b,
#' theta <- seq(-2,2,1), return_additional = TRUE)
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

  Ps_d = Ps*Qs  # first derivatives of Ps with respect to theta (under Rasch model)
  Qs_d = -Ps_d  # first derivatives of Qs with respect to theta (under Rasch model)

  # Three lines below prepare a 5-dimensional array for LW calculation [assertion+2 by assertion+2 by assertion by nodes by person]
  PPs = rbind(cbind(rep(0,n.ass+1),diag(n.ass+1)),rep(0,n.ass+2)) %o% Ps
  QQs = rbind(rep(0,n.ass+2),cbind(rep(0,n.ass+1),diag(n.ass+1))) %o% Qs
  PQs = PPs + QQs
  # Same as the three lines above but stores the fist derivatives
  PPs_d = rbind(cbind(rep(0,n.ass+1),diag(n.ass+1)),rep(0,n.ass+2)) %o% Ps_d
  QQs_d = rbind(rep(0,n.ass+2),cbind(rep(0,n.ass+1),diag(n.ass+1))) %o% Qs_d
  PQs_d = PPs_d + QQs_d

  # prk: conditional (on nodes) results of the following quantities for a cluster k
  #   W0: a matrix containing the probability of raw scores after seeing assertion j, where j = 1, 2, 3,..., n_k
  #   W0_d: a matrix containing the derivatives of the probability of raw scores after seeing assertion j, where j = 1, 2, 3,..., n_k
  #   W1: a matrix containing a special set of quantities after seeing assertion j, where j = 1, 2, 3,..., n_k.
  #   W2: a matrix containing another special set of quantities similar to W1
  # W1 and W2 are necessary components for computing the derivative of the expectation of the loglikelihood for a cluster
  # W0_d is necessary components for computing the variance of loglikelihood for a cluster
  prk = lapply(1:dim(PQs)[5], function(z) { # sapply over the fifth dimension (the person dimension)
    prk.cond.u = matrix(c(c(0,1),rep(0,n.ass)), dim(PQs)[1], dim(PQs)[4])
    prk.cond.u_previous = matrix(c(c(0,1),rep(0,n.ass)), dim(PQs)[1], dim(PQs)[4])
    prk.cond.u_W1 = matrix(0, dim(PQs)[1], dim(PQs)[4])
    prk.cond.u_W1_previous = matrix(0, dim(PQs)[1], dim(PQs)[4])
    prk.cond.u_W2 = matrix(0, dim(PQs)[1], dim(PQs)[4])
    prk.cond.u_d = matrix(0, dim(PQs)[1], dim(PQs)[4])
    for (j in 1:n.ass) {

      # These two components are necessary components for W1 and W2. They are previous step's
      # W0 and W1 results with the last row removed and with an all 0 row added as the first row
      prk.cond.u_previous = rbind(0, prk.cond.u_previous[-nrow(prk.cond.u_previous),])
      prk.cond.u_W1_previous = rbind(0, prk.cond.u_W1_previous[-nrow(prk.cond.u_W1_previous),])

      # W0 & w0_d
      prk.cond.u = matrix(c(c(rbind(prk.cond.u[,1:(dim(PQs)[4]-1)],matrix(0, dim(PQs)[1]*dim(PQs)[4], dim(PQs)[4]-1))), prk.cond.u[,dim(PQs)[4]]), dim(PQs)[1]*dim(PQs)[4], dim(PQs)[4])
      PQs_u = PQs[,,j,1:dim(PQs)[4],z]
      PQs_u = matrix(aperm(PQs_u, c(1,3,2)), dim(PQs)[1]*dim(PQs)[4], dim(PQs)[1])

      prk.cond.u_d = matrix(c(c(rbind(prk.cond.u_d[,1:(dim(PQs_d)[4]-1)],matrix(0, dim(PQs_d)[1]*dim(PQs_d)[4], dim(PQs_d)[4]-1))), prk.cond.u_d[,dim(PQs_d)[4]]), dim(PQs_d)[1]*dim(PQs_d)[4], dim(PQs_d)[4])
      PQs_u_d = PQs_d[,,j,1:dim(PQs_d)[4],z]
      PQs_u_d = matrix(aperm(PQs_u_d, c(1,3,2)), dim(PQs_d)[1]*dim(PQs_d)[4], dim(PQs_d)[1])

      prk.cond.u_d = t(PQs_u_d) %*% prk.cond.u + t(PQs_u) %*% prk.cond.u_d
      prk.cond.u = t(PQs_u) %*% prk.cond.u

      # W1
      prk.cond.u_W1 = matrix(c(c(rbind(prk.cond.u_W1[,1:(dim(PQs)[4]-1)],matrix(0, dim(PQs)[1]*dim(PQs)[4], dim(PQs)[4]-1))), prk.cond.u_W1[,dim(PQs)[4]]), dim(PQs)[1]*dim(PQs)[4], dim(PQs)[4])
      # t(PQs_u) %*% prk.cond.u_W1 ---> This is p_j* W1(j-1,r-1) + q*W1(j-1,r). Needed component 1 for W1
      # Ps[j,,z] * (theta[z]- b[j])^2 ---> This p_j*(theta-b_j). Needed component 2 for W1
      prk.cond.u_W1 = t(PQs_u) %*% prk.cond.u_W1 + sweep(prk.cond.u_previous, 2, (Ps[j,,z] * (theta[z] - b[j])), "*")

      # W2
      prk.cond.u_W2 = matrix(c(c(rbind(prk.cond.u_W2[,1:(dim(PQs)[4]-1)],matrix(0, dim(PQs)[1]*dim(PQs)[4], dim(PQs)[4]-1))), prk.cond.u_W2[,dim(PQs)[4]]), dim(PQs)[1]*dim(PQs)[4], dim(PQs)[4])
      # t(PQs_u) %*% prk.cond.u_W2 ---> This is p_j* W2(j-1,r-1) + q*W2(j-1,r). Needed component 1 for W1
      # Ps[j,,z] * (theta[z]- b[j])^2 ---> This p_j*(theta-b_j)^2 Needed component 2 for W2
      # 2 * Ps[j,,z] * (theta[z] - b[j]) ---> This 2*p_j*(theta-b_j) Needed component 3 for W2
      prk.cond.u_W2 = t(PQs_u) %*% prk.cond.u_W2 +
        sweep(prk.cond.u_previous, 2, (Ps[j,,z] * (theta[z] - b[j])^2), "*") +
        sweep(prk.cond.u_W1_previous, 2, (2*(Ps[j,,z] * (theta[z] - b[j]))), "*")

      # Save for the next iteration
      prk.cond.u_previous = prk.cond.u
      prk.cond.u_W1_previous = prk.cond.u_W1
    }
    list(W0=prk.cond.u, W0_d=prk.cond.u_d, W1=prk.cond.u_W1, W2=prk.cond.u_W2)
  })
  # return additional by-product or not?
  if(return_additional == TRUE) {
    return(list(prk=prk, probs=probs, qrobs=qrobs, nodes=nodes, whts=whts))
  } else {
    return(prk)
  }
}
