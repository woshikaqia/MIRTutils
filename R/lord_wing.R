####################  Lord and Wingersky algorithm with matrix operations ################################
## Description:
##  For N persons (theta values), this function computes the marginal probability of the raw scores for
##  a cluster item (with J dichotomously scored assertions) modeled by the testlet model, where "marginal" means
##  the nuisance dimension being integrated out

## INPUT:
## 1. cluster_var:  a vector of  of legnth J with repeated values of the cluster variance(e.g., rep(0.8, 9) for a 9-assertion item)
##                  Alternativaly, a scalar value of the cluster variance for the item.
## 2. a: a vector of length J for the a parameters
## 3. b: a vector of of length J for the b parameters
## 4. theta: a vector of thetas of length N
## 5. n.node: number of nodes used when integrating out the specific dimension
## 6. return_additional: if TRUE, return a list containing the mraginal probablity as well as some additional by-product of the function
                        # such as the conditional probabilty tables (currently support: prob, qrob)

## Return values:
##  When return_additional=FALSE: returns the marginal probabilty of raws scores, which is
##                                a J+1 by N matrix, where J+1 is the number of possible raw scores
##  When return_additional=TRUE: a list containing the mraginal probablity (prk.marginal) as well as prob and qrob, where
##                               prob is a N by n.nodes by J array containing the conditional probablity of correct
##                               response for each theta at each node of the nuisance dimension, for each assertion, and
##                               qrob = 1 - prob

## Author: Zhongtian Lin
## Last Updated: 09/21/2020
#'@export

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
  Ps = aperm(probs, c(3,2,1)) # permutate the conditional probabilities array so it's [assertion by nodes by person]
  Qs = aperm(qrobs, c(3,2,1)) # 1-p

  # Three lines below prepare a 5-dimenisonal array for LW calculation [assertion+2 by assertion+2 by assertion by nodes by person]
  PPs = rbind(cbind(rep(0,n.ass+1),diag(n.ass+1)),rep(0,n.ass+2)) %o% Ps
  QQs = rbind(rep(0,n.ass+2),cbind(rep(0,n.ass+1),diag(n.ass+1))) %o% Qs
  PQs = PPs + QQs

  # prk.marginal: marginal probablity of raw scores for the item "k" (notation k stems from upper level function notation)
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

  # retrun additional by-product or not?
  if(return_additional == TRUE) {
    return(list(prk.marginal=prk.marginal, probs=probs, qrobs=qrobs))
  } else {
    return(prk.marginal)
  }
}
