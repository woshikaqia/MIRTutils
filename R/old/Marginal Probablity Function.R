####################  Marginal probability of correct response ##########################################
## INPUT: 
## 1. theta: a vector of thetas or a scalar value of theta
## 2. SA_parm: a matrix or dataframe of a,b,g parameters (column must be in this order), ItemID, and Assertion_ID for SA items
## 3. Cluster_parm: a matrix or dataframe of a,b and variance parameters for each assertion, a column of cluster position, a column of cluster ItemID, and a column of Assertion_ID for Cluster items
## 4. Dv: scaling factor for IRT model [1 or 1.7]
## 5. n.nodes: number of nodes used when integrating out the specific dimension 
# *** Its okay to treat SA item as clusters. To do so, simply store them in the "Cluster_parm" argument with 0 variances 
##########################################################################################################
marginal.prob = function(theta, SA_parm=NULL, Cluster_parm=NULL, Dv=1, n.nodes = 21) {
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")} 
  
  if(is.null(SA_parm)) {
    probs.SA = NULL
  } else {
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
    colnames(probs.SA) = SA_parm$Assertion_ID
  }
  
  if(is.null(Cluster_parm)) {
    marginal_prob = NULL
  } else {
    gq = statmod::gauss.quad.prob(n.nodes, dist = 'normal', sigma = 1)
    nodes = gq$nodes
    whts = gq$weights
    
    names(Cluster_parm) = c("a","b","cluster_var","position", "ItemID", "Assertion_ID")
    cluster_var = Cluster_parm$cluster_var
    rescaled.nodes = nodes %o% sqrt(cluster_var) # rescaled nodes for these assertions
    All_thetas = outer(theta, rescaled.nodes, "+")
    ma = Cluster_parm$a # a parameters for these assertions
    mb = Cluster_parm$b # b parameters for these assertions
    a_long = rep(ma, each = dim(All_thetas)[1] * dim(All_thetas)[2])  # each value of 5-element b parameter vector repeated [person by nodes] times
    b_long = rep(mb, each = dim(All_thetas)[1] * dim(All_thetas)[2])
    lin_pred_cluster = Dv * a_long * (All_thetas - b_long) # a*(theta + u - b) for all theta, all nodes of u, and all assertions [person by nodes by assertion]
    probs = plogis(lin_pred_cluster)
    marginal_prob = matrix(apply(probs,3, function(x) x %*% whts), nrow = length(theta))
    colnames(marginal_prob) = Cluster_parm$Assertion_ID
  }
  
  # -------------------------------------------------
  # -------------------- Output ---------------------
  # -------------------------------------------------
  output = data.frame(theta=theta, probs.SA, marginal_prob)
  return(output)
}
