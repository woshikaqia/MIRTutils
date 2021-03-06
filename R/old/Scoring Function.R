####################  Person Scoring ##########################################
## INPUT: run this function for one person at a time
## 1. SA_dat: a vector, or a row matrix/dataframe of responses to standalone items
## 2. Cluster_dat:
# In general: a matrix or dataframe of responses to cluster items. One assertion per column. Column order must match row order in Cluster_parm
# For one student: row matrix/dataframe preferred but it can simply be a vector of responses
## 3. SA_parm: a matrix or dataframe of a,b,g parameters (column must be in this order), ItemID, and Assertion_ID for SA items
## 4. Cluster_parm: a matrix or dataframe of a,b and variance parameters for each assertion, a column of cluster position, a column of cluster ItemID, and a column of Assertion_ID for Cluster items
## 5. Dv: scaling factor for IRT model [1 or 1.7]
## 6. n.nodes: number of nodes used when integrating out the specific dimension
## 7. SE: if TRUE, return standard error
# *** Its okay to treat SA item as clusters. To do so, simply store them in the "Cluster_parm" argument with 0 variances, and store all student responses in "Cluster_dat"

# *** DONOTUSE:
  # SA_dat and Cluster_dat can be multi-row matrix/dataframe where it is one person per row. In that case, it simultanousely solves
  # all thetas but it assumes thetas are not independent. So it shouldn't be used for scoring purpose.
##########################################################################################################
scoring = function(SA_dat=NULL, Cluster_dat=NULL, SA_parm=NULL, Cluster_parm=NULL, Dv=1, n.nodes = 21, censor=c(-4,4), SE=F) {
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")}
  if(is.null(SA_dat) & is.null(Cluster_dat)) {stop("No data found!!!")}

  if(!is.null(SA_dat)) {
    if(is.vector(SA_dat)) SA_dat = matrix(SA_dat, nrow = 1) else SA_dat = as.matrix(SA_dat)
  }
  if(!is.null(Cluster_dat)) {
    if(is.vector(Cluster_dat)) Cluster_dat = matrix(Cluster_dat, nrow = 1) else Cluster_dat = as.matrix(Cluster_dat)
  }
  combined_dat = cbind(SA_dat, Cluster_dat)

  # starting values
  pvalue = rowMeans(combined_dat, na.rm = T)
  start_val = log(pvalue/(1-pvalue))

  #objective function
  obj_fn = function(theta) {
    ll = MIRTutils::obs.data.loglik(theta, SA_dat, Cluster_dat, SA_parm, Cluster_parm, Dv, n.nodes, return_additional=F)
    return(-sum(ll$loglik))
  }

  if(prod(rowSums(combined_dat, na.rm = T) == 0)) {
    out = nlminb(censor[1], obj_fn)
  } else if (prod(rowSums(combined_dat, na.rm = T) == ncol(combined_dat))) {
    out = nlminb(censor[2], obj_fn)
  } else {
    out = nlminb(start_val, obj_fn)
  }

  ### Return standard errors
  if(out$par < censor[1]){
    se = 1/sqrt(optimHess(censor[1], obj_fn))
  } else if (out$par > censor[2]){
    se = 1/sqrt(optimHess(censor[2], obj_fn))
  } else {
    se = 1/sqrt(optimHess(out$par, obj_fn))
  }

  if(SE==T) out$SE = se
  out
}
