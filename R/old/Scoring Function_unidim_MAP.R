####################  Person Scoring with MAP method ##########################################
## INPUT: run this function for one person at a time
## 1. SA_dat: a vector, or a row matrix/dataframe of responses to standalone items
## 2. Cluster_dat: NULL; a placeholder 
## 3. SA_parm: a matrix or dataframe of a,b,g parameters (column must be in this order), ItemID, and Assertion_ID for SA items
## 4. Cluster_parm: NULL; a placeholder 
## 5. prior_mu: mean of the prior distribution
## 6. prior_var: variance of the prior distribution
## 7. Dv: scaling factor for IRT model [1 or 1.7]
## 8. n.nodes: use default; a placehoder
## 9. SE: if TRUE, return standard error

# *** CATUION:
# *** This function currently only does unidimensional MAP scoring
##########################################################################################################
scoring.MAP = function(SA_dat=NULL, Cluster_dat=NULL, SA_parm=NULL, Cluster_parm=NULL, prior_mu=0, prior_var=1, Dv=1, n.nodes = 21, censor=c(-4,4), SE=F) {
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
    ll = obs.data.loglik(theta, SA_dat, Cluster_dat=NULL, SA_parm, Cluster_parm=NULL, Dv, n.nodes, return_additional=F)
    ll$prior.ll = log(dnorm(theta, prior_mu, sqrt(prior_var)))
    ll$post.ll = ll$loglik + ll$prior.ll
    return(-sum(ll$post.ll))
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