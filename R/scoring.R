####################  Person Scoring ##########################################
## Description: 
## Compute IRT latent score. Model currently supported are 1-3 PL, GPCM, Rasch Testlet, 
## and a mix of these models. see dat_loglik() for model details. 
## Scoring method currently supported are:
##   MMLE: marginal maximum likelihood estimation. When no Rasch testlet model item is invovled, 
##         it reduce to regular MLE
##
## Dependency: data_loglik()
##
## INPUT: run this function for one person at a time
## 1. SA_dat: a vector, or a row matrix/dataframe of responses to standalone items
## 2. Cluster_dat: 
# In general: a matrix or dataframe of responses to cluster items. One assertion per column. Column order must match row order in Cluster_parm
# For one student: row matrix/dataframe preferred but it can simply be a vector of responses
## 3. SA_parm: a matrix or dataframe of a,b,g parameters (column must be in this order), ItemID, and Assertion_ID for SA items
## 4. Cluster_parm: a matrix or dataframe of a,b and variance parameters for each assertion, a column of cluster position, a column of cluster ItemID, and a column of Assertion_ID for Cluster items
## 5. Dv: scaling factor for IRT model [1 or 1.7]
## 6. n.nodes: number of nodes used when integrating out the specific dimension 
## 7. censor: when there's perfect or all wrong score, this is a 2-element vector of upper and lower limit of 
##      1): starting value in point estimate optimization
##      2): pseudo parameter estimate used in SE computation to avoid extremely large SE estimate
## 8. correction_val: a value to add or subtract when there's perfect or all wrong score to avoid extremely large theta estimate
## 9. SE: if TRUE, return standard error
##
# *** Its okay to treat SA item as clusters. To do so, simply store them in the "Cluster_parm" argument with 0 variances, and store all student responses in "Cluster_dat"

# *** DONOTUSE:
  # SA_dat and Cluster_dat can be multi-row matrix/dataframe where it is one person per row. In that case, it simultanousely solves
  # all thetas but it assumes thetas are not independent. So it shouldn't be used for scoring purpose.
##########################################################################################################
scoring = function(SA_dat=NULL, Cluster_dat=NULL, SA_parm=NULL, Cluster_parm=NULL, Dv=1, n.nodes = 21, censor=c(-4,4), correction_val = 0.5, SE=FALSE) {
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
  possible = ncol(combined_dat) + sum(!is.na(SA_parm[,c(-1,-2,-(ncol(SA_parm)-2):-ncol(SA_parm))])) # 
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
  } else if (prod(rowSums(combined_dat, na.rm = TRUE) == ncol(combined_dat))) {
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