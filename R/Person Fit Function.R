####################  PERSON FIT FUNCTION ##########################################
## INPUT: 
## 1. theta: a vector of thetas or a scalar value of theta
## 2. SA_dat: (use NA for missing reponses)
# For one student, a vector of response
# For more than one student: a matrix or dataframe of response to standalone items. One assertion per column. Column order must match row order in SA_parm
## 3. Cluster_dat: (use NA for missing reponses)
# For one student, a vector of response
# For more than one student: a matrix or dataframe of student response cluster items. One assertion per column. Column order must match row order in Cluster_parm
## 4. SA_parm: a matrix or dataframe of a,b,g parameters (column must be in this order), ItemID, and Assertion_ID for SA items
## 5. Cluster_parm: a matrix or dataframe of a,b and variance parameters for each assertion, a column of cluster position, a column of cluster ItemID, and a column of Assertion_ID for Cluster items
## 6. Dv: scaling factor for IRT model [1 or 1.7]
## 7. n.nodes: number of nodes used when integrating out the specific dimension 
# *** Its okay to treat SA item as clusters. To do so, simply store them in the "Cluster_parm" argument with 0 variances, and store all student responses in "Cluster_dat"
# *** This function can be run for all students in one shot. Not recommended when number of students is large. It's faster to run one student at a time using parallel processing 
##########################################################################################################
source("J:\\NGSS\\Zhongtian\\R code\\Person fit\\function\\Observed Data Marginal Likelihood Function.R")
source("J:\\NGSS\\Zhongtian\\R code\\Item info\\LW Function.R")
person.fit = function(theta, SA_dat=NULL, Cluster_dat=NULL, SA_parm=NULL, Cluster_parm=NULL, Dv=1, n.nodes = 21) {
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")} 
  if(is.null(SA_dat) & is.null(Cluster_dat)) {stop("No data found!!!")} 
  if((is.null(SA_parm) | is.null(SA_dat)) && (is.null(Cluster_parm) | is.null(Cluster_dat)))  {stop("Data or Parameter file missing!!!")} 
  # -------------------------------------------------
  # ------------ Standalone Part --------------------
  # -------------------------------------------------
  if(is.null(SA_parm) | is.null(SA_dat)) {
    loglik_SA = EXP_LL_SA = VAR_LL_SA = RAW_SA = 0
  } else {
    SA_dat = matrix(as.matrix(SA_dat), nrow = length(theta))
    SA_parm = as.data.frame(SA_parm)
    names(SA_parm) = c("a","b","g", "ItemID", "Assertion_ID")
    # run observed data loglikelihood function (with return additional = T)
    loglik_output_SA = obs.data.loglik(theta, SA_dat=SA_dat, Cluster_dat=NULL, SA_parm=SA_parm, Cluster_parm=NULL, Dv=Dv, n.nodes = n.nodes, return_additional = T)
    loglik_SA = loglik_output_SA$loglik$loglik
    probs.SA = loglik_output_SA$probs.SA
    EXP_LL_SA = rowSums(probs.SA * log(probs.SA) + (1 - probs.SA) * log(1 - probs.SA))
    VAR_LL_SA = rowSums(probs.SA * (1 - probs.SA) * (log(probs.SA/(1 - probs.SA)))^2)
    RAW_SA = rowSums(SA_dat, na.rm = T)
  }
  
  # -------------------------------------------------
  # --------------- Cluster Part --------------------
  # -------------------------------------------------
  if(is.null(Cluster_parm) | is.null(Cluster_dat)) {
    loglik_cluster = EXP_LL_cluster = VAR_LL_cluster = RAW_cluster = 0
  } else {
    gq = statmod::gauss.quad.prob(n.nodes, dist = 'normal', sigma = 1)
    nodes = gq$nodes
    whts = gq$weights
    
    Cluster_parm = as.data.frame(Cluster_parm)
    colnames(Cluster_parm) = c("a","b","cluster_var","position", "ItemID", "Assertion_ID")
    Cluster_parm$position = match(Cluster_parm$position, sort(unique(Cluster_parm$position))) # reorder the position so it starts from 1
    
    Cluster_dat = matrix(as.matrix(Cluster_dat), nrow = length(theta))
    colnames(Cluster_dat) = paste0(Cluster_parm$ItemID, "_", Cluster_parm$Assertion_ID)
    loglik_cluster = EXP_LL_cluster = VAR_LL_cluster = list()
    for (k in 1:length(unique(Cluster_parm$position))) {
      one_cluster_parm = filter(Cluster_parm, position == k)
      one_cluster_dat = Cluster_dat[,grep(paste0("^",unique(one_cluster_parm$ItemID)), colnames(Cluster_dat))]
      loglik_output_one_cluster = obs.data.loglik(theta, SA_dat=NULL, Cluster_dat=one_cluster_dat, SA_parm=NULL, Cluster_parm=one_cluster_parm, Dv=Dv, n.nodes = n.nodes, return_additional = F)
      loglik_cluster[[k]] = loglik_output_one_cluster$loglik
      
      mvars = one_cluster_parm$cluster_var
      rescaled.nodes = nodes %o% sqrt(mvars) # rescaled nodes for these assertions
      ma = one_cluster_parm$a # a parameters for these assertions
      mb = one_cluster_parm$b # b parameters for these assertions
      mtheta = theta  # students' theta in a vector
      n.ass = length(mb)
      # Execute the Lord and Wingersky function
      LW_results = LW(cluster_var = mvars, a = ma, b = mb, theta = mtheta, n.nodes = n.nodes, return_additional = T, Dv=Dv)
      prk.marginal = LW_results$prk.marginal
      probs = LW_results$probs
      qrobs = LW_results$qrobs
      
      EXP_LL_cluster_term1 = probs * outer(matrix(replicate(n.nodes,mtheta), nrow = length(mtheta)), mb, "-")
      EXP_LL_cluster_term1 = apply(EXP_LL_cluster_term1, c(1,2), sum) %*% whts
      
      rawscores = 0:n.ass
      rawscoresXnodes = 0:n.ass %o% rescaled.nodes[,1]  # polish this later
      
      temp1 = apply(log(qrobs), c(1,2), sum)
      temp2 = lapply(1:nrow(rawscoresXnodes), function(x) sweep(temp1, 2, rawscoresXnodes[x,], "+"))
      temp3 = do.call(cbind, lapply(temp2, function(x) log(exp(x) %*% whts)))
      EXP_LL_cluster_term2 = rowSums(temp3 * t(prk.marginal))
      
      EXP_LL_cluster[[k]] = as.vector(EXP_LL_cluster_term1 + EXP_LL_cluster_term2)
      
      # LL variance calculation: treat clusters as SAs and run observed data loglikelihood function (with return additional = T)
      one_cluster_parm_temp = select(one_cluster_parm, -position) %>% rename(g=cluster_var) %>% mutate(g=0)
      loglik_output_one_cluster_as_SA = obs.data.loglik(theta, SA_dat=one_cluster_dat, Cluster_dat=NULL, SA_parm=one_cluster_parm_temp, Cluster_parm=NULL, Dv=Dv, n.nodes = n.nodes, return_additional = T)
      probs.one_cluster_as_SA = loglik_output_one_cluster_as_SA$probs.SA
      VAR_LL_cluster[[k]] = rowSums(probs.one_cluster_as_SA * (1 - probs.one_cluster_as_SA) * (log(probs.one_cluster_as_SA/(1 - probs.one_cluster_as_SA)))^2)
    }
    loglik_cluster = rowSums(matrix(simplify2array(loglik_cluster), nrow = length(mtheta)))
    EXP_LL_cluster = rowSums(matrix(simplify2array(EXP_LL_cluster), nrow = length(mtheta)))
    VAR_LL_cluster = rowSums(matrix(simplify2array(VAR_LL_cluster), nrow = length(mtheta)))
    RAW_cluster = rowSums(Cluster_dat, na.rm = T)
  }
  
  # -------------------------------------------------
  # ------------ Combine two parts ------------------
  # -------------------------------------------------
  loglik = as.vector(loglik_SA) + loglik_cluster
  EXP_LL = EXP_LL_SA + EXP_LL_cluster
  VAR_LL = VAR_LL_SA + VAR_LL_cluster
  z_ll = (loglik - EXP_LL)/sqrt(VAR_LL)
  pfit_flag = ifelse(abs(z_ll)>3, 1, 0)
  RAW = RAW_SA + RAW_cluster
  
  # -------------------------------------------------
  # -------------------- Output ---------------------
  # -------------------------------------------------
  output = data.frame(theta=theta, raw = RAW, ll=loglik, exp_ll=EXP_LL, var_ll=VAR_LL, z_ll=z_ll, pfit_flag=pfit_flag)
  return(output)
}
