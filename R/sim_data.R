####################  Data simulation function ##########################################
## Description: 
## Simulate IRT Model data. Models currently supported are 1-3 PL, GPCM, Rasch Testlet, 
## and a mix of these models. When Rasch Testlet Model is mixed with other unidimensonal models,
## the model is similar to a Rasch Testlet model but allows for standalone items that loads only on the overall dimension
##
## INPUT: 
## 1. thetas: 
##      When tests only contain standalone items: 
##          a vector of theta when tests only contain standalone item, where the length of thetas equals to number
##          of students. Works for a scalar theta value.
##      When tests contain cluster items:
##          In general: a matrix or dataframe of theta values when tests contain cluster item, where first column is 
##                      theta for the general dimension and the rest of columns are thetas for the specific (nuisance) dimensions
##          Special case: when simulating for only one student and the test has cluster item, it can also be a vector of 
##                        theta values, where the first element is general dimension theta and the rest are
##                        thetas for the specific (nuisance) dimensions. This is just for convenience. 
## 2. SA_parm: a matrix or dataframe of item parameters, where
##  a: slope
##  b1,b2...b_k: for 3PL items, b1 is the diffiuclty parameter. All other b parameters should be set to NA
##               for GPCM items, b1,b2,...,b_k are step difficulty parameters where k is the maximum possible score
##                            of item that has the largest number of score categories in the test. score starts from 0.
##                            E.g., for two items with 3 and 4 categories (0,1,2 & 0,1,2,3), there are b1, b2, b3, where
##                            b3 for the first item should be set to NA.
##  g: gussing parameters for 3PL items. Set to NA for GPCM items.
##  ItemID: Item IDs for SA items 
##  AssertionID: Assertion IDs for SA items
##  ***Note that columns must be in the above order
## 3. Cluster_parm: a matrix or dataframe of a,b and variance parameters for each assertion, 
##    a column of cluster position, a column of cluster ItemID, and a column of AssertionID for Cluster items.
##    Columns must be in the above order
## 4. Dv: scaling factor for IRT model [1 or 1.7]
# *** If the test does not have SA items or Cluster items, use default (NULL) for the corresponding parameter argument
# *** Rasch SA items can be treated as clusters. To do so, store b parameters in "Cluster_parm" and set variances to 0
##########################################################################################################
sim_data = function(thetas, SA_parm=NULL, Cluster_parm=NULL, Dv=1) {
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")}
  
  if (!is.null(Cluster_parm)) Cluster_parm = as.data.frame(Cluster_parm)
  if (!is.null(SA_parm)) SA_parm = as.data.frame(SA_parm)
  
  if(is.vector(thetas)) {
    if(!is.null(Cluster_parm)) {
      if (length(thetas) != (length(unique(Cluster_parm[,4]))+1)) {
        stop("Number of dimension does not match between theta and Cluster_parm")
      } else {
        thetas = matrix(thetas, nrow = 1)
      }
    } else {
      thetas = as.matrix(thetas)
    }
  } else {
    thetas = as.matrix(thetas)
    if (ncol(thetas) != (length(unique(Cluster_parm[,4]))+1)) {
      stop("Number of dimension does not match between theta and Cluster_parm")
    }
  }
  if(is.vector(thetas) & !is.null(Cluster_parm) & length(thetas)==length(unique(Cluster_parm[,4]))) {
    thetas = matrix(thetas, nrow = 1)
  } else {
    thetas = as.matrix(thetas)
  }
  theta = thetas[,1] #overall dimension theta
  if(is.null(SA_parm)) {
    data.SA = NULL
  } else {
    SA_parm = as.data.frame(SA_parm)
    names(SA_parm)[1] = c("a")
    names(SA_parm)[(ncol(SA_parm)-2):ncol(SA_parm)] = c("g", "ItemID", "AssertionID")
    names(SA_parm)[2:(ncol(SA_parm)-3)] = paste0("b",2:(ncol(SA_parm)-3)-1)
    # +++ Add SA_parm validity checks later
    # --- Separate out 3PL and GPC items ---
    dich.pos = which(apply(is.na(SA_parm[grepl("^b",names(SA_parm))][,-1]), 1, prod) == 1)
    poly.pos = which(apply(is.na(SA_parm[grepl("^b",names(SA_parm))][,-1]), 1, prod) == 0)
    SA_parm_3pl = SA_parm[dich.pos,]
    SA_parm_gpc = SA_parm[poly.pos,]
    if (nrow(SA_parm_3pl) != 0) {
      a = SA_parm_3pl$a
      b = SA_parm_3pl$b1
      g = rep(1,length(theta)) %o% SA_parm_3pl$g
      lin_pred_SA = Dv * t(apply(outer(theta, b, "-"), 1, function(x) a*x))
      probs.SA.3pl = g + (1 - g) * plogis(lin_pred_SA)
      colnames(probs.SA.3pl) = SA_parm_3pl$AssertionID
      UNI = replicate(length(b),runif(length(theta)))
      data.SA.3pl = ifelse(UNI<probs.SA.3pl, 1, 0)
    } else {data.SA=NULL}
    
    if (nrow(SA_parm_gpc) != 0) {
      a.gpc = SA_parm_gpc[,1]
      b.gpc = SA_parm_gpc[,grepl("^b",names(SA_parm_gpc))]
      maxscr = rowSums(!is.na(b.gpc))
      b.gpc.list = split(b.gpc, seq(nrow(b.gpc)))
      b.gpc.list = lapply(b.gpc.list,function(x) x[!is.na(x)])
      # theta is a vector of student ability
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
      probs.SA.gpc = mapply(gpcm, rep(list(theta),nrow(SA_parm_gpc)), a.gpc, b.gpc.list, maxscr)
      data.SA.gpc = sapply(probs.SA.gpc, function(x) {
        temp = sign(t(apply(x, 1, cumsum)) - runif(length(theta)))
        apply(temp, 1, function(y) which(y==1)[1]) - 1
      })
      data.SA.gpc = matrix(data.SA.gpc, nrow = length(theta))
      colnames(data.SA.gpc) = SA_parm_gpc$AssertionID
    } else {data.SA.gpc = NULL}
    data.SA = cbind(data.SA.3pl, data.SA.gpc)
  }
  
  if(is.null(Cluster_parm)) {
    data.cluster = NULL
  } else {
    colnames(Cluster_parm) = c("a","b","cluster_var","position", "ItemID", "AssertionID")
    Cluster_parm$position = match(Cluster_parm$position, sort(unique(Cluster_parm$position))) # reorder the position so it starts from 1
    data.cluster = list()
    for (k in 1:length(unique(Cluster_parm$position))) {
      one_cluster_parm = filter(Cluster_parm, position == k)
      mvars = unique(one_cluster_parm$cluster_var)
      All_thetas = theta + thetas[,1+k]
      ma = one_cluster_parm$a
      mb = one_cluster_parm$b
      lin_pred_cluster = Dv * t(apply(outer(All_thetas, mb, "-"), 1, function(x) ma*x))
      prob.cluster = plogis(lin_pred_cluster)
      colnames(prob.cluster) = one_cluster_parm$AssertionID
      UNI = replicate(length(mb),runif(length(theta)))
      data.cluster[[k]] = ifelse(UNI<prob.cluster, 1, 0)
    }
    data.cluster = do.call(cbind, data.cluster)
  }
  
  
  # -------------------------------------------------
  # -------------------- Output ---------------------
  # -------------------------------------------------
  output = cbind(data.SA, data.cluster)
  return(output)
}
