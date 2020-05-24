####################  Marginal probability of correct response ##########################################
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
## 8. missing_as_incorrect: by default, missings (NAs) are treated as missing; if TRUE, missings are treated as incorrect 
# *** Its okay to treat SA item as clusters. To do so, simply store them in the "Cluster_parm" argument with 0 variances, and store all student responses in "Cluster_dat"
##########################################################################################################
item.score.residual = function(theta, SA_dat=NULL, Cluster_dat=NULL, SA_parm=NULL, Cluster_parm=NULL, Dv=1, n.nodes = 21, missing_as_incorrect = F) {
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")} 
  if(is.null(SA_dat) & is.null(Cluster_dat)) {stop("No data found!!!")} 
  
  if(!is.null(SA_dat)) SA_dat = matrix(as.matrix(SA_dat), nrow = length(theta))
  if(!is.null(Cluster_dat)) Cluster_dat = matrix(as.matrix(Cluster_dat), nrow = length(theta))
  combined_dat = cbind(SA_dat, Cluster_dat)
  if(any(rowSums(is.na(combined_dat)) == ncol(combined_dat))) {stop("one or more students did not respond to any item!!!")}
  if(missing_as_incorrect == T) combined_dat[is.na(combined_dat)] = 0  # recode if missing is treated as incorrect
  
  marginal_prob = marginal.prob(theta=theta, SA_parm=SA_parm, Cluster_parm=Cluster_parm, Dv=Dv, n.nodes = n.nodes)
  
  cbind(theta = marginal_prob[,1], combined_dat - marginal_prob[,-1])
}
