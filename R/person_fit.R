#' Person fit index computation
#' @description Calculates lz* person fit index. Model currently supported are 1-3 PL, GPCM, Rasch Testlet, and a mix of these models.
#'              For traditional test (standalone items only), this is lz* described in Snijders (2001).
#'              For Cluster items, this is a extension of lz* calculated by a method
#'              that extends the lz* to Rasch testlet model when theta is estimated by MMLE.
#'
#' @section Dependencies:
#' \describe{
#'  \item{statmod}
#' }
#'
#' @param theta a scalar or a vector of student ability
#' @param SA_dat For one student, a vector of response to standalone items.
#' For more than one student, a matrix or dataframe of response to standalone items. One assertion per column.
#' Column order must match row order in \code{SA_parm}. Use NA for missing responses
#' @param Cluster_dat For one student, a vector of response to cluster items.
#' For more than one student, a matrix or dataframe of response to cluster items. One assertion per column.
#' Column order must match row order in \code{Cluster_parm}. Use NA for missing responses.
#' @param SA_parm a matrix or dataframe of item parameters for standalone items, where columns are
#'  a (slope), b1, b2, ..., b_k (difficulty or step difficulty), g (guessing), ItemID, and AssertionID.
#'  Columns must follow the above order.
#'  See \code{example_SA_parm} for an example. Use \code{?example_SA_parm} for detailed column descriptions
#' @param Cluster_parm a matrix or dataframe of item parameters for cluster items, where columns are
#'  a (slope), b (difficulty), cluster variance, cluster position, ItemID, and AssertionID.
#'  Columns must follow the above order.
#'  See \code{example_Cluster_parm} for an example. Use \code{?example_Cluster_parm} for detailed column descriptions
#' @param Dv scaling factor for IRT model (usually 1 or 1.7)
#' @param Alpha a numeric vector for one or more nominal type I error rates for flagging aberrant responses.
#' For example, \code{c(0.01,0.05)}
#' @param n.nodes number of nodes used when integrating out the nuisance dimension
#' @param n.ass.limit the number of assertions above which when the simulated assertion score would be used
#' instead of exact response pattern when calculating the variance of cluster item loglikelihood
#' @param LW.nodes same as n.nodes, but specifically used in the Lord-Wingersky algorithm for computing marginal probability
#' of raw scores
#'
#' @return a one-row dataframe with the following columns: \cr
#' \code{theta}\cr
#' \code{raw}, raw score \cr
#' \code{ll}, observed data loglikelihood based on MMLE \cr
#' \code{exp_ll}, expected value of the data loglikelihood \cr
#' \code{var_ll_org}, variance of the data loglikelihood before correction \cr
#' \code{correction}, correction of the variance of the data loglikelihood \cr
#' \code{var_ll}, variance of the data loglikelihood after correction \cr
#' \code{lz}, lz* person fit index \cr
#' \code{flag_[alpha level]}, one or more columns of flags for aberrant responses based on lz and the
#' type I error rates specified in \code{Alpha}
#'
#' @references Snijders, T. A. (2001). Asymptotic null distribution of person fit statistics with estimated person parameter.
#' \emph{Psychometrika}, 66(3), 331-342.
#'
#' @note If the test does not have SA items or Cluster items, use default (NULL) for the corresponding data and parameter arguments \cr\cr
#' Rasch SA items can be treated as clusters. To do so, store SA item parameters in the \code{Cluster_parm}
#' argument with 0 variances, and store all student responses in \code{Cluster_dat} \cr\cr
#' This function should be run for one student at a time. Use parallel processing for higher processing speed. \cr
#' !!!!! With correction version 3: dll based correction using numerical integration for the cluster correction term \cr
#' !!!!! Lord and Wingersky algorithm uses 21 nodes by default \cr
#' !!!!! When computing the original vll, 11 nodes were used for the term that involves all possible response patterns \cr
#' !!!!! If number of assertion > n.ass.limit (default 10), use simulation approach for original vll
#'
#' @author Zhongtian Lin lzt713@gmail.com
#' @examples
#' data(example_SA_parm)
#' data(example_Cluster_parm)
#' sigma <- diag(c(1, sqrt(unique(example_Cluster_parm$cluster_var))))
#' mu <- rep(0, nrow(sigma))
#' thetas <- MASS::mvrnorm(7,mu,sigma)
#' thetas[,1] <- seq(-3,3,1) #overall dimension theta values
#' itmDat <- sim_data(thetas = thetas, SA_parm = example_SA_parm, Cluster_parm = example_Cluster_parm)
#' SA_dat <- itmDat[,1:20]
#' Cluster_dat <- itmDat[,-1:-20]
#' # Scoring for the first student
#' out_scoring <- scoring(SA_dat[1,], Cluster_dat[1,], example_SA_parm, example_Cluster_parm, n.nodes = 11, SE=TRUE)
#' est_theta <- out_scoring$par
#' # person fit for the first student
#' rst <- person_fit(est_theta, SA_dat[1,], Cluster_dat[1,], example_SA_parm, example_Cluster_parm, Alpha=c(0.01,0.05))
#' @export
#' @useDynLib MIRTutils
person_fit = function(theta, SA_dat=NULL, Cluster_dat=NULL, SA_parm=NULL, Cluster_parm=NULL, Dv=1, Alpha=0.05, n.nodes = 100, n.ass.limit = 10, LW.nodes=21) {
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")}
  if(is.null(SA_dat) & is.null(Cluster_dat)) {stop("No data found!!!")}
  if((is.null(SA_parm) | is.null(SA_dat)) && (is.null(Cluster_parm) | is.null(Cluster_dat)))  {stop("Data or Parameter file missing!!!")}
  # -------------------------------------------------
  # ------------ Standalone Part --------------------
  # -------------------------------------------------
  if(is.null(SA_parm) | is.null(SA_dat)) {
    loglik_SA = EXP_LL_SA = VAR_LL_SA = RAW_SA = correction_SA = SA_info = 0
  } else {
    SA_dat = matrix(as.matrix(SA_dat), nrow = length(theta))
    loglik_output_SA = data_loglik(theta, SA_dat=SA_dat, Cluster_dat=NULL, SA_parm=SA_parm, Cluster_parm=NULL, Dv=Dv, n.nodes = n.nodes, return_additional = TRUE)
    loglik_SA = loglik_output_SA$loglik$loglik
    probs.SA.3pl = loglik_output_SA$probs.SA$probs.SA.3pl
    probs.SA.gpc = loglik_output_SA$probs.SA$probs.SA.gpc
    SA_parm_3pl = loglik_output_SA$parms$SA_parm_3pl
    SA_parm_gpc = loglik_output_SA$parms$SA_parm_gpc

    if (nrow(SA_parm_3pl) != 0) {
      a = SA_parm_3pl$a; b = SA_parm_3pl$b1; g = rep(1,length(theta)) %o% SA_parm_3pl$g
      EXP_LL_SA_3pl = rowSums(probs.SA.3pl * log(probs.SA.3pl) + (1 - probs.SA.3pl) * log(1 - probs.SA.3pl))
      VAR_LL_SA_3pl = rowSums(probs.SA.3pl * (1 - probs.SA.3pl) * (log(probs.SA.3pl/(1 - probs.SA.3pl)))^2)
      correction_SA_3pl = -rowSums( Dv*a *(probs.SA.3pl - g)/(1-g) * (1 - probs.SA.3pl) * log(probs.SA.3pl/(1 - probs.SA.3pl)))
      SA_info_3pl = sum((Dv*a)^2 * (1 - probs.SA.3pl)/probs.SA.3pl * ((probs.SA.3pl - g)/(1 - g))^2)
    } else {
      EXP_LL_SA_3pl = VAR_LL_SA_3pl = correction_SA_3pl = SA_info_3pl = 0
    }

    if (nrow(SA_parm_gpc) != 0) {
      EXP_LL_SA_gpc = Reduce("+", lapply(probs.SA.gpc, function(x) sum(x * log(x))))
      VAR_LL_SA_gpc = Reduce("+", lapply(probs.SA.gpc, function(x) sum((log(x) - sum(x * log(x)))^2 * x)))
      # VAR_LL_SA_gpc = Reduce("+", lapply(probs.SA.gpc, function(x) sum(sapply(x, function(y) y*log(y)*x*log(y/x)))))
      a.gpc = SA_parm_gpc$a
      SA_info_gpc = sum((Dv*a.gpc)^2 * sapply(probs.SA.gpc,function(x) (sum((seq_along(x)-1)^2*x) - sum((seq_along(x)-1)*x)^2)))
      # SA_info_gpc = sum((Dv*a.gpc)^2 * sapply(probs.SA.gpc,function(x) sum(((seq_along(x)-1) - sum((seq_along(x)-1)*x))^2 * x)))
      probs.SA.gpc.deriv0 = lapply(probs.SA.gpc,function(x)  x * ((seq_along(x)-1) - sum((seq_along(x)-1)*x)) )
      probs.SA.gpc.deriv = mapply("*", Dv*a.gpc, probs.SA.gpc.deriv0)
      correction_SA_gpc0 =  mapply(function(x,y) x * (log(y) + 1), probs.SA.gpc.deriv, probs.SA.gpc)
      correction_SA_gpc = -sum(unlist(correction_SA_gpc0))
    } else {
      EXP_LL_SA_gpc = VAR_LL_SA_gpc = SA_info_gpc = correction_SA_gpc = 0
    }

    EXP_LL_SA = EXP_LL_SA_3pl + EXP_LL_SA_gpc
    VAR_LL_SA = VAR_LL_SA_3pl + VAR_LL_SA_gpc
    RAW_SA = rowSums(SA_dat, na.rm = TRUE)
    correction_SA = correction_SA_3pl + correction_SA_gpc
    SA_info = SA_info_3pl + SA_info_gpc
  }

  # -------------------------------------------------
  # --------------- Cluster Part --------------------
  # -------------------------------------------------
  if(is.null(Cluster_parm) | is.null(Cluster_dat)) {
    loglik_cluster = EXP_LL_cluster = VAR_LL_cluster = RAW_cluster = correction_CL = cluster_info = 0
  } else {
    gq = statmod::gauss.quad.prob(n.nodes, dist = 'normal', sigma = 1)
    nodes = gq$nodes
    whts = gq$weights

    gq_small = statmod::gauss.quad.prob(11, dist = 'normal', sigma = 1)
    whts_small = gq_small$weights

    Cluster_parm = as.data.frame(Cluster_parm)
    colnames(Cluster_parm) = c("a","b","cluster_var","position", "ItemID", "Assertion_ID")
    Cluster_parm$position = match(Cluster_parm$position, sort(unique(Cluster_parm$position))) # reorder the position so it starts from 1

    Cluster_dat = matrix(as.matrix(Cluster_dat), nrow = length(theta))
    colnames(Cluster_dat) = paste0(Cluster_parm$ItemID, "_", Cluster_parm$Assertion_ID)
    loglik_cluster = EXP_LL_cluster = VAR_LL_cluster = correction_CL = cluster_info = EXP_LL_cluster_new = list()
    for (k in 1:length(unique(Cluster_parm$position))) {
      one_cluster_parm = Cluster_parm[Cluster_parm$position == k,]
      one_cluster_dat = Cluster_dat[,grep(paste0("^",unique(one_cluster_parm$ItemID)), colnames(Cluster_dat))]
      loglik_output_one_cluster = data_loglik(theta, SA_dat=NULL, Cluster_dat=one_cluster_dat, SA_parm=NULL, Cluster_parm=one_cluster_parm, Dv=Dv, n.nodes = n.nodes, return_additional = TRUE)
      loglik_cluster[[k]] = loglik_output_one_cluster$loglik$loglik
      probs = loglik_output_one_cluster$probs.cluster[[1]]
      qrobs = 1 - probs
      qrobs[qrobs==0] = .Machine$double.xmin

      mvars = one_cluster_parm$cluster_var
      rescaled.nodes = nodes %o% sqrt(mvars) # rescaled nodes for these assertions
      ma = one_cluster_parm$a # a parameters for these assertions
      mb = one_cluster_parm$b # b parameters for these assertions
      mtheta = theta  # students' theta in a vector
      n.ass = length(mb)
      # Execute the Lord and Wingersky function
      LW_results = lord_wing(cluster_var = mvars, a = ma, b = mb, theta = mtheta, n.nodes = LW.nodes, return_additional = TRUE, Dv=Dv)
      prk.marginal = LW_results$prk.marginal
      # probs = LW_results$probs
      # qrobs = LW_results$qrobs

      EXP_LL_cluster_term1_temp = apply(probs * outer(matrix(replicate(n.nodes,mtheta), nrow = length(mtheta)), mb, "-"), c(1,2), sum)
      EXP_LL_cluster_term1 = EXP_LL_cluster_term1_temp %*% whts

      rawscores = 0:n.ass
      rawscoresXnodes = rawscores %o% rescaled.nodes[,1]  # polish this later

      temp1 = apply(log(qrobs), c(1,2), sum)
      temp2 = lapply(1:nrow(rawscoresXnodes), function(x) sweep(temp1, 2, rawscoresXnodes[x,], "+"))
      temp3 = do.call(cbind, lapply(temp2, function(x) log(exp(x) %*% whts)))
      EXP_LL_cluster_term2 = rowSums(temp3 * t(prk.marginal))

      EXP_LL_cluster[[k]] = as.vector(EXP_LL_cluster_term1 + EXP_LL_cluster_term2)

      h = 1e-9
      mtheta_new = mtheta + h
      # Execute the Lord and Wingersky function again for theta + small value
      LW_results_new = lord_wing(cluster_var = mvars, a = ma, b = mb, theta = mtheta_new, n.nodes = LW.nodes, return_additional = T, Dv=Dv)
      prk.marginal_new = LW_results_new$prk.marginal

      loglik_output_one_cluster_new = data_loglik(mtheta_new, SA_dat=NULL, Cluster_dat=one_cluster_dat, SA_parm=NULL, Cluster_parm=one_cluster_parm, Dv=Dv, n.nodes = n.nodes, return_additional = TRUE)
      probs_new = loglik_output_one_cluster_new$probs.cluster[[1]]
      qrobs_new = 1 - probs_new
      qrobs_new[qrobs_new==0] = .Machine$double.xmin

      EXP_LL_cluster_term1_temp_new = apply(probs_new * outer(matrix(replicate(n.nodes,mtheta_new), nrow = length(mtheta_new)), mb, "-"), c(1,2), sum)
      EXP_LL_cluster_term1_new = EXP_LL_cluster_term1_temp_new %*% whts

      temp1_new = apply(log(qrobs_new), c(1,2), sum)
      temp2_new = lapply(1:nrow(rawscoresXnodes), function(x) sweep(temp1_new, 2, rawscoresXnodes[x,], "+"))
      temp3_new = do.call(cbind, lapply(temp2_new, function(x) log(exp(x) %*% whts)))
      EXP_LL_cluster_term2_new = rowSums(temp3_new * t(prk.marginal_new))

      EXP_LL_cluster_new[[k]] = as.vector(EXP_LL_cluster_term1_new + EXP_LL_cluster_term2_new)

      # Compute numerical differantiation of the correction term (from cluster items)
      correction_CL[[k]] = -(EXP_LL_cluster_new[[k]] - EXP_LL_cluster[[k]]) / h

      # Compute Expected CSEM
      ECSEM_term1 = t(sapply(temp2, exp))
      ECSEM_term2 = outer(rawscores, as.vector(apply(probs, c(1,2), sum)), "-")
      cluster_info[[k]] = -(t(-(((ECSEM_term1 * ECSEM_term2) %*% whts) / (ECSEM_term1 %*% whts))^2) %*% prk.marginal)


      # Equation 2 in Tao's write up
      EXP_squared_LL_term1 = EXP_LL_cluster_term1_temp^2 +
        apply(probs * qrobs* (outer(matrix(replicate(n.nodes,mtheta), nrow = length(mtheta)), mb, "-"))^2, c(1,2), sum)
      EXP_squared_LL_term1 = EXP_squared_LL_term1 %*% whts
      EXP_squared_LL_term3 = temp3^2 %*% prk.marginal

      if (n.ass<=n.ass.limit) { #if n.ass is small, use exact response patterns
        loglik_output_one_cluster_small = data_loglik(mtheta, SA_dat=NULL, Cluster_dat=one_cluster_dat, SA_parm=NULL, Cluster_parm=one_cluster_parm, Dv=Dv, n.nodes = 11, return_additional = TRUE)
        probs_small = loglik_output_one_cluster_small$probs.cluster[[1]]
        qrobs_small = 1 - probs_small
        qrobs_small[qrobs_small==0] = .Machine$double.xmin
        # +++ pattern score based computation for EXP_squared_LL (in Tao's 08182020 document)
        all_patterns = expand.grid(replicate(n.ass, 0:1, simplify = FALSE))
        probs_matrix = aperm(probs_small, c(3,2,1))[,,1]
        qrobs_matrix = aperm(qrobs_small, c(3,2,1))[,,1]
        all_patterns_all_nodes_prob = lapply(1:ncol(probs_matrix), function(x) mmult(as.matrix(all_patterns), probs_matrix[,x]))
        all_patterns_all_nodes_qrob = lapply(1:ncol(qrobs_matrix), function(x) mmult(as.matrix(1-all_patterns), qrobs_matrix[,x]))
        # all_patterns_all_nodes_prob = lapply(1:ncol(probs_matrix), function(x) t(t(all_patterns) * probs_matrix[,x]))
        # all_patterns_all_nodes_qrob = lapply(1:ncol(qrobs_matrix), function(x) t(t(1-all_patterns) * qrobs_matrix[,x]))
        # all_patterns_all_nodes_prob = lapply(1:ncol(probs_matrix), function(x) sweep(all_patterns, 2, probs_matrix[,x], "*"))
        # all_patterns_all_nodes_qrob = lapply(1:ncol(qrobs_matrix), function(x) sweep(1-all_patterns, 2, qrobs_matrix[,x], "*"))
        all_patterns_all_nodes_pq = Map("+", all_patterns_all_nodes_prob, all_patterns_all_nodes_qrob)
        all_patterns_pq_marginal = sapply(all_patterns_all_nodes_pq, Rfast::rowprods) %*% whts_small
        # all_patterns_pq_marginal = sapply(all_patterns_all_nodes_pq, function(x) apply(x, 1, prod)) %*% whts
        EXP_squared_LL_term2_multiplier2 = as.matrix(all_patterns) %*% (mtheta-mb)
        EXP_squared_LL_term2_multiplier3 = as.matrix(as.vector(temp3)[rowSums(all_patterns)+1])
        EXP_squared_LL_term2 = 2 * sum(all_patterns_pq_marginal * EXP_squared_LL_term2_multiplier2 * EXP_squared_LL_term2_multiplier3)
      } else { #if n.ass is large, use simulated response patterns
        EXP_squared_LL_term2_addend1 = EXP_LL_cluster_term1 * (temp3 %*% prk.marginal)
        set.seed(1234)
        thetas = cbind(mtheta,rnorm(5000, 0, sqrt(one_cluster_parm$cluster_var[1])))
        sim_scores = sim_data(thetas, SA_parm=NULL, Cluster_parm=one_cluster_parm)
        EXP_squared_LL_term2_addend2_p1_cor1 = as.vector(sim_scores %*% (mtheta-mb))
        EXP_squared_LL_term2_addend2_p1_cor2 = as.vector(temp3)[rowSums(sim_scores)+1]
        EXP_squared_LL_term2_addend2_p1 = cor(EXP_squared_LL_term2_addend2_p1_cor1, EXP_squared_LL_term2_addend2_p1_cor2)
        EXP_squared_LL_term2_addend2_p2 = sqrt(EXP_squared_LL_term1 - EXP_LL_cluster_term1^2)
        EXP_squared_LL_term2_addend2_p3 = sqrt(EXP_squared_LL_term3 - (temp3 %*% prk.marginal)^2)
        EXP_squared_LL_term2_addend2 = EXP_squared_LL_term2_addend2_p1 * EXP_squared_LL_term2_addend2_p2 *EXP_squared_LL_term2_addend2_p3
        EXP_squared_LL_term2 = 2 * (EXP_squared_LL_term2_addend1 + EXP_squared_LL_term2_addend2)
      }
      EXP_squared_LL = EXP_squared_LL_term1 + EXP_squared_LL_term2 + EXP_squared_LL_term3

      EXP_LL_squared = EXP_LL_cluster[[k]]^2
      VAR_LL_cluster[[k]] = EXP_squared_LL - EXP_LL_squared
    }
    loglik_cluster = rowSums(matrix(simplify2array(loglik_cluster), nrow = length(mtheta)))
    EXP_LL_cluster = rowSums(matrix(simplify2array(EXP_LL_cluster), nrow = length(mtheta)))
    VAR_LL_cluster = rowSums(matrix(simplify2array(VAR_LL_cluster), nrow = length(mtheta)))
    RAW_cluster = rowSums(Cluster_dat, na.rm = T)
    cluster_info = do.call(sum, cluster_info)
  }

  # -------------------------------------------------
  # ------------ Combine two parts ------------------
  # -------------------------------------------------
  test_info = SA_info + cluster_info
  loglik = as.vector(loglik_SA) + loglik_cluster
  EXP_LL = EXP_LL_SA + EXP_LL_cluster
  VAR_LL_original = VAR_LL_SA + VAR_LL_cluster
  VAR_LL_correction = (do.call(sum, as.list(correction_CL)) + correction_SA)^2 * (1/test_info)
  # +++
  VAR_LL = VAR_LL_original - VAR_LL_correction
  lz = (loglik - EXP_LL)/sqrt(VAR_LL)
  flag = data.frame(matrix(NA, nrow = 1, ncol = length(Alpha)))
  names(flag) = paste0("flag_", Alpha)
  for (i in seq_along(Alpha)) {
    flag[1,i] = ifelse(lz < qnorm(Alpha[i]), 1, 0)
  }
  RAW = RAW_SA + RAW_cluster

  # -------------------------------------------------
  # -------------------- Output ---------------------
  # -------------------------------------------------
  output = data.frame(theta=theta, raw = RAW, ll=loglik, exp_ll=EXP_LL, var_ll_org=VAR_LL_original, correction=VAR_LL_correction, var_ll=VAR_LL, lz=lz)
  output = cbind(output,flag)
  return(output)
}
