#' Item fit statistics and fit plot
#' @description Computing item fit statistics including Pearson's Chi-squared test,
#' G2 logliklihood-ratio test, Phi coefficient difference, as well as fit plots
#' Model currently supported are 1-3 PL, and Rasch Testlet, and a mix of these models.
#' Use \code{what} argument to specify one or more results to return. See \strong{Value} section for details.
#' @section Dependencies:
#' \describe{
#'  \item{statmod, dplyr}
#' }
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
#' @param n.nodes number of nodes used when integrating out the nuisance dimension
#' @param est_theta a vector of estimated theta. The length of the vector should be equal to the number of students
#' @param what a character vector specifying what to compute. Possible values are "X2", "G2", "phi", and "fitplot".
#' See \strong{Value} section for details.
#' @param sd_overall variance of the (overall) theta
#' @param missing_code code for missing response
#' @note If the test does not have SA items or Cluster items, use default (NULL) for the corresponding data and parameter arguments. \cr\cr
#' Rasch SA items can be treated as clusters. To do so, store SA item parameters in the \code{Cluster_parm} argument with 0 variances.
#'
#' @return A list of results specified in \code{what}. Possible results in the list are: \cr
#'  \code{X2} A list of 2. The first and second element are lower triangle matrices of Chi-square values and
#'  standardized Chi-square values, respectively, for every pair of assertions \cr
#'  \code{G2} A lower triangle matrix of G2 values for every pair of assertions \cr
#'  \code{phi} A lower triangle matrix of the phi coefficient difference between
#'  the observed and expected responses, for every pair of assertions \cr
#'  \code{fitplot} Observed and expected ICC plot for every assertion
#'
#' @author Zhongtian Lin lzt713@gmail.com
#' @examples
#' data(example_SA_parm)
#' data(example_Cluster_parm)
#' theta <- seq(-3,3,1)
#' rst <- utility(theta, example_SA_parm, example_Cluster_parm, n.nodes = 11)
#' @export
item_fit <- function(SA_dat=NULL,Cluster_dat=NULL,SA_parm=NULL,Cluster_parm=NULL,Dv=1,n.nodes=21,
                     est_theta = NULL, what=c("X2","G2","phi","fitplot"), sd_overall=1,missing_code=NULL,
                     fitplot_arg = list(binsize=50)) {
  # -------------------------------------------------
  # ------------ Data Cleansing ---------------------
  # -------------------------------------------------
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")}
  if(is.null(SA_dat) & is.null(Cluster_dat)) {stop("No data found!!!")}
  if((is.null(SA_parm) | is.null(SA_dat)) && (is.null(Cluster_parm) | is.null(Cluster_dat)))  {stop("Data or Parameter file missing!!!")}

  if(!is.null(SA_dat)) {
    if(is.vector(SA_dat)) SA_dat = matrix(SA_dat, nrow = 1) else SA_dat = as.matrix(SA_dat)
    if (ncol(SA_dat) != nrow(SA_parm)) {stop("Number of items does not match between data and parameter file for standalone items")}
  }
  if(!is.null(Cluster_dat)) {
    if(is.vector(Cluster_dat)) Cluster_dat = matrix(Cluster_dat, nrow = 1) else Cluster_dat = as.matrix(Cluster_dat)
    if (ncol(Cluster_dat) != nrow(Cluster_parm)) {stop("Number of items does not match between data and parameter file for cluster items")}
  }
  if (nrow(combined_dat)!=1 & nrow(combined_dat)<30) warning("Number of students less than 30. Results may not be accurate!")

  combined_dat = cbind(SA_dat, Cluster_dat)
  if (any(rowSums(is.na(combined_dat)) == ncol(combined_dat))) {stop("one or more students did not respond to any item!!!")}
  if (nrow(combined_dat)<30) warning("Number of students less than 30. Results may not be accurate!")
  combined_dat[combined_dat == missing_code] <- NA

  if(!is.null(SA_parm)) {
    SA_parm = as.data.frame(SA_parm)
    names(SA_parm)[1] = c("a")
    names(SA_parm)[(ncol(SA_parm)-2):ncol(SA_parm)] = c("g", "ItemID", "AssertionID")
    names(SA_parm)[2:(ncol(SA_parm)-3)] = paste0("b",2:(ncol(SA_parm)-3)-1)
    colnames(SA_dat) = paste0(SA_parm$ItemID, "_", SA_parm$AssertionID)
  }
  if(!is.null(Cluster_parm)) {
    Cluster_parm = as.data.frame(Cluster_parm)
    colnames(Cluster_parm) = c("a","b","cluster_var","position", "ItemID", "AssertionID")
    Cluster_parm$position = match(Cluster_parm$position, sort(unique(Cluster_parm$position))) # reorder the position so it starts from 1
    colnames(Cluster_dat) = paste0(Cluster_parm$ItemID, "_", Cluster_parm$AssertionID)
  }

  # -------------------------------------------------
  # ------ Pairwise Item Fit Statistics -------------
  # -------------------------------------------------
  # Get combination of all pairs of assertions in data, and identify which ones are from the same item
  n.ass = ncol(combined_dat)
  temp_parm=dplyr::bind_rows(SA_parm,Cluster_parm)
  All_AS_Combo = combn(1:n.ass,2, simplify = FALSE)
  All_AS_Combo_names = sapply(All_AS_Combo, function(x) paste0(colnames(combined_dat)[x], collapse = "|"))
  All_AS_Combo_is_same = sapply(All_AS_Combo_names, function(x) {
    item.id = temp_parm$ItemID[grepl(x,temp_parm$AssertionID)]
    item.id[1]==item.id[2]
  })
  All_AS_Combo_final = data.frame(AS_Combo_names = All_AS_Combo_names, is_same = All_AS_Combo_is_same)
  All_AS_Combo_final$col_num1 = sapply(All_AS_Combo, function(x) x[1])
  All_AS_Combo_final$col_num2 = sapply(All_AS_Combo, function(x) x[2])

  # Observed probs to populate the contingency table from a pair of assertions
  combnXcell_obs <- lapply(All_AS_Combo,
                           function(x) {
                             temp <- apply(combined_dat[,x], 1, paste0, collapse="")
                             temp <- temp[!grepl("NA", temp)]
                             c(table(temp)/length(temp), n.stu=length(temp))
                           })
  combnXcell_obs <- as.data.frame(do.call(dplyr::bind_rows, combnXcell_obs))
  combnXcell_obs[is.na(combnXcell_obs)] = 0
  rownames(combnXcell_obs) = All_AS_Combo_names
  obs.list=list()
  for (i in 1:(ncol(combnXcell_obs)-1)) {
    obs.list[[i]] = matrix(NA, n.ass, n.ass)
    obs.list[[i]][lower.tri(obs.list[[i]])] <- combnXcell_obs[,i]
  }
  n.stu_matrix = matrix(NA, n.ass, n.ass)
  n.stu_matrix[lower.tri(n.stu_matrix)] <- combnXcell_obs$n.stu

  # Expected probs to populate the contingency table cells from a pair of assertions
  # First marginalized over nodes of the nuisance dimension
  # Then marginalized over nodes of the overall dimension
  gq = statmod::gauss.quad.prob(n.nodes, dist = 'normal', sigma = 1)
  theta = gq$nodes * sd_overall
  whts = gq$weights
  rst <- utility(theta, SA_parm, Cluster_parm, Dv=Dv, n.nodes = n.nodes, what = "prob")

  # Convert SA item probs to the same format as Cluster item probs
  if(!is.null(nrow(rst$prob$probs_SA_3pl))) {
    cond_probs_SA=lapply(1:nrow(rst$prob$probs_SA_3pl), function(x) {
      xx = rep(rst$prob$probs_SA_3pl[x,], each=n.nodes)
      xx = matrix(xx, ncol = n.nodes, byrow = TRUE)
      xx
    })
  } else (cond_probs_SA=list(NULL))
  # Cluster item conditional probs
  if(!is.null(rst$prob$probs_CL)) {
    cond_probs_CL = lapply(1:n.nodes, function (k) do.call(rbind, lapply(rst$prob$probs_CL, function(x) x[[k]])))
  } else (cond_probs_CL=list(NULL))
  # Combined conditional probs
  cond_probs = Map(rbind, cond_probs_SA, cond_probs_CL)

  foo <- function(x, indat, whts) {
    is_same_cluster = All_AS_Combo_final$is_same[All_AS_Combo_final$col_num1==x[1] & All_AS_Combo_final$col_num2==x[2] ]
    if (is_same_cluster) {
      t(data.frame(cell_00 = (1 - indat[,x[1]])* (1 - indat[,x[2]]),
                   cell_01 = (1 - indat[,x[1]])* indat[,x[2]],
                   cell_10 = indat[,x[1]]* (1 - indat[,x[2]]),
                   cell_11 = indat[,x[1]] * indat[,x[2]])) %*% whts
    } else {
      mprob1 = indat[,x[1]] %*% whts
      mprob2 = indat[,x[2]] %*% whts
      t(data.frame(cell_00 = (1 - mprob1) * (1 - mprob2),
                   cell_01 = (1 - mprob1) * mprob2,
                   cell_10 = mprob1 * (1 - mprob2),
                   cell_11 = mprob1 * mprob2))
    }
  }
  combnXcell_exp = lapply(cond_probs,
                          function(k) t(do.call(cbind, lapply(All_AS_Combo, foo, indat=t(k), whts=whts))))
  combnXcell_exp = sapply(1:4, function(i) do.call(cbind, lapply(combnXcell_exp, function(x) x[,i])) %*% whts)
  combnXcell_exp = matrix(combnXcell_exp, ncol = 4)
  exp.list=list()
  for (i in 1:ncol(combnXcell_exp)) {
    exp.list[[i]] = matrix(NA, n.ass, n.ass)
    exp.list[[i]][lower.tri(exp.list[[i]])] <- combnXcell_exp[,i]
  }

  # ===================================================
  # Pearson's Chi-squared and standardized Chi-square
  if ("X2" %in% what) {
    X2 <- n.stu_matrix * Reduce("+", lapply(1:4, function(x) (obs.list[[x]] - exp.list[[x]])^2 / exp.list[[x]]))
    X2_std <- (X2 - 1)/sqrt(2) # std = sqrt(2)*df (df=1)
    X2 = list(X2=X2, X2_std=X2_std)
  } else {X2=NULL}
  # ===================================================
  # G-test (G2 logliklihood-ratio test)
  if ("G2" %in% what) {
    G2 <- 2 * n.stu_matrix * Reduce("+", lapply(1:4, function(x){
      outdat = obs.list[[x]] * log(obs.list[[x]]/exp.list[[x]])
      replace(outdat, is.nan(outdat), 0) }))
  } else {G2=NULL}
  # ===================================================
  # Phi coefficient and phi difference. See wiki page
  if ("phi" %in% what) {
    phi1 <- sqrt(X2$X2/n.stu_matrix)
    phi_exp <- (exp.list[[1]] * exp.list[[4]] - exp.list[[2]] * exp.list[[3]])/sqrt((exp.list[[1]]+exp.list[[2]])*(exp.list[[3]]+exp.list[[4]])*(exp.list[[1]]+exp.list[[3]])*(exp.list[[2]]+exp.list[[4]]))
    phi_obs <- (obs.list[[1]] * obs.list[[4]] - obs.list[[2]] * obs.list[[3]])/sqrt((obs.list[[1]]+obs.list[[2]])*(obs.list[[3]]+obs.list[[4]])*(obs.list[[1]]+obs.list[[3]])*(obs.list[[2]]+obs.list[[4]]))
    phi_diff <- phi_obs - phi_exp
  } else {phi_diff=NULL}

  # -------------------------------------------------
  # ------------- Fit plot --------------------------
  # -------------------------------------------------
  if ("fitplot" %in% what) {
    combined_dat2 = as.data.frame(combined_dat)
    if(!is.null(est_theta)) {
      combined_dat2$est_theta = est_theta
    } else {
      # scoring(SA_dat=SA_dat,Cluster_dat=Cluster_dat,SA_parm=SA_parm,Cluster_parm=Cluster_parm,Dv=1,n.nodes=21)
    }

    p_list = list()
    binsize = fitplot_arg$binsize
    for (i in 1:(ncol(combined_dat2)-1)) {
      temp_dat = combined_dat2[,c(i, ncol(combined_dat2))]
      temp_dat = temp_dat[!is.na(temp_dat[,1]),]
      temp_dat = transform(temp_dat, bin = cut(temp_dat[,"est_theta"], binsize))
      pd = by(temp_dat, list(temp_dat$bin), function (x)
        c(meanVal = mean(x[,1], na.rm = TRUE),
          count = nrow(x)))
      pd=as.data.frame(do.call(rbind, pd))
      pd$bin = rownames(pd); rownames(pd) <- NULL
      pd$xval = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", pd$bin)),
                               as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", pd$bin))))

      # collapse bins with counts smaller than binsize
      if (collapsed) {
        pd0 = pd
        pd_temp = pd
        pd = data.frame(stringsAsFactors = FALSE)

        # the last collapsed bin sometime would have count < binsize
        count_temp = pd_temp$count
        while (sum(count_temp)>binsize | (sum(pd_temp$count) < binsize && sum(pd_temp$count) != 0 )) {
          if (sum(pd_temp$count) < binsize) {
            pd_chuck = pd_temp
            ind = nrow(pd_chuck)
          } else {
            ind = which(cumsum(count_temp) >50)[1]
            pd_chuck = pd_temp[1:ind,]
          }

          pt1 = as.numeric(sub("\\((.+),.*", "\\1", pd_chuck$bin[1]))
          pt2 = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", pd_chuck$bin[nrow(pd_chuck)]))
          nn = sum(pd_chuck$count)
          mm = t(pd_chuck$meanVal) %*% pd_chuck$count / nn
          bb = paste0("(",pt1,",",pt2,"]")
          xx = (pt1+pt2)/2

          pd = rbind(pd, c(mm,nn,bb,xx), stringsAsFactors = FALSE)
          count_temp = count_temp[-1:-ind]
          pd_temp = pd_temp[-1:-ind,]
        }
        names(pd) = names(pd0)
        pd$meanVal = as.numeric(pd$meanVal)
        pd$count = as.numeric(pd$count)
        pd$xval = as.numeric(pd$xval)
        pd$SE = sqrt(pd$meanVal*(1 - pd$meanVal)/pd$count)
      }

      # prepare data for expected ICC curve
      fixed_theta = seq(LOT, HOT, 0.1)
      if (names(temp_dat)[1] %in% SA_parm$AssertionID) {
        SA_parm_one = SA_parm[SA_parm$AssertionID==names(temp_dat)[1],]
        rst2 = utility(fixed_theta, SA_parm_one, NULL, Dv=Dv, n.nodes = n.nodes, what = "mprob")
        exp_icc = cbind(fixed_theta, rst2$mprob$probs_SA_3pl)
      }
      if (names(temp_dat)[1] %in% Cluster_parm$AssertionID) {
        Cluster_parm_one = Cluster_parm[Cluster_parm$AssertionID==names(temp_dat)[1],]
        rst2 = utility(fixed_theta, NULL, Cluster_parm_one, Dv=Dv, n.nodes = n.nodes, what = "mprob")
        exp_icc = cbind(fixed_theta, rst2$mprob$mprobs_CL)
      }
      exp_icc = as.data.frame(exp_icc)
      var_name = names(exp_icc)[2]
      names(exp_icc)[2] = "prob"

      # make plots
      p_list[[i]] = ggplot(pd, aes(x = xval, y = meanVal)) + geom_point() + geom_errorbar(aes(ymin=meanVal-1.96*SE, ymax=meanVal+1.96*SE)) +
        geom_line(data = exp_icc, aes(x = fixed_theta, y = prob), col = "red", lwd = 1) + xlab("Theta") + ylab("Probability/Proportion Correct") +
        ggtitle(var_name) + theme(plot.title = element_text(hjust = 0.5)) + xlim(c(-4, 4))

    }
  } else {fitplot=NULL}

  return(list(X2=X2, G2=G2, phi_diff=phi_diff, fitplot=p_list))

}

out = item_fit(SA_dat=SA_dat,Cluster_dat=Cluster_dat,SA_parm=SA_parm,Cluster_parm=Cluster_parm,Dv=1,n.nodes=21,
               est_theta = rnorm(nrow(SA_dat)), sd_overall=1,missing_code="-1")

data("example_SA_parm")
sigma <- diag(c(1, sqrt(unique(example_Cluster_parm$cluster_var))))
mu <- rep(0, nrow(sigma))
thetas <- MASS::mvrnorm(7,mu,sigma)
thetas[,1] <- seq(-3,3,1) #overall dimension theta values
itemdata <- sim_data(thetas = thetas, SA_parm = example_SA_parm, Cluster_parm = example_Cluster_parm)
SA_dat <- itemdata[,1:17]
Cluster_dat <- itemdata[,-1:-20]
example_SA_parm = example_SA_parm[1:17,]
out = item_fit(SA_dat=SA_dat,Cluster_dat=Cluster_dat,SA_parm=example_SA_parm,Cluster_parm=example_Cluster_parm,Dv=1,n.nodes=21,
               est_theta = rnorm(nrow(SA_dat)), sd_overall=1,missing_code="-1")
