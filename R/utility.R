#' A utility function computing some commonly used values in IRT
#' @description Compute some commonly used values in IRT such as
#' probability of correct response, expected item scores, and item information.
#' Use \code{what} argument to specify one or more results to return. See \strong{Value} section for details. \cr\cr
#' Run '?MIRTutils-package' for more details, such as the context of the current package and models supported.
#'
#' @param theta a scalar or a vector of examinee ability
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
#' @param what a character vector specifying what to compute. Possible values are "prob", "escore", "iteminfo", and "mprob".
#' See \strong{Value} section for details.
#'
#' @note If the test does not have SA items or Cluster items, use default (NULL) for the corresponding data and parameter arguments. \cr\cr
#'
#' @return A list of results specified in \code{what}: \cr
#' * \code{prob}: probability of correct response. When Rasch testlet model is involved, this is conditional probability given both overall and each node of the nuisance dimension
#' * \code{escore}: expected item score
#' * \code{iteminfo}: item information. When Rasch testlet model is involved, this is item info based on marginal prob of raw scores calculated using Lord-Wingerskey algorithm
#' * \code{mprob}: When Rasch testlet model is involved, this is the marginal probability of correct response for each assertion. For unidimensional models, this is the same as prob
#'
#' @author Zhongtian Lin <lzt713@gmail.com>
#' @examples
#' data(example_SA_parm)
#' data(example_Cluster_parm)
#' theta <- seq(-3,3,1)
#' rst <- utility(theta, example_SA_parm, example_Cluster_parm, n.nodes = 11)
#' @export
utility <- function(theta, SA_parm=NULL, Cluster_parm=NULL, Dv=1, n.nodes = 50, what = c("prob","escore","iteminfo","mprob")){
  if(is.null(SA_parm) & is.null(Cluster_parm)) {stop("No item found!!!")}
  if ("mprob" %in% what && is.null(Cluster_parm)) {warning("No cluster item found. mprob=prob for unidimensional models")}
  # -------------------------------------------------
  # ------------ Standalone Part --------------------
  # -------------------------------------------------
  if(is.null(SA_parm)) {
    SA_parm_3pl = SA_parm_gpc =  NA
    probs_SA = escore_SA = info_SA = mprobs_SA = NULL
  } else {
    # +++ Add SA_parm validity checks later
    SA_parm = as.data.frame(SA_parm)
    names(SA_parm)[1] = c("a")
    names(SA_parm)[(ncol(SA_parm)-2):ncol(SA_parm)] = c("g", "ItemID", "AssertionID")
    names(SA_parm)[2:(ncol(SA_parm)-3)] = paste0("b",2:(ncol(SA_parm)-3)-1)

    # Separate out 3PL and GPC item par
    dich.pos = which(apply(as.matrix(is.na(SA_parm[grepl("^b",names(SA_parm))][,-1])), 1, prod) == 1)
    poly.pos = which(apply(as.matrix(is.na(SA_parm[grepl("^b",names(SA_parm))][,-1])), 1, prod) == 0)
    SA_parm_3pl = SA_parm[dich.pos,]
    SA_parm_gpc = SA_parm[poly.pos,]

    if (nrow(SA_parm_3pl) != 0) {
      a = SA_parm_3pl$a
      b = SA_parm_3pl$b1
      g = rep(1,length(theta)) %o% SA_parm_3pl$g
      lin_pred_SA = Dv * t(apply(outer(theta, b, "-"), 1, function(x) a*x))
      lin_pred_SA = matrix(lin_pred_SA, ncol = length(a))
      # conditional prob
      probs_SA_3pl = g + (1 - g) * plogis(lin_pred_SA)
      colnames(probs_SA_3pl) = SA_parm_3pl$AssertionID
      # item info
      if ("iteminfo" %in% what) {
        info_SA_3pl = t((Dv * a)^2 * t(1 - probs_SA_3pl))/probs_SA_3pl * ((probs_SA_3pl - g)/(1 - g))^2
      } else {
        info_SA_3pl = NULL
      }
      # Expected score
      if ("escore" %in% what) {
        escore_SA_3pl = probs_SA_3pl
      } else {
        escore_SA_3pl = NULL
      }
    } else {
      probs_SA_3pl = NA
      escore_SA_3pl = info_SA_3pl = NULL
    }
    if (nrow(SA_parm_gpc) != 0) {
      a.gpc = SA_parm_gpc[,1]
      b.gpc = SA_parm_gpc[,grepl("^b",names(SA_parm_gpc))]
      maxscr = rowSums(!is.na(b.gpc))
      b.gpc.list = split(b.gpc, seq(nrow(b.gpc)))
      b.gpc.list = lapply(b.gpc.list,function(x) x[!is.na(x)])
      # theta is a vector of examinee ability
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
      # conditional prob
      probs_SA_gpc = mapply(gpcm, rep(list(theta),nrow(SA_parm_gpc)), a.gpc, b.gpc.list, maxscr, Dv, SIMPLIFY = FALSE)
      names(probs_SA_gpc) = SA_parm_gpc$AssertionID
      # item info
      if ("iteminfo" %in% what) {
        info_SA_gpc = sapply(probs_SA_gpc,
                             function(y) apply(y, 1,
                                               function(x) (sum((seq_along(x)-1)^2*x) - sum((seq_along(x)-1)*x)^2)))
        info_SA_gpc_names = names(probs_SA_gpc)
        info_SA_gpc = matrix(info_SA_gpc,nrow = length(theta))
        colnames(info_SA_gpc) = info_SA_gpc_names
        info_SA_gpc = sweep(info_SA_gpc, 2, (Dv * a.gpc)^2, "*")
      } else {
        info_SA_gpc= NULL
      }
      # Expected score
      if ("escore" %in% what) {
        escore_SA_gpc = sapply(probs_SA_gpc, function(x) rowSums(sweep(x,2,0:(ncol(x)-1),"*")))
        escore_SA_gpc_names = names(probs_SA_gpc)
        escore_SA_gpc = matrix(escore_SA_gpc,nrow = length(theta))
        colnames(escore_SA_gpc) = escore_SA_gpc_names
      } else {
        escore_SA_gpc = NULL
      }
    } else {
      probs_SA_gpc = NA
      escore_SA_gpc = info_SA_gpc = NULL
    }
    probs_SA = list(probs_SA_3pl = probs_SA_3pl, probs_SA_gpc = probs_SA_gpc)
    escore_SA = cbind(escore_SA_3pl, escore_SA_gpc)
    info_SA = cbind(info_SA_3pl, info_SA_gpc)
    mprobs_SA = probs_SA
  }
  # -------------------------------------------------
  # --------------- Cluster Part --------------------
  # -------------------------------------------------
  if(is.null(Cluster_parm)) {
    probs_CL = info_CL = escore_CL = mprobs_CL = NULL
  } else {
    gq = statmod::gauss.quad.prob(n.nodes, dist = 'normal', sigma = 1)
    nodes = gq$nodes
    whts = gq$weights

    Cluster_parm = as.data.frame(Cluster_parm)
    colnames(Cluster_parm) = c("a","b","cluster_var","position", "ItemID", "AssertionID")
    Cluster_parm$position = match(Cluster_parm$position, sort(unique(Cluster_parm$position))) # reorder the position so it starts from 1

    probs_CL = info_CL = escore_CL = mprobs_CL = list()
    for (k in 1:length(unique(Cluster_parm$position))) {
      one_cluster_parm = Cluster_parm[Cluster_parm$position == k,]
      mvars = one_cluster_parm$cluster_var
      ma = one_cluster_parm$a # a parameters for these assertions
      mb = one_cluster_parm$b # b parameters for these assertions
      n.ass = length(mb)
      rescaled.nodes = nodes %o% sqrt(mvars)  # rescaled nodes for these assertions
      probs_CL[[k]] = lapply(theta, function(x) plogis(Dv * (ma * (x + t(rescaled.nodes) - mb))))
      # line above: list[[cluster]], list[person], row(assertion), col(node)

      # item info & escore
      if ("iteminfo" %in% what | "escore" %in% what) {
        LW_results = lord_wing(cluster_var = mvars, a = ma, b = mb, theta = theta, n.nodes = n.nodes, return_additional = T, Dv=Dv)
        # probs_CL[[k]] = lapply(seq(dim(LW_results$probs)[1]), function(x) LW_results$probs[ x, , ])
        # # line above: list[[cluster]], list[person], row(node), col(assertion)
        prk.marginal = sapply(1:length(LW_results$prk), function(x) as.vector(LW_results$prk[[x]]$W0[-1,] %*% LW_results$whts))
        probs = LW_results$probs
        qrobs = LW_results$qrobs
        rawscores = 0:n.ass
        rawscoresXnodes = 0:n.ass %o% rescaled.nodes[,1]
        temp1 = apply(log(qrobs), c(1,2), sum)
        temp2 = lapply(1:nrow(rawscoresXnodes), function(x) sweep(temp1, 2, rawscoresXnodes[x,], "+"))
        temp3 = exp(array(as.numeric(unlist(temp2)), dim=c(dim(temp2[[1]]), length(temp2))))  # every person, every node, every possible raw score
        temp4 = temp3 * outer(-apply(probs, c(1,2), sum), rawscores, "+")
        numer = apply(temp4, 3, function(x) x %*% whts)
        denom = apply(temp3, 3, function(x) x %*% whts)
        info_CL[[k]] = -rowSums(t(prk.marginal) * -(numer/denom)^2)
        escore_CL[[k]] = colSums(prk.marginal * 0:n.ass)
      }
      # marginal prob
      if ("mprob" %in% what) {
        mprobs_CL[[k]] = t(do.call(cbind, lapply(probs_CL[[k]], function(x) x %*% whts)))
      }
    }
    info_CL = do.call(cbind, info_CL)
    escore_CL = do.call(cbind, escore_CL)
    mprobs_CL = do.call(cbind, mprobs_CL)
    names(probs_CL) = unique(Cluster_parm$ItemID)
    if ("iteminfo" %in% what) colnames(info_CL) = unique(Cluster_parm$ItemID)
    if ("escore" %in% what) colnames(escore_CL) = unique(Cluster_parm$ItemID)
    if ("mprob" %in% what) colnames(mprobs_CL) = Cluster_parm$AssertionID
  }

  # Outputing Results
  prob = list(probs_SA_3pl=probs_SA$probs_SA_3pl, probs_SA_gpc=probs_SA$probs_SA_gpc, probs_CL = probs_CL)
  iteminfo = as.data.frame(cbind(theta, info_SA, info_CL))
  escore = as.data.frame(cbind(theta, escore_SA, escore_CL))
  mprob = list(probs_SA_3pl=mprobs_SA$probs_SA_3pl, probs_SA_gpc=mprobs_SA$probs_SA_gpc, mprobs_CL = mprobs_CL)

  return_val = list()
  if ("prob" %in% what) return_val$prob = prob
  if ("escore" %in% what) return_val$escore = escore
  if ("iteminfo" %in% what) return_val$iteminfo = iteminfo
  if ("mprob" %in% what) return_val$mprob = mprob
  return_val
}



