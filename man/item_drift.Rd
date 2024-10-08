% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/item_drift.R
\name{item_drift}
\alias{item_drift}
\title{Item drift analysis}
\usage{
item_drift(
  theta,
  SA_dat = NULL,
  Cluster_dat = NULL,
  SA_parm = NULL,
  Cluster_parm = NULL,
  drift = 0.3,
  Alpha = 0.05,
  Dv = 1,
  n.nodes = 21,
  missing_as_incorrect = FALSE
)
}
\arguments{
\item{theta}{a vector of examinee ability. The length of the vector should be equal to the number of examinees}

\item{SA_dat}{For one examinee, a vector of response to standalone items.
For more than one examinee, a matrix or dataframe of response to standalone items. One assertion per column.
Column order must match row order in \code{SA_parm}. Use NA for missing responses}

\item{Cluster_dat}{For one examinee, a vector of response to cluster items.
For more than one examinee, a matrix or dataframe of response to cluster items. One assertion per column.
Column order must match row order in \code{Cluster_parm}. Use NA for missing responses.}

\item{SA_parm}{a matrix or dataframe of item parameters for standalone items, where columns are
a (slope), b1, b2, ..., b_k (difficulty or step difficulty), g (guessing), ItemID, and AssertionID.
Columns must follow the above order.
See \code{example_SA_parm} for an example. Use \code{?example_SA_parm} for detailed column descriptions}

\item{Cluster_parm}{a matrix or dataframe of item parameters for cluster items, where columns are
a (slope), b (difficulty), cluster variance, cluster position, ItemID, and AssertionID.
Columns must follow the above order.
See \code{example_Cluster_parm} for an example. Use \code{?example_Cluster_parm} for detailed column descriptions}

\item{drift}{a numeric scalar for the amount of item difficulty parameter drift tested}

\item{Alpha}{a numeric vector for one or more nominal type I error rates for flagging aberrant responses.}

\item{Dv}{scaling factor for IRT model (usually 1 or 1.7)}

\item{n.nodes}{number of nodes used when integrating out the nuisance dimension}

\item{missing_as_incorrect}{by default, missings (NAs) are treated as missing; if TRUE, missings are treated as incorrect}
}
\value{
a list of
\itemize{
\item \code{RR} = N (number of examinees) by n (number of assertions) dataframe of assertion residuals \cr
\item \code{LRR} = N by n dataframe of the lower bound assertion residuals\cr
\item \code{URR} = N by n dataframe of the upper bound assertion residuals\cr
If run for more than one examinees, returns an additional \code{drift_output} table in the list \cr
}
}
\description{
Conduct item difficulty drift analysis and flag items with potential drift
Model currently supported are 1-3 PL, GPCM, and Rasch Testlet, and a mix of these models.
}
\note{
If the test does not have SA items or Cluster items, use default (NULL) for the corresponding data and parameter arguments. \cr\cr
}
\examples{
data(example_SA_parm)
data(example_Cluster_parm)
sigma <- diag(c(1, sqrt(unique(example_Cluster_parm$cluster_var))))
mu <- rep(0, nrow(sigma))
thetas <- MASS::mvrnorm(7,mu,sigma)
thetas[,1] <- seq(-3,3,1) #overall dimension theta values
itmDat <- sim_data(thetas = thetas, SA_parm = example_SA_parm, Cluster_parm = example_Cluster_parm)
SA_dat <- itmDat[,1:20]
Cluster_dat <- itmDat[,-1:-20]
out_scoring <- list()
for (i in 1:nrow(itmDat)) {
  out_scoring[[i]] <- scoring(SA_dat[i,], Cluster_dat[i,], example_SA_parm, example_Cluster_parm, n.nodes = 11, SE=TRUE)
}
est_theta <- sapply(out_scoring, function(x) x$par)
rst <- item_drift(est_theta, SA_dat=SA_dat, Cluster_dat=Cluster_dat, SA_parm=example_SA_parm, Cluster_parm=example_Cluster_parm,
                  drift = 0.3, Alpha=0.05, Dv=1, n.nodes = 21, missing_as_incorrect = FALSE)
}
\author{
Zhongtian Lin \href{mailto:lzt713@gmail.com}{lzt713@gmail.com}
}
