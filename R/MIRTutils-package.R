#' MIRT utility functions
#' @description A set of functions to calculate some commonly used values when doing IRT analysis
#' especially with new item types such as item bundles (testlets).
#' IRT model currently supported are 1-3 PL, GPCM, Rasch Testlet Model, and a mix of these models.
#' Values that can be computed include the probability of correct response, item information,
#' item score residuals, person scores, person fit index, etc.
#' It also provides a function for multidimensional IRT data simulation.
#' @details
#' In the context of the current package, test items can be either standalone items (SA) or cluster items (Cluster).
#' Both SA and cluster items consist of assertions, where an assertion is the same as a traditional item that was
#' typically modeled by a unidimensional IRT model. An SA item typically includes 1 to a handful of assertions,
#' whereas a cluster item typically includes 6 or more assertions. An SA item is more like a small set
#' of traditional items assuming local dependence among assertions, whereas a cluster is essentially a testlet.
#'
#'
#' The IRT model used is a special model similar to a Rasch Testlet model but allows for SA items to load only on the
#' overall dimension. When there are only standalone items, the model reduces to a regular unidimensional IRT model, and
#' all the values computed in this package fits to the traditional unidimensional IRT paradigm.
#' When there are only cluster items, the model reduces to a Rasch testlet model.
#' Unlike the unidimensional model where theta is usually estimated
#' by the univariate MLE, EAP or MAP, and the multidimensional model where theta is usually estimated by the multivarite
#' MLE or EAP, most of the values computed by this package involving testlets is based on the marginal maximum
#' likelihood estimation (MMLE) of theta. MMLE is a hybrid of MLE and EAP, where the nuisance dimension of a testlet (cluster)
#' is first integrated out, and the resulting marginal likelihood is then maximized to find the theta estimate.
#' @author Zhongtian Lin <lzt713@gmail.com>
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
