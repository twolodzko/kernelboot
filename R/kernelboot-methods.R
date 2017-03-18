

#' kernelboot class object
#'
#' @param x          \code{kernelboot} class object.
#' @param quantiles  returned quantiles (see \code{\link{quantile}}).
#' @param \dots      further arguments passed to or from other methods.
#'
#' @details
#'
#' Object of class \code{"kernelboot"}, is a list with components including
#'
#' \tabular{ll}{
#' \code{orig.stat}          \tab  estimates from \code{statistic} on the original data, \cr
#' \code{boot.samples}       \tab  samples drawn, \cr
#' \code{call}               \tab  function call, \cr
#' \code{statistic}          \tab  actual \code{statistic} function that was used, \cr
#' \code{orig.data}          \tab  original data used for bootstrapping, \cr
#' \code{variables}          \tab  used variables: it is \code{NULL} for univariate data and
#'                                 for multivaiate data it contains two lists of \code{smoothed}
#'                                 and \code{ignored} variables (names or column indexes) during
#'                                 the smoothing phase. \cr
#' \code{type}               \tab  type of kernel density that was used, \cr
#' \code{param}              \tab  list of parameters that were used.
#' }
#'
#' \code{param} section contains:
#'
#' \tabular{ll}{
#' \code{R}                  \tab  number of bootstrap iterations, \cr
#' \code{bw}                 \tab  the bandwidth that was used, \cr
#' \code{weights}            \tab  vector of the weights that were applied, \cr
#' \code{kernel}             \tab  name of the kernel that was used, \cr
#' \code{shrinked}           \tab  value of the \code{shrinked} parameter, \cr
#' \code{parallel}           \tab  indicates if parallel computation was used.
#' }
#'
#' @seealso \code{\link{kernelboot}}
#'
#' @name kernelboot-class
NULL


#' @rdname kernelboot-class
#' @export

summary.kernelboot <- function(x, quantiles = c(0.025, 0.5, 0.975), ...) {
  samp <- getSamples(x, simplify = TRUE)
  if (is.data.frame(samp) || is.matrix(samp)) {
    t(apply(samp, 2, function(x) {
      c(mean = mean(x),
        sd = sd(x),
        quantile(x, quantiles))
    }))
  } else summary(samp)
}


#' @export

print.kernelboot <- function(x, ...) {

  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "")
  cat("\n\n")
  cat(length(x$boot.samples), " samples were generated", sep = "")
  if (x$type == "none") {
    cat(" using standard bootstrap. ", sep = "")
  } else {
    cat(" from ", x$type, " kernel density with ",
        x$param$kernel, " kernel. ", sep = "")
  }

  if (!is.null(length(x$variables))) {
    if (length(x$variables$smoothed) > 0L && !is.numeric(x$variables$smoothed)) {
      cat("The following columns were smoothed: ",
          paste(sprintf("'%s'", x$variables$smoothed), collapse = ", "), ". ", sep = "")
    }
    k <- length(x$variables$ignored)
    if (k > 0L) {
      if (k == 1L)
        cat("1 column was ignored during the smoothing phase. ")
      else
        cat(k, "columns were ignored during the smoothing phase. ")
    }
  }

}


#' @param simplify  if \code{TRUE} object is simplified to data.frame.
#'
#' @rdname kernelboot-class
#' @export

getSamples <- function(x, simplify = TRUE) {
  samp <- x$boot.samples
  if (simplify)
    return(do.call(rbind, samp))
  samp
}


#' @rdname kernelboot-class
#' @export

as.data.frame.kernelboot <- function(x, ...) {
  as.data.frame(getSamples(x, simplify = TRUE))
}

