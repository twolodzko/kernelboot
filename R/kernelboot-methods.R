

#' kernelboot class object
#'
#' @param x,object   \code{kernelboot} class object.
#' @param probs      quantiles staturned by \code{summary} (see \code{\link{quantile}}).
#' @param \dots      further arguments passed to or from other methods.
#'
#' @details
#'
#' Function \code{getSamples} can be used to extract the bootstrap samples from
#' \code{"kernelboot"} class object. When calling \code{getSamples} with
#' \code{simplify = TRUE} (default) the samples are converted to a matrix, or
#' data.frame, otherwise they are returned as a list. Method \code{as.data.frame}
#' extracts samples from \code{"kernelboot"} class object converting them to
#' data.frame. \code{summary} prints numeric summary of the results if feasable.
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
#' @importFrom stats sd quantile
#' @name kernelboot-class
NULL


#' @rdname kernelboot-class
#' @export

summary.kernelboot <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  samp <- getSamples(object, simplify = TRUE)
  res <- lapply(1:ncol(samp), function(i) {
    x <- samp[, i]
    if (is.numeric(x)) {
      c(mean = mean(x), sd = sd(x), quantile(x, probs = probs))
    } else {
      warning("skipping non-numeric variable")
      NA
    }
  })
  names(res) <- colnames(samp)
  do.call(rbind, res)
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

  invisible(x)
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
  as.data.frame(getSamples(x, simplify = TRUE), ...)
}

