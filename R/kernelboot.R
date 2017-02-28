
#' Kernel density bootstrap
#'
#' @param data         Data.
#' @param statistic    A function which when applied to data returns a vector containing
#'                     the statistic(s) of interest. The first argument passed will always
#'                     be the original data. Any further arguments can be passed to
#'                     \code{statistic} through the \code{...} argument.
#' @param R            The number of bootstrap replicates.
#' @param bw           the smoothing bandwidth to be used. The kernels are scaled such that
#'                     this is the standard deviation, or covariance matrix of the smoothing kernel.
#'                     If missing, by default \code{\link[stats]{bw.nrd0}} is used for univariate data,
#'                     and \code{\link{bw.scott}} is used for multivariate data.
#' @param kernel       a character string giving the smoothing kernel to be used.
#' @param preserve.var logical, if \code{TRUE}, then the bootstrap samples preserve sample variance.
#' @param adjust       the bandwidth used is actually \code{adjust*bw}. This makes it easy to specify
#'                     values like 'half the default' bandwidth.
#' @param weights      Vector of importance weights. It should have as many
#'                     elements as there are observations in \code{data}.
#' @param parallel     if \code{TRUE} uses parallel processing (see \code{\link[parallel]{mclapply}}).
#' @param mc.cores     number of cores used for parallel computing (see \code{\link[parallel]{mclapply}}).
#' @param \dots        further arguments passed to \code{statistic}.
#'
#' @details
#'
#' Samples are drawn from kernel density using the following procedure (Silverman, 1986):
#'
#' \emph{Step 1} Sample \eqn{i} uniformly with replacement from \eqn{1,\dots,n}.
#'
#' \emph{Step 2} Generate \eqn{\varepsilon}{\epsilon} to have probability density \eqn{K}.
#'
#' \emph{Step 3} Set \eqn{Y = X_i + h\varepsilon}{Y = X[i] + h\epsilon}.
#'
#' If samples are required to have the same variance as \code{data}
#'  (\code{preserve.var} is set to \code{TRUE}), then \emph{Step 3} is modified
#' as following:
#'
#' \emph{Step 3'} \eqn{Y = \hat X + (X_i - \hat X + h\varepsilon)/(1 + h^2 \sigma^2_K/\sigma^2_X)^{1/2}}{Y = m + (X[i] - m + h\epsilon)/(1 + h^2 var(K)/var(X))^(1/2)}
#'
#' @references
#' Silverman, B. W. (1986). Density estimation for statistics and data analysis.
#' Chapman and Hall/CRC.
#'
#' @references
#' Wand, M. P. and Jones, M. C. (1995). Kernel Smoothing. Chapman and Hall/CRC.
#'
#' @references
#' Scott, D. W. (1992). Multivariate density estimation: theory, practice,
#' and visualization. John Wiley & Sons.
#'
#' @seealso \code{\link{bw.scott}}, \code{\link[stats]{density}},
#'          \code{\link[stats]{bw.nrd}}, \code{\link[boot]{boot}}
#'
#' @importFrom stats bw.nrd0
#' @importFrom parallel mclapply
#' @export

kernelboot <- function(data, statistic, R = 500L, bw,
                       kernel = c("gaussian", "epanechnikov", "rectangular",
                                  "triangular", "biweight", "triweight",
                                  "cosine", "optcosine"), preserve.var = TRUE,
                       adjust = 1, weights = NULL,
                       parallel = FALSE, mc.cores = getOption("mc.cores", 2L),
                       ...) {

  call <- match.call()
  kernel <- match.arg(kernel)
  n <- NROW(data)
  k <- NCOL(data)

  if (!(is.vector(data) || is.data.frame(data) || is.matrix(data)))
    stop("data must be a vector, data.frame, or matrix.")

  if (missing(bw)) {
    if (is.vector(data)) {
      if (n < 2)
        stop("need at least 2 points to select a bandwidth automatically")
      bw <- bw.nrd0(data)
    } else {
      bw <- bw.scott(data)
    }
  }
  if (!is.numeric(bw))
    stop("non-numeric bw value")

  bw <- bw * adjust

  if (is.vector(bw)) {
    if (length(bw) > k) {
      bw <- bw[1:k]
      warning("bw has length > number of dimensions of the data")
    }
    if (any(bw <= 0))
      stop("bw is not positive.")
  }

  if (!all(is.finite(bw)))
    stop("non-finite bw")

  if (!is.null(weights)) {
    if (length(weights) != n)
      stop("data and weights have unequal sizes")
    if (!all(is.finite(weights)))
      stop("weights must all be finite")
    if (any(weights < 0))
      stop("weights must not be negative")
  }

  tryCatch(
    orig.stat <- statistic(data, ...),
    error = function(e) {
      message("Applying the statistic on the original data resulted in an error")
      stop(e)
    }
  )

  if (mc.cores > 1 && parallel) {
    repeatFun <- function(i, FUN, mc.cores) mclapply(i, FUN, mc.cores = mc.cores)
  } else {
    repeatFun <- function(i, FUN, mc.cores) lapply(i, FUN)
  }

  if (is.data.frame(data) || is.matrix(data)) {

    num_cols <- is_numeric(data)

    if (length(num_cols) == 0) {

      res <- repeatFun(1:R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        boot.data <- data[idx, ]
        statistic(boot.data, ...)

      }, mc.cores = mc.cores)

    } else {

      if (is.null(weights))
        weights <- rep(1/n, n)

      if (kernel != "gaussian") {
        kernel <- "gaussian"
        warning("for multivariate data only Gaussian kernel is supported; defaulting to Gaussian")
      }

      data_mtx <- as.matrix(data[, num_cols])
      if (mtxrank(data_mtx) < min(n, k))
        warning("x is rank deficient")

      bw <- bw[num_cols, num_cols]
      bw_chol <- chol(bw)

      res <- repeatFun(1:R, function(i) {

        samp <- cpp_rmvkde(n, data_mtx, bw_chol, weights, is_chol = TRUE)
        idx <- samp$boot_index
        boot.data <- data[idx, ]
        boot.data[, num_cols] <- samp$sample
        statistic(boot.data, ...)

      }, mc.cores = mc.cores)

    }

  } else if (is.vector(data)) {

    if (is.numeric(data)) {

      res <- repeatFun(1:R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        boot.data <- data[idx]
        statistic(boot.data, ...)

      }, mc.cores = mc.cores)

    } else {

      if (is.null(weights))
        weights <- rep(1/n, n)

      res <- repeatFun(1:R, function(i) {

        samp <- cpp_ruvkde(n, data, bw, weights, kernel, preserve.var)
        boot.data <- drop(samp$sample)
        statistic(boot.data, ...)

      }, mc.cores = mc.cores)

    }

  } else {

    stop("unsupported data type")

  }

  structure(list(
    orig.stat   = orig.stat,
    boot.sample = do.call(rbind, res),
    call        = call,
    statistic   = statistic,
    param = list(
      R            = R,
      bw           = bw,
      adjust       = adjust,
      kernel       = kernel,
      preserve.var = preserve.var,
      weights      = weights,
      parallel     = parallel
    )
  ), class = "kernelboot")

}

