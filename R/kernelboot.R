
#' Kernel density bootstrap
#'
#' @param data
#' @param statistic   A function which when applied to data returns a vector containing
#'                    the statistic(s) of interest. The first argument passed will always
#'                    be the original data. Any further arguments can be passed to
#'                    \code{statistic} through the \code{...} argument.
#' @param R           The number of bootstrap replicates.
#' @param bw          the smoothing bandwidth to be used. The kernels are scaled such that
#'                    this is the standard deviation, or covariance matrix of the smoothing kernel.
#'                    If provided as numeric value it must be: single numeric value for univariate data;
#'                    vector of numbers for multivariate data (length equal to \code{ncol(data)},
#'                    not counting the excluded columns declared in \code{exclude});
#'                    or matrix of numbers for \code{mvn} kernel (number of columns and rows
#'                    equal to \code{ncol(data)}, not counting the excluded columns declared
#'                    in \code{exclude}).
#' @param kernel      a character string giving the smoothing kernel to be used.
#' @param adjust      the bandwidth used is actually \code{adjust*bw}. This makes it easy to specify
#'                    values like ‘half the default’ bandwidth.
#' @param weights     Vector or matrix of importance weights. It should have as many
#'                    elements as there are observations in \code{data}.
#' @param exclude     a character vector naming columns \emph{not} to be sampled from the kernel
#'                    density. Regular bootstrap is applied to those columns.
#'
#' @export

kernelboot <- function(data, statistic, R = 500,
                       bw = c("default", "hns", "silverman",
                              "bw.nrd0", "bw.nrd", "bw.ucv",
                              "bw.bcv", "bw.SJ"),
                       kernel = c("mvn", "gaussian"),
                       #covariance = c("uniform", "variance", "covariance"),
                       adjust = 1, weights, exclude, ...) {

  call <- match.call()
  bw <- match.arg(bw)
  kernel <- match.arg(kernel)

  tryCatch(
    orig.stat <- statistic(data, ...),
    error = function(e) {
      message("Applying the statistic on the original data resulted in an error")
      stop(e)
    }
  )

  if (is.matrix(data) || is.data.frame(data)) {

    if (bw == "default") {
      message("switching bandwidth estimation method to 'hns'")
      bw <- "hns"
    } else if (!is.numeric(bw) && bw %in% c("bw.nrd0", "bw.nrd", "bw.ucv",
                                     "bw.bcv", "bw.SJ")) {
      stop("only the 'hns' and 'silvernam' methods for estimating bandwidth are supported for multivariate data")
    }

    n <- nrow(data)
    k <- ncol(data)

    if (!missing(exclude))
      kernel.cols <- !(colnames(data) %in% exclude)
    else
      kernel.cols <- rep(TRUE, k)

    mx <- colMeans(data)
    sx <- cov(data)
    H <- NULL

  } else {

    if (bw == "default") {
      message("switching bandwidth estimation method to 'nrd0'")
      bw <- "nrd0"
    } else if (!is.numeric(bw) && bw %in% c("hns", "silverman")) {
      stop("only the 'bw.nrd0', 'bw.nrd', 'bw.ucv', 'bw.bcv', 'bw.SJ' methods for estimating bandwidth are supported for multivariate data")
    }

    n <- length(data)

    mx <- mean(x)
    sx <- var(x)
    H <- NULL

    if (kernel == "mvn") {
      message("'mvn' kernel is not available for univariate data, 'gaussian' kernel is going to be used instead")
      kernel <- "gaussian"
    }


  }

  kdeFun <- switch (kernel,
    "gaussian" = function() rnorm(n*k, sd = H),
    function() mvtnorm::rmvnorm(n, sigma = H)
  )

  res <- replicate(R, {

    idx <- sample.int(n, n, replace = TRUE)
    data.boot <- data[idx, ]

    if (any(kernel.cols)) {
      # KDE sim
    }

    statistic(data.boot, ...)

  })

  structure(list(
    orig.stat   = orig.stat,
    boot.sample = res,
    call = call,
    statistic = statistic,
    param = list(
      R           = R,
      bw          = bw,
      adjust      = adjust,
      H           = H,
      kernel      = kernel,
      weights     = if (missing(weights)) "uniform" else weights,
      exclude     = if (missing(exclude)) NA else exclude
    )
  ), class = "kernelboot")

}

