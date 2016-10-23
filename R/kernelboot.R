
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
#'                     and \code{\link{bw.silv86}} is used for multivariate data.
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
#' @seealso \code{\link{bandwidth}}, \code{\link[stats]{density}},
#'          \code{\link[stats]{bw.nrd}}, \code{\link[boot]{boot}}
#'
#' @export

kernelboot <- function(data, statistic, R = 500, bw,
                       kernel = c("epanechnikov", "gaussian", "rectangular",
                                  "triangular", "biweight", "triweight",
                                  "cosine", "optcosine"), preserve.var = TRUE,
                       adjust = 1, weights = NULL,
                       parallel = FALSE, mc.cores = getOption("mc.cores", 2L),
                       ...) {

  call <- match.call()
  kernel <- match.arg(kernel)
  n <- NROW(data)
  d <- NCOL(data)

  if (!(is.vector(data) || is.data.frame(data) || is.matrix(data)))
    stop("'data' must be a vector, data.frame, or matrix.")

  if (missing(bw)) {
    if (is.vector(data)) {
      if (n < 2)
        stop("need at least 2 points to select a bandwidth automatically")
      bw <- bw.nrd0(data)
    } else {
      bw <- bw.silv86(data)
    }
  }
  if (!is.numeric(bw))
    stop("non-numeric 'bw' value")

  if (length(bw) > d) {
    bw <- bw[1:d]
    warning("'bw' has length > number of dimensions of the data")
  }

  if (!all(is.finite(bw)))
    stop("non-finite 'bw'")
  if (any(bw <= 0))
    stop("'bw' is not positive.")

  bw <- bw * adjust

  if (!is.null(weights)) {
    if (length(weights) != n)
      stop("'data' and 'weights' have unequal sizes")
    if (!all(is.finite(weights)))
      stop("'weights' must all be finite")
    if (any(weights < 0))
      stop("'weights' must not be negative")
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

    num_cols <- apply(data, 2, is.numeric)

    if (preserve.var) {

      mx <- colMeans(data)
      sx <- diag(cov(data))

      res <- repeatFun(1:R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        boot.data <- data[idx, ]
        boot.data[, num_cols] <- add_noise(as.matrix(boot.data[, num_cols]),
                                           kernel, bw, mean = mx, var = sx,
                                           preserve_var = preserve.var)

        statistic(boot.data, ...)

      }, mc.cores = mc.cores)

    } else {

      res <- repeatFun(1:R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        boot.data <- data[idx, ]
        boot.data[, num_cols] <- add_noise(as.matrix(boot.data[, num_cols]),
                                           kernel, bw, preserve_var = FALSE)

        statistic(boot.data, ...)

      }, mc.cores = mc.cores)

    }

  } else {

    rng_kern <- switch(kernel,
                       epanechnikov = rempan,
                       rectangular = rrect,
                       triangular = rtriang,
                       biweight = rbiweight,
                       triweight = rtriweight,
                       cosine = rcosine,
                       optcosine = roptcos,
                       rnorm)

    if (preserve.var) {

      mx <- mean(data)
      sx <- var(data)

      res <- repeatFun(1:R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        eps <- rng_kern(n) * bw
        boot.data <- mx + (data[idx]-mx+eps)/sqrt(1 + bw^2/sx)

        statistic(boot.data, ...)

      }, mc.cores = mc.cores)

    } else {

      res <- repeatFun(1:R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        eps <- rng_kern(n) * bw
        boot.data <- data[idx] + eps

        statistic(boot.data, ...)

      }, mc.cores = mc.cores)

    }

  }

  res <- do.call(rbind, res)

  structure(list(
    orig.stat   = orig.stat,
    boot.sample = res,
    call        = call,
    statistic   = statistic,
    param = list(
      R            = R,
      bw           = bw,
      adjust       = adjust,
      kernel       = kernel,
      preserve.var = preserve.var,
      weights      = if (missing(weights)) "uniform" else weights,
      parallel     = parallel
    )
  ), class = "kernelboot")

}

