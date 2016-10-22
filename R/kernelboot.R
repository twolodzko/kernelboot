
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
#'                     If provided as numeric value it must be: single numeric value for univariate data;
#'                     vector of numbers for multivariate data (length equal to \code{ncol(data)}.
#' @param kernel       a character string giving the smoothing kernel to be used.
#' @param preserve.var logical, if \code{TRUE}, then the bootstrap samples preserve sample variance.
#' @param adjust       the bandwidth used is actually \code{adjust*bw}. This makes it easy to specify
#'                     values like ‘half the default’ bandwidth.
#' @param weights      Vector or matrix of importance weights. It should have as many
#'                     elements as there are observations in \code{data}.
#'
#' @references
#' Silverman, B. W. (1986). Density estimation for statistics and data analysis. Chapman and Hall/CRC.
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

kernelboot <- function(data, statistic, R = 500,
                       bw = bw.silv86(data),
                       kernel = c("epanechnikov", "gaussian", "rectangular",
                                  "triangular", "biweight", "triweight",
                                  "cosine", "optcosine"), preserve.var = TRUE,
                       adjust = 1, weights, ...) {

  call <- match.call()
  kernel <- match.arg(kernel)
  bw <- bw*adjust

  data <- as.matrix(data)
  n <- nrow(data)

  tryCatch(
    orig.stat <- statistic(data, ...),
    error = function(e) {
      message("Applying the statistic on the original data resulted in an error")
      stop(e)
    }
  )

  if (preserve.var) {

    means <- apply(data, 2, mean)
    vars <- apply(data, 2, var)

    res <- replicate(R, {

      idx <- sample.int(n, n, replace = TRUE)
      data.new <- add_noise(data[idx, ], kernel, bw,
                            mean = means, var = vars,
                            preserve_var = preserve.var)

      statistic(data.new, ...)

    })

  } else {

    res <- replicate(R, {

      idx <- sample.int(n, n, replace = TRUE)
      data.new <- add_noise(data[idx, ], kernel, bw,
                            preserve_var = FALSE)

      statistic(data.new, ...)

    })

  }

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
      weights      = if (missing(weights)) "uniform" else weights
    )
  ), class = "kernelboot")

}

