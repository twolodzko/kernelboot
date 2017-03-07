
#' Smoothed bootstrap
#'
#' Smoothed bootstrap is an extension of standard bootstrap using kernel densities.
#'
#' @param data         vector, matrix, or data.frame. For non-numeric values standard bootstrap
#'                     is applied (see below).
#' @param statistic    a function that is applied to the \code{data}. The first argument of
#'                     the function will always be the original data. Any further arguments
#'                     can be passed to \code{statistic} through \dots argument.
#' @param R            the number of bootstrap replicates.
#' @param bw           the smoothing bandwidth to be used (see \code{\link[stats]{density}}).
#'                     The kernels are scaled such that this is the standard deviation,
#'                     or covariance matrix of the smoothing kernel. By default
#'                     \code{\link[stats]{bw.nrd0}} is used for univariate data,
#'                     and \code{\link{bw.silv}} is used for multivariate data.
#' @param kernel       a character string giving the smoothing kernel to be used.
#'                     This must partially match one of "gaussian", "rectangular",
#'                     "triangular", "epanechnikov", "biweight", "triweight", "cosine"
#'                     or "optcosine", with default "gaussian", and may be abbreviated.
#' @param adjust       scalar; the bandwidth used is actually \code{adjust*bw}. This makes
#'                     it easy to specify values like 'half the default' bandwidth.
#' @param weights      vector of importance weights. It should have as many elements
#'                     as there are observations in \code{data}. It defaults to uniform
#'                     weights.
#' @param preserve.var logical; if \code{TRUE} random generation algorithm preserves
#'                     variance of the variables (see \code{\link{ruvk}} for details).
#'                     This parameter is used only for univariate and product kernels.
#' @param ignore       vector of names of columns to be ignored during the smoothing phase of
#'                     bootstrap procedure (their values are not altered using random noise).
#' @param parallel     if \code{TRUE} uses parallel processing (see \code{\link[parallel]{mclapply}}).
#' @param mc.cores     number of cores used for parallel computing (see \code{\link[parallel]{mclapply}}).
#' @param \dots        optional arguments passed to \code{statistic}.
#'
#'
#' @details
#'
#' \emph{Smoothed bootstrap} (Efron, 1981; Silverman, 1986; Scott, 1992; Hall, DiCiccio
#' and Romano, 1989) is an extension of standard bootstrap procedure, where instead
#' of drawing samples with replacement from the empirical distribution, they are drawn
#' from kernel density estimate of the distribution.
#'
#' For smoothed bootstrap, points (in univariate case), or rows (in multivariate case), are drawn with
#' replacement, to obtain samples of size \eqn{n} from the initial dataset of size \eqn{n}, as with
#' standard bootstrap. Next, random noise from kernel density \eqn{K} is added to each of the drawn
#' values. The proceure is repeated \eqn{R} times and \code{statistic} is evaluated on each of the
#' samples.
#'
#' The noise is added \emph{only} to the numeric columns, while non-numeric columns (i.e.
#' \code{character}, \code{factor}, \code{logical}) are not altered. What follows, to the
#' non-numeric columns and columns listed in \code{ignore} standard bootstrap procedure
#' is applied.
#'
#' With multivariate data, when using \code{kernel = "gaussian"} and \code{bw} is a non-diagonal
#' matrix, multivariate Gaussian kernel is applied (see \code{\link{rmvn}} and \code{\link{rmvk}}).
#' When \code{kernel = "gaussian"} and \code{bw} is a diagonal matrix, or a vector, product kernel
#' is used (see \code{\link{rmvpk}}). In other cases, depending on the data, univariate, or product
#' kernels, are used.
#'
#' @return
#' An object of class \code{"kernelboot"}, i.e., a list with components including
#'
#' \tabular{ll}{
#' \code{orig.stat}          \tab  estimates from \code{statistic} on the original data, \cr
#' \code{boot.sample}        \tab  samples drawn, \cr
#' \code{call}               \tab  function call, \cr
#' \code{statistic}          \tab  actual \code{statistic} function that was used, \cr
#' \code{orig.data}          \tab  original data used for bootstrapping, \cr
#' \code{smoothed.variables} \tab  names of variables that were included in the smoothing phase
#'                                 (\code{NULL} by default); those are the numeric columns and
#'                                 the variables not mentioned in \code{ignore} parameter, \cr
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
#' \code{preserve.var}       \tab  logical; value of \code{preserve.var} parameter, \cr
#' \code{parallel}           \tab  logical; states if parallel computation was used.
#' }
#'
#'
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
#' @references
#' Efron, B. (1981). Nonparametric estimates of standard error: the jackknife,
#' the bootstrap and other methods. Biometrika, 589-599.
#'
#' @references
#' Hall, P., DiCiccio, T.J., and Romano, J.P. (1989). On smoothing and the bootstrap.
#' The Annals of Statistics, 692-704. \url{http://projecteuclid.org/euclid.aos/1176347135}
#'
#'
#' @seealso \code{\link{bandwidth}}, \code{\link[stats]{density}},
#'          \code{\link[stats]{bandwidth}}, \code{\link{dmvn}},
#'          \code{\link{duvk}}, \code{\link{dmvk}}, \code{\link{dmvpk}}
#'
#'
#' @examples
#'
#' kernelboot(mtcars, function(data) coef(lm(mpg ~ ., data = data)) , R = 250)
#' kernelboot(mtcars, function(data) median(data$mpg) , R = 250)
#'
#'
#' @importFrom stats bw.SJ bw.bcv bw.nrd bw.nrd0 bw.ucv
#' @importFrom parallel mclapply
#'
#' @export

kernelboot <- function(data, statistic, R = 500L, bw = "default", ...,
                       kernel = c("gaussian", "epanechnikov", "rectangular",
                                  "triangular", "biweight", "triweight",
                                  "cosine", "optcosine"),
                       weights = NULL, adjust = 1,
                       preserve.var = TRUE, ignore = NULL,
                       parallel = FALSE, mc.cores = getOption("mc.cores", 2L)) {

  call <- match.call()
  kernel <- match.arg(kernel)
  n <- NROW(data)
  m <- NCOL(data)
  smoothed_variables <- NULL

  if (!(is.vector(data) || is.data.frame(data) || is.matrix(data)))
    stop("data is not vector, data.frame, or matrix")

  if (is.character(bw)) {
    bw <- tolower(bw)
    if (bw == "default") {
      if (is.vector(data)) {
        bw <- bw.nrd0(data)
      } else {
        bw <- bw.silv(data)
      }
    } else {
      bw <- switch(bw, nrd0 = bw.nrd0(data), nrd = bw.nrd(data),
                   ucv = bw.ucv(data), bcv = bw.bcv(data), sj = ,
                   `sj-ste` = bw.SJ(data, method = "ste"),
                   `sj-dpi` = bw.SJ(data, method = "dpi"),
                   silv = bw.silv(data), scott = bw.scott(data),
                   stop("unknown bandwidth rule"))
    }
  }
  if (!is.numeric(bw))
    stop("non-numeric bw value")

  if (!is.vector(adjust))
    stop("adjust is not a scalar")

  bw <- bw * adjust[1L]

  if (!all(is.finite(bw)))
    stop("inappropriate values of bw")

  if (!is.null(weights)) {
    if (!all(is.finite(weights)))
      stop("inappropriate values of weights")
  }

  # try evaluating statistic() on the original data

  tryCatch(
    orig.stat <- statistic(data, ...),
    error = function(e) {
      message("applying the statistic on the original data resulted in an error")
      stop(e)
    }
  )

  # looping functions - paralell or lapply

  if (parallel && mc.cores > 1L) {
    repeatFun <- function(i, FUN, mc.cores) mclapply(i, FUN, mc.cores = mc.cores)
  } else {
    repeatFun <- function(i, FUN, mc.cores) lapply(i, FUN)
  }

  if (is.data.frame(data) || is.matrix(data)) {

    # data is data.frame or matrix

    num_cols <- is_numeric(data)
    ignored_cols <- FALSE

    if (!is.null(ignore)) {
      ignored_cols <- colnames(data) %in% ignore
      if (any(ignored_cols)) {
        msg <- paste(colnames(data)[ignored_cols | !num_cols], collapse = ", ")
        message(paste0("the following variables are ignored during smoothing phase: ", msg))
      }
    }

    incl_cols <- num_cols & !ignored_cols
    smoothed_variables <- colnames(data)[incl_cols]

    if (!any(incl_cols)) {

      # standard bootstrap

      res <- repeatFun(1:R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        boot.data <- data[idx, ]
        statistic(boot.data, ...)

      }, mc.cores = mc.cores)

    } else {

      # smoothed bootstrap

      data_mtx <- as.matrix(data[, incl_cols])

      if (is.null(weights))
        weights <- rep(1/n, n)

      if (kernel != "gaussian" || is.vector(bw) || is.diag(bw)) {

        # product kernel

        if (is.vector(bw)) {
          if (length(bw) == 1L)
            bw <- rep(bw, m)
        } else {
          if (!is.square(bw))
            stop("bw is not a square matrix")
          bw <- diag(bw)
        }

        bw <- bw[incl_cols]

        res <- repeatFun(1:R, function(i) {

          samp <- cpp_rmvpk(n, data_mtx, bw, weights, kernel)
          idx <- samp$boot_index
          boot.data <- data[idx, ]
          boot.data[, incl_cols] <- samp$samples
          statistic(boot.data, ...)

        }, mc.cores = mc.cores)

      } else {

        # MVN kernel

        if (qr(data_mtx)$rank < min(dim(data_mtx)))
          warning("data matrix is rank deficient")

        if (kernel != "gaussian") {
          kernel <- "gaussian"
          message("for multivariate data only Gaussian kernel is supported; defaulting to Gaussian")
        }

        if (is.vector(bw)) {
          if (length(bw) == 1L)
            bw <- diag(bw, nrow = ncol(data))
          else
            bw <- diag(bw)
        }

        bw <- as.matrix(bw)
        bw <- bw[incl_cols, incl_cols]
        bw_chol <- chol(bw)

        res <- repeatFun(1:R, function(i) {

          samp <- cpp_rmvk(n, data_mtx, bw_chol, weights, is_chol = TRUE)
          idx <- samp$boot_index
          boot.data <- data[idx, ]
          boot.data[, incl_cols] <- samp$samples
          statistic(boot.data, ...)

        }, mc.cores = mc.cores)

      }
    }

  } else if (is.vector(data)) {

    # data is a vector

    if (!is.numeric(data)) {

      # standard bootstrap

      res <- repeatFun(1:R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        boot.data <- data[idx]
        statistic(boot.data, ...)

      }, mc.cores = mc.cores)

    } else {

      # smoothed bootstrap

      if (!is.vector(bw))
        stop("bw is not a scalar")
      if (length(bw) > 1L) {
        bw <- bw[1L]
        message("bw has length > 1 and only the first element will be used")
      }

      if (is.null(weights))
        weights <- rep(1/n, n)

      res <- repeatFun(1:R, function(i) {

        samp <- cpp_ruvk(n, data, bw, weights, kernel, preserve.var)
        boot.data <- drop(samp$samples)
        statistic(boot.data, ...)

      }, mc.cores = mc.cores)

    }

  } else {

    stop("unsupported data type")

  }

  structure(list(
    orig.stat          = orig.stat,
    boot.sample        = do.call(rbind, res),
    call               = call,
    statistic          = statistic,
    orig.data          = data,
    smoothed.variables = smoothed_variables,
    param = list(
      R            = R,
      bw           = bw,
      adjust       = adjust,
      weights      = weights,
      kernel       = kernel,
      preserve.var = preserve.var,
      parallel     = parallel
    )
  ), class = "kernelboot")

}

