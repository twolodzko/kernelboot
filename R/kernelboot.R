
#' Smoothed bootstrap
#'
#' Smoothed bootstrap is an extension of standard bootstrap using kernel densities.
#'
#' @param data       vector, matrix, or data.frame. For non-numeric values standard bootstrap
#'                   is applied (see below).
#' @param statistic  a function that is applied to the \code{data}. The first argument of
#'                   the function will always be the original data. Any further arguments
#'                   can be passed to \code{statistic} through \dots argument.
#' @param R          the number of bootstrap replicates.
#' @param bw         the smoothing bandwidth to be used (see \code{\link[stats]{density}}).
#'                   The kernels are scaled such that this is the standard deviation,
#'                   or covariance matrix of the smoothing kernel. By default
#'                   \code{\link[stats]{bw.nrd0}} is used for univariate data,
#'                   and \code{\link{bw.silv}} is used for multivariate data. For using
#'                   multivariate gaussian kernel this parameter should be a \emph{covariance
#'                   matrix}.
#' @param kernel     a character string giving the smoothing kernel to be used.
#'                   This must partially match one of "gaussian", "rectangular",
#'                   "triangular", "epanechnikov", "biweight", "cosine"
#'                   or "optcosine", with default "gaussian", and may be abbreviated.
#' @param adjust     scalar; the bandwidth used is actually \code{adjust*bw}. This makes
#'                   it easy to specify values like 'half the default' bandwidth.
#' @param weights    vector of importance weights. It should have as many elements
#'                   as there are observations in \code{data}. It defaults to uniform
#'                   weights.
#' @param shrinked   logical; if \code{TRUE} random generation algorithm preserves
#'                   means and variances of the variables (see \code{\link{ruvk}} for details).
#'                   This parameter is used only for univariate and product kernels.
#' @param ignore     vector of names of columns to be ignored during the smoothing phase of
#'                   bootstrap procedure (their values are not altered using random noise).
#' @param parallel   if \code{TRUE} uses parallel processing.
#' @param workers    number of workers used for parallel computing.
#'
#'
#' @details
#'
#' \emph{Smoothed bootstrap} is an extension of standard bootstrap procedure, where instead
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
#' non-numeric columns and columns listed in \code{ignore} parmeter standard bootstrap procedure
#' is applied.
#'
#' With multivariate data, when using \code{kernel = "gaussian"} and \code{bw} is a covariance
#' matrix, multivariate Gaussian kernel is applied. With multivariate data, when \code{bw} is
#' a vector, or \code{kernel} is other then \code{"gaussian"}, product kernel is used.
#'
#'
#' \strong{Univariate kernels}
#'
#' Univariate kernel density estimator is defined as
#'
#' \deqn{
#' \hat{f_h}(x) = \sum_{i=1}^n w_i \, K_h\left(\frac{x-y_i}{h}\right)
#' }{
#' f(x) = sum[i](w[i] * Kh((x-y[i])/h))
#' }
#'
#' where \eqn{w} is a vector of weights such that \eqn{\sum_i w_i = 1}{sum(w) = 1}
#' (by default \eqn{w_i=1/n}{w[i]=1/n} for all \eqn{i}), \eqn{K_h = K(x/h)/h}{Kh = K(x/h)/h} is
#' kernel \eqn{K} parametrized by bandwidth \eqn{h} and \eqn{y} is a vector of
#' data points used for estimating the kernel density.
#'
#' To draw samples from univariate kernel density, the following procedure can be applied (Silverman, 1986):
#'
#' \emph{Step 1} Sample \eqn{i} uniformly with replacement from \eqn{1,\dots,n}.
#'
#' \emph{Step 2} Generate \eqn{\varepsilon}{\epsilon} to have probability density \eqn{K}.
#'
#' \emph{Step 3} Set \eqn{x = y_i + h\varepsilon}{x = y[i] + h\epsilon}.
#'
#' If samples are required to have the same variance as \code{data}
#' (i.e. \code{shrinked = TRUE}), then \emph{Step 3} is modified
#' as following:
#'
#' \emph{Step 3'} \eqn{
#' x = \bar y + (y_i - \bar y + h\varepsilon)/(1 + h^2 \sigma^2_K/\sigma^2_Y)^{1/2}
#' }{
#' x = m + (y[i] - m + h\epsilon)/sqrt(1 + h^2 var(K)/var(y))
#' }
#'
#' where \eqn{\sigma_K^2}{sK} is variance of the kernel (fixed to 1 for kernels used in this package).
#'
#' When shrinkage described in \emph{Step 3'} is applied, the smoothed bootstrap density function changes it's form to
#'
#' \deqn{
#' \hat{f}_{h,b}(x) = (1 + r) \hat{f_h}(x + r(x - \bar{y}))
#' }{
#' fb(x) = (1+r) f(x + r (x-mean(y)))
#' }
#'
#' where \eqn{r = \left(1 + h^2 \sigma_K^2 / \sigma_y^2 \right)^{1/2}-1}{r = sqrt(1 + h^2 sK/var(y)) - 1}.
#'
#' This package offers the following univariate kernels:
#'
#' \tabular{ll}{
#' \emph{Gaussian}     \tab \eqn{\frac{1}{\sqrt{2\pi}} e^{-{u^2}/2}}{1/sqrt(2\pi) exp(-(u^2)/2)} \cr
#' \emph{Rectangular}  \tab \eqn{\frac{1}{2} \ \mathbf{1}_{(|u|\leq1)}}{1/2} \cr
#' \emph{Triangular}   \tab \eqn{(1-|u|) \ \mathbf{1}_{(|u|\leq1)}}{1 - |u|} \cr
#' \emph{Epanchenikov} \tab \eqn{\frac{3}{4}(1-u^2) \ \mathbf{1}_{(|u|\leq1)}}{3/4 (1 - u^2)} \cr
#' \emph{Biweight}     \tab \eqn{\frac{15}{16}(1-u^2)^2 \ \mathbf{1}_{(|u|\leq1)}}{15/16 (1 - u^2)^2} \cr
#' \emph{Cosine}       \tab \eqn{\frac{1}{2} \left(1 + \cos(\pi u)\right) \ \mathbf{1}_{(|u|\leq1)}}{1/2 (1 + cos(\pi u))} \cr
#' \emph{Optcosine}    \tab \eqn{\frac{\pi}{4}\cos\left(\frac{\pi}{2}u\right) \ \mathbf{1}_{(|u|\leq1)}}{\pi/4 cos(\pi/2 u)}
#' }
#'
#' All the kernels are re-scalled so that their standard deviations are equal to 1,
#' so that bandwidth parameter controls their standard deviations.
#'
#' Random generation from Epachenikov kernel is done using algorithm
#' described by Devoye (1986). For optcosine kernel inverse transform
#' sampling is used. For biweight kernel random values are drawn from
#' \eqn{\mathrm{Beta}(3, 3)}{Beta(3, 3)} distribution and
#' \eqn{\mathrm{Beta}(3.3575, 3.3575)}{Beta(3.3575, 3.3575)}
#' distribution serves as a close approximation of cosine kernel.
#' Random generation for triangular kernel is done by taking difference
#' of two i.i.d. uniform random variates. To sample from rectangular
#' and Gaussian kernels standard random generation algorithms are used
#' (see \code{\link[stats]{runif}} and \code{\link[stats]{rnorm}}).
#'
#'
#' \strong{Product kernels}
#'
#' Univariate kernels may easily be extended to multiple dimensions by
#' using product kernel
#'
#' \deqn{
#' \hat{f_H}(x_1,\dots,x_n) = \sum_{i=1}^n w_i \prod_{j=1}^m
#' K_{h_j} \left( \frac{x_i - y_{ij}}{h_j} \right)
#' }{
#' f(x) = sum[i](w[i] * prod[j]( Kh[j]((x[i]-y[i,j])/h[j]) ))
#' }
#'
#' where \eqn{w} is a vector of weights such that \eqn{\sum_i w_i = 1}{sum(w) = 1}
#' (by default \eqn{w_i=1/n}{w[i]=1/n} for all \eqn{i}), and \eqn{K_{h_j}}{Kh[j]}
#' are univariate kernels \eqn{K} parametrized by bandwidth \eqn{h_j}{h[j]}, where
#' \eqn{\boldsymbol{y}}{y} is a matrix of data points used for estimating the
#' kernel density.
#'
#' Random generation from product kernel is done by drawing with replacement
#' rows of \eqn{y}, and then adding random noise from univariate kernel \eqn{K},
#' parametrized by corresponding bandwidth parameter \eqn{h}, to the sampled values.
#'
#'
#' \strong{Multivariate kernels}
#'
#' Multivariate kernel density estimator may also be defined in terms of multivariate kernels
#' (e.g. multivariate normal distribution, as in this package)
#'
#' \deqn{
#' \hat{f_H}(x_1,\dots,x_n) = \sum_{i=1}^n w_i \, K_H \left( \mathbf{x}-\boldsymbol{y}_i \right)
#' }{
#' f(x) = sum[i](w[i] * KH(x-y[i]))
#' }
#'
#' where \eqn{w} is a vector of weights such that \eqn{\sum_i w_i = 1}{sum(w) = 1}
#' (by default \eqn{w_i=1/n}{w[i]=1/n} for all \eqn{i}), \eqn{K_H}{KH} is
#' kernel \eqn{K} parametrized by bandwidth matrix \eqn{H} and \eqn{\boldsymbol{y}}{y}
#' is a matrix of data points used for estimating the kernel density.
#'
#' \emph{Notice:} When using multivariate normal (Gaussian) distribution as a kernel \eqn{K}, the
#' bandwidth parameter \eqn{H} is a \emph{covariance matrix} as compared to standard deviations
#' used in univariate and product kernels.
#'
#' Random generation from multivariate kernel is done by drawing with replacement
#' rows of \eqn{y}, and then adding random noise from multivariate kernel \eqn{K},
#' parametrized by corresponding bandwidth matrix \eqn{H}, to the sampled values.
#'
#'
#' @references
#' Silverman, B. W. (1986). Density estimation for statistics and data analysis.
#' Chapman and Hall/CRC.
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
#' Hall, P., DiCiccio, T.J. and Romano, J.P. (1989). On smoothing and the bootstrap.
#' The Annals of Statistics, 692-704.
#'
#' @references
#' Silverman, B.W. and Young, G.A. (1987). The bootstrap: To smooth or not to smooth?
#' Biometrika, 469-479.
#'
#' @references
#' Scott, D.W. (1992). Multivariate density estimation: theory, practice,
#' and visualization. John Wiley & Sons.
#'
#' @references
#' Wang, S. (1995). Optimizing the smoothed bootstrap. Annals of the Institute of
#' Statistical Mathematics, 47(1), 65-80.
#'
#' @references
#' Young, G.A. (1990). Alternative smoothed bootstraps. Journal of the Royal
#' Statistical Society. Series B (Methodological), 477-484.
#'
#' @references
#' De Angelis, D. and Young, G.A. (1992). Smoothing the bootstrap.
#' International Statistical Review/Revue Internationale de Statistique, 45-56.
#'
#' @references
#' Polansky, A.M. and Schucany, W. (1997). Kernel smoothing to improve bootstrap
#' confidence intervals. Journal of the Royal Statistical Society: Series B
#' (Statistical Methodology), 59(4), 821-838.
#'
#' @references
#' Devroye, L. (1986). Non-uniform random variate generation. New York: Springer-Verlag.
#'
#' @references
#' Parzen, E. (1962). On estimation of a probability density function and mode.
#' The annals of mathematical statistics, 33(3), 1065-1076.
#'
#' @references
#' Silverman, B.W. and Young, G.A. (1987). The bootstrap: To smooth or not to smooth?
#' Biometrika, 469-479.
#'
#' @references
#' Jones, M.C. (1991). On correcting for variance inflation in kernel density estimation.
#' Computational Statistics & Data Analysis, 11, 3-15.
#'
#'
#' @seealso \code{\link{bw.silv}}, \code{\link[stats]{density}},
#'          \code{\link[stats]{bandwidth}}, \code{\link{kernelboot-class}}
#'
#'
#' @examples
#'
#' set.seed(1)
#'
#' # smooth bootstrap of parameters of linear regression
#'
#' b1 <- kernelboot(mtcars, function(data) coef(lm(mpg ~ drat + wt, data = data)) , R = 250)
#' b1
#' summary(b1)
#'
#' b2 <- kernelboot(mtcars, function(data) coef(lm(mpg ~ drat + wt, data = data)) , R = 250,
#'                  kernel = "epanechnikov")
#' b2
#' summary(b2)
#'
#' # smooth bootstrap of parameters of linear regression
#' # smoothing phase is not applied to "am" and "cyl" variables
#'
#' b3 <- kernelboot(mtcars, function(data) coef(lm(mpg ~ drat + wt + am + cyl, data = data)) , R = 250,
#'                  ignore = c("am", "cyl"))
#' b3
#' summary(b3)
#'
#' # standard bootstrap (without kernel smoothing)
#'
#' b4 <- kernelboot(mtcars, function(data) coef(lm(mpg ~ drat + wt + am + cyl, data = data)) , R = 250,
#'                  ignore = colnames(mtcars))
#' b4
#' summary(b4)
#'
#' # smooth bootstrap for median of univariate data
#'
#' b5 <- kernelboot(mtcars$mpg, function(data) median(data) , R = 250)
#' b5
#' summary(b5)
#'
#' # draw samples from different kernels
#'
#' \dontrun{
#'
#' kernels <- c("gaussian", "epanechnikov", "rectangular", "triangular",
#'              "biweight", "cosine", "optcosine")
#'
#' data <- mtcars[, c(1, 3)]
#' bw <- bw.silv(data)
#' R <- 250
#'
#' partmp <- par(mfrow = c(2, 4), mar = c(3, 3, 3, 3))
#' for (k in kernels) {
#'   plot(kernelboot(data, identity, R = R, kernel = k, bw = sqrt(diag(bw)))$boot.samples,
#'        xlim = c(-10, 50), ylim = c(-100, 600), col = "#ADD8E640")
#'   points(data, pch = 2, lwd = 2, col = "red")
#'   title(k)
#' }
#' plot(kernelboot(data, identity, R = R, kernel = "g", bw = bw)$boot.samples,
#'      xlim = c(-10, 50), ylim = c(-100, 600), col = "#ADD8E640")
#' points(data, pch = 2, lwd = 2, col = "red")
#' title("multivariate gaussian")
#' par(partmp)
#'
#' }
#'
#'
#' @importFrom stats rnorm bw.SJ bw.bcv bw.nrd bw.nrd0 bw.ucv
#' @importFrom future plan multiprocess future_lapply availableCores
#'
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG registerDoRNG %dorng%
#'
#' @export

kernelboot <- function(data, statistic, R = 500L, bw = "default",
                       kernel = c("gaussian", "epanechnikov", "rectangular",
                                  "triangular", "biweight", "cosine", "optcosine"),
                       weights = NULL, adjust = 1,
                       shrinked = TRUE, ignore = NULL,
                       parallel = FALSE, workers = getOption("mc.cores", 2L)) {

  call <- match.call()
  kernel <- match.arg(kernel)
  n <- NROW(data)
  m <- NCOL(data)
  vars <- NULL
  kd_type <- NULL

  if (!(is.vector(data) || is.data.frame(data) || is.matrix(data)))
    stop("data is not a vector, data.frame, or matrix")

  if (is.character(bw)) {
    bw <- tolower(bw)
    if (bw == "default") {
      if (is.vector(data)) {
        bw <- bw.nrd0(data)
      } else {
        bw <- bw.silv(data)
        if (kernel != "gaussian")
          bw <- sqrt(diag(bw))
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

  if (!is.vector(adjust))
    stop("adjust is not a scalar")

  bw <- bw * adjust[1L]

  # check for non-numeric, NAs, NaNs, infinite values
  if (!all(is.finite(bw)))
    stop("inappropriate values of bw")

  if (!is.null(weights)) {
    if (!all(is.finite(weights)))
      stop("inappropriate values of weights")
  }

  # equally weighted
  if (length(weights) == 1L)
    weights <- NULL

  # try evaluating statistic() on the original data

  tryCatch(
    orig.stat <- statistic(data),
    error = function(e) {
      message("applying the statistic on the original data resulted in an error")
      stop(e)
    }
  )

  # looping functions - paralell or lapply

  if (parallel && workers > 1L) {

    # plan(multiprocess, workers = workers)
    # repeatFun <- function(n, FUN, workers) future_lapply(1:n, FUN, future.seed = TRUE)

    cl <- makeCluster(workers)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
    repeatFun <- function(n, FUN, workers) foreach(i = 1:n) %dorng% FUN(i)

  } else {

    parallel <- FALSE
    repeatFun <- function(n, FUN, workers) lapply(1:n, FUN)

  }

  if (is.data.frame(data) || is.matrix(data)) {

    # data is data.frame or matrix

    num_cols <- is_numeric(data)             # find numeric columns
    ignr_cols <- colnames(data) %in% ignore
    incl_cols <- num_cols & !ignr_cols

    if (!is.null(colnames(data))) {
      vars <- list(
        smoothed = colnames(data)[incl_cols],
        ignored  = colnames(data)[!incl_cols]
      )
    } else {
      vars <- list(
        smoothed = which(incl_cols),
        ignored  = which(!incl_cols)
      )
    }

    if (!any(incl_cols)) {

      # standard bootstrap

      kd_type <- "none"

      res <- repeatFun(R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        boot.data <- data[idx, ]
        statistic(boot.data)

      }, workers = workers)

    } else {

      # smoothed bootstrap

      data_mtx <- as.matrix(data[, incl_cols])

      if (is.null(weights))
        weights <- rep(1/n, n)

      if (kernel != "gaussian" || is.vector(bw)) {

        # product kernel

        kd_type <- "product"

        if (is.vector(bw)) {
          if (length(bw) == 1L)
            bw <- rep(bw, m)
        } else {
          if (!is.square(bw))
            stop("bw is not a square matrix")
          bw <- diag(bw)
        }

        bw <- bw[incl_cols]

        res <- repeatFun(R, function(i) {

          samp <- cpp_rmvk(n, data_mtx, bw, weights, kernel)
          idx <- samp$boot_index
          boot.data <- data[idx, ]
          boot.data[, incl_cols] <- samp$samples
          statistic(boot.data)

        }, workers = workers)

      } else {

        # MVN kernel

        kd_type <- "multivariate"

        # is this check really needed?
        # if (qr(data_mtx)$rank < min(dim(data_mtx)))
        #   warning("data matrix is rank deficient")

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

        if (ncol(bw) != m)
          stop("dimensions of data and bw do not match")

        if (!is.square(bw))
          stop("bw is not a square matrix")

        bw <- as.matrix(bw)
        bw <- bw[incl_cols, incl_cols]
        bw_chol <- chol(bw)
        mm <- sum(incl_cols)

        res <- repeatFun(R, function(i) {

          idx <- sample.int(n, n, replace = TRUE, prob = weights)
          boot.data <- data[idx, ]
          samp <- matrix(rnorm(n*mm), n, mm) %*% bw_chol
          boot.data[, incl_cols] <- boot.data[, incl_cols] + samp
          statistic(boot.data)

        }, workers = workers)

      }
    }

  } else if (is.vector(data)) {

    # data is a vector

    if (!is.numeric(data)) {

      # standard bootstrap

      kd_type <- "none"

      res <- repeatFun(R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        boot.data <- data[idx]
        statistic(boot.data)

      }, workers = workers)

    } else {

      # smoothed bootstrap

      kd_type <- "univariate"

      if (!is.vector(bw))
        stop("bw is not a scalar")
      if (length(bw) > 1L) {
        bw <- bw[1L]
        message("bw has length > 1 and only the first element will be used")
      }

      if (is.null(weights))
        weights <- rep(1/n, n)

      res <- repeatFun(R, function(i) {

        samp <- cpp_ruvk(n, data, bw, weights, kernel, shrinked)
        boot.data <- drop(samp$samples)
        statistic(boot.data)

      }, workers = workers)

    }

  } else {

    stop("unsupported data type")

  }

  # simplify the results to data.frame
  samples <- do.call(rbind, res)

  structure(list(
    orig.stat     = orig.stat,
    boot.samples  = samples,
    call          = call,
    statistic     = statistic,
    orig.data     = data,
    variables     = vars,
    type          = kd_type,
    param = list(
      R         = R,
      bw        = bw,
      adjust    = adjust,
      weights   = weights,
      kernel    = kernel,
      shrinked  = shrinked,
      parallel  = parallel
    )
  ), class = "kernelboot")

}

