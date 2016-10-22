

#' Draw samples from kernel density
#'
#' @param n            number of observations to be drawn. If \code{length(n) > 1}, the length is taken to be the number required.
#'                     If \code{n} is not provided, it defaults to \code{length(x)}.
#' @param x            the data from which the estimate is to be computed. \emph{Warning:} infinite values are omiitted.
#'                     If not provided, draws samples from given \code{kernel} and uses \code{bw} only if
#'                     it was provided as numeric value.
#' @param bw           the smoothing bandwidth to be used. The kernels are scaled such that this is
#'                     the standard deviation of the smoothing kernel. (Note this differs from the
#'                     reference books cited below, and from S-PLUS.)
#'                     \code{bw} can also be a character string giving a rule to choose the bandwidth. See
#'                     \code{\link[stats]{bw.nrd}}. The default, \code{"nrd0"}, has remained the default for historical
#'                     and compatibility reasons, rather than as a general recommendation, where e.g., \code{"SJ"}
#'                     would rather fit, see also Venables and Ripley (2002).
#' @param kernel       a character string giving the smoothing kernel to be used.
#' @param preserve.var logical, if \code{TRUE}, then the bootstrap samples preserve sample variance.
#' @param adjust       the bandwidth used is actually \code{adjust*bw}. This makes it easy to specify values like
#'                     ‘half the default’ bandwidth.
#' @param weights      numeric vector of non-negative observation weights, hence of same length as \code{x}.
#'                     The default NULL is equivalent to \code{weights = rep(1/nx, nx)} where \code{nx}
#'                     is the length of (the finite entries of) \code{x[]}.
#' @param na.rm        logical; if \code{TRUE}, missing values are removed from \code{x}.
#'                     If \code{FALSE} any missing values cause an error.
#'
#' @seealso \code{\link[stats]{density}}
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
#' @references
#' Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S.
#' New York: Springer.
#'
#' @examples
#'
#' hist(rkde(1e5, mtcars$disp), 100, freq = FALSE)
#' lines(density(mtcars$disp), col = "red")
#' rug(mtcars$disp, lwd = 2)
#'
#' @export

rkernel <- function(n, x, bw = "nrd0",
                    kernel = c("epanechnikov", "gaussian", "rectangular",
                               "triangular", "biweight", "triweight",
                               "cosine", "optcosine"), preserve.var = TRUE,
                    adjust = 1, weights = NULL, na.rm = FALSE) {

  kernel <- match.arg(kernel)

  rng_kern <- switch(kernel,
                     epanechnikov = rempan,
                     rectangular = rrect,
                     triangular = rtriang,
                     biweight = rbiweight,
                     triweight = rtriweight,
                     cosine = rcosine,
                     optcosine = roptcos,
                     rnorm)

  if (missing(x)) {
    if (!is.numeric(bw)) {
      warning(paste0("missing x, not using '", bw, "' bandwidth"))
      bw <- 1
    }
    return(rng_kern(n) * bw*adjust)
  }

  x <- as.vector(x)

  if (length(x) == 1 || isTRUE(all.equal(var(x), 0))) {
    preserve.var <- FALSE
    warning("var(x) is 0: preserve.var was changed to FALSE")
  }

  x.na <- is.na(x)
  if (any(x.na)) {
    if (na.rm)
      x <- x[!x.na]
    else stop("'x' contains missing values")
  }

  nx <- N <- length(x)

  x.finite <- is.finite(x)
  if (any(!x.finite)) {
    x <- x[x.finite]
    nx <- length(x)
  }

  if (is.character(bw)) {
    if (nx < 2)
      stop("need at least 2 points to select a bandwidth automatically")
    bw <- switch(tolower(bw),
                 nrd0 = bw.nrd0(x),
                 nrd = bw.nrd(x),
                 ucv = bw.ucv(x),
                 bcv = bw.bcv(x), sj = ,
                 `sj-ste` = bw.SJ(x, method = "ste"),
                 `sj-dpi` = bw.SJ(x, method = "dpi"),
                 stop("unknown bandwidth rule"))
  }
  if (length(bw) > 1) {
    bw <- bw[1]
    warning("'bw' has length > 1 and only the first element will be used")
  }
  if (!is.finite(bw))
    stop("non-finite 'bw'")
  if (bw <= 0)
    stop("'bw' is not positive.")

  if (missing(n)) {
    n <- length(x)
  } else if (length(n) > 1) {
    n <- length(n)
  }

  if (!is.null(weights)) {
    if (length(weights) != N)
      stop("'x' and 'weights' have unequal length")
    if (!all(is.finite(weights)))
      stop("'weights' must all be finite")
    if (any(weights < 0))
      stop("'weights' must not be negative")
    if (any(!x.finite))
      weights <- weights[x.finite]
  }

  idx <- sample.int(nx, n, replace = TRUE, prob = weights)

  bw <- adjust * bw
  eps <- rng_kern(n) * bw

  if (preserve.var) {
    mx <- mean(x)
    sx <- var(x)
    mx + (x[idx]-mx+eps)/sqrt(1 + bw^2/sx)
  } else {
    x[idx]+eps
  }

}

