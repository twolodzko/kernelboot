
#' @importFrom stats quantile
#' @export

print.kernelboot <- function(x, ...) {

  if (!is.null(x$call))
    cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n")

  cat("\nSummary Statistics:\n")

  prc <- c(0.025, 0.5, 0.975)

  stat <- cbind(
    colMeans(x$boot.sample),
    apply(x$boot.sample, 2, sd),
    t(apply(x$boot.sample, 2, quantile, prc))
  )
  colnames(stat) <- c("mean", "sd", paste0(prc*100, "%"))

  print(stat)

}
