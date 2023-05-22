#' Butcher methods for tidyclust
#'
#' These methods allow you to use the butcher package to reduce the size of
#' a tidyclust cluster object. After calling `butcher::butcher()` on an object,
#' the only guarantee is that you will still be able to `predict()` from that
#' object. Other functions may not work as expected.
#'
#' @param x A tidyclust cluster object.
#' @param verbose Should information be printed about how much memory is freed
#'   from butchering?
#' @param ... Extra arguments possibly used by underlying methods.
#'
#' @name tidyclust-butcher

# @export - onLoad
#' @rdname tidyclust-butcher
axe_env.k_means <- function(x, verbose = FALSE, ...) {
  x$args <- x$args |> lapply(axe_env, verbose = verbose, ...)
  x$method$fit$args <- x$method$fit$args |> lapply(axe_env, verbose = verbose, ...)

  return(x)
}

# @export - onLoad
#' @rdname tidyclust-butcher
axe_env.gmm <- function(x, verbose = FALSE, ...) {
  x$args <- x$args |> lapply(axe_env, verbose = verbose, ...)
  x$method$fit$args <- x$method$fit$args |> lapply(axe_env, verbose = verbose, ...)

  return(x)
}

# @export - onLoad
#' @rdname tidyclust-butcher
axe_env.hier_clust <- function(x, verbose = FALSE, ...) {
  x$args <- x$args |> lapply(axe_env, verbose = verbose, ...)
  x$method$fit$args <- x$method$fit$args |> lapply(axe_env,
                                                   verbose = verbose, ...)

  return(x)
}

# @export - onLoad
#' @rdname tidyclust-butcher
axe_env.cluster_fit <- function(x, verbose = FALSE, ...) {
  x$spec <- axe_env(x$spec, verbose = verbose, ...)
  x$fit <- axe_env(x$fit, verbose = verbose, ...)

  return(x)
}

# @export - onLoad
#' @rdname tidyclust-butcher
axe_env.workflow <- function(x, verbose = FALSE, ...) {
  x <- workflows:::axe_env.workflow(x)
  x$fit$actions$model$spec <- extract_fit_parsnip(x)$spec

  return(x)
}
