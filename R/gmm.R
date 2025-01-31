#' Gaussian Mixture Model (GMM) Clustering
#'
#' @description
#'
#' `gmm()` defines a model that fits clusters based on distances to a number
#' of centers.
#'
#' @param mode A single character string for the type of model. The only
#'   possible value for this model is "partition".
#' @param engine A single character string specifying what computational engine
#'   to use for fitting. Possible engines are listed below. The default for this
#'   model is `"ClusterR"`.
#' @param num_clusters Positive integer, number of clusters in model.
#'
#' @return A `gmm` cluster specification.
#'
#' @examples
#' # Show all engines
#' modelenv::get_from_env("gmm")
#'
#' gmm()
#' @export
gmm <-
  function(mode = "partition",
           engine = "ClusterR",
           num_clusters = NULL) {
    args <- list(
      num_clusters = enquo(num_clusters)
    )

    new_cluster_spec(
      "gmm",
      args = args,
      eng_args = NULL,
      mode = mode,
      method = NULL,
      engine = engine
    )
  }

#' @export
print.gmm <- function(x, ...) {
  cat("Gaussian Mixture Model (GMM) Clustering Specification (", x$mode, ")\n\n", sep = "")
  model_printer(x, ...)

  if (!is.null(x$method$fit$args)) {
    cat("Model fit template:\n")
    print(show_call(x))
  }

  invisible(x)
}

#' @export
translate_tidyclust.gmm <- function(x, engine = x$engine, ...) {
  x <- translate_tidyclust.default(x, engine, ...)
  x
}

# ------------------------------------------------------------------------------

#' @method update gmm
#' @rdname tidyclust_update
#' @export
update.gmm <- function(object,
                           parameters = NULL,
                           num_clusters = NULL,
                           fresh = FALSE, ...) {
  eng_args <- parsnip::update_engine_parameters(
    object$eng_args, fresh = fresh, ...
  )

  if (!is.null(parameters)) {
    parameters <- parsnip::check_final_param(parameters)
  }
  args <- list(
    num_clusters = enquo(num_clusters)
  )

  args <- parsnip::update_main_parameters(args, parameters)

  if (fresh) {
    object$args <- args
    object$eng_args <- eng_args
  } else {
    null_args <- map_lgl(args, null_value)
    if (any(null_args)) {
      args <- args[!null_args]
    }
    if (length(args) > 0) {
      object$args[names(args)] <- args
    }
    if (length(eng_args) > 0) {
      object$eng_args[names(eng_args)] <- eng_args
    }
  }

  new_cluster_spec(
    "gmm",
    args = object$args,
    eng_args = object$eng_args,
    mode = object$mode,
    method = NULL,
    engine = object$engine
  )
}

# ------------------------------------------------------------------------------

#' @export
check_args.gmm <- function(object) {
  args <- lapply(object$args, rlang::eval_tidy)

  if (all(is.numeric(args$num_clusters)) && any(args$num_clusters < 0)) {
    rlang::abort("The number of centers should be >= 0.")
  }

  invisible(object)
}

# ------------------------------------------------------------------------------

#' Simple Wrapper around ClusterR Gaussian Mixture Model (GMM) Clustering
#'
#' This wrapper runs `ClusterR::GMM`, creates `clusters` and `probabilities` fields,
#' adds column names to the `centroids` and `covariance_matrices` fields, and
#' add summary stats such as degrees of freedom, log-likelihood, AIC, BIC, and SSEs.
#'
#' @param data matrix or data frame
#' @param clusters the number of gaussian mixture components
#' @param dist_mode the distance used during the seeding of initial means and k-means clustering. One of, \emph{eucl_dist}, \emph{maha_dist}.
#' @param seed_mode how the initial means are seeded prior to running k-means and/or EM algorithms. One of, \emph{static_subset}, \emph{random_subset}, \emph{static_spread}, \emph{random_spread}.
#' @param km_iter the number of iterations of the k-means algorithm
#' @param em_iter the number of iterations of the EM algorithm
#' @param verbose either TRUE or FALSE; enable or disable printing of progress during the k-means and EM algorithms
#' @param var_floor the variance floor (smallest allowed value) for the diagonal covariances
#' @param seed integer value for random number generator (RNG)
#' @param full_covariance_matrices a boolean. If FALSE "diagonal" covariance matrices (i.e. in each covariance matrix, all entries outside the main diagonal are assumed to be zero) otherwise "full" covariance matrices will be returned. Be aware in case of "full" covariance matrices a cube (3-dimensional) rather than a matrix for the output "covariance_matrices" value will be returned.
#' @return a list consisting of the centroids, covariance matrix ( where each row of the matrix represents a diagonal covariance matrix), weights and the log-likelihoods for each gaussian component. In case of Error it returns the error message and the possible causes.
#' @details
#' This function is an R implementation of the 'gmm_diag' class of the Armadillo library. The only exception is that user defined parameter settings are not supported, such as seed_mode = 'keep_existing'.
#' For probabilistic applications, better model parameters are typically learned with dist_mode set to maha_dist.
#' For vector quantisation applications, model parameters should be learned with dist_mode set to eucl_dist, and the number of EM iterations set to zero.
#' In general, a sufficient number of k-means and EM iterations is typically about 10.
#' The number of training samples should be much larger than the number of Gaussians.
#' Seeding the initial means with static_spread and random_spread can be much more time consuming than with static_subset and random_subset.
#' The k-means and EM algorithms will run faster on multi-core machines when OpenMP is enabled in your compiler (eg. -fopenmp in GCC)
#' @references
#' http://arma.sourceforge.net/docs.html
#' @export
ClusterR_gmm_fit <- function(data, clusters,
                             dist_mode = 'eucl_dist',
                             seed_mode = 'random_subset',
                             km_iter = 10,
                             em_iter = 5,
                             verbose = FALSE,
                             var_floor = 1e-10,
                             seed = 1,
                             full_covariance_matrices = FALSE) {
  res <- ClusterR::GMM(data, clusters,
                       dist_mode = dist_mode,
                       seed_mode = seed_mode,
                       km_iter = km_iter,
                       em_iter = em_iter,
                       verbose = verbose,
                       var_floor = var_floor,
                       seed = seed,
                       full_covariance_matrices = full_covariance_matrices)

  # Dimensions
  m <- ncol(res$centroids)
  k <- nrow(res$centroids)
  n <- nrow(data)

  # Degrees of freedom
  if (length(dim(res$covariance_matrices)) == 3) {
    # Full covariance matrices
    n_params <- k * m * (m + 1) / 2
  } else if (length(dim(res$covariance_matrices)) == 2) {
    # Diagonal covariance matrices
    n_params <- k * m
  } else {
    stop("Unknown size of the covariance matrices!")
  }

  names <- paste0("Cluster_", seq_len(k))

  preds <- .gmm_predict_ClusterR(res, data, fuzzy = TRUE, prefix = "Cluster_")

  res$clusters <- as.integer(preds$.pred_cluster)
  res$probabilities <- preds$.pred_prob

  if (!is.null(res$Log_likelihood)) {
    res$log_likelihood <- res$Log_likelihood
  } else if (!is.null(preds$log_likelihood)) {
    res$log_likelihood <- preds$log_likelihood
  }
  res$Log_likelihood <- NULL
  colnames(res$log_likelihood) <- names
  res$log_likelihood <- tibble::as_tibble(res$log_likelihood)

  rownames(res$centroids) <- names
  colnames(res$centroids) <- colnames(data)
  res$centroids <- tibble::as_tibble(res$centroids)

  names(res$weights) <- names

  if (is.matrix(res$covariance_matrices)) {
    rownames(res$covariance_matrices) <- names
    colnames(res$covariance_matrices) <- colnames(data)
    res$covariance_matrices <- tibble::as_tibble(res$covariance_matrices)
  }

  # Log-Likelihood
  loglik <- 0
  for (i in unique(res$clusters)) {
    loglik <- loglik + sum(res$log_likelihood[which(res$clusters == i), i])
  }

  res$obs_per_cluster <- table(res$clusters)
  names(res$obs_per_cluster) <- names

  res$n <- n
  res$m <- m
  res$k <- k
  res$df <- n_params

  res$aic <- -2 * loglik + 2.0 * n_params
  res$bic <- -2 * loglik + log(n) * n_params
  res$loglik <- loglik

  by_clust <- data %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      .cluster = factor(names[res$clusters], levels = names)
    ) %>%
    dplyr::group_by(.cluster) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.cluster)

  res$se <- matrix(NA, k, k)
  res$me <- matrix(NA, k, k)
  res$sse <- matrix(NA, k, k)
  res$mse <- matrix(NA, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      distances <- Rfast::dista(res$centroids[i,], by_clust$data[[j]])
      res$se[i,j] <- sum(distances, na.rm = TRUE)
      res$me[i,j] <- mean(distances, na.rm = TRUE)
      res$sse[i,j] <- sum(distances^2, na.rm = TRUE)
      res$mse[i,j] <- mean(distances^2, na.rm = TRUE)
    }
  }
  colnames(res$se) <- names
  rownames(res$se) <- names
  colnames(res$me) <- names
  rownames(res$me) <- names

  colnames(res$sse) <- names
  rownames(res$sse) <- names
  colnames(res$mse) <- names
  rownames(res$mse) <- names

  res$se_within <- diag(res$se)
  res$me_within <- diag(res$me)

  res$sse_within <- diag(res$sse)
  res$mse_within <- diag(res$mse)

  res$se_between <- res$se
  diag(res$se_between) <- NA

  res$me_between <- res$me
  diag(res$me_between) <- NA

  res$sse_between <- res$sse
  diag(res$sse_between) <- NA

  res$mse_between <- res$mse
  diag(res$mse_between) <- NA

  distances_total <- Rfast::dista(t(colMeans(data)), data)
  res$se_total <- sum(distances_total, na.rm = TRUE)
  res$me_total <- mean(distances_total, na.rm = TRUE)
  res$sse_total <- sum(distances_total^2, na.rm = TRUE)
  res$mse_total <- mean(distances_total^2, na.rm = TRUE)

  res$se_ratio <- sum(res$se_within, na.rm = TRUE) / res$se_total
  res$me_ratio <- mean(res$me_within, na.rm = TRUE) / res$me_total

  res$sse_ratio <- sum(res$sse_within, na.rm = TRUE) / res$sse_total
  res$mse_ratio <- mean(res$mse_within, na.rm = TRUE) / res$mse_total

  structure(res[order(names(res))],
            class = c("GMMCluster", 'Gaussian Mixture Models'))
}
