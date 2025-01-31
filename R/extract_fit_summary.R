#' S3 method to get fitted model summary info depending on engine
#'
#' @param object a fitted [`cluster_spec`] object
#' @param ... other arguments passed to methods
#'
#' @return A list with various summary elements
#'
#' @details
#'
#' The elements `cluster_names` and `cluster_assignments` will be factors.
#'
#' @examples
#' kmeans_spec <- k_means(num_clusters = 5) %>%
#'   set_engine("stats")
#'
#' kmeans_fit <- fit(kmeans_spec, ~., mtcars)
#'
#' kmeans_fit %>%
#'   extract_fit_summary()
#' @export
extract_fit_summary <- function(object, ...) {
  UseMethod("extract_fit_summary")
}

#' @export
extract_fit_summary.cluster_spec <- function(
  object,
  ...,
  call = rlang::caller_env(n = 0)
) {
  cli::cli_abort(
    c(
      "This function requires a fitted model.",
      "i" = "Please use {.fn fit} on your cluster specification."
    ),
    call = call
  )
}

#' @export
extract_fit_summary.cluster_fit <- function(object, ...) {
  extract_fit_summary(object$fit, ...)
}

#' @export
extract_fit_summary.workflow <- function(object, ...) {
  extract_fit_summary(object$fit$fit$fit, ...)
}

#' @export
extract_fit_summary.kmeans <- function(object, ..., prefix = "Cluster_") {
  names <- paste0(prefix, seq_along(object$size))
  names <- factor(names)

  cluster_asignments <- factor(
    names[object$cluster],
    levels = levels(names)
  )

  centroids <- tibble::as_tibble(object$centers)

  list(
    cluster_names = names,
    centroids = centroids,
    n_members = object$size,
    sse_within_total_total = object$withinss,
    sse_total = object$totss,
    orig_labels = unname(object$cluster),
    cluster_assignments = cluster_asignments
  )
}

#' @export
extract_fit_summary.KMeansCluster <- function(
  object,
  ...,
  prefix = "Cluster_"
) {
  names <- paste0(prefix, seq_len(nrow(object$centroids)))
  names <- factor(names)

  cluster_asignments <- factor(
    names[object$clusters],
    levels = levels(names)
  )

  centroids <- tibble::as_tibble(object$centroids)

  list(
    cluster_names = names,
    centroids = centroids,
    n_members = as.integer(object$obs_per_cluster),
    sse_within_total_total = as.numeric(object$WCSS_per_cluster),
    sse_total = object$total_SSE,
    orig_labels = object$clusters,
    cluster_assignments = cluster_asignments
  )
}

#' @export
extract_fit_summary.kproto <- function(object, ..., prefix = "Cluster_") {
  names <- paste0(prefix, seq_len(nrow(object$centers)))
  names <- factor(names)

  cluster_asignments <- factor(
    names[object$cluster],
    levels = levels(names)
  )

  centroids <- tibble::as_tibble(object$centers)

  list(
    cluster_names = names,
    centroids = centroids,
    n_members = as.integer(object$size),
    sse_within_total_total = object$withinss,
    sse_total = object$tot.withinss,
    orig_labels = seq_len(length(table(object$cluster))),
    cluster_assignments = cluster_asignments
  )
}

#' @export
extract_fit_summary.kmodes <- function(object, ..., prefix = "Cluster_") {
  names <- paste0(prefix, seq_len(nrow(object$modes)))
  names <- factor(names)

  cluster_asignments <- factor(
    names[object$cluster],
    levels = levels(names)
  )

  centroids <- tibble::as_tibble(object$modes)

  list(
    cluster_names = names,
    centroids = centroids,
    n_members = as.integer(object$size),
    sse_within_total_total = object$withinss,
    sse_total = object$tot.withinss,
    orig_labels = seq_len(length(table(object$cluster))),
    cluster_assignments = cluster_asignments
  )
}

#' @export
extract_fit_summary.GMMCluster <- function(object, ...) {
  names <- paste0("Cluster_", seq_len(nrow(object$centroids)))

  cluster_asignments <- factor(
    names[object$clusters],
    levels = names
  )

  centroids <- tibble::as_tibble(object$centroids)
  covariance_matrices <- tibble::as_tibble(object$covariance_matrices)

  list(
    cluster_names = names,
    centroids = centroids,
    n_members = object$obs_per_cluster,
    n = object$n,
    m = object$m,
    k = object$k,
    df = object$df,
    aic = object$aic,
    bic = object$bic,
    loglik = object$loglik,
    log_likelihood = object$log_likelihood,
    probabilities = object$probabilities,
    covariance_matrices = object$covariance_matrices,
    weights = object$weights,
    me = object$me,
    me_between = object$me_between,
    me_within = object$me_within,
    me_total = object$me_total,
    me_ratio = object$me_ratio,
    mse = object$mse,
    mse_between = object$mse_between,
    mse_within = object$mse_within,
    mse_total = object$mse_total,
    mse_ratio = object$mse_ratio,
    se = object$se,
    se_within = object$se_within,
    se_between = object$se_between,
    se_total = object$se_total,
    se_ratio = object$se_ratio,
    sse = object$sse,
    sse_between = object$sse_between,
    sse_within = object$sse_within,
    sse_within_total_total = object$sse_within,
    sse_total = object$sse_total,
    sse_ratio = object$sse_ratio,
    orig_labels = object$clusters,
    cluster_assignments = cluster_asignments
  )
}

#' @export
extract_fit_summary.hclust <- function(object, ...) {
  clusts <- extract_cluster_assignment(object, ...)$.cluster
  n_clust <- dplyr::n_distinct(clusts)

  training_data <- attr(object, "training_data")

  overall_centroid <- colMeans(training_data)

  by_clust <- training_data %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      .cluster = clusts
    ) %>%
    dplyr::group_by(.cluster) %>%
    tidyr::nest()

  centroids <- by_clust$data %>%
    map(dplyr::summarize_all, mean) %>%
    dplyr::bind_rows()

  sse_within_total_total <- map2_dbl(
    by_clust$data,
    seq_len(n_clust),
    ~sum(
      philentropy::dist_many_many(
        as.matrix(centroids[.y, ]),
        as.matrix(.x),
        method = "euclidean"
      )
    )
  )

  list(
    cluster_names = unique(clusts),
    centroids = centroids,
    n_members = unname(as.integer(table(clusts))),
    sse_within_total_total = sse_within_total_total,
    sse_total = sum(
      philentropy::dist_many_many(
        t(overall_centroid),
        as.matrix(training_data),
        method = "euclidean"
      )
    ),
    orig_labels = NULL,
    cluster_assignments = clusts
  )
}
