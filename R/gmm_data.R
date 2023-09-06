# nocov start

make_gmm <- function() {
  modelenv::set_new_model("gmm")

  modelenv::set_model_mode("gmm", "partition")

  # ----------------------------------------------------------------------------

  modelenv::set_model_engine("gmm", "partition", "ClusterR")
  modelenv::set_dependency(
    model = "gmm",
    mode = "partition",
    eng = "ClusterR",
    pkg = "ClusterR"
  )
  modelenv::set_dependency(
    model = "gmm",
    mode = "partition",
    eng = "ClusterR",
    pkg = "tidyclust"
  )

  modelenv::set_fit(
    model = "gmm",
    eng = "ClusterR",
    mode = "partition",
    value = list(
      interface = "matrix",
      data = c(x = "data"),
      protect = c("data", "clusters"),
      func = c(pkg = "tidyclust", fun = "ClusterR_gmm_fit"),
      defaults = list()
    )
  )

  modelenv::set_encoding(
    model = "gmm",
    eng = "ClusterR",
    mode = "partition",
    options = list(
      predictor_indicators = "traditional",
      compute_intercept = TRUE,
      remove_intercept = TRUE,
      allow_sparse_x = FALSE
    )
  )

  modelenv::set_model_arg(
    model = "gmm",
    eng = "ClusterR",
    exposed = "num_clusters",
    original = "clusters",
    func = list(pkg = "dials", fun = "num_clusters"),
    has_submodel = TRUE
  )

  modelenv::set_pred(
    model = "gmm",
    eng = "ClusterR",
    mode = "partition",
    type = "cluster",
    value = list(
      pre = NULL,
      post = NULL,
      func = c(fun = ".gmm_predict_ClusterR"),
      args =
        list(
          object = rlang::expr(object$fit),
          new_data = rlang::expr(new_data)
        )
    )
  )

  modelenv::set_pred(
    model = "gmm",
    eng = "ClusterR",
    mode = "partition",
    type = "raw",
    value = list(
      pre = NULL,
      post = NULL,
      func = c(fun = ".gmm_predict_ClusterR"),
      args =
        list(
          object = rlang::expr(object$fit),
          new_data = rlang::expr(new_data),
          fuzzy = TRUE
        )
    )
  )
}

# nocov end
