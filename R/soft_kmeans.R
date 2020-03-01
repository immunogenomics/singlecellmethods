## colums are observations
#' @export
soft_kmeans <- function(X, k, w, max_iter=20, sigma=0.1) {
    message('WARNING: soft_kmeans fxn uses cosine distance only')
    Z <- cosine_normalize_cpp(X, 2)
    if (missing(w))
    Y <- stats::kmeans(t(Z), centers = k, iter.max = 25, nstart = 10)$centers %>% t() ## D x K
    res <- soft_kmeans_cpp(Y, Z, max_iter, sigma)
    return(res)
}

