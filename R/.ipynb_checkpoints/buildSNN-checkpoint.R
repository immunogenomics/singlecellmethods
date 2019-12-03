#' @export
buildSNN_fromFeatures <- function(V, prune_snn, nn_k = 30, nn_eps = 0) {
    nn_idx <- RANN::nn2(data = V, k = nn_k, searchtype = "standard", eps = nn_eps)$nn.idx
    buildSNN_fromKNN(nn_idx, prune_snn, nn_k)
}

#' @export
buildSNN_fromKNN <- function(nn_idx, prune_snn, nn_k = 30) {
    adj <- Matrix::sparseMatrix(i = rep(1:nrow(nn_idx), each = ncol(nn_idx)), 
                                j = c(t(nn_idx)), 
                                x = rep(1, prod(dim(nn_idx))))
    snn <- Matrix::tcrossprod(adj)
    snn@x <- snn@x / (2 * nn_k - snn@x)
    snn@x[snn@x < (prune_snn)] <- 0
    snn
}