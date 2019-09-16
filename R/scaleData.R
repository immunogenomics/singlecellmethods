scaleData <- function(A, margin = 1, thresh = 10) {
    if (!"dgCMatrix" %in% class(A))
        A <- as(A, "dgCMatrix")
    
    if (margin != 1) A <- t(A)
    
    res <- scaleRows_dgc(A@x, A@p, A@i, ncol(A), nrow(A), thresh)
    if (margin != 1) res <- t(res)
    row.names(res) <- row.names(A)
    colnames(res) <- colnames(A)
    return(res)
}

rowSDs <- function(A) {
    if (!"dgCMatrix" %in% class(A))
        A <- as(A, "dgCMatrix")

    if (margin != 1) A <- t(A)
    row_means <- Matrix::rowMeans(A)
    res <- rowSDs_dgc(A@x, A@p, A@i, row_means, ncol(A), nrow(A))
    names(res) <- row.names(A)
    return(res)
}
