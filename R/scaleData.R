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

scaleDataWithStats <- function(A, margin = 1, thresh = 10) {
    if (!"dgCMatrix" %in% class(A))
        A <- as(A, "dgCMatrix")
    
    if (margin != 1) A <- t(A)
    
    mean_vec <- Matrix::rowMeans(A)
    sd_vec <- rowSDs_dgc(A@x, A@p, A@i, mean_vec, ncol(A), nrow(A))
    res <- scaleRowsWithStats_dgc(A@x, A@p, A@i, mean_vec, sd_ved, 
                                  ncol(A), nrow(A), thresh)
    if (margin != 1) res <- t(res)
    row.names(res) <- row.names(A)
    colnames(res) <- colnames(A)
    return(res)
}


rowSDs <- function(A, row_means=NULL) {
    if (!is.null(row_means)) {
        row_means <- Matrix::rowMeans(A)
    }
    res <- rowSDs_dgc(A@x, A@p, A@i, row_means, ncol(A), nrow(A))
    names(res) <- row.names(A)
    return(res)
}
