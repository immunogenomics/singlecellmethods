#' @export
scaleData <- function(A, margin = 1, thresh = 10) {
    A <- as(A, "dgCMatrix")
    
    if (margin != 1) A <- t(A)
    
    res <- scaleRows_dgc(A@x, A@p, A@i, ncol(A), nrow(A), thresh)
    if (margin != 1) res <- t(res)
    row.names(res) <- row.names(A)
    colnames(res) <- colnames(A)
    return(res)
}

#' @export
scaleDataWithStats <- function(A, mean_vec, sd_vec, margin = 1, thresh = 10) {
    if (!"dgCMatrix" %in% class(A))
        A <- as(A, "dgCMatrix")
    
    if (margin != 1) A <- t(A)
    
    res <- scaleRowsWithStats_dgc(A@x, A@p, A@i, mean_vec, sd_vec, 
                                  ncol(A), nrow(A), thresh)
    if (margin != 1) res <- t(res)
    row.names(res) <- row.names(A)
    colnames(res) <- colnames(A)
    return(res)
}


#' @export
rowSDs <- function(A, row_means=NULL, weights=NULL) {
    if (is.null(row_means)) {
#         row_means <- Matrix::rowMeans(A)
        row_means <- singlecellmethods::rowMeans(A, weights)
    }
    if (is.null(weights)) {
        res <- as.numeric(rowSDs_dgc(A@x, A@p, A@i, row_means, ncol(A), nrow(A)))
    } else {
        res <- as.numeric(rowSDsWeighted_dgc(A@x, A@p, A@i, row_means, weights, ncol(A), nrow(A)))
    }
    names(res) <- row.names(A)
    return(res)
}

#' @export
rowMeans <- function(A, weights=NULL) {
    if (is.null(weights)) {
        res <- Matrix::rowMeans(A)
    } else {
        res <- as.numeric(rowMeansWeighted_dgc(A@x, A@p, A@i, weights, ncol(A), nrow(A)))
    }
    names(res) <- row.names(A)
    return(res)
}