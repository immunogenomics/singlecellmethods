normalizeData <- function(A, scaling_factor = 1e4, method) {
    if (method == "log") {
        A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
        A@x <- scaling_factor * A@x
        A@x <- log(1 + A@x)
    } else if (method == "fft") {
        A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
        A@x <- scaling_factor * A@x
        A@x <- sqrt(A@x) + sqrt(1 + A@x)
    } else if (method == "geneCLR") {
        A@x <- as.numeric(normalizeCLR_dgc(A@x, A@p, A@i, ncol(A), nrow(A), 1))        
    } else if (method == "cellCLR") {
        A@x <- as.numeric(normalizeCLR_dgc(A@x, A@p, A@i, ncol(A), nrow(A), 2))
    } else {
        stop(sprintf("ERROR: method %s not implemented", method))
    }

	return(A)
}