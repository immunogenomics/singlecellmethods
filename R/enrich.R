meanScaledRows <- function(X, genes_use, sigma, mu) {
    if (missing(genes_use)) {
        genes_use <- rownames(X)
    } else {
        if (!all(genes_use %in% rownames(X))) {
            l_orig <- length(genes_use)
            genes_use <- intersect(genes_use, rownames(X))
            warning(sprintf('removing some genes, keeping %d/%d', length(genes_use), l_orig))
        }
    }
    X <- X[genes_use, ]
    ## TODO: check that mu and sigma are named vectors 
    if (missing(sigma)) {
        sigma <- rowSDs(X)
    }
    if (missing(mu)) {
        mu <- rowMeans(X)
    }
    
    if ('dgCMatrix' %in% class(X)) {
#         geneset_mask <- genes_use %in% rownames(X)
#         res <- divRows_dgc(X@x, X@p, X@i, sigma[genes_use], geneset_mask, ncol(X))[, 1]
        res <- enrich_dgc(X@x, X@p, X@i, sigma[genes_use], ncol(X))[, 1]
                
    } else {
        res <- Matrix::colSums(X / sigma[genes_use])
    }

    res <- (1 / length(genes_use)) * (res - sum(mu[genes_use] / sigma[genes_use]))        
    return(res)
}


do_enrich <- function(exprs_norm, genesets, min_size=0, max_size=Inf) {
    ## TODO: what if only 1 geneset? 
    ## TODO: add parallel
    ## TODO: what if only 1 gene?
    ## TODO: check for SD=0
    
    ## filter genes and genesets
    genesets <- lapply(genesets, function(x) intersect(x, rownames(exprs_norm)))
    gsl <- lapply(genesets, length)
    genesets <- genesets[which(gsl >= min_size & gsl <= max_size)]

    ## precompute some statistics
    genes_all <- Reduce(union, genesets)
    exprs_norm <- exprs_norm[genes_all, ]
    sigma <- rowSDs(exprs_norm)
    mu <- rowMeans(exprs_norm)

    X <- genesets %>% 
        lapply(function(genes_use) {
            meanScaledRows(exprs_norm, genes_use, sigma, mu)#[, 1]
        }) %>% 
        purrr::reduce(cbind)
    colnames(X) <- names(genesets)
    return(X)
}
                       
                       