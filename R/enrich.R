sqrtMeanScaledRows <- function(X, genes_use, sigma, mu) {
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

    res <- (1 / sqrt(length(genes_use))) * (res - sum(mu[genes_use] / sigma[genes_use]))        
    return(res)
}


enrich_cells <- function(
    exprs_norm, genesets, min_size=0, max_size=Inf, 
    mode=c('unweighted', 'weighted')[1],
    do_par=FALSE
) {
    ## TODO: what if only 1 geneset? 
    ## TODO: what if only 1 gene?
    ## TODO: check for SD=0
    
    ## filter genes and genesets
    genesets <- lapply(genesets, function(x) intersect(x, rownames(exprs_norm)))
    gsl <- lapply(genesets, length)
    genesets <- genesets[which(gsl >= min_size & gsl <= max_size)]

    ## Keep only the relevant genes 
    genes_all <- Reduce(union, genesets)
    exprs_norm <- exprs_norm[genes_all, ]
    
    ## Setup stuff for parallel execution
    if (do_par) {
        plan(multicore)
        iter_fxn <- furrr::future_map
    } else {
        iter_fxn <- purrr::map
    }

    ## Do the appropriate method 
    if (mode == 'unweighted') {
        res <- enrich_cells.unweighted(exprs_norm, genesets, iter_fxn)
    } else if (mode == 'weighted') {
        res <- enrich_cells.weighted(exprs_norm, genesets, iter_fxn)
    }
    return(res)
}

enrich_cells.unweighted <- function(exprs_norm, genesets, iter_fxn) {
    ## precompute some statistics
    sigma <- rowSDs(exprs_norm)
    mu <- rowMeans(exprs_norm)
    
    X <- genesets %>% 
        iter_fxn(function(genes_use) {
            sqrtMeanScaledRows(exprs_norm, genes_use, sigma, mu)
        }) %>% 
        bind_cols() %>% 
        data.frame()
    
    colnames(X) <- names(genesets)
    rownames(X) <- colnames(exprs_norm)
    return(X)
}

enrich_cells.weighted <- function(exprs_norm, genesets, do_par=TRUE) {
    ## Get PC1 embeddings and loadings for each pathway
    pca_res <- genesets %>% 
        map(intersect, rownames(exprs_norm)) %>% 
        future_map(function(genes_use) {
            X <- Matrix::t(exprs_norm[genes_use, ])
            ## RSpectra does implicit scaling, which is more memory efficient
            res <- RSpectra::svds(X, k=1, opts=list(center=TRUE, scale=TRUE))
            loadings <- as.numeric(res$v)
            names(loadings) <- genes_use
            list(scores = res$u * res$d * sqrt(max(dim(X)) - 1), loadings = loadings)
        }) 
    
    ## bind the cell-specific scores for each pathway into table
    scores_df <- pca_res %>% map('scores') %>% bind_cols() %>% data.frame()
    rownames(scores_df) <- colnames(exprs_norm)
    
    ## number of genes is different for each pathway, so we can't bind into table
    loadings <- pca_res %>% map('loadings')
    
    ## because PCA is unique up to a sign flip, let's make sure that each pathway
    ## has mostly positive loadings, making it easier to interpret the scores
    flip_sign <- loadings %>% map(function(x) as.numeric(sum(x > 0)) / length(x))
    flip_sign <- c(-1, 1)[1 + as.integer(flip_sign > .5)]
    names(flip_sign) <- names(loadings)
    
    loadings <- map2(loadings, flip_sign, `*`)
    scores_df <- data.frame(as.matrix(scores_df) %*% diag(x = flip_sign))
    colnames(scores_df) <- names(loadings)
    # return(scores_df)
    return(list(scores = scores_df, loadings = loadings))
}






