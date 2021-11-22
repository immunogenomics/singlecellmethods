#' @export
weighted_pca <- function(X, weights, genes_use=NULL, npc=20, do_corr=TRUE, scale_thresh=10) {
    if (!identical(length(weights), ncol(X))) {
        stop('Columns in X must match length of weights')
    }
    
#     y <- factor(y)
#     weights <- as.numeric((1 / prop.table(table(y)))[y]) / nlevels(y)
    if (any(is.na(weights))) {
        idx_keep <- which(is.na(weights))
#         y <- y[idx_keep]
        weights <- weights[idx_keep]
        X <- X[, idx_keep]
    }
    if (is.null(genes_use)) {
        genes_use <- row.names(X)
    } else if (length(genes_use) < nrow(X)) {
        if (any(!genes_use %in% row.names(X))) {
            stop('genes_use not in rownames of X')
        }
        X <- X[genes_use, ]
    }
    
    ## weighted z-scores
    mu <- rowMeans(X, weights = weights)
    sig <- rowSDs(X, weights = weights)
    
    vargenes_means_sds <- tibble(
        symbol = genes_use,
        mean = mu
    )
    vargenes_means_sds$stddev <- sig
    
    X <- scaleDataWithStats(X, mu, sig) 
    X <- X[which(is.na(rowSums(X)) == 0), ]
    if (do_corr) {
        X <- X %>% scale() %>% pmin(scale_thresh) %>% pmax(-scale_thresh)
    }
    
    ## weighted SVD
    pres <- RSpectra::svds(X %*% Matrix::Diagonal(x = sqrt(weights)), npc)
    V <- (Matrix::Diagonal(x = 1 / sqrt(weights)) %*% pres$v) %*% diag(pres$d)
    V <- as.matrix(V)
    colnames(V) <- paste0('PC', 1:npc)
    row.names(V) <- colnames(X)
    colnames(pres$u) <- paste0('PC', 1:npc)
    row.names(pres$u) <- row.names(X)
    return(list(loadings = pres$u, embeddings = V, vargenes = vargenes_means_sds))
}


## Wrapper to Seurat objects
#' @export
RunBalancedPCA <- function(obj, weight.by='orig.ident', npcs=20, assay.use='RNA', reduction.name = "pca", reduction.key = "PC_") {
    if (!weight.by %in% colnames(obj@meta.data)) 
        stop(glue('Variable {weight.by} not defined in this object'))
    if (length(unique(obj@meta.data[[weight.by]])) == 1) 
        stop(glue('Only 1 level found for {weight.by}. Balanced PCA requires at least 2.'))
    if (any(is.na(obj@meta.data[[weight.by]])))
        stop(glue('Balanced PCA must have non-NA values in {weight.by} variable for all cells'))
    
    y <- obj@meta.data[[weight.by]]
    weights <- as.numeric(((1 / table(y))[y]) * (length(y) / length(unique(y))))    
    pca_res <- weighted_pca(
        obj@assays[[assay.use]]@data, 
        weights, 
        obj@assays[[assay.use]]@var.features,
        npcs,
        scale_thresh=10
    )

    ## Put it back into Seurat 
    obj@assays[[assay.use]]@meta.features <- obj@assays[[assay.use]]@meta.features %>% 
        tibble::rownames_to_column('symbol') %>% 
        left_join(dplyr::rename(pca_res$vargenes, balanced_pca_mean = mean, balanced_pca_stddev = stddev), by = 'symbol') %>% 
        tibble::column_to_rownames('symbol')
    obj[[reduction.name]] <- Seurat::CreateDimReducObject(
        embeddings = pca_res$embeddings,
        loadings = pca_res$loadings, 
        stdev = as.numeric(apply(pca_res$embeddings, 2, stats::sd)),
        assay = assay.use,
        key = reduction.key
    )
    return(obj)
}

