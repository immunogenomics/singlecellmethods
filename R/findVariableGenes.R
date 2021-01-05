#' @export
findVariableGenes <- function(X, groups, min_expr = .1, max_expr = Inf, 
                               min_dispersion = 0, max_dispersion = Inf, 
                               num.bin = 20, binning.method = "equal_width", return_top_n = 0) {
    ## TODO: check that groups are 0 indexed
    groups <- factor(groups)
    groups_int <- as.integer(factor(groups)) - 1
    groups_table <- table(groups_int)
    
    ## initially compute means in non-log space, to use in vmr function below
    means_nonlog <- singlecellmethods:::exp_mean(X@x, X@p, X@i, ncol(X), nrow(X), groups_int, groups_table)
    colnames(means_nonlog) <- levels(groups)
    
    vmr <- singlecellmethods:::log_vmr(X@x, X@p, X@i, ncol(X), nrow(X), means_nonlog, groups_int, groups_table)    
    colnames(vmr) <- levels(groups)    

    ## transform means to logspace and join means and VMR  
    vargenes_df <- dplyr::inner_join(
        means_nonlog %>% log1p %>% as_tibble() %>% 
            cbind(symbol = row.names(X)) %>% 
            tidyr::gather(group, gene_mean, -symbol),
        vmr %>% as_tibble() %>% 
            cbind(symbol = row.names(X)) %>% 
            tidyr::gather(group, gene_dispersion, -symbol), 
        by = c("symbol", "group")
    )
    
    if (num.bin > 0) {
        if (binning.method == "equal_width") {
            .breaks <- num.bin
        }
        else if (binning.method == "equal_frequency") {
            .breaks <- c(-1, quantile(vargenes_df$gene_mean[vargenes_df$gene_mean > 0], probs = seq(0, 1, length.out = num.bin)))
        }
        else {
            stop(paste0("Invalid selection: '", binning.method, "' for 'binning.method'."))
        }
        
        vargenes_df <- data.table(vargenes_df)[
            , the_bin := cut(gene_mean, .breaks), by = group
        ][]
        
        vargenes_df <- data.table(vargenes_df)[
            , gene_dispersion_scaled := scale(gene_dispersion), by = .(the_bin, group)
        ][]
        
        vargenes_df <- data.table(vargenes_df)[, the_bin := NULL][]
    }
    
    vargenes_df <- vargenes_df %>% 
        dplyr::arrange(-gene_dispersion) %>% 
        subset(gene_mean >= min_expr & gene_mean <= max_expr) %>%
        subset(gene_dispersion >= min_dispersion & gene_dispersion <= max_dispersion) 

    return(vargenes_df)
#     if (return_top_n > 0) {
#         vargenes_union <- unique(data.table(vargenes_df)[, head(.SD, return_top_n), by = group][, symbol])
#         return(vargenes_union)
#     } else {
#         return(vargenes_df)
#     }
    
}


#' @export
vargenes_vst <- function(object, groups, topn, loess.span = 0.3) {
    clip.max <- sqrt(ncol(object))

    N <- ncol(object)
    if (missing(groups)) {
        groups <- rep('A', N)
    }
    
    res <- split(seq_len(N), groups) %>% lapply(function(idx) {
        object_group <- object[, idx]
        ## row means
        hvf.info <- data.frame(
          symbol = rownames(object_group), 
          mean = Matrix::rowMeans(object_group)
        )

        ## row vars
        hvf.info$variance <- rowVars(object_group, hvf.info$mean)

        ## initialize
        hvf.info$variance.expected <- 0
        hvf.info$variance.standardized <- 0

        not.const <- hvf.info$variance > 0

        ## loess curve fit 
        suppressWarnings({
            fit <- loess(formula = log10(variance) ~ log10(mean), 
                data = hvf.info[not.const, ], span = loess.span)            
        })

        ## extract fitted variance 
        hvf.info$variance.expected[not.const] <- 10^fit$fitted

        ## get row standard deviations after clipping
        hvf.info$variance.standardized <- rowVarsStd(
            object_group, 
            hvf.info$mean, 
            sqrt(hvf.info$variance.expected), 
            clip.max
        )

        hvf.info <- hvf.info %>% 
#             tibble::rownames_to_column('symbol') %>% 
            arrange(-variance.standardized) %>% 
            tibble::rowid_to_column('rank') %>% 
            transform(group = unique(groups[idx]))

        return(hvf.info)        
    })
    
    
    if (missing(topn)) {
        ## MODE 1: return table 
        res <- Reduce(rbind, res) %>% 
            dplyr::select(group, symbol, rank, everything())

        if (length(unique(res$group)) == 1) {
            res$group <- NULL
        }
    } else {
        ## MODE 2: return genes
        res <- Reduce(union, lapply(res, function(x) head(x, topn)$symbol))
    }
    return(res)
}
