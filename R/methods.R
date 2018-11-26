findVariableGenes <- function(X, groups, min_expr = .1, max_expr = Inf, 
                              min_dispersion = 0, max_dispersion = Inf, return_top_n = 0) {
    ## TODO: check that groups are 0 indexed
    groups <- factor(groups)
    groups_int <- as.integer(factor(groups)) - 1
    groups_table <- table(groups_int)
    
    ## initially compute means in non-log space, to use in vmr function below
    means_nonlog <- exp_mean(X@x, X@p, X@i, ncol(X), nrow(X), groups_int, groups_table)
    colnames(means_nonlog) <- levels(groups)
    row.names(means_nonlog) <- row.names(X)
    
    vmr <- log_vmr(X@x, X@p, X@i, ncol(X), nrow(X), means_nonlog, groups_int, groups_table)    
    colnames(vmr) <- levels(groups)
    row.names(vmr) <- row.names(X)
    
    ## transform means to logspace and join means and VMR  
    vargenes_df <- dplyr::inner_join(
        means_nonlog %>% log1p %>% as.data.frame() %>% 
            tibble::rownames_to_column("symbol") %>% 
            tidyr::gather(group, gene_mean, -symbol),
        vmr %>% as.data.frame() %>% 
            tibble::rownames_to_column("symbol") %>% 
            tidyr::gather(group, gene_dispersion, -symbol), 
        by = c("symbol", "group")
    ) %>% 
        dplyr::arrange(-gene_dispersion) %>% 
        subset(gene_mean >= min_expr & gene_mean <= max_expr) %>%
        subset(gene_dispersion >= min_dispersion & gene_dispersion <= max_dispersion)    
    
    if (return_top_n > 0) {
        vargenes_union <- unique(data.table(vargenes_df)[, head(.SD, return_top_n), by = group][, symbol])
        return(vargenes_union)
    } else {
        return(vargenes_df)
    }
    
}

