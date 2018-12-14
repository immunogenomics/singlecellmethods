## TODO: re-cluster NA points
leidenCluster <- function(snn, resolution_list, min_cluster_size = 20, verbose = TRUE, pythondir) {
    library(reticulate)
    use_python(pythondir)
    leiden <- import("leidenalg")
    ig <- import("igraph")
    
    tmpfile <- tempfile()
    snn <- as(snn, "dgTMatrix")
    data.table(i = snn@i, j = snn@j) %>% fwrite(tmpfile, sep = " ", col.names = FALSE)
    g <- ig$Graph(directed = FALSE)
    g <- g$Read_Edgelist(tmpfile)
    
    .res <- Reduce(cbind, lapply(resolution_list, function(resolution) {
        cres <- leiden$find_partition(g, leiden$CPMVertexPartition, resolution_parameter = resolution)$membership + 1
        clusters_keep <- names(table(cres))[which(table(cres) >= min_cluster_size)]
        cres[which(!cres %in% clusters_keep)] <- NA
        n_na <- sum(is.na(cres))
        message(sprintf("Resolution %f yielded %d clusters", resolution, length(clusters_keep)))
        if (verbose & n_na > 0) message(sprintf("WARNING: for resolution %f, %d/%d points removed from small clusters", resolution, n_na, length(cres)))
        return(cres)
    }))
    
#     if (max(1 + ceiling(-log10(resolution_list))) > 4) {
#         message("WARNING: colnames may not have enough precision")
#     }

    colnames(.res) <- sprintf("res_%0.4e", resolution_list)
    return(as.data.frame(.res))    
}