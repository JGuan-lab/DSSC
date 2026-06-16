cran_packages <- c("dplyr", "Matrix", "gtools", "reshape2", "data.table",
                   "scales", "pheatmap", "BiocManager")

bioc_packages <- c("limma", "edgeR", "scater")

cran_install <- sapply(cran_packages, function(pck) {
    if (!require(pck, quietly = TRUE))
        install.packages(pck)
})

bioc_install <- sapply(bioc_packages, function(pck) {
    if (!require(pck, quietly = TRUE))
        BiocManager::install(pck)
})