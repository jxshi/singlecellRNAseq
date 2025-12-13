run_soupx <- function(toc,tod,rho=NULL) {

    #i make sure gene name are consistent
    tod <- tod[rownames(toc),]

    all <- toc
    all <- CreateSeuratObject(all)
    all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
    all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
    all.genes <- rownames(all)
    all <- ScaleData(all, features = all.genes)
    all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
    all <- FindNeighbors(all, dims = 1:30)
    all <- FindClusters(all, resolution = 0.5)
    all <- RunUMAP(all, dims = 1:30)
    
    matx <- all@meta.data

    sc = SoupChannel(tod, toc)
    sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
    # automated analysis of contamination rates
    if (is.null(rho)) {
        tryCatch(
            {sc = autoEstCont(sc)}, 
            error=function(e) { 
            # If encounters error, set rho rate to 0.2
                print("autoEstCont Error !")
            },
	    sc = setContaminationFraction(sc, 0.2)
        )
    }else{
    # set contamination fraction
        sc = setContaminationFraction(sc, rho)
    }
    # adjust matrix file
    out = adjustCounts(sc)
    # save data matrix
    #saveRDS(sc,"sc.rds")
    # Save adjusted matrix and export as 10X format
    DropletUtils:::write10xCounts(paste0("../datamatrix/", samples[i],"/soupX_filtered"), out, version="3")
}
