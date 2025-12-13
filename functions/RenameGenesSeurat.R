# RenameGenesSeurat_Spatial  ------------------------------------------------------------------------------------
RenameGenesSeurat_Spatial <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial

  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}
# RenameGenesSeurat(obj = SeuratObj, newnames = HGNC.updated.genes)

# RenameGenesSeurat  ------------------------------------------------------------------------------------
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
# RenameGenesSeurat(obj = SeuratObj, newnames = HGNC.updated.genes)
