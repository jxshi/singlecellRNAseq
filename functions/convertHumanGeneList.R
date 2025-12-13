# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){

require("biomaRt")
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",verbose = TRUE, host = "https://asia.ensembl.org")
#mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",verbose = TRUE, host = "https://asia.ensembl.org")

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",verbose = TRUE, host = "https://dec2021.archive.ensembl.org")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",verbose = TRUE, host = "https://dec2021.archive.ensembl.org")

genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

humanx <- unique(genesV2[, 2])

# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}
