
############ GO and KEGG enrichment analysis for each sub-cluster #################
library(clusterProfiler)
library(KEGG.db)
library(AnnotationDbi)
library(org.Mm.eg.db)
org = "org.Mm.eg.db"
library(stringr)
library(dplyr)

# keggpathid2name <- toTable(KEGG.db::KEGGPATHID2NAME)
# saveRDS(keggpathid2name, "~/software/functions/keggpathid2name.Rds")

Idents(sce.all) <- sce.all$subcelltype

sce.all.list <- SplitObject(sce.all, split.by = "subcelltype")
sce.all.list
cname <- names(sce.all.list)
fname <- gsub(' ', '', names(sce.all.list))
fname
saveRDS(cname, "subcelltype_name.Rds")
saveRDS(fname, "subcelltype_file_folder.Rds")

### Read keggpathid2name 
keggpathid2name  <- readRDS("~/software/functions/keggpathid2name.Rds")

# Separate up- and down-regulated genes and run GO and KEGG enrichment analysis
# Revised on Jan 13, 2024
top = 200

upTopMarkers <- sce.markers %>% filter(avg_log2FC > 0.25 & p_val_adj < 0.01) %>% group_by(cluster) %>% 
  slice_max(n = top, order_by = avg_log2FC)
write.csv(upTopMarkers,'top200UPmarkers.csv', row.names=F)

downTopMarkers <- sce.markers %>% filter(avg_log2FC < -0.25 & p_val_adj < 0.01) %>% group_by(cluster) %>% 
  slice_max(n = top, order_by = -avg_log2FC)
write.csv(downTopMarkers,'top200DOWNmarkers.csv', row.names=F)

TopMarkers <- rbind(upTopMarkers, downTopMarkers)
write.csv(TopMarkers,'top400markers.csv', row.names=F)
##TopMarkers为前面获得的每个亚群的top200的高表达的基因，ann为自己手动整理的注释及基因转换id的文件，将TopMarkers的geneid为标准，
# 进行取交集，获得TopMarkers里面基因的注释结果和geneid号

org = "org.Mm.eg.db"

all.genes <- rownames(sce.all)

ann <- bitr(all.genes, fromType = "SYMBOL", toType = c("ENTREZID",'ENSEMBL'), OrgDb = org, drop = T)
colnames(ann) <- c("gene","ENTREZID","ENSEMBL")

alllocus <- left_join(TopMarkers, ann, by = "gene", relationship = "many-to-many")
alllocus <- alllocus[which(alllocus$ENTREZID != "NA"),]

uplocus <- left_join(upTopMarkers, ann, by = "gene", relationship = "many-to-many")
uplocus <- uplocus[which(uplocus$ENTREZID != "NA"),]

downlocus <- left_join(downTopMarkers, ann, by = "gene", relationship = "many-to-many")
downlocus <- downlocus[which(downlocus$ENTREZID != "NA"),]


# using the following loop to do GO and KEGG enrichment analysis
for (j in seq(length(unique(alllocus$cluster)))){
  # GO and KEGG enrichmeng analysis for both up- and down-regulated genes
  allcluster <- alllocus[which(alllocus$cluster == cname[j]),]
  
  allego <- enrichGO(allcluster$ENTREZID, OrgDb = org, keyType = "ENTREZID", pvalueCutoff = 1, qvalueCutoff = 1, ont = "all")
  allego <- setReadable(allego, OrgDb = org, keyType = "ENTREZID")
  write.table(allego@result, file = paste0("cluster_",fname[j],"_all_GO_enrichment.tsv"), quote = F, sep="\t", col.names = NA)
  # dotplot(allego, split="ONTOLOGY",showCategory = 20,label_format=100) + facet_grid(ONTOLOGY~., scale="free")
  barplot(allego, split="ONTOLOGY",showCategory = 20,label_format=100) + facet_grid(ONTOLOGY~., scale="free")
  ggsave(paste0("cluster_", fname[j], "_all_GO_barplot.pdf"), device=cairo_pdf, width =10, height = 16)
  
  # KEGG enrichment analysis
  allcompKEGG = enrichKEGG(allcluster$ENTREZID, organism="mmu", pvalueCutoff=1, keyType='kegg', use_internal_data = T)
  allcompKEGG <- setReadable(allcompKEGG, OrgDb = org, keyType = "ENTREZID")
  # Merge the two datasets based on the ID and path_id columns
  merged_data <- merge(allcompKEGG@result, keggpathid2name, by.x = "ID", by.y = "path_id", all.x = TRUE)
  merged_data <- dplyr::left_join(allcompKEGG@result, merged_data)
  
  # Assign the path_name to Description where there is a match
  allcompKEGG@result$Description <- paste0(merged_data$path_name)
  
  write.table(allcompKEGG@result, file = paste0("cluster_",fname[j],"_all_KEGG_enrichment.tsv"), quote = F, sep="\t", col.names = NA)
  barplot(allcompKEGG, showCategory = 15,label_format=100, title = paste0("Cluster ",cname[j]," KEGG Pathway Enrichment Analysis"))
  ggsave(paste0("cluster_", fname[j], "_all_KEGG_barplot.pdf"), device=cairo_pdf, width =10, height = 16)
  
  # GO and KEGG enrichmeng analysis for up-regulated genes
  upcluster <- uplocus[which(uplocus$cluster == cname[j]),]
  
  upego <- enrichGO(upcluster$ENTREZID, OrgDb = org, keyType = "ENTREZID", pvalueCutoff = 1, qvalueCutoff = 1, ont = "all")
  upego <- setReadable(upego, OrgDb = org, keyType = "ENTREZID")
  write.table(upego@result, file = paste0("cluster_",fname[j],"_up_GO_enrichment.tsv"), quote = F, sep="\t", col.names = NA)
  # dotplot(upego, split="ONTOLOGY",showCategory = 20,label_format=100) + facet_grid(ONTOLOGY~., scale="free")
  barplot(upego, split="ONTOLOGY",showCategory = 20,label_format=100) + facet_grid(ONTOLOGY~., scale="free")
  ggsave(paste0("cluster_", fname[j], "_up_GO_barplot.pdf"), device=cairo_pdf, width =10, height = 16)
  
  # KEGG enrichment analysis
  upcompKEGG = enrichKEGG(upcluster$ENTREZID, organism="mmu", pvalueCutoff=1, keyType='kegg', use_internal_data = T)
  upcompKEGG <- setReadable(upcompKEGG, OrgDb = org, keyType = "ENTREZID")
  # Merge the two datasets based on the ID and path_id columns
  merged_data <- merge(upcompKEGG@result, keggpathid2name, by.x = "ID", by.y = "path_id", all.x = TRUE)
  merged_data <- dplyr::left_join(upcompKEGG@result, merged_data)
  
  # Assign the path_name to Description where there is a match
  upcompKEGG@result$Description <- paste0(merged_data$path_name)
  write.table(upcompKEGG@result, file = paste0("cluster_", fname[j],"_up_KEGG_enrichment.tsv"), quote = F, sep="\t", col.names = NA)
  barplot(upcompKEGG, showCategory = 15,label_format=100, title = paste0("Cluster ", cname[j]," KEGG Pathway Enrichment Analysis"))
  ggsave(paste0("cluster_", fname[j], "_up_KEGG_barplot.pdf"), device=cairo_pdf, width =10, height = 16)
  
  # GO and KEGG enrichmeng analysis for down-regulated genes
  downcluster <- downlocus[which(downlocus$cluster == cname[j]),]
  
  downego <- enrichGO(downcluster$ENTREZID, OrgDb = org, keyType = "ENTREZID", pvalueCutoff = 1, qvalueCutoff = 1, ont = "all")
  downego <- setReadable(downego, OrgDb = org, keyType = "ENTREZID")
  write.table(downego@result, file = paste0("cluster_",fname[j],"_down_GO_enrichment.tsv"), quote = F, sep="\t", col.names = NA)
  # dotplot(downego, split="ONTOLOGY",showCategory = 20,label_format=100) + facet_grid(ONTOLOGY~., scale="free")
  barplot(downego, split="ONTOLOGY",showCategory = 20,label_format=100) + facet_grid(ONTOLOGY~., scale="free")
  ggsave(paste0("cluster_", fname[j], "_down_GO_barplot.pdf"), device=cairo_pdf, width =10, height = 16)
  
  # KEGG enrichment analysis
  downcompKEGG = enrichKEGG(downcluster$ENTREZID, organism="mmu", pvalueCutoff=1, keyType='kegg',use_internal_data = T)
  downcompKEGG <- setReadable(downcompKEGG, OrgDb = org, keyType = "ENTREZID")
  # Merge the two datasets based on the ID and path_id columns
  merged_data <- merge(downcompKEGG@result, keggpathid2name, by.x = "ID", by.y = "path_id", all.x = TRUE)
  merged_data <- dplyr::left_join(downcompKEGG@result, merged_data)
  
  # Assign the path_name to Description where there is a match
  downcompKEGG@result$Description <- paste0(merged_data$path_name)
  
  write.table(downcompKEGG@result, file = paste0("cluster_",fname[j],"_down_KEGG_enrichment.tsv"), quote = F, sep="\t", col.names = NA)
  barplot(downcompKEGG, showCategory = 15, label_format=100, title = paste0("Cluster ",cname[j]," KEGG Pathway Enrichment Analysis"))
  ggsave(paste0("cluster_", fname[j], "_down_KEGG_barplot.pdf"), device=cairo_pdf, width =10, height = 16)
}