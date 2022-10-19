library(readr)
library(msigdbr)
library(Seurat)
library(patchwork)
library(tidyr)
library(DESeq2)
library(readxl)
library(magrittr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(forcats)
library(ggstance)
library(ggplot2)
library(gprofiler2)

dataset1_topgenes <- read_csv("NTU_BMDS/BS6202_Techniques_in_Biomedical_Data_Mining/Group project/dataset1_topgenes.csv")
dataset1_gfs <- read_csv("NTU_BMDS/BS6202_Techniques_in_Biomedical_Data_Mining/Group project/GSE76275GFS_DEGs.csv")
### 25 intersected genes in dataset 2 ###
length(intersect(dataset1_gfs$Probe, dataset1_topgenes$PROBEID))

index_gfs <- as.vector((intersect(dataset1_gfs$Probe, dataset2_gfs$Probe)))
gfs_inter <- dataset1_gfs[dataset1_gfs$Probe %in% index_gfs,]

index_limma <- as.vector((intersect(dataset1_topgenes$PROBEID, dataset2_topgenes$ID)))
limma_inter <- dataset1_topgenes[dataset1_topgenes$PROBEID %in% index_limma,]

dataset2_gfs <- read_csv("NTU_BMDS/BS6202_Techniques_in_Biomedical_Data_Mining/Group project/GSE43358_DEG.csv")
dataset2_topgenes <- read_csv("NTU_BMDS/BS6202_Techniques_in_Biomedical_Data_Mining/Group project/dataset2_topgenes.csv")
length(intersect(dataset2_gfs$Probe, dataset2_topgenes$ID))

tmp <- "gp__vIDE_1rtb_4YY"
go_df <- read.gmt("./gprofiler_full_hsapiens.name.gmt")
GSE43358_kru <- read_csv("GSE43358_DEG.csv")

# upload GMT file( "x.gmt") 
# tmp <- 'gp__4unl_zQBZ_rTg' #ensg ids

ORA_slim <- function(query, threshold = 1, ordered_query = FALSE, background = dge_analysis,
                     exclude_iea = FALSE, drop_parents = FALSE, intersection_min = 2, significant = F,
                     term_size_min = 15, term_size_max = 500, correction_method = "fdr") {
  
  gp_res <- gost(query = as.character(query$Row.names),
                 organism = tmp,#"hsapiens",
                 multi_query = FALSE,
                 significant = significant, # Only return significant results  
                 exclude_iea = exclude_iea, # Whether to exclude electronically inferred annotations.
                 measure_underrepresentation = FALSE,
                 evcodes = TRUE, # Get evidence codes and intersection column
                 user_threshold = threshold,
                 correction_method = correction_method,  
                 domain_scope = "custom", #  background is our custom set of genes.
                 custom_bg = as.character(background$Row.names)
  )
  gp_res_terms <- gp_res$result
  
  # Filter by term size between 15 and 500.
  gp_res_terms <- gp_res_terms[gp_res_terms$term_size >= term_size_min & gp_res_terms$term_size <= term_size_max,]
  gp_res_terms <- gp_res_terms[gp_res_terms$intersection_size >= intersection_min,]
  
  # Create -log(FDR)
  #gp_res_terms$logp <- -log10(gp_res_terms$p_value)
  
  return(gp_res_terms)
}

limma_dataset1 <- as.data.frame(cbind(dataset1_topgenes$SYMBOL, dataset1_topgenes$logFC, dataset1_topgenes$adj.P.Val))

colnames(limma_dataset2) <- c("Row.names","log2FC","adj.P")
selected_genes_1_limma <-subset(limma_dataset1, (log2FC > 1.2 | log2FC < -1.5))   
write_csv(limma_dataset1, "limma_dataset1.csv")

sel_limma <- read_csv("sel_limma.csv")
sel_limma_dataset1 <-read_csv("limma_dataset1.csv")


limma_dataset2 <- read_csv("NTU_BMDS/BS6202_Techniques_in_Biomedical_Data_Mining/dataset2_back.csv")
sel_limma_dataset2 <- read_csv("NTU_BMDS/BS6202_Techniques_in_Biomedical_Data_Mining/sel_limma2.csv")


sel_gfs1 <- read_csv("NTU_BMDS/BS6202_Techniques_in_Biomedical_Data_Mining/sel_gfs1.csv")
back_gfs1 <- read_csv("NTU_BMDS/BS6202_Techniques_in_Biomedical_Data_Mining/back_gfs1.csv")

sel_gfs2 <- read_csv("NTU_BMDS/BS6202_Techniques_in_Biomedical_Data_Mining/sel_gfs2.csv")
back_gfs2 <- read_csv("NTU_BMDS/BS6202_Techniques_in_Biomedical_Data_Mining/back_gfs2.csv")

ora_limma_1 <- ORA_slim(query=sel_limma_dataset1, background = limma_dataset1)
ora_limma_2<- ORA_slim(query=sel_limma_dataset2, background = limma_dataset2)
ora_gfs_1<- ORA_slim(query=sel_gfs1, background = back_gfs1)
ora_gfs_2<- ORA_slim(query=sel_gfs2, background = back_gfs2)

fm_df <- go_df %>% filter(term %in% ora_limma$term_id)
fm_df <- merge(fm_df, ora_limma, by.x = "term", by.y="term_id") %>%
  select(term_name,gene)
geneList_limma_1 <- sel_limma$log2FC
sel_limma <- sel_limma[order(-sel_limma$log2FC),]
names(geneList_limma_1) <- sel_limma$Row.names

# names(geneList_tnbc_normal_epithelial) <- resdata_epi$Row.names
goo_all_tnbc_normal_epithelial <- GSEA(geneList_limma_1, 
                                       TERM2GENE = fm_df, pvalueCutoff = 1)
goo_all_tnbc_normal_epithelial <- goo_all_tnbc_normal_epithelial@result

library(forcats)
library(ggstance)
y <- arrange(goo_all_tnbc_normal_epithelial, abs(NES)) %>% group_by(sign(NES))
ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalues),showCatagory) + ggtitle("Epithelial cells (TNBC vs Normal) [After ORA]") + theme(plot.title = element_text(hjust = 0.5))  +
  geom_barh(stat='identity') + scale_fill_continuous(low = 'red', high = 'blue') + theme_minimal() + ylab(NULL)   + theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()) + theme(axis.text = element_text(size = 15))

##############################
limma_all <-read_csv("limma_all.csv")

### GSEA ###
go_slim_GeneNames <- read_csv("./scRNA_TNBC/go_slim_GeneNames.csv")
deseq2_tnbc_normal_genes <- read_csv("Desktop/deseq2_tnbc_normal_genes.csv")
goslim68 <- read.csv("./scRNA_TNBC/genes_658_goslim.txt", header=F)
go_slim_GOTERMS <- read.csv("./scRNA_TNBC/go_slim_GOTERMS.csv")

m_df <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  filter(gs_subcat == "GO:BP")#%>% select(gs_name, human_gene_symbol)
m_df <- go_slim_GOTERMS %>% 
  filter(term %in% goslim68$V1) %>%
  select(GO, name)

limma_all <- limma_all[order(-limma_all$log2FC),]
limma_geneList <- limma_all$log2FC
names(limma_geneList) <- limma_all$Row.names


deg_GSE43358 <- deg_GSE43358[order(-deg_GSE43358$logFC),]
geneList <- deg_GSE43358$logFC
names(geneList) <- deg_GSE43358$gene_name
goo_all <- GSEA(limma_geneList, TERM2GENE = m_df, pvalueCutoff = 1)
goo_all <- goo_all@result
y <- arrange(goo_all, abs(NES)) %>% group_by(sign(NES))
ggplot(y, aes(NES, fct_reorder(Description, NES), fill=pvalue),showCatagory) + ggtitle("Limma intersected genes enrichment") + theme(plot.title = element_text(hjust = 0.5))  +
  geom_barh(stat='identity') + scale_fill_continuous(low = 'red', high = 'blue') + theme_minimal() + ylab(NULL)   + theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()) + theme(axis.text = element_text(size = 15))

dataset2_topgenes <- read_csv("dataset2_topgenes.csv")
dataset2_topgenes <- dataset2_topgenes[order(-dataset2_topgenes$logFC),]
geneList2 <- dataset2_topgenes$logFC
names(geneList2) <- dataset2_topgenes$`Gene Symbol`
goo_all <- GSEA(geneList2, TERM2GENE = m_df, pvalueCutoff = 1)
goo_all <- goo_all@result
y <- arrange(goo_all, abs(NES)) %>% group_by(sign(NES))
ggplot(y, aes(NES, fct_reorder(Description, NES), fill=pvalue),showCatagory) + ggtitle("Dataset2") + theme(plot.title = element_text(hjust = 0.5))  +
  geom_barh(stat='identity') + scale_fill_continuous(low = 'red', high = 'blue') + theme_minimal() + ylab(NULL)   + theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()) + theme(axis.text = element_text(size = 15))

ora_limma_1 <- ora_limma_1[,-14]
ora_limma_2 <- ora_limma_2[,-14]
ora_gfs_1 <- ora_gfs_1[,-14]
ora_gfs_2 <- ora_gfs_2[,-14]

write.csv(ora_gfs_2, "ora_gfs_dataset2.csv")
write.csv(ora_limma_2, "ora_limma_dataset2.csv")
#ora_limma_1 <- cbind(ora_limma$query, ora_limma$significant, ora_limma$p_value, ora_limma$term_size, ora_limma$query_size, ora_limma$intersection_size, ora_limma$precision, ora_limma$recall, ora_limma$term_id, ora_limma$)
#colnames(ora_limma_1) <- c("query", "significant", "p_value", "term_size", "query_size", "intersection_size", "precision", "recall")

gfs_dataset_1_DEG_entrez <- read_csv("gfs_dataset_1_DEG_entrez (1).csv")
gfs_dataset_2_DEG_entrez_1_ <- read_csv("gfs_dataset_2_DEG_entrez (1).csv")
limma_dataset_1_DEG_entrez <- read_csv("limma_dataset_1_DEG_entrez (1).csv")
limma_dataset_2_DEG_entrez <- read_csv("limma_dataset_2_DEG_entrez (1).csv")
library(ggvenn)

vene_list <- list(
  limma_dataset1 = limma_dataset_1_DEG_entrez$ENTREZID,
  limma_dataset2 = limma_dataset_2_DEG_entrez$ENTREZID
)
intersect(ora_gfs_1$term_name, ora_gfs_2$term_name)
intersect(ora_limma_1$term_name,ora_limma_2$term_name)
ggvenn(
  vene_list, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 5
)

ggvenn(
  vene_list, 
  fill_color = c("#8491B4B2", "#DC0000B2"),
  stroke_size = 0.5, set_name_size = 5
)

rf = list(1:10)
dt = list(3)
########
vene_list <- list(
  Random_forest = 5:15, 
  SVM = 2:11, 
  Decision_tree = 8:17,
  KNN = 9:18
  
)


ggvenn(
  vene_list, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF","#FFFF00"),
  stroke_size = 0.5, set_name_size = 5
)

