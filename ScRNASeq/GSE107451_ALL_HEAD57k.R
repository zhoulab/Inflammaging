#############################################set working directory and library package################################################################################################

setwd("/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/")
#load("GSE107451_head_ALL_57k.RData")


library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(data.table)
library(tidyverse)
library(clustree)
library(openxlsx)
library(DESeq2)
library(presto)
library(EnhancedVolcano)
library(ggrepel)
library(gplots)
#devtools::install_github("immunogenomics/presto")


#############################################################################################################################################




#############################################Import data and create Seurat object################################################################################################

ALL_HEAD_57k <- Seurat::Read10X("./GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv")

ALL_HEAD_57k <- CreateSeuratObject(counts = ALL_HEAD_57k,
                               project = "ALL_HEAD_57k_GSE107451_Summer2024", 
                               min.cells = 0,  
                               min.features = 0)
#look at the meta data
View(ALL_HEAD_57k@meta.data)

#import metadata of author's and add that to Seurat OBJ
meta_author <- as.data.frame(fread("./GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv",header = T))

table(rownames(ALL_HEAD_57k@meta.data)==meta_author$new_barcode)
#TRUE 
#56902 
# 

ALL_HEAD_57k@meta.data$"new_barcode" <- rownames(ALL_HEAD_57k@meta.data)

ALL_HEAD_57k@meta.data <- cbind(ALL_HEAD_57k@meta.data,meta_author)

#look at the meta.data after adding author's meta data
View(ALL_HEAD_57k@meta.data)
#############################################################################################################################################




#####################################################This is author's 57k data, no need to filter########################################################################################
# n.Feature.min <- 200  #least gene # >200
# n.Feature.max <- 4000  #max gene # <4000
# n.Count.min <-500 #least count >500
# n.Mt <- 10  #max mt <20%
# n.Rb <- 10   #max ribosome <10% *** changeable
# 
# cat("Before filter :",nrow(ALL_HEAD_57k@meta.data),"cells\n")
# 
# 
# ALL_HEAD_57k <- subset(ALL_HEAD_57k, 
#                    subset = 
#                      nFeature_RNA > n.Feature.min &
#                      nFeature_RNA<n.Feature.max &
#                      nCount_RNA >n.Count.min&
#                      percent.mt < n.Mt&
#                      percent.rb < n.Rb)
# 
#cat("After filter :",nrow(ALL_HEAD_57k@meta.data),"cells\n")
#############################################################################################################################################




####################################################Normalization and find clusters under different resolutions#########################################################################################
# run standard anlaysis workflow (82 PCs is what autho's parameter)
ALL_HEAD_57k2 <- NormalizeData(ALL_HEAD_57k)
ALL_HEAD_57k2 <- FindVariableFeatures(ALL_HEAD_57k2)
ALL_HEAD_57k2 <- ScaleData(ALL_HEAD_57k2)
min(nrow(ALL_HEAD_57k2), ncol(ALL_HEAD_57k2))
ALL_HEAD_57k2 <- RunPCA(ALL_HEAD_57k2,npcs = 82,verbose = F)
ALL_HEAD_57k2 <- FindNeighbors(ALL_HEAD_57k2, dims = 1:82, reduction = "pca", verbose=TRUE)

#find clusters under different resolutions
for (i in seq(0.80,0.95,by=0.01)) {
  ALL_HEAD_57k2 <- FindClusters(ALL_HEAD_57k2, resolution = i, cluster.name = paste0("ALL_HEAD_57k_clusters_",i))
}

#save the resolution adn clustering results with clustree
ggsave(clustree(ALL_HEAD_57k2,prefix = "ALL_HEAD_57k_clusters_"), device = "pdf",filename = "./cluster_tree.pdf",width = 30,height = 20)

#change active.ident to 87 clusters for next step (using author's annotation)
ALL_HEAD_57k2@active.ident <- as.factor(ALL_HEAD_57k2@meta.data$res.2)
#############################################################################################################################################





#############################################################################################################################################
#############################################################################################################################################
######################################################cell level analysis#######################################################################################
#############################################################################################################################################
#############################################################################################################################################


#######################################################Run tSNE reductons######################################################################################
#Run tSNE
ALL_HEAD_57k3 <- RunTSNE(ALL_HEAD_57k2, dims = 1:82, reduction = "pca", reduction.name = "tSNE_ALL_HEAD_57k_Author_cluster",verbose=F,perplexity=30)


#change tSNE coordinates to authors' values (so the topology looks if not the same, very similar)
table(rownames(ALL_HEAD_57k3@meta.data)%in%rownames(ALL_HEAD_57k3@reductions[["tSNE_ALL_HEAD_57k_Author_cluster"]]@cell.embeddings))
# TRUE 
# 56902 
ALL_HEAD_57k3@reductions[["tSNE_ALL_HEAD_57k_Author_cluster"]]@cell.embeddings[,1] <- ALL_HEAD_57k3@meta.data$seurat_tsne1
ALL_HEAD_57k3@reductions[["tSNE_ALL_HEAD_57k_Author_cluster"]]@cell.embeddings[,2] <- ALL_HEAD_57k3@meta.data$seurat_tsne2



#tSNE plot of cell clusters based on author's annotation
pdf(paste0("./ALL_HEAD_57k_tSNE.pdf"),width = 30,height = 20)
DimPlot(ALL_HEAD_57k3, reduction = "tSNE_ALL_HEAD_57k_Author_cluster",group.by = "annotation",
        label = TRUE,
        label.size = 5,
        repel = TRUE)
dev.off()

#tSNE plot of different ages based on author's annotation
pdf(paste0("./ALL_HEAD_57k_tSNE_days.pdf"),width = 15,height = 10)
DimPlot(ALL_HEAD_57k3, reduction = "tSNE_ALL_HEAD_57k_Author_cluster", group.by = "Age")
dev.off()
#############################################################################################################################################





#########################################################No need to find markers, there is author's annotation####################################################################################
#find markers

# 
# all.markers <- FindAllMarkers(ALL_HEAD_57k3, test.use = 'wilcox', only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25,group.by ="ALL_HEAD_57k_clusters_0.86")
# 
# marker.sig <- all.markers %>%
#   filter(p_val_adj < 0.05)
# 
# topmarker <- marker.sig %>% group_by(cluster) %>% 
#   top_n(n = 2, wt = avg_log2FC)
# 
# 
# 
# 
# 
# 
# pdf(paste0("./ALL_HEAD_57k_DotPlot.pdf"),width = 30,height = 20)
# DotPlot(ALL_HEAD_57k3, features = unique(marker.sig$gene),group.by = "ALL_HEAD_57k_clusters_0.86")+
#   theme(axis.text.x = element_text(angle = 90, size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# dev.off()
# 
# 
# 
# pdf(paste0("./ALL_HEAD_57k_Heatmap.pdf"),width = 30,height = 20)
# DoHeatmap(ALL_HEAD_57k3, features = topmarker$gene,size = 2,group.by = "ALL_HEAD_57k_clusters_0.86") +
#   theme(legend.position = "none", 
#         axis.text.y = element_text(size = 6))
# dev.off()
# 
# 
# all_markers_dotplot <- DotPlot(ALL_HEAD_57k3, features = unique(marker.sig$gene))+
#   theme(axis.text.x = element_text(angle = 90, size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# View(all_markers_dotplot$data)
# 
#############################################################################################################################################






#################################################Creat a Gene List of inflammaging markers and corresponding feature plots############################################################################################
#creating inflammaging marker gene list
III_Marker_Genes <- c("SPH93",
                      "CecA1",
                      "CecA2",
                      "AttA",
                      "AttB",
                      "AttC",
                      "DptA",
                      "DptB",
                      "Dro",
                      "Mtk",
                      #"IM18",
                      "edin")

#This is an easy loop to save feature plots of all the inflammaging markers
for (i in 1:length(III_Marker_Genes)) {
  tmp <- Seurat::FeaturePlot(ALL_HEAD_57k3, features = III_Marker_Genes[i],  cols = c("grey","darkred"), reduction = "tSNE_ALL_HEAD_57k_Author_cluster",keep.scale = "feature")
  ggsave(tmp,device = "pdf",file= paste0("./Featureplot/Featureplot_HEAD57k_",III_Marker_Genes[i]),width = 15,height = 10)
}

#these two additional feature plots of noe and pros are for comparision of plot colors when there are many cells expressing certain genes
ggsave(Seurat::FeaturePlot(ALL_HEAD_57k3, features = "noe",  cols = c("grey","darkred"), reduction = "tSNE_ALL_HEAD_57k_Author_cluster",keep.scale = "feature"),
       device = "pdf",file= paste0("./Featureplot/Featureplot_HEAD57k_","noe"),width = 15,height = 10)
ggsave(Seurat::FeaturePlot(ALL_HEAD_57k3, features = "pros",  cols = c("grey","darkred"), reduction = "tSNE_ALL_HEAD_57k_Author_cluster",keep.scale = "feature"),
       device = "pdf",file= paste0("./Featureplot/Featureplot_HEAD57k_","pros"),width = 15,height = 10)
#############################################################################################################################################



########################################################violin plot of Inflammaging markers####################################################################################

#This is an easy loop to save voilin plots of all the inflammaging markers based on Genotype
for (i in 1:length(III_Marker_Genes)) {
  tmp2 <- Seurat::VlnPlot(ALL_HEAD_57k3, features = III_Marker_Genes[i],group.by = "Genotype")
  ggsave(tmp2,device = "pdf",file= paste0("./Vlnplot/Vlnplot_HEAD57k_Genotype_",III_Marker_Genes[i]),width = 15,height = 10)
}

#This vlnplot of noe is used to demonstrate what vlnplot looks like when there is enough cells expressing a gene
ggsave(Seurat::VlnPlot(ALL_HEAD_57k3, features = "noe",group.by = "Genotype")
       ,device = "pdf",file= paste0("./Vlnplot/Vlnplot_HEAD57k_Genotype_","noe"),width = 15,height = 10)

#save voilin plots of all the inflammaging markers based on days
for (i in 1:length(III_Marker_Genes)) {
  tmp3 <- Seurat::VlnPlot(ALL_HEAD_57k3, features = III_Marker_Genes[i],group.by = "Age")
  ggsave(tmp3,device = "pdf",file= paste0("./Vlnplot/Vlnplot_HEAD57k_Age_",III_Marker_Genes[i]),width = 15,height = 10)
}
#This vlnplot of noe is used to demonstrate what vlnplot looks like when there is enough cells expressing a gene
ggsave(Seurat::VlnPlot(ALL_HEAD_57k3, features = "noe",group.by = "Age")
       ,device = "pdf",file= paste0("./Vlnplot/Vlnplot_HEAD57k_Age_","noe"),width = 15,height = 10)

#############################################################################################################################################





##################################################make patial feature plot of III markers in plasmatocyte###########################################################################################

#tSNE of the plasmatocytes cluster
pdf(file = "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Featureplot/partial_Plasmatocytes/Plasmatocytes_tSNE.pdf",width = 15,height = 10)
DimPlot(subset(ALL_HEAD_57k3,annotation=="Plasmatocytes"),
        reduction = "tSNE_ALL_HEAD_57k_Author_cluster",
        group.by = "annotation",
        #label = TRUE,
        #label.size = 2,
        repel = TRUE)
dev.off()


#partial feature plot of III markers in plasmatocytes cluster
for (i in 1:length(III_Marker_Genes)) {
  tmp <- Seurat::FeaturePlot(subset(ALL_HEAD_57k3,annotation=="Plasmatocytes"),
                             features = III_Marker_Genes[i], 
                             cols = c("grey","darkred"),
                             reduction = "tSNE_ALL_HEAD_57k_Author_cluster",
                             keep.scale = "feature",
                             pt.size = 1
                             )
  ggsave(tmp,
         device = "pdf",
         file= paste0(
           "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Featureplot/partial_Plasmatocytes/",
           III_Marker_Genes[i],
           "_Plasmatocytes.pdf"),
         width = 15,
         height = 10
  )
  }


#partial feature plot of III markers in plasmatocytes cluster and split by Age
for (i in 1:length(III_Marker_Genes)) {
  tmp <- Seurat::FeaturePlot(subset(ALL_HEAD_57k3,annotation=="Plasmatocytes"),
                             features = III_Marker_Genes[i],
                             cols = c("grey","darkred"),
                             reduction = "tSNE_ALL_HEAD_57k_Author_cluster",
                             keep.scale = "all",
                             split.by = "Age",
                             pt.size = 1
  )+
    patchwork::plot_layout(ncol = 3, nrow = 3)& theme(legend.position = "right")
  
  ggsave(tmp,
         device = "pdf",
         file= paste0(
           "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Featureplot/partial_Plasmatocytes/",
           III_Marker_Genes[i],
           "_Plasmatocytes_Age.pdf"),
         width = 15,
         height = 10,
         limitsize = FALSE
  )
}
#############################################################################################################################################








################################################***test, find all markers, lowest standard#############################################################################################
#change active idents in ALL_HEAD_57k3
Idents(ALL_HEAD_57k3) <- ALL_HEAD_57k3@meta.data$annotation
#find marker genes with the lowest standard (so it will definitely include III markers)
lowest_standard_markers <-
  FindAllMarkers(ALL_HEAD_57k3,
                 test.use = 'wilcox',
                 only.pos = F,
                 min.pct = 0,
                 logfc.threshold = 0,
                 group.by ="annotation")
#dotplot of III markers
lowest_standard_markers_dotplot <- DotPlot(ALL_HEAD_57k3, features = III_Marker_Genes)+
    theme(axis.text.x = element_text(angle = 90, size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12))+
  coord_flip()
View(lowest_standard_markers_dotplot$data)

#export dotplot of III markers
pdf("/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Dotplot/III_Markers_dotplot.pdf",
    width = 40,
    height = 15)
lowest_standard_markers_dotplot
dev.off()

#export dotplot data () of III markers
openxlsx::write.xlsx(lowest_standard_markers_dotplot$data[lowest_standard_markers_dotplot$data$features.plot%in%III_Marker_Genes,],
                     "lowest_standard_markers_III_markers.xlsx",
                     rowNames=F)

### Try another method of finding percentage of gene expression
#find the percentage of expression of III markers in cell clusters.
#https://samuel-marsh.github.io/scCustomize/reference/Percent_Expressing.html
#https://samuel-marsh.github.io/scCustomize/articles/Statistics.html
percent_EXPR_III_markers <- scCustomize::Percent_Expressing(seurat_object = ALL_HEAD_57k3,
                                                            features = III_Marker_Genes,
                                                            group_by = "annotation")
View(t(percent_EXPR_III_markers))
###
#############################################################################################################################################








#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################












#############################################################################################################################################
#############################################################################################################################################
#############################################################Pseudo_bulk level analysis################################################################################
#############################################################################################################################################
#############################################################################################################################################












####################################################Make Psuedo bulk count to compare with our BulkRNAseq***#########################################################################################
#***Comment: In Seurat v5, we encourage the use of the AggregateExpression function to perform pseudobulk analysis.
#https://satijalab.org/seurat/articles/announcements.html


#extract Psudo count based on Genotype, sex, and Age then make subset of it tom compare with our BulkRNAseq
EXPR_Genotype_sex_Age <- 
  as.data.frame(AggregateExpression(object = ALL_HEAD_57k3, group.by = c("Genotype","sex","Age"))$RNA)

#DGRP-551_female
count_DGRP551_female <- 
  EXPR_Genotype_sex_Age[,
                        intersect(
                          which(str_split(colnames(EXPR_Genotype_sex_Age),pattern = "_",3,simplify =T)[,1]=="DGRP-551"),
                          which(str_split(colnames(EXPR_Genotype_sex_Age),pattern = "_",3,simplify =T)[,2]=="female")
                          )
                        ]
#DGRP-551_male
count_DGRP551_male <- 
  EXPR_Genotype_sex_Age[,
                        intersect(
                          which(str_split(colnames(EXPR_Genotype_sex_Age),pattern = "_",3,simplify =T)[,1]=="DGRP-551"),
                          which(str_split(colnames(EXPR_Genotype_sex_Age),pattern = "_",3,simplify =T)[,2]=="male")
                        )
                        ]
#w1118_female
count_w1118_female <- 
  EXPR_Genotype_sex_Age[,
                        intersect(
                          which(str_split(colnames(EXPR_Genotype_sex_Age),pattern = "_",3,simplify =T)[,1]=="w1118"),
                          which(str_split(colnames(EXPR_Genotype_sex_Age),pattern = "_",3,simplify =T)[,2]=="female")
                        )
                        ]
#w1118_male
count_w1118_male <- 
  EXPR_Genotype_sex_Age[,
                        intersect(
                          which(str_split(colnames(EXPR_Genotype_sex_Age),pattern = "_",3,simplify =T)[,1]=="w1118"),
                          which(str_split(colnames(EXPR_Genotype_sex_Age),pattern = "_",3,simplify =T)[,2]=="male")
                        )
                        ]

#export these psudo count for further differential expression analysis to compare with our Bulk_RNA data
write.xlsx(count_DGRP551_female,"count_DGRP551_female_HEAD57k.xlsx",rowNames=T)
write.xlsx(count_DGRP551_male,"count_DGRP551_male_HEAD57k.xlsx",rowNames=T)
write.xlsx(count_w1118_female,"count_w1118_female_HEAD57k.xlsx",rowNames=T)
write.xlsx(count_w1118_male,"count_w1118_male_HEAD57k.xlsx",rowNames=T)
#############################################################################################################################################




####################################################heatmap of pseudo bulk expression#########################################################################################

make_heatmap <- function(EXPR_OBJ,GENELIST,TITLE) {
  tmp_hm <- EXPR_OBJ[rownames(EXPR_OBJ)%in%GENELIST,];
  tmp_hm <- as.matrix(tmp_hm);
  heatmap.2(x=tmp_hm,
            Colv=F,
            Rowv = F,
            dendrogram="none",
            scale="row",
            col="bluered",
            main=TITLE,
            margins = c(15,15),
            cexRow = 1.0,
            cexCol = 1.0,
            trace = "none"
            )
  }

#count_DGRP551_female
pdf(file = "./Heatmap/DGRP551_female_heatmap.pdf",
       width = 15,
       height = 10)
make_heatmap(EXPR_OBJ = count_DGRP551_female[,-c(1:2)],GENELIST = III_Marker_Genes,TITLE = "DGRP551_female")
dev.off()

#count_DGRP551_male
pdf(file = "./Heatmap/DGRP551_male_heatmap.pdf",
    width = 15,
    height = 10)
make_heatmap(EXPR_OBJ = count_DGRP551_male[,-c(1:2)],GENELIST = III_Marker_Genes,TITLE = "DGRP551_male")
dev.off()

#count_w1118_female
pdf(file = "./Heatmap/w1118_female_heatmap.pdf",
    width = 15,
    height = 10)
make_heatmap(EXPR_OBJ = count_w1118_female[,-c(1:2)],GENELIST = III_Marker_Genes,TITLE = "w1118_female")
dev.off()

#count_w1118_male
pdf(file = "./Heatmap/w1118_male_heatmap.pdf",
    width = 15,
    height = 10)
make_heatmap(EXPR_OBJ = count_w1118_male[,-c(1:2)],GENELIST = III_Marker_Genes,TITLE = "w1118_male")
dev.off()
#############################################################################################################################################






#####################################################Perform DEseq on Pseudo count to compare with our data########################################################################################

#DEseq function

Do_DEseq <- 
  function(Counts_matrix,
           n_TR,                # 
           n_CR,                #
           p_thresh,            #
           log2FC_thresh,       #
           File_Title){        #
    # make meta for DEseq (make a variable that tells DEseq which columns are controls and which are treatment ones)
    
    meta <- colnames(Counts_matrix)
    
    meta <- as.data.frame(meta)
    
    rownames(meta) <- meta[,1]
    
    meta$"Condition" <- c(rep("TR",n_TR),rep("CR",n_CR))
    
    
    meta$Condition <- factor(meta$Condition,levels = c("CR","TR"))
    
    # make the dds object for DEseq
    
    dds <- DESeqDataSetFromMatrix(countData = Counts_matrix,
                                  colData = meta,
                                  design= ~ Condition)
    dds <- DESeq(dds)
    res <- results(dds, pAdjustMethod = "BH")
    
    #extracting result
    res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
    
    res1 <- na.omit(res1)
    res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
    
    
    res1[which(res1$log2FoldChange >= log2FC_thresh & res1$pvalue < p_thresh),'sig'] <- 'up'
    res1[which(res1$log2FoldChange <= -log2FC_thresh & res1$pvalue < p_thresh),'sig'] <- 'down'
    res1[which(abs(res1$log2FoldChange) <= log2FC_thresh | res1$pvalue >= p_thresh),'sig'] <- 'none'
    
    
    write.xlsx(res1, paste0(File_Title,".xlsx"),rowNames=T)
    
    write.xlsx(res1[res1$'sig'=='up',],paste0(File_Title,"_UP",p_thresh,".xlsx"),rowNames=T)
    
    write.xlsx(res1[res1$'sig'=='down',],paste0(File_Title,"_DOWN",p_thresh,".xlsx"),rowNames=T)
  }



##DGRP-551_female
colnames(count_DGRP551_female)
#30vs1
Do_DEseq(Counts_matrix=count_DGRP551_female[,c(7,2)], 
         #select the columns in readcount matrix and select treatment sample first then control samples
         n_TR=1,   #number of treatment sample             
         n_CR=1,    #number of control sample            
         p_thresh=0.01,         #p value threshold   
         log2FC_thresh=1,       #log2 fold change threshold
         File_Title="DGRP551_female_30vs1")  #xlsx file title
#50vs1
Do_DEseq(Counts_matrix=count_DGRP551_female[,c(8,2)], 
         #select the columns in readcount matrix and select treatment sample first then control samples
         n_TR=1,   #number of treatment sample             
         n_CR=1,    #number of control sample            
         p_thresh=0.01,         #p value threshold   
         log2FC_thresh=1,       #log2 fold change threshold
         File_Title="DGRP551_female_50vs1")  #xlsx file title

##DGRP-551_male
#30vs1

#50vs1

##w1118_female
#30vs1

#50vs1

##w1118_male
#30vs1

#50vs1



#!!! DEseq cannot compare only two samples (1 replicate), have to use another tool for differential expression


#############################################################################################################################################





#####################################################Differential expression with presto package all comparison########################################################################################
#https://satijalab.org/seurat/articles/announcements.html (recommended by Seurat team)
#https://github.com/immunogenomics/presto

#for example: wilcoxauc(ALL_HEAD_57k3,"annotation") gives differentially expressed genes in each clusters

#Find the DEGs at different age
##DGRP-551_female
DE_DGRP551_female <- wilcoxauc(
  subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="female"),
  group_by ="Age")

##DGRP-551_male
DE_DGRP551_male <- wilcoxauc(
  subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="male"),
  group_by ="Age")

##w1118_female
DE_w1118_female <- wilcoxauc(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="female"),
  group_by ="Age")

##w1118_male
DE_w1118_male <- wilcoxauc(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="male"),
  group_by ="Age")
#############################################################################################################################################








########################################################Differential expression with presto package, 1vs1 comparison#####################################################################################
#This area of code demonstrates the case where Sometimes, you don't want to test all groups in the dataset against all other groups. 
#For instance, I want to compare only observations in group A to those in group B. This is achieved with the groups_use argument.


#### These two code chunks are different in terms of the order of Age put into the function
View(
  wilcoxauc(
  subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="female"),
  group_by = "Age",
  groups_use = c("50","3"))
  )


View(
  wilcoxauc(
    subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="female"),
    group_by = "Age",
    groups_use = c("3","50"))
)
#### The order of age put into function does not influence comparison result.



#find the DEGs by comparing one age to another age

#DGRP551_female day50 vs day3 and day30 vs day3
DGRP551_female_50vs3 <- 
  wilcoxauc(
  subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="female"),
  group_by = "Age",
  groups_use = c("50","3"))
DGRP551_female_30vs3 <- 
  wilcoxauc(
  subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="female"),
  group_by = "Age",
  groups_use = c("30","3"))

#DGRP551_male day50 vs day3 and day30 vs day3
DGRP551_male_50vs3 <- 
  wilcoxauc(
    subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="male"),
    group_by = "Age",
    groups_use = c("50","3"))
DGRP551_male_30vs3 <- 
  wilcoxauc(
    subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="male"),
    group_by = "Age",
    groups_use = c("30","3"))

#w1118_female day30 vs day3
w1118_female_30vs3 <- 
  wilcoxauc(
    subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="female"),
    group_by = "Age",
    groups_use = c("30","3"))

#w1118_male day30 vs day3
w1118_male_30vs3 <- 
  wilcoxauc(
    subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="male"),
    group_by = "Age",
    groups_use = c("30","3"))
#############################################################################################################################################





############################################################visualization of presto Differential expression result#################################################################################

####The first part is volcano plot of all comparisons(i.e. not one Age vs another Age; it's a holistic comparison)
#DE_DGRP551_female
ggsave(EnhancedVolcano(DE_DGRP551_female,
                lab = paste0(DE_DGRP551_female$feature,"_day",DE_DGRP551_female$group),
                x = 'logFC',
                y = 'padj',
                title = "DE_DGRP551_female",
                pCutoff = 0.01,
                FCcutoff = 0.6,
                pointSize = 3.0,
                labSize = 3.0,
                max.overlaps = 100,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = T,
                widthConnectors = 0.3,
                lengthConnectors = unit(0.01, "npc"),
                directionConnectors = "both",
                arrowheads = F,
                maxoverlapsConnectors =50,
                typeConnectors = "open",
                endsConnectors = "first"),
       device = "pdf",
       filename = "./DE_VolcanoPlot/DE_DGRP551_female.pdf",
       width = 15,
       height = 10)
#DE_DGRP551_male
ggsave(EnhancedVolcano(DE_DGRP551_male,
                       lab = paste0(DE_DGRP551_male$feature,"_day",DE_DGRP551_male$group),
                       x = 'logFC',
                       y = 'padj',
                       title = "DE_DGRP551_male",
                       pCutoff = 0.01,
                       FCcutoff = 0.6,
                       pointSize = 3.0,
                       labSize = 3.0,
                       max.overlaps = 100,
                       colAlpha = 1,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = T,
                       widthConnectors = 0.3,
                       lengthConnectors = unit(0.01, "npc"),
                       directionConnectors = "both",
                       arrowheads = F,
                       maxoverlapsConnectors =50,
                       typeConnectors = "open",
                       endsConnectors = "first"),
       device = "pdf",
       filename = "./DE_VolcanoPlot/DE_DGRP551_male.pdf",
       width = 15,
       height = 10)
#DE_w1118_female
ggsave(EnhancedVolcano(DE_w1118_female,
                       lab = paste0(DE_w1118_female$feature,"_day",DE_w1118_female$group),
                       x = 'logFC',
                       y = 'padj',
                       title = "DE_w1118_female",
                       pCutoff = 0.01,
                       FCcutoff = 0.6,
                       pointSize = 3.0,
                       labSize = 3.0,
                       max.overlaps = 100,
                       colAlpha = 1,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = T,
                       widthConnectors = 0.3,
                       lengthConnectors = unit(0.01, "npc"),
                       directionConnectors = "both",
                       arrowheads = F,
                       maxoverlapsConnectors =50,
                       typeConnectors = "open",
                       endsConnectors = "first"),
       device = "pdf",
       filename = "./DE_VolcanoPlot/DE_w1118_female.pdf",
       width = 15,
       height = 10)
#DE_w1118_male
ggsave(EnhancedVolcano(DE_w1118_male,
                       lab = paste0(DE_w1118_male$feature,"_day",DE_w1118_male$group),
                       x = 'logFC',
                       y = 'padj',
                       title = "DE_w1118_male",
                       pCutoff = 0.01,
                       FCcutoff = 0.6,
                       pointSize = 3.0,
                       labSize = 3.0,
                       max.overlaps = 100,
                       colAlpha = 1,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = T,
                       widthConnectors = 0.3,
                       lengthConnectors = unit(0.01, "npc"),
                       directionConnectors = "both",
                       arrowheads = F,
                       maxoverlapsConnectors =50,
                       typeConnectors = "open",
                       endsConnectors = "first"),
       device = "pdf",
       filename = "./DE_VolcanoPlot/DE_w1118_male.pdf",
       width = 15,
       height = 10)
####




#####The second part is the volcano plot of the second type of comparison (i.e. one Age vs another Age)

#DGRP551_female_50vs3
ggsave(EnhancedVolcano(subset(DGRP551_female_50vs3,group=="50"),
                       lab = subset(DGRP551_female_50vs3,group=="50")$feature,
                       x = 'logFC',
                       y = 'padj',
                       title = "DGRP551_female_50vs3",
                       pCutoff = 0.01,
                       FCcutoff = 0.5,
                       pointSize = 4.0,
                       labSize = 3.0,
                       max.overlaps = 100,
                       colAlpha = 1,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = T,
                       widthConnectors = 0.3,
                       lengthConnectors = unit(0.01, "npc"),
                       directionConnectors = "both",
                       arrowheads = F,
                       maxoverlapsConnectors =50,
                       typeConnectors = "open",
                       endsConnectors = "first"),
       device = "pdf",
       filename = "./DE_VolcanoPlot/DGRP551_female_50vs3.pdf",
       width = 15,
       height = 10)

#DGRP551_female_30vs3
ggsave(EnhancedVolcano(subset(DGRP551_female_30vs3,group=="30"),
                       lab = subset(DGRP551_female_30vs3,group=="30")$feature,
                       x = 'logFC',
                       y = 'padj',
                       title = "DGRP551_female_30vs3",
                       pCutoff = 0.01,
                       FCcutoff = 0.5,
                       pointSize = 4.0,
                       labSize = 3.0,
                       max.overlaps = 100,
                       colAlpha = 1,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = T,
                       widthConnectors = 0.3,
                       lengthConnectors = unit(0.01, "npc"),
                       directionConnectors = "both",
                       arrowheads = F,
                       maxoverlapsConnectors =50,
                       typeConnectors = "open",
                       endsConnectors = "first"),
       device = "pdf",
       filename = "./DE_VolcanoPlot/DGRP551_female_30vs3.pdf",
       width = 15,
       height = 10)

#DGRP551_male_50vs3
ggsave(EnhancedVolcano(subset(DGRP551_male_50vs3,group=="50"),
                       lab = subset(DGRP551_male_50vs3,group=="50")$feature,
                       x = 'logFC',
                       y = 'padj',
                       title = "DGRP551_male_50vs3",
                       pCutoff = 0.01,
                       FCcutoff = 0.5,
                       pointSize = 4.0,
                       labSize = 3.0,
                       max.overlaps = 100,
                       colAlpha = 1,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = T,
                       widthConnectors = 0.3,
                       lengthConnectors = unit(0.01, "npc"),
                       directionConnectors = "both",
                       arrowheads = F,
                       maxoverlapsConnectors =50,
                       typeConnectors = "open",
                       endsConnectors = "first"),
       device = "pdf",
       filename = "./DE_VolcanoPlot/DGRP551_male_50vs3.pdf",
       width = 15,
       height = 10)

#DGRP551_male_30vs3
ggsave(EnhancedVolcano(subset(DGRP551_male_30vs3,group=="30"),
                       lab = subset(DGRP551_male_30vs3,group=="30")$feature,
                       x = 'logFC',
                       y = 'padj',
                       title = "DGRP551_male_30vs3",
                       pCutoff = 0.01,
                       FCcutoff = 0.5,
                       pointSize = 4.0,
                       labSize = 3.0,
                       max.overlaps = 100,
                       colAlpha = 1,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = T,
                       widthConnectors = 0.3,
                       lengthConnectors = unit(0.01, "npc"),
                       directionConnectors = "both",
                       arrowheads = F,
                       maxoverlapsConnectors =50,
                       typeConnectors = "open",
                       endsConnectors = "first"),
       device = "pdf",
       filename = "./DE_VolcanoPlot/DGRP551_male_30vs3.pdf",
       width = 15,
       height = 10)

#w1118_female_30vs3
ggsave(EnhancedVolcano(subset(w1118_female_30vs3,group=="30"),
                       lab = subset(w1118_female_30vs3,group=="30")$feature,
                       x = 'logFC',
                       y = 'padj',
                       title = "w1118_female_30vs3",
                       pCutoff = 0.01,
                       FCcutoff = 0.25,
                       pointSize = 4.0,
                       labSize = 3.0,
                       max.overlaps = 100,
                       colAlpha = 1,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = T,
                       widthConnectors = 0.3,
                       lengthConnectors = unit(0.01, "npc"),
                       directionConnectors = "both",
                       arrowheads = F,
                       maxoverlapsConnectors =50,
                       typeConnectors = "open",
                       endsConnectors = "first"),
       device = "pdf",
       filename = "./DE_VolcanoPlot/w1118_female_30vs3.pdf",
       width = 15,
       height = 10)

#w1118_male_30vs3
ggsave(EnhancedVolcano(subset(w1118_male_30vs3,group=="30"),
                       lab = subset(w1118_male_30vs3,group=="30")$feature,
                       x = 'logFC',
                       y = 'padj',
                       title = "w1118_male_30vs3",
                       pCutoff = 0.01,
                       FCcutoff = 0.25,
                       pointSize = 4.0,
                       labSize = 3.0,
                       max.overlaps = 100,
                       colAlpha = 1,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = T,
                       widthConnectors = 0.3,
                       lengthConnectors = unit(0.01, "npc"),
                       directionConnectors = "both",
                       arrowheads = F,
                       maxoverlapsConnectors =50,
                       typeConnectors = "open",
                       endsConnectors = "first"),
       device = "pdf",
       filename = "./DE_VolcanoPlot/w1118_male_30vs3.pdf",
       width = 15,
       height = 10)

#############################################################################################################################################










#####################################################Find which cell type III Markers increased significantly in########################################################################################
#based on presto package that provides statistical significance on a holistic level of comparison for all cell clusters

DE_Cell_Cluster <- wilcoxauc(ALL_HEAD_57k3,group_by = "annotation")

#look at the result for III markers
View(subset(DE_Cell_Cluster,DE_Cell_Cluster$feature%in%III_Marker_Genes))
#export the statistical result fo III markers for further analysis
openxlsx::write.xlsx(subset(DE_Cell_Cluster,DE_Cell_Cluster$feature%in%III_Marker_Genes),
                     "DE_cell_cluster_result_III_Markers.xlsx",
                     rowNames=F)
#############################################################################################################################################




############################################generate the pseudo-bulk expression of heat map of III markers in cell types#################################################################################################
#obtain the pseudo-bulk expression of all genes in different cell clusters
count_celltype <- as.data.frame(
  AggregateExpression(ALL_HEAD_57k3,group.by = "annotation")$"RNA"
  )
#check how many cell types there are
ncol(count_celltype)
#[1] 116

#make a heat map of the pseudo-bulk expression of III markers in different cell types.
pdf(file = "./Heatmap/III_Markers_Expression_cell_types.pdf",
    width = 30,
    height = 15)
make_heatmap(EXPR_OBJ = count_celltype, GENELIST = III_Marker_Genes, TITLE = "III Markers Expression in different cell types")
dev.off()

#make a heat map that include noe and pros for to debug
pdf(file = "./Heatmap/III_Markers_Expression_cell_types_test.pdf",
    width = 50,
    height = 50)
make_heatmap(EXPR_OBJ = count_celltype, GENELIST = c("pros","noe"), TITLE = "III Markers Expression in different cell types")
dev.off()

#See the data for making the above heat map and determine what is wrong that casued no blue in heatmap
View(count_celltype[rownames(count_celltype)%in%append(III_Marker_Genes,c("noe","pros")),])
#too many zero values


#############################################################################################################################################









######################################Obtain Pseudo-bulk expression based on Strain, Sex, Cell cluster, and Age#################################################################################
#extract the Pseudo-bulk expression based on Strain, Sex, Cell cluster, and Age
count_Strain_Sex_Cellcluster_Age <- as.data.frame(
  AggregateExpression(ALL_HEAD_57k3,
                      group.by = c("Genotype",
                                   "sex",
                                   "annotation",
                                   "Age")
                      )$"RNA"
  )

#check if str_split is successful by looking at how many Ages there are.
table(str_split(colnames(count_Strain_Sex_Cellcluster_Age),pattern = "_",4,simplify =T)[,4])
#remove day 0 and day 1 from the count_Strain_Sex_Cellcluster_Age
count_Strain_Sex_Cellcluster_Age <- 
  count_Strain_Sex_Cellcluster_Age[,str_split(colnames(count_Strain_Sex_Cellcluster_Age),pattern = "_",4,simplify =T)[,4]!="0"&
                                   str_split(colnames(count_Strain_Sex_Cellcluster_Age),pattern = "_",4,simplify =T)[,4]!="1"]
#check if there are still day 0 and day 1.
table(str_split(colnames(count_Strain_Sex_Cellcluster_Age),pattern = "_",4,simplify =T)[,4])
#There is no day 0 and day 1 now
#############################################################################################################################################




###########################make a function to plot heat map of genes at different Age in a specific cell cluster####################################################################################

heatmap_specific_cluster <- function(EXPR_COUNT,GENELIST,STRAIN,SEX,CLUSTER_NAME,TITLE) {
  tmp_count <- EXPR_COUNT[,str_split(colnames(EXPR_COUNT),pattern = "_",4,simplify =T)[,1]==STRAIN];
  tmp_count <- tmp_count[,str_split(colnames(tmp_count),pattern = "_",4,simplify =T)[,2]==SEX];
  tmp_count <- tmp_count[,str_split(colnames(tmp_count),pattern = "_",4,simplify =T)[,3]==CLUSTER_NAME];
  tmp_count2 <-tmp_count[rownames(tmp_count)%in%GENELIST,];
  
  tmp_count2 <- as.matrix(tmp_count2);
  heatmap.2(x=tmp_count2,
            Colv=F,
            Rowv = F,
            dendrogram="none",
            scale="row",
            col="bluered",
            main=TITLE,
            margins = c(20,20),
            cexRow = 1.0,
            cexCol = 1.0,
            trace = "none"
  )
}

#EXAMPLE of this function by ploting expression of III markers in plasmotocytes where most III markers increased 
#on a holistic level of all clusters

###
#Plasmatocytes_DGRP-551_female
pdf("./Heatmap/Plasmatocytes_DGRP-551_female.pdf",width = 15,height = 10)
heatmap_specific_cluster(EXPR_COUNT=count_Strain_Sex_Cellcluster_Age,
                         GENELIST=III_Marker_Genes,
                         STRAIN = "DGRP-551",
                         SEX = "female",
                         CLUSTER_NAME="Plasmatocytes",
                         TITLE = "III_Markers in Plasmatocytes DGRP-551 female")
dev.off()

#Plasmatocytes_DGRP-551_male
pdf("./Heatmap/Plasmatocytes_DGRP-551_male.pdf",width = 15,height = 10)
heatmap_specific_cluster(EXPR_COUNT=count_Strain_Sex_Cellcluster_Age,
                         GENELIST=III_Marker_Genes,
                         STRAIN = "DGRP-551",
                         SEX = "male",
                         CLUSTER_NAME="Plasmatocytes",
                         TITLE = "III_Markers in Plasmatocytes DGRP-551 male")
dev.off()

#Plasmatocytes_w1118_female
pdf("./Heatmap/Plasmatocytes_w1118_female.pdf",width = 15,height = 10)
heatmap_specific_cluster(EXPR_COUNT=count_Strain_Sex_Cellcluster_Age,
                         GENELIST=III_Marker_Genes,
                         STRAIN = "w1118",
                         SEX = "female",
                         CLUSTER_NAME="Plasmatocytes",
                         TITLE = "III_Markers in Plasmatocytes w1118 female")
dev.off()

#Plasmatocytes_w1118_male
pdf("./Heatmap/Plasmatocytes_w1118_male.pdf",width = 15,height = 10)
heatmap_specific_cluster(EXPR_COUNT=count_Strain_Sex_Cellcluster_Age,
                         GENELIST=III_Marker_Genes,
                         STRAIN = "w1118",
                         SEX = "male",
                         CLUSTER_NAME="Plasmatocytes",
                         TITLE = "III_Markers in Plasmatocytes w1118 male")
dev.off()
###



####make similar heat map in other cluster (if needed), TBA



####

#############################################################################################################################################






#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

save.image("GSE107451_head_ALL_57k.RData")
