## this is the R script for the analysis of the head ScRNA-Seq data of the GSE107451 


## Objectives: 
### 1. Monitor the expression pattern of the inflammaging markers across different cell types (clusters) and ages.
### 2. Identify which cell type(s) account for the inflammaging pattern observed in whole head sampesl.
#### 2.1 Calculate aggregated expression for individual cluster - function to extract expression of a gene in a particular cluster.
#### 2.2 Identify which cell cluster demonstrate inflammaging pattern.
#### 2.3 Evaluate how each cell type contributes to the inflammaging pattern observed in whole head samples.
### 3. Identify if the percentage of cells expressing a certain III marker is increased across ages
### 4. Identify if there is a change in the percentage of plasmatocytes in different strains across ages


## Dataset
### GSE107451 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107451 
### we download the GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv.tar.gz 
### we used Author's annotation as metadata of this analysis: GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv.gz
### which was based on analysis down by the authors of the paper 
### A Single-Cell Transcriptome Atlas of the Aging Drosophila Brain. 
###  - Cell 2018 Aug 9;174(4):982-998.e20. PMID: 29909982
###  - https://pubmed.ncbi.nlm.nih.gov/29909982/
### For analysis of scRNAseq data, Seurat V5 is used. For more information, check https://satijalab.org/seurat/, and https://satijalab.org/seurat/articles/install.html for specific instruction on installation.


## Functions written, primarily used for analyzing this dataset
###count_cells_subset
####a function to see number of cells in different Strains, Sex at different Ages
#### Example usage:
#### To get the number of cells for a specific strain, sex, age, and cell cluster:
count_cells_subset(seurat_obj=ALL_HEAD_57k3,
                                  strain = "DGRP-551",
                                  sex = "female",
                                  age = "50",
                                  cell_cluster = "Plasmatocytes")
###calculate_gene_expression_percentage
#### calculate the persentage of cells expressing different genes within a subset of a seurat object of a certain Strain, Sex, Cluster, and Age 
#### Example usage:
calculate_gene_expression_percentage(
  seurat_obj=ALL_HEAD_57k3,
  genes=III_Marker_Genes, 
  strain = "DGRP-551", 
  Sex = "female", 
  cluster = "Plasmatocytes", 
  ages = c(3,6,9,15,30,50)
  )
###calculate_cell_percentage
#### A function to determine the persentage of plasmatocytes of a specific Strain, sex at different ages
#### Example usage
calculate_cell_percentage(seurat_obj = ALL_HEAD_57k3,
                          cell_types = c("Plasmatocytes"),
                          strains =c("DGRP-551", "w1118"),
                          sexes = c("female", "male")
                          )
###plot_percentage
####Barplot of percentage of Plasmatocytes across different ages
#### Example usage
#### Barplot of DGRP-551_female_Plasmatocytes
plot_percentage(calculate_cell_percentage(seurat_obj = ALL_HEAD_57k3,
                                          cell_types = c("Plasmatocytes"),
                                          strains =c("DGRP-551", "w1118"),
                                          sexes = c("female", "male")
                                          ),
                Strain = "DGRP-551",
                Cell_type = "Plasmatocytes", 
                Sex = "female"
                )
###make_heatmap
####heatmap of pseudo bulk expression to compare with our BulkRNAseq
####Example usage
####count_DGRP551_female
make_heatmap(EXPR_OBJ = count_DGRP551_female[,-c(1:2)],
             GENELIST = III_Marker_Genes,
             TITLE = "DGRP551_female"
             )
###heatmap_specific_cluster
####a function to plot heat map of genes at different Age in a specific cell cluster
####Example Usage
####Plasmatocytes_DGRP-551_female
heatmap_specific_cluster(EXPR_COUNT=count_Strain_Sex_Cellcluster_Age,
                         GENELIST=III_Marker_Genes,
                         STRAIN = "DGRP-551",
                         SEX = "female",
                         CLUSTER_NAME="Plasmatocytes",
                         AGES = c("3","6","9","15","30","50"),
                         TITLE = "III_Markers in Plasmatocytes DGRP-551 female"
                         )

#############################################set working directory and library package################################################################################################

## run at HiPerGator (HPG)
setwd("/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/")  

#This is the work.image saved on HPG
#load("GSE107451_head_ALL_57k.RData")


library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(data.table)
library(tidyverse)
library(clustree)
library(openxlsx)
library(presto)
library(EnhancedVolcano)
library(ggrepel)
library(gplots)
#devtools::install_github("immunogenomics/presto")
#############################################################################################################################################




#############################################Import data and create Seurat object################################################################################################

#import barcode, matrix, and feature files obtained by cellranger under folder GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv
ALL_HEAD_57k <- Seurat::Read10X("./GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv")
#create seurat object, setting minimum standard (57k is filtered dataset by author, no need to filter)
ALL_HEAD_57k <- CreateSeuratObject(counts = ALL_HEAD_57k,
                               project = "ALL_HEAD_57k_GSE107451_Summer2024", 
                               min.cells = 0,  
                               min.features = 0)
#look at the meta data
View(ALL_HEAD_57k@meta.data)

#import metadata of author's and add that to Seurat OBJ
meta_author <- as.data.frame(fread("./GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv",header = T))
#check number of cells in the 57k dataset
table(rownames(ALL_HEAD_57k@meta.data)==meta_author$new_barcode)
#TRUE 
#56902 

#adding new_barcode column in meta.data slot of Seurat OBJ from author's metadata
ALL_HEAD_57k@meta.data$"new_barcode" <- rownames(ALL_HEAD_57k@meta.data)
#merge author's metadata to the meta.data slot of Seurat OBJ
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
###NOTE: Running the standard workflow from normalization, pca, etc.. to tSNE is simply for adding the "reduction" slot in the Seurat OBJ
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


ALL_HEAD_57k3 <- JoinLayers(ALL_HEAD_57k3)
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
                             pt.size = 0.5
                             )
  ggsave(tmp,
         device = "pdf",
         file= paste0(
           "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Featureplot/partial_Plasmatocytes/",
           III_Marker_Genes[i],
           "_Plasmatocytes.pdf"),
         width = 9,
         height = 6
  )
  }


#partial feature plot of III markers in plasmatocytes cluster and split by Age
for (i in 1:length(III_Marker_Genes)) {
  tmp <- Seurat::FeaturePlot(subset(ALL_HEAD_57k3,annotation=="Plasmatocytes"),
                             features = III_Marker_Genes[i],
                             cols = c("grey","darkred"),
                             reduction = "tSNE_ALL_HEAD_57k_Author_cluster",
                             keep.scale = NULL,#set to null so it is easier to see if an III marker is expressed at a certain age
                             split.by = "Age",
                             pt.size = 0.5
  )+
    patchwork::plot_layout(ncol = 3, nrow = 3)& theme(legend.position = "right")
  
  ggsave(tmp,
         device = "pdf",
         file= paste0(
           "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Featureplot/partial_Plasmatocytes/",
           III_Marker_Genes[i],
           "_Plasmatocytes_Age.pdf"),
         width = 9,
         height = 6,
         limitsize = FALSE
  )
}

#partial feature plot of III markers in plasmatocytes cluster DGRP-551 female and split by Age
for (i in 1:length(III_Marker_Genes)) {
  tmp <- Seurat::FeaturePlot(subset(ALL_HEAD_57k3,annotation=="Plasmatocytes"&Genotype=="DGRP-551"&sex=="female"),
                             features = III_Marker_Genes[i],
                             cols = c("grey","darkred"),
                             reduction = "tSNE_ALL_HEAD_57k_Author_cluster",
                             keep.scale = NULL,#set to null so it is easier to see if an III marker is expressed at a certain age
                             split.by = "Age",
                             pt.size = 0.5
  )+
    patchwork::plot_layout(ncol = 3, nrow = 3)& theme(legend.position = "right")
  
  ggsave(tmp,
         device = "pdf",
         file= paste0(
           "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Featureplot/partial_Plasmatocytes/",
           III_Marker_Genes[i],
           "_Plasmatocytes_Age_DGRP551_female.pdf"),
         width = 9,
         height = 6,
         limitsize = FALSE
  )
}


#partial feature plot of III markers in plasmatocytes cluster DGRP-551 male and split by Age
for (i in 1:length(III_Marker_Genes)) {
  tmp <- Seurat::FeaturePlot(subset(ALL_HEAD_57k3,annotation=="Plasmatocytes"&Genotype=="DGRP-551"&sex=="male"),
                             features = III_Marker_Genes[i],
                             cols = c("grey","darkred"),
                             reduction = "tSNE_ALL_HEAD_57k_Author_cluster",
                             keep.scale = NULL,#set to null so it is easier to see if an III marker is expressed at a certain age
                             split.by = "Age",
                             pt.size = 0.5
  )+
    patchwork::plot_layout(ncol = 3, nrow = 3)& theme(legend.position = "right")
  
  ggsave(tmp,
         device = "pdf",
         file= paste0(
           "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Featureplot/partial_Plasmatocytes/",
           III_Marker_Genes[i],
           "_Plasmatocytes_Age_DGRP551_male.pdf"),
         width = 9,
         height = 6,
         limitsize = FALSE
  )
}


#partial feature plot of III markers in plasmatocytes cluster w1118 female and split by Age
for (i in 1:length(III_Marker_Genes)) {
  tmp <- Seurat::FeaturePlot(subset(ALL_HEAD_57k3,annotation=="Plasmatocytes"&Genotype=="w1118"&sex=="female"),
                             features = III_Marker_Genes[i],
                             cols = c("grey","darkred"),
                             reduction = "tSNE_ALL_HEAD_57k_Author_cluster",
                             keep.scale = NULL,#set to null so it is easier to see if an III marker is expressed at a certain age
                             split.by = "Age",
                             pt.size = 0.5
  )+
    patchwork::plot_layout(ncol = 3, nrow = 3)& theme(legend.position = "right")
  
  ggsave(tmp,
         device = "pdf",
         file= paste0(
           "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Featureplot/partial_Plasmatocytes/",
           III_Marker_Genes[i],
           "_Plasmatocytes_Age_w1118_female.pdf"),
         width = 9,
         height = 6,
         limitsize = FALSE
  )
}


#partial feature plot of III markers in plasmatocytes cluster w1118 male and split by Age
for (i in 1:length(III_Marker_Genes)) {
  tmp <- Seurat::FeaturePlot(subset(ALL_HEAD_57k3,annotation=="Plasmatocytes"&Genotype=="w1118"&sex=="male"),
                             features = III_Marker_Genes[i],
                             cols = c("grey","darkred"),
                             reduction = "tSNE_ALL_HEAD_57k_Author_cluster",
                             keep.scale = NULL,#set to null so it is easier to see if an III marker is expressed at a certain age
                             split.by = "Age",
                             pt.size = 0.5
  )+
    patchwork::plot_layout(ncol = 3, nrow = 3)& theme(legend.position = "right")
  
  ggsave(tmp,
         device = "pdf",
         file= paste0(
           "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Featureplot/partial_Plasmatocytes/",
           III_Marker_Genes[i],
           "_Plasmatocytes_Age_w1118_male.pdf"),
         width = 9,
         height = 6,
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






############################make a function to see number of cells in different Strains, Sex at different Ages####################################################################################################
#This is for making the heatmap later, because for example, if there is no day 50 cell in DGRP-551 male, 
#then the heatmap will not show expression as its pseudo-bulk expression will not be calculated


count_cells_subset <- function(seurat_obj, strain = NULL, sex = NULL, age = NULL, cell_cluster = NULL) {
  # Extract metadata from the Seurat object
  meta_data <- seurat_obj@meta.data
  
  # Apply filters based on the provided parameters
  if (!is.null(strain)) {
    meta_data <- meta_data %>% filter(Genotype == strain)
  }
  if (!is.null(sex)) {
    meta_data <- meta_data %>% filter(sex == sex)
  }
  if (!is.null(age)) {
    meta_data <- meta_data %>% filter(Age == age)
  }
  if (!is.null(cell_cluster)) {
    meta_data <- meta_data %>% filter(annotation == cell_cluster)
  }
  
  # Count the number of cells for each combination of Strain, Sex, Age, and cell cluster
  cell_counts <- meta_data %>%
    group_by(Genotype, sex, Age, annotation) %>%
    summarise(CellCount = n())
  
  return(cell_counts)
  print(cell_counts,n=Inf)
}

# Example usage:
# To get the number of cells for a specific strain, sex, age, and cell cluster:
count_cells_subset(seurat_obj=ALL_HEAD_57k3,
                                  strain = "DGRP-551",
                                  sex = "female",
                                  age = "50",
                                  cell_cluster = "Plasmatocytes")


# To get the number of cells for all combinations:
print(
  count_cells_subset(seurat_obj=ALL_HEAD_57k3),
  n=Inf)

cell_numbers <- print(
  count_cells_subset(seurat_obj=ALL_HEAD_57k3),
  n=Inf)

#to see if the function generates correct number of cells
#in plasmatocytes cluster, there are 35 cells in DGRP-551 female day 50; 

nrow(
  subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="female"&Age=="50"&annotation=="Plasmatocytes")@meta.data
)
#35

#in plasmatocytes cluster, there are 14 cells in DGRP-551 male day 50; 
nrow(
  subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="male"&Age=="50"&annotation=="Plasmatocytes")@meta.data
)
#14

## THe function provides correct number of cells based on different requirements

#export the cell number into xlsx
openxlsx::write.xlsx(cell_numbers,"cell_numbers.xlsx",rowNames=F)
#############################################################################################################################################




##################################percentage expressing function ###########################################################################################################
#calculate the persentage of cells expressing different genes within a subset of a seurat object of a certain Strain, Sex, Cluster, and Age 
calculate_gene_expression_percentage <- function(seurat_obj, genes, strain, Sex, cluster, ages) {
  # Create a data frame to store the results
  results <- data.frame(
    Gene = character(),
    Age = character(),
    Percentage = numeric(),
    Cells=numeric(),
    TotalCells=numeric()
  )
  
  # Iterate over each gene
  for (gene_i in genes) {
    # Iterate over each age
    for (age_i in ages) {
      # Subset the Seurat object based on the specified metadata
      subset_cells <- subset(seurat_obj,Genotype == strain & sex == Sex & annotation == cluster & Age == age_i)
      
      
      # Calculate the number of cells expressing the gene (expression > 0)
      
      expressing_cells <-
        sum(
          as.data.frame(
            FetchData(subset_cells, vars = gene_i,layer = "counts")
            ) > 0
          )
      
      
      # Calculate the total number of cells in the subset
      total_cells <- nrow(subset_cells@meta.data)
      
      # Calculate the percentage of cells expressing the gene
      #expression_percentage <- expressing_cells/total_cells*100
      expression_percentage <- if (total_cells > 0) (expressing_cells / total_cells) * 100 else NA
      
      # Add the result to the data frame
      results <- rbind(
        results, data.frame(Gene = gene_i,
                            Age = age_i, 
                            Percentage = expression_percentage,
                            Cells=expressing_cells,
                            TotalCells=total_cells
                            )
        )
    }
  }
  
  return(results)
}

# Example usage:
View(
calculate_gene_expression_percentage(
  seurat_obj=ALL_HEAD_57k3,
  genes=III_Marker_Genes, 
  strain = "DGRP-551", 
  Sex = "female", 
  cluster = "Plasmatocytes", 
  ages = c(3,6,9,15,30,50)
  )
)


#check with simple code
#  gene   day  percentage at day 50 in DGRP-551 female in Plasmatocytes
#SPH93  50 0.05714286

#check cells expressing SPH93 at day 50 in DGRP-551 female in Plasmatocytes
nrow(
  subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="female"&Age=="50"&annotation=="Plasmatocytes")@'meta.data'
)
#35 cells in total

#check number of these cell that express SPH93
sum(
  FetchData(subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="female"&Age=="50"&annotation=="Plasmatocytes"),
            vars = 'SPH93',layer = "counts") > 0
  )

#2 cells express SPH93 at day 50 in DGRP-551 female in Plasmatocytes

nrow(FetchData(subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&sex=="female"&Age=="50"&annotation=="Plasmatocytes"), vars = 'SPH93'))
#35 cells in total


2/35*100
#5.714286
####

#export PercentageExpression_III_Markers_Ages_DGRP551_female_Plasmatocytes
openxlsx::write.xlsx(calculate_gene_expression_percentage(
  seurat_obj=ALL_HEAD_57k3,
  genes=III_Marker_Genes, 
  strain = "DGRP-551", 
  Sex = "female", 
  cluster = "Plasmatocytes", 
  ages = c(3,6,9,15,30,50)
),
"PercentageExpression_III_Markers_Ages_DGRP551_female_Plasmatocytes.xlsx",
rowNames=F
)

#export PercentageExpression_III_Markers_Ages_DGRP551_male_Plasmatocytes
openxlsx::write.xlsx(calculate_gene_expression_percentage(
  seurat_obj=ALL_HEAD_57k3,
  genes=III_Marker_Genes, 
  strain = "DGRP-551", 
  Sex = "male", 
  cluster = "Plasmatocytes", 
  ages = c(3,6,15,30,50)
),
"PercentageExpression_III_Markers_Ages_DGRP551_male_Plasmatocytes.xlsx",
rowNames=F
)

#export PercentageExpression_III_Markers_Ages_w1118_female_Plasmatocytes
openxlsx::write.xlsx(calculate_gene_expression_percentage(
  seurat_obj=ALL_HEAD_57k3,
  genes=III_Marker_Genes, 
  strain = "w1118", 
  Sex = "female", 
  cluster = "Plasmatocytes", 
  ages = c(3,6,9,30)
),
"PercentageExpression_III_Markers_Ages_w1118_female_Plasmatocytes.xlsx",
rowNames=F
)
#check what is wrong with day 15
nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="female"&Age=="15"&annotation=="Plasmatocytes")@meta.data
  )
#there is only one cell... in w1118 female day 15 Plasmatocytes
sum(
  FetchData(subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="female"&Age=="15"&annotation=="Plasmatocytes"),
            vars = 'SPH93') > 0
)

View(subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="female"&Age=="15"&annotation=="Plasmatocytes")@meta.data)
#check all cells in all ages
nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="female"&Age=="3"&annotation=="Plasmatocytes")@meta.data
)#12 cells
nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="female"&Age=="6"&annotation=="Plasmatocytes")@meta.data
)#35 cells
nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="female"&Age=="9"&annotation=="Plasmatocytes")@meta.data
)#20 cells
nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="female"&Age=="15"&annotation=="Plasmatocytes")@meta.data
)#1 cells...
nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="female"&Age=="30"&annotation=="Plasmatocytes")@meta.data
)# 10 cells


#export PercentageExpression_III_Markers_Ages_w1118_male_Plasmatocytes
openxlsx::write.xlsx(calculate_gene_expression_percentage(
  seurat_obj=ALL_HEAD_57k3,
  genes=III_Marker_Genes, 
  strain = "w1118", 
  Sex = "male", 
  cluster = "Plasmatocytes", 
  ages = c(3,6,9,30)
),
"PercentageExpression_III_Markers_Ages_w1118_male_Plasmatocytes.xlsx",
rowNames=F
)

#check number of cells in all ages, because if there is only one cell, result cannot be generated (Seurat design)
nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="male"&Age=="3"&annotation=="Plasmatocytes")@meta.data
)#11 cells
nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="male"&Age=="6"&annotation=="Plasmatocytes")@meta.data
)#19 cells
nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="male"&Age=="9"&annotation=="Plasmatocytes")@meta.data
)#8 cells
nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="male"&Age=="15"&annotation=="Plasmatocytes")@meta.data
)# 1 cell ...
nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&sex=="male"&Age=="30"&annotation=="Plasmatocytes")@meta.data
)# 3 cells
#############################################################################################################################################



############################Look at number of cells in DGRP-551 and w1118#################################################################################################################

nrow(
  subset(ALL_HEAD_57k3,Genotype=="DGRP-551")@meta.data
)
#29137 cells in DGRP-551


nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118")@meta.data
)
#27765 cells in w1118

###let's check in Plasmatocytes

nrow(
  subset(ALL_HEAD_57k3,Genotype=="DGRP-551"&annotation=="Plasmatocytes")@meta.data
)
#200 cells in DGRP-551 in Plasmatocytes


nrow(
  subset(ALL_HEAD_57k3,Genotype=="w1118"&annotation=="Plasmatocytes")@meta.data
)
#129 cells in w1118 in Plasmatocytes

nrow(
  subset(ALL_HEAD_57k3,annotation=="Plasmatocytes")@meta.data
)
#329 total cells in plasmatocytes


#############################################################################################################################################



#############################################################################################################################################
#A function to determine the persentage of plasmatocytes of a specific Strain, sex at different ages

calculate_cell_percentage <- function(seurat_obj, cell_types, strains, sexes) {
  # Initialize an empty dataframe to store results
  result_df <- data.frame()
  
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Iterate over each combination of strain, sex, and cell type
  for (strain in strains) {
    for (sex_i in sexes) {
      for (cell_type in cell_types) {
        # Filter metadata for the specific strain and sex
        filtered_metadata <- metadata %>%
          filter(Genotype == strain, sex == sex_i)
        
        # Group by Age and calculate the percentage of the specified cell type
        percentage_df <- filtered_metadata %>%
          group_by(Age) %>%
          summarise(
            TotalCells = n(),
            CellTypeCells = sum(annotation == cell_type),
            Percentage = (CellTypeCells / TotalCells) * 100
          ) %>%
          mutate(Genotype = strain, sex = sex_i, annotation = cell_type)
        
        # Append the results to the result_df
        result_df <- bind_rows(result_df, percentage_df)
      }
    }
  }
  
  return(result_df)
}

# Example usage
calculate_cell_percentage(seurat_obj = ALL_HEAD_57k3,
                          cell_types = c("Plasmatocytes"),
                          strains =c("DGRP-551", "w1118"),
                          sexes = c("female", "male")
                          )
#export the result into xlsx

openxlsx::write.xlsx(
  calculate_cell_percentage(seurat_obj = ALL_HEAD_57k3,
                            cell_types = c("Plasmatocytes"),
                            strains =c("DGRP-551", "w1118"),
                            sexes = c("female", "male")
  ),
  "Plasmatocytes_percentage_across_ages.xlsx",
  rowNames=F
)
#############################################################################################################################################








############################################Barplot of percentage of Plasmatocytes across different ages#################################################################################################

plot_percentage <- function(df, Strain, Cell_type, Sex) {
  # Filter the dataframe based on the input parameters
  df_filtered <- df[df$Genotype == Strain & df$annotation == Cell_type & df$sex == Sex, ]
  
  # Load necessary library
  library(ggplot2)
  
  # Create the barplot
  ggplot(df_filtered, aes(x = as.factor(Age), y = Percentage, fill = as.factor(Age))) +
    geom_bar(stat = "identity") +
    labs(title = paste("Percentage of", Cell_type, "for", Strain, "(", Sex, ")", "across Ages"),
         x = "Age",
         y = "Percentage") +
    theme_minimal() +
    theme(legend.position = "none")
}

# Example of how to use the function

#DGRP-551_female_Plasmatocytes
pdf(file = "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Barplot/DGRP-551_female_Plasmatocytes.pdf",width = 6,height = 4)
plot_percentage(calculate_cell_percentage(seurat_obj = ALL_HEAD_57k3,
                                          cell_types = c("Plasmatocytes"),
                                          strains =c("DGRP-551", "w1118"),
                                          sexes = c("female", "male")
                                          ),
                Strain = "DGRP-551",
                Cell_type = "Plasmatocytes", 
                Sex = "female"
                )
dev.off()

#DGRP-551_male_Plasmatocytes
pdf(file = "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Barplot/DGRP-551_male_Plasmatocytes.pdf",width = 6,height = 4)
plot_percentage(calculate_cell_percentage(seurat_obj = ALL_HEAD_57k3,
                                          cell_types = c("Plasmatocytes"),
                                          strains =c("DGRP-551", "w1118"),
                                          sexes = c("female", "male")
),
Strain = "DGRP-551",
Cell_type = "Plasmatocytes", 
Sex = "male"
)
dev.off()

#w1118_female_Plasmatocytes
pdf(file = "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Barplot/w1118_female_Plasmatocytes.pdf",width = 6,height = 4)
plot_percentage(calculate_cell_percentage(seurat_obj = ALL_HEAD_57k3,
                                          cell_types = c("Plasmatocytes"),
                                          strains =c("DGRP-551", "w1118"),
                                          sexes = c("female", "male")
),
Strain = "w1118",
Cell_type = "Plasmatocytes", 
Sex = "female"
)
dev.off()

#w1118_male_Plasmatocytes
pdf(file = "/orange/zhou/projects/II_Cancer/GSE107451_HEAD_57K_filtered/Barplot/w1118_male_Plasmatocytes.pdf",width = 6,height = 4)
plot_percentage(calculate_cell_percentage(seurat_obj = ALL_HEAD_57k3,
                                          cell_types = c("Plasmatocytes"),
                                          strains =c("DGRP-551", "w1118"),
                                          sexes = c("female", "male")
),
Strain = "w1118",
Cell_type = "Plasmatocytes", 
Sex = "male"
)
dev.off()









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




####################################################heatmap of pseudo bulk expression to compare with our BulkRNAseq#########################################################################################

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
heatmap_specific_cluster <- function(EXPR_COUNT,GENELIST,STRAIN,SEX,CLUSTER_NAME,AGES,TITLE) {
  tmp_count <- EXPR_COUNT[,str_split(colnames(EXPR_COUNT),pattern = "_",4,simplify =T)[,1]==STRAIN];
  tmp_count <- tmp_count[,str_split(colnames(tmp_count),pattern = "_",4,simplify =T)[,2]==SEX];
  tmp_count <- tmp_count[,str_split(colnames(tmp_count),pattern = "_",4,simplify =T)[,3]==CLUSTER_NAME];
  tmp_count <- tmp_count[,str_split(colnames(tmp_count),pattern = "_",4,simplify =T)[,4]%in%AGES];
  tmp_count2 <-tmp_count[rownames(tmp_count)%in%GENELIST,];
  
  #remove rows where genes have all zero values across columns
  tmp_count2 <- tmp_count2[rowSums(tmp_count2)!=0,];
  #change column names to days
  colnames(tmp_count2) <- paste0("Day ", AGES);
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
                         AGES = c("3","6","9","15","30","50"),
                         TITLE = "III_Markers in Plasmatocytes DGRP-551 female")
dev.off()

#Plasmatocytes_DGRP-551_male
pdf("./Heatmap/Plasmatocytes_DGRP-551_male.pdf",width = 15,height = 10)
heatmap_specific_cluster(EXPR_COUNT=count_Strain_Sex_Cellcluster_Age,
                         GENELIST=III_Marker_Genes,
                         STRAIN = "DGRP-551",
                         SEX = "male",
                         CLUSTER_NAME="Plasmatocytes",
                         AGES = c("3","6","15","30","50"),
                         TITLE = "III_Markers in Plasmatocytes DGRP-551 male")
dev.off()

#Plasmatocytes_w1118_female
pdf("./Heatmap/Plasmatocytes_w1118_female.pdf",width = 15,height = 10)
heatmap_specific_cluster(EXPR_COUNT=count_Strain_Sex_Cellcluster_Age,
                         GENELIST=III_Marker_Genes,
                         STRAIN = "w1118",
                         SEX = "female",
                         CLUSTER_NAME="Plasmatocytes",
                         AGES = c("3","6","9","15","30"),
                         TITLE = "III_Markers in Plasmatocytes w1118 female")
dev.off()

#Plasmatocytes_w1118_male
pdf("./Heatmap/Plasmatocytes_w1118_male.pdf",width = 15,height = 10)
heatmap_specific_cluster(EXPR_COUNT=count_Strain_Sex_Cellcluster_Age,
                         GENELIST=III_Marker_Genes,
                         STRAIN = "w1118",
                         SEX = "male",
                         CLUSTER_NAME="Plasmatocytes",
                         AGES = c("3","6","9","15","30"),
                         TITLE = "III_Markers in Plasmatocytes w1118 male")
dev.off()
###



####make similar heat map in other cluster (if needed), TBA



####

#############################################################################################################################################





##################################Based on copilot###########################################################################################################
#NOTE: This function generated by copilot is only for an attempt; It is not used in the real analysis
#make a function to plot heat map of genes at different Age in a specific cell cluster based on the pseudo-bulk expression of different Strain, Age, Sex, and Cell cluster

# Function to obtain pseudo-bulk expression and plot heatmap
analyze_seurat_object <- function(seurat_object, gene_list, strain, sex, cell_cluster) {
  # Filter the Seurat object based on specified conditions
  filtered_seurat <- subset(seurat_object,
                            Strain == strain & Sex == sex & `Cell cluster` == cell_cluster)
  
  # Calculate pseudo-bulk expression
  pseudo_bulk_expression <- aggregateexpression(filtered_seurat)
  
  # Subset the expression matrix to include only selected genes
  expression_matrix <- t(as.matrix(pseudo_bulk_expression))
  expression_matrix <- expression_matrix[gene_list, ]
  
  # Create a heatmap
  heatmap_data <- as.data.frame(expression_matrix)
  colnames(heatmap_data) <- rownames(expression_matrix)
  heatmap_data$Age <- filtered_seurat$Age
  
  # Plot the heatmap
  library(ggplot2)
  ggplot(heatmap_data, aes(x = Age, y = rownames(heatmap_data), fill = .)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = paste("Heatmap of Gene Expression (", strain, ", ", sex, ", ", cell_cluster, ")", sep = ""),
         x = "Age", y = "Genes") +
    theme_minimal()
}

# Example usage:
seurat_object <- Read10X(data.dir = "path/to/your/10x_data")  # Load your Seurat object
gene_list_of_interest <- c("Gene1", "Gene2", "Gene3")  # Specify your gene list
analyze_seurat_object(seurat_object, gene_list_of_interest,
                      strain = "StrainA", sex = "Male", cell_cluster = "Cluster1")
#############################################################################################################################################






##################################Based on chatgpt###########################################################################################################
#NOTE: This function generated by chatgpt is only for an attempt; It is not used in the real analysis
#make a function to plot heat map of genes at different Age in a specific cell cluster based on the pseudo-bulk expression of different Strain, Age, Sex, and Cell cluster

library(Seurat)
library(pheatmap)
library(dplyr)
library(tidyr)

plot_heatmap_genes <- function(seurat_obj, genelist, specific_age, cell_cluster, specific_strain, specific_sex) {
  # Step 1: Obtain the pseudo-bulk expression
  get_pseudo_bulk <- function(seurat_obj) {
    # Aggregate expression based on metadata
    pseudo_bulk <- AggregateExpression(seurat_obj, 
                                       group.by = c("Strain", "Sex", "Age", "cell_cluster"), 
                                       assays = "RNA")
    return(pseudo_bulk$RNA)
  }
  
  # Step 2: Subset pseudo-bulk expression matrix for the specified cell cluster, strain, and sex
  subset_pseudo_bulk <- function(pseudo_bulk, specific_age, cell_cluster, specific_strain, specific_sex) {
    # Convert the Seurat aggregated matrix to a data frame for easier subsetting
    pseudo_bulk_df <- as.data.frame(t(pseudo_bulk))
    pseudo_bulk_df <- rownames_to_column(pseudo_bulk_df, var = "CellCluster")
    
    # Filter based on specific cell cluster, strain, sex, and age
    filtered_df <- pseudo_bulk_df %>%
      filter(grepl(cell_cluster, CellCluster),
             grepl(specific_strain, CellCluster),
             grepl(specific_sex, CellCluster),
             grepl(paste0("_", specific_age, "_"), CellCluster))
    
    # Retain only the genes of interest
    filtered_df <- filtered_df[, c("CellCluster", genelist)]
    
    # Set the row names as the CellCluster
    rownames(filtered_df) <- filtered_df$CellCluster
    filtered_df <- filtered_df[, -1] # Remove CellCluster column
    
    return(filtered_df)
  }
  
  # Step 3: Plot heatmap for the filtered pseudo-bulk data
  plot_heatmap <- function(filtered_df) {
    pheatmap(as.matrix(filtered_df), cluster_rows = TRUE, cluster_cols = TRUE,
             main = paste("Heatmap of genes in cell cluster", cell_cluster, 
                          "for Strain:", specific_strain, 
                          "Sex:", specific_sex, 
                          "Ages:", paste(specific_age, collapse=", ")))
  }
  
  # Obtain the pseudo-bulk expression
  pseudo_bulk <- get_pseudo_bulk(seurat_obj)
  
  # Subset the pseudo-bulk expression matrix
  filtered_df <- subset_pseudo_bulk(pseudo_bulk, specific_age, cell_cluster, specific_strain, specific_sex)
  
  # Plot the heatmap
  plot_heatmap(filtered_df)
}

# Example usage:
# Assuming `seurat_obj` is your Seurat object and you are interested in genes "Gene1", "Gene2", and "Gene3"
# and you want to visualize data for cell cluster "Cluster1", strain "Strain1", sex "Male", and ages 3 and 6

plot_heatmap_genes(seurat_obj = ALL_HEAD_57k3,
                   genelist = III_Marker_Genes,
                   specific_age = c(3, 6, 9, 15, 30, 50),
                   cell_cluster = "Plasmatocytes",
                   specific_strain = "DGRP-551",
                   specific_sex = "female"
                   )

#############################################################################################################################################




###################################I asked chatgpt to optmize my function##########################################################################################################
#NOTE: This function generated by chatgpt is only for an attempt; It is not used in the real analysis
library(dplyr)
library(pheatmap)

# heatmap_specific_cluster <- function(EXPR_COUNT, GENELIST, STRAIN, SEX, CLUSTER_NAME, TITLE) {
#   # Extract metadata from column names
#   metadata <- str_split(colnames(EXPR_COUNT),pattern = "_",4,simplify =T)
#   colnames(metadata) <- c("Strain", "Sex", "Cluster")
#   
#   # Combine EXPR_COUNT and metadata into a single data frame
#   expr_data <- as.data.frame(EXPR_COUNT)
#   expr_data <- bind_cols(metadata, expr_data)
#   
#   # Filter data based on the specified STRAIN, SEX, and CLUSTER_NAME
#   filtered_data <- expr_data %>%
#     filter(Strain == STRAIN & Sex == SEX & Cluster == CLUSTER_NAME)
#   
#   # Subset the expression counts for the genes of interest
#   gene_data <- filtered_data %>%
#     select(starts_with("V")) %>%
#     select(one_of(GENELIST))
#   
#   # Convert to matrix
#   gene_matrix <- as.matrix(gene_data)
#   rownames(gene_matrix) <- filtered_data$Cluster
#   
#   # Plot heatmap
#   pheatmap(gene_matrix, 
#            scale = "row",
#            color = colorRampPalette(c("blue", "white", "red"))(100),
#            main = TITLE,
#            cluster_rows = FALSE,
#            cluster_cols = FALSE,
#            display_numbers = TRUE,
#            fontsize_row = 10,
#            fontsize_col = 10,
#            border_color = NA)
# }

# Example usage:
# Assuming `EXPR_COUNT` is your expression matrix
# and you are interested in genes "Gene1", "Gene2", and "Gene3"
# and you want to visualize data for strain "Strain1", sex "Male", and cluster "Cluster1"
heatmap_specific_cluster(EXPR_COUNT=count_Strain_Sex_Cellcluster_Age,
                         GENELIST = III_Marker_Genes, 
                         STRAIN = "DGRP-551",
                         SEX = "female", 
                         CLUSTER_NAME = "Plasmatocytes", 
                         TITLE = "Heatmap of Plasmatocytes")




#############################################################################################################################################








#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

save.image("GSE107451_head_ALL_57k.RData")
