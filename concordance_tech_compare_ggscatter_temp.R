library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(SCpubr)
library(sctransform)
sample1.data=Read10X(data.dir ='/path/outs/filtered_feature_bc_matrix/')
sample1=CreateSeuratObject(counts=sample1.data, project ='sample1')
sample1
sample2.data= Read10X(data.dir = '/path/outs/filtered_feature_bc_matrix/')
sample2=CreateSeuratObject(counts = sample2.data, project = 'sample2')
sample2
  #add identification by tech
sample1$tech="tech1"
sample2$tech='tech2'
  
  ##SCT normalize
S1=SCTransform(sample1, vst.flavor = "v2", ncells = 10000, verbose = FALSE) %>%
    RunPCA(npcs = 35, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:35, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:35, verbose = FALSE) %>%
    FindClusters(resolution = 0.8, verbose = FALSE)
p1 <- DimPlot(S1, label = T, repel = T) + ggtitle("sample1")
p1
  
  ##Perform integration using pearson residuals
S2=SCTransform(sample2, vst.flavor = "v2", ncells = 10000, verbose = FALSE) %>%
    RunPCA(npcs = 35, verbose = FALSE)
Scomb.list = list(S1=S1, S2=S2)
S_features <- SelectIntegrationFeatures(object.list = Scomb.list, nfeatures = 5000)
Scomb.list=PrepSCTIntegration(object.list = Scomb.list, anchor.features = S_features)
Scomb.anchors=FindIntegrationAnchors(object.list = Scomb.list, normalization.method = "SCT",
                                       anchor.features = S_features)
Scomb.sct=IntegrateData(anchorset = Scomb.anchors, normalization.method = 'SCT')

  
  ##Perform an integrated analysis
Scomb.sct=RunPCA(Scomb.sct, verbose = FALSE)
Scomb.sct=RunUMAP(Scomb.sct, reduction = "pca", dims = 1:35, verbose = FALSE)
Scomb.sct=FindNeighbors(Scomb.sct, reduction = "pca", dims = 1:35)
Scomb.sct=FindClusters(Scomb.sct)
  
DimPlot(Scomb.sct, reduction = 'umap', split.by = 'tech')
p2=DimPlot(Scomb.sct, reduction = 'umap', group.by = 'tech')
p2
  
  ###Identify differential expressed genes across conditions
##Corrected counts are obtained by setting the sequencing depth for all the cells to a fixed value and reversing the learned regularized negative-binomial regression model.
###Prior to performing differential expression, we first run PrepSCTFindMarkers, which ensures that the fixed value is set properly.
Scomb.sct=PrepSCTFindMarkers(Scomb.sct)
  
  
  ###assign idents by tech, default was by cluster.
Idents(Scomb.sct) = Scomb.sct$tech
  
  
  ##############
  
cell.markers=FindMarkers(Scomb.sct, ident.1 = 'tech1', ident.2 = 'tech2', test.use = 'wilcox')

  ###average genes expression in each sample
  
avg.cells <- as.data.frame(log1p(AverageExpression(Scomb.sct, verbose = FALSE, group.by = 'tech')$RNA))

  
avg.cells$gene <- rownames(avg.cells)

  
  
  ##order() will sort the given numbers according to its index in the ascending order.
  ###x[order(x)] Here the indexing of order is done where the actual values are printed in the ascending order.
  ###24h-fresh are acceding order, means negative will be first, means fresh high expression
  
topdif = c(rownames(head(avg.cells[order(-avg.cells$tech1 + avg.cells$'tech2'),],10)),
             rownames(head(avg.cells[order(+avg.cells$tech1 - avg.cells$'tech2'),],10)))
topdif
  ###add column named color, can equal to anything, for label purpose, set it empty
avg.cells$color='none'
  
  ##top different genes in color column named 'topdif'
avg.cells[avg.cells$gene %in% topdif, ]$color="topdif"
  #
plot_gg= ggscatter(avg.cells, x = "tech1", y ="tech2", fill =  'color',shape =21,
                     palette = c(none = "#F4F8FB", topdif="#BC1119"),
                     size = 1.2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray")) +
    stat_cor(method=c("pearson", "kendall", "spearman"), label.x = 0.5, label.y = 5)+
    ggtitle("tech1 vs tech2")
  
plot_label=LabelPoints(plot=plot_gg, points = topdif, color = "red", repel = TRUE )
plot_label
  
  
saveRDS(avg.cells, file = "/volumes/...../concordance.rds")
  
