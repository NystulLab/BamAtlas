####################################
# data preprocessing and analysis  #
####################################

# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

library(future)
options(future.globals.maxSize = 100 * 1024 * 1024^2)
plan("multiprocess", workers = 4)

### Creat Seurat Object
message( "==>Reading 10x data<==" )
##rawdata/*sample  contains : barcodes.tsv,genes.tsv,matrix.mtx
samples <- c("ck","Bam")
for ( i in samples){
        mat <- Read10X(data.dir = paste0("rawdata/",i), gene.column = 1)
        object.list[[i]] <- CreateSeuratObject(counts = mat, project = parameter$data$name[i], assay = "RNA" )
}
obj <- merge(x = object.list[[1]], y = unlist(object.list[-1]), add.cell.ids = samples)

color.sample <- scales::hue_pal()(length(levels(obj@meta.data$orig.ident)))
names(color.sample) <- levels(obj@meta.data$orig.ident)
obj@misc[["color.sample"]] <- color.sample

##add Feature Data
ref_name_file <- "name_list.xls"
fdata <- read.table(ref_name_file, row.names = 1, stringsAsFactors = F, sep = "\t")
colnames(fdata) <- c("merge_name", "name", "type")
fdata$merge_name <- fdata$name
fdata$merge_name[fdata$merge_name == "-"] <- rownames(fdata)[fdata$merge_name == "-"]
index <- c(which(duplicated(fdata$merge_name, fromLast=T)), which(duplicated(fdata$merge_name, fromLast=F)))
fdata$merge_name[index] <- paste0(fdata$merge_name[index], " (", rownames(fdata)[index], ")")
fdata <- AddUnderscore(fdata)
obj@misc[["fdata"]] <- fdata

##add p data
obj@misc[["pdata"]] <- FetchData(obj, c("orig.ident", "nFeature_RNA", "nCount_RNA"))

##add mito.percent
mito_list <- c(  "FBgn0013672",  #ATP6
  "FBgn0013673",  #ATP8
  "FBgn0013674",  #COX1
  "FBgn0013675",  #COX2
  "FBgn0013676",  #COX3
  "FBgn0013678",  #CYTB
  "FBgn0013679",  #ND1
  "FBgn0013680",  #ND2
  "FBgn0013681",  #ND3
  "FBgn0013683",  #ND4L
  "FBgn0013684",  #ND5
  "FBgn0013685",  #ND6
  "FBgn0262952"  #ND4
)

metadata <- Matrix::colSums(x = GetAssayData(object = obj, slot = "counts", assay = "RNA")[mito_list, , drop = FALSE])
metadata <- metadata / obj@meta.data[["nCount_rna"]] * 100
obj <- AddMetaData(object = obj, metadata = metadata, col.name = "percent.mito")


###Filter Cells
cells.use <- Cells(obj)
#nCount_RNA: [ -.inf,25000 ]
standard <- c(-.inf,25000)
cells.use <- obj@meta.data %>% tibble::rownames_to_column(var = "cells") %>% filter(.data[["nCount_RNA"]] >= standard[1] & .data[["nCount_RNA"]] <= standard[2] & cells %in% cells.use) %>% select(cells) %>% unlist()
#nFeature_RNA: [ 100,3000 ]
standard <- c(100,3000)
cells.use <- obj@meta.data %>% tibble::rownames_to_column(var = "cells") %>% filter(.data[["nFeature_RNA"]] >= standard[1] & .data[["nFeature_RNA"]] <= standard[2] & cells %in% cells.use) %>% select(cells) %>% unlist()
#percent.mito: [ -.inf,30 ]
standard <- c(-.inf,30)
cells.use <- obj@meta.data %>% tibble::rownames_to_column(var = "cells") %>% filter(.data[["percent.mito"]] >= standard[1] & .data[["percent.mito"]] <= standard[2] & cells %in% cells.use) %>% select(cells) %>% unlist()

##Do Filt
obj <- obj[,cells.use]

### Normalization Data
obj <- NormalizeData(obj, normalization.method = normalization.method, scale.factor = scale.factor)

### Find Variable Genes
pdf("Variable_gene.pdf")
obj <- FindVariableGenes( object = obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.01, y.cutoff = 0.5, col.use = "blue", contour.col = "red", do.recalc = TRUE, contour.lwd = 1, contour.lty = 1 )
dev.off()

###Confound factors and Scale Data
vars.regress <- c( "nUMI", "percent.mito" )
obj <- ScaleData(object = obj, vars.to.regress = obj@misc$vars.regress, features = rownames(obj))

### 进行CCA 校正
old.assay <- DefaultAssay(obj)
obj.list <- SplitObject(obj, split.by = 'orig.ident')
anchor.features <- 3000
normalization.method <- "LogNormalize"
k.filter <- min(200, ceiling(min(sapply(obj.list, ncol))/2))
anchors <- FindIntegrationAnchors(obj.list = obj.list, dims = 1:50, normalization.method = normalization.method, anchor.features = anchor.features, k.filter = k.filter)
integrated <- IntegrateData(anchorset = anchors, dims = 1:50, normalization.method = normalization.method)
integrated <- ScaleData(integrated, verbose = FALSE)
integrated@misc <- obj@misc
integrated[[old.assay]] <- obj[[old.assay]]
integrated@reductions <- obj@reductions

### PCA
message( "==>PCA<==" )
pc.num <- 51
obj <- RunPCA( object = obj, pcs.compute = pc.num, pc.genes = obj@var.genes, do.print = FALSE )
pdf("pcaPlot.pdf")
PCAPlot( object = obj, dim.1 = 1, dim.2 = 2, group.by = "orig.ident" )
if ( "CC.Difference" %in% vars.regress ) PCAPlot( object = obj, dim.1 = 1, dim.2 = 2, group.by = "Phase" )
dev.off()
pdf("pcHeatmap.pdf")
#PCHeatmap( object = obj, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE )
PCHeatmap( object = obj, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE )
dev.off()

### JackStraw for find significant PCs
obj <- JackStraw( object = obj, num.pc = pc.num, num.replicate = 100 )
pdf("Jackstraw.pdf")
JackStrawPlot( object = obj, PCs = draw.pcs )
dev.off()

##根据JackStraw图选取 1:50 PCs
sig.PCs <- 1:50

### tSNE 
obj <- RunTSNE( object = obj, dims.use = sig.PCs, do.fast = TRUE, perplexity = 50, max_iter = 2000, verbose = T )
### UMAP
obj <- RunUMAP(obj, dims = sig.PCs, umap.method = "uwot")


### Find clusters
clustering_resolution <- 0.6   # Resolution parameter for Seurat clustering     
obj <- FindClusters( object = obj, reduction.type = "pca", dims.use = sig.PCs, resolution = cluster_resolution, print.output = T, save.SNN = TRUE, force.recalc = TRUE, temp.file.location = getwd() )
color.cluster <- rainbow(length(levels(obj@ident)))
names(color.cluster) <- levels(obj@ident)


## Draw t-SNE plot (对应part6)
message( "==>Draw t-SNE plot<==" )
p1 <- TSNEPlot( object = obj, do.return = T, group.by = "orig.ident", colors.use = color.sample, pt.size = 0.5 )
p2 <- TSNEPlot( object = obj, do.return = T, do.label = TRUE, colors.use = color.cluster, pt.size = 0.5 )
p  <- plot_grid(p1, p2)
ggsave( p, file = "tSNE.pdf", width = 12, height = 6 )

## Draw UMAP plot
p1 <- DimPlot(obj, group.by = "orig.ident", colors.use = color.sample, pt.size = 0.5, reduction = 'umap')
p2 <- DimPlot(obj, colors.use = color.sample, pt.size = 0.5, reduction = 'umap')
p  <- plot_grid(p1, p2)
ggsave( p, file = "UMAP.pdf", width = 12, height = 6 )

### Save data object 
message( "==>Output obj.Rda<==" )
obj <- AddMetaData( object = obj, metadata = obj@ident, col.name = "cluster" )
save( obj, file = "obj.Rda" )

message( "==>All Done!<==" )

