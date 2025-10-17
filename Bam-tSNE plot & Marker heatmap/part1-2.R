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
sample <- c("Bam")
mat <- Read10X(data.dir = paste0("rawdata/",sample), gene.column = 1)
obj <- CreateSeuratObject(counts = mat, project = sample, assay = "RNA" )

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
#nCount_RNA: [ -.inf,20000 ]
standard <- c(-.inf,20000)
cells.use <- obj@meta.data %>% tibble::rownames_to_column(var = "cells") %>% filter(.data[["nCount_RNA"]] >= standard[1] & .data[["nCount_RNA"]] <= standard[2] & cells %in% cells.use) %>% select(cells) %>% unlist()
#nFeature_RNA: [ 200,5000 ]
standard <- c(200,5000)
cells.use <- obj@meta.data %>% tibble::rownames_to_column(var = "cells") %>% filter(.data[["nFeature_RNA"]] >= standard[1] & .data[["nFeature_RNA"]] <= standard[2] & cells %in% cells.use) %>% select(cells) %>% unlist()
#percent.mito: [ -.inf,10 ]
standard <- c(-.inf,10)
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

##根据JackStraw图选取 1:30 PCs
sig.PCs <- 1:30

### tSNE 
obj <- RunTSNE( object = obj, dims.use = sig.PCs, do.fast = TRUE, perplexity = 50, max_iter = 2000, verbose = T )


## Draw t-SNE plot (对应part1)
message( "==>Draw t-SNE plot<==" )
p1 <- TSNEPlot( object = obj, do.return = T, group.by = "orig.ident", colors.use = color.sample, pt.size = 0.5 )
p2 <- TSNEPlot( object = obj, do.return = T, do.label = TRUE, colors.use = color.cluster, pt.size = 0.5 )
p  <- plot_grid(p1, p2)
ggsave( p, file = "tSNE.pdf", width = 12, height = 6 )

### Find clusters
clustering_resolution <- 0.6   # Resolution parameter for Seurat clustering     
obj <- FindClusters( object = obj, reduction.type = "pca", dims.use = sig.PCs, resolution = cluster_resolution, print.output = T, save.SNN = TRUE, force.recalc = TRUE, temp.file.location = getwd() )
color.cluster <- rainbow(length(levels(obj@ident)))
names(color.cluster) <- levels(obj@ident)

### Save data object 
message( "==>Output obj.Rda<==" )
obj <- AddMetaData( object = obj, metadata = obj@ident, col.name = "cluster" )
save( obj, file = "obj.Rda" )


###Find maker genes
min.pct <- 0.25    #minimum percent of gene-expressed cells 
logfc <- 0.25      #log Fold change
pvalue <- 0.01     #pvalue threshold

Idents(object) <- "seurat_clusters"
obj.markers <- FindAllMarkers(object = obj, only.pos = TRUE,min.pct = min.pct, logfc.threshold = logfc,return.thresh = pvalue, pseudocount.use = 0.00000001 )
save( obj.markers, file = "markers.Rda" )

#Find top20 marker
top_num <- 20
top <- object.markers %>% group_by( cluster ) %>%
		arrange(desc(avg_logFC), p_val, p_val_adj, .by_group = TRUE) %>% filter(1:n() <= top_num)

###plot heatmap (对应part2)
p <- DoHeatmap( object = obj, genes.use = unique(top$gene), slim.col.label = TRUE, remove.key = FALSE, do.plot = FALSE )
ggsave(p,filename="Top.Heatmap.pdf",width = length(unique(top$gene)) * 0.11, height = length(unique(top$gene)) * 0.11, limitsize = FALSE)

message( "==>All Done!<==" )

