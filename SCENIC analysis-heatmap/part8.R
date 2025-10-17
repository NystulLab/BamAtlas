library(ComplexHeatmap)
library(SCENIC)

scenicOptions <- readRDS(scenicOptions_file)

aucMat = AUCell::getAUC(loadInt(scenicOptions, 'aucell_regulonAUC'))
binaryMat = loadInt(scenicOptions, 'aucell_binary_full')
cellInfo = loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"))

#draw AUC heatmap
avg_mat <- do.call(cbind, by(rownames(cellInfo), cellInfo$group, 
        function(x) {
            Matrix::rowMeans(aucMat[, x, drop = F])
        }, simplify = FALSE))

hmp <- Heatmap(avg_mat, 
        col = NULL,
        name = "Regulon activity",
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T
        )

pdf("RegulonHeatmap.all.Avg.bygroup.pdf",width = 7, height = 8)
draw(hmp)
dev.off()

#draw binary heatmap
avg_mat <- do.call(cbind, by(rownames(cellInfo), cellInfo$group,   
        function(x) {
            Matrix::rowMeans(binaryMat[, x, drop = F])
        }, simplify = FALSE))

col <- colorRamp2(c(0, 1), c("white", "black"))
hmp <- Heatmap(avg_mat,
        col = col,
        name = "Binary activity",
        cluster_rows = T, 
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T
        )
pdf("BinaryHeatmap.all.Avg.bygroup.pdf",width = 7, height = 8)
draw(hmp)
dev.off()
