library(pheatmap)

data <- read.table("Cluster.sample.stat.xls",header = T, sep = "\t", check.names = F, row.names = 1)

#计算ROE
ROE <- function(data, samples = NULL) {
        r <- chisq.test(data)
        roe <- r$observed/r$expected
        result <- cbind(Cluster = rownames(roe), roe)
        return(result)
}

roe_result <- data.frame(ROE(data))
roe_row <- roe_result$Cluster
roe_result$Cluster <- NULL

#绘制热图
roe_result <- apply(roe_result,2,as.numeric)
rownames(roe_result) <- roe_row

pheatmap(roe_result,filename="ROE.heatmap.pdf",cluster_rows = F, cluster_cols = F, scale = 'row')
