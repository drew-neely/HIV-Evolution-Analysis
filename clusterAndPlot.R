install.packages("tidyverse")
install.packages("cluster")
install.packages("factoextra")
install.packages("dendextend")
install.packages("ape")

library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
library(ape)
library(data.table)

# matrix <- as.matrix(read.csv(file='data/dissimilarityMatrix.csv', header=F))
matrix<-read.table("disimilarity_matrix.txt", sep=" ", na.strings = "" ,fill=TRUE)
colnames(matrix) <- c("1-25","1-33","1-40","1-53","1-83","1-112","1-140","1-196","1-243","2-7","2-14","2-21","2-44","2-68","2-97","2-127","2-181","3-3","3-12","3-15","3-72","3-99","3-129","3-197","4-6","4-12","4-21","4-42","4-48","4-57","4-69","4-98","4-154","4-210","4-273","4-350","5-0","5-21","5-35","5-70","5-98","5-455","6-0","6-21","6-36","6-63","6-140","6-805","7-0","7-21","7-36","7-70","7-168")
rownames(matrix) <- c("1-25","1-33","1-40","1-53","1-83","1-112","1-140","1-196","1-243","2-7","2-14","2-21","2-44","2-68","2-97","2-127","2-181","3-3","3-12","3-15","3-72","3-99","3-129","3-197","4-6","4-12","4-21","4-42","4-48","4-57","4-69","4-98","4-154","4-210","4-273","4-350","5-0","5-21","5-35","5-70","5-98","5-455","6-0","6-21","6-36","6-63","6-140","6-805","7-0","7-21","7-36","7-70","7-168")

patientMatrices = list()
patientClusterings = list()
patientPhylo = list()
for (i in 1:7) {
        basePatientMatrix = matrix[rownames(matrix) %like% paste0(i, "-"), colnames(matrix) %like% paste0(i, "-")]
        patientMatrices[[i]] = as.dist(basePatientMatrix)
        patientClusterings[[i]] = hclust(patientMatrices[[i]], method = "complete")
        patientPhylo[[i]] = as.phylo(patientClusterings[[i]])
}

dissimilarityMatrix <- as.dist(matrix)

# Hierarchical clustering using Complete Linkage
hc <- hclust(dissimilarityMatrix, method = "complete" )

# Plot the obtained dendrogram
# Note: Can specify `hang = -1` to put all labels at the bottom of the chart
plot(hc)



plot(as.phylo(hc), cex = 0.6, label.offset = 0.5)

plot(as.phylo(hc), type = "cladogram", cex = 0.6, 
     label.offset = 0.5)

plot(as.phylo(hc), type = "unrooted", cex = 0.6,
     no.margin = TRUE)

plot(as.phylo(hc), type = "fan")

plot(as.phylo(hc), type = "radial")

# cluster




colors = c("red", "blue", "green", "purple", "orange", "deeppink", "brown", "orange")
clus4 = cutree(hc, 7)
clus4 = cutree(patientClusterings[[1]], 4)

ph = as.phylo(hc)
ph[["edge.length"]] = log2((ph[["edge.length"]]) * 1000000 + 1.2) ** 2

plot(patientPhylo[[1]], type = "phylogram",
     label.offset = 0.001, cex = 1.3, use.edge.length=TRUE)

plot(ph, type = "unrooted", tip.color = colors[clus4],
     label.offset = 0.1, cex = 1, use.edge.length=TRUE)





plot(as.phylo(hc), type = "phylogram", tip.color = colors[clus4],
     label.offset = 0.001, cex = 1, use.edge.length=TRUE)

plot(ph, type = "cladogram", tip.color = colors[clus4],
     label.offset = 0.001, cex = 1, use.edge.length=TRUE)

plot(ph, type = "fan", tip.color = colors[clus4],
     label.offset = 0.001, cex = 1, use.edge.length=TRUE)

plot(ph, type = "radial", tip.color = colors[clus4],
     label.offset = 0.001, cex = 1, use.edge.length=TRUE)

plot(ph, type = "unrooted", tip.color = colors[clus4],
     label.offset = 0.001, cex = 1, use.edge.length=TRUE)





plot(as.phylo(hc), type = "fan", tip.color = colors[clus4],
     label.offset = 0.001, cex = 0.7)

plot(as.phylo(hc), type = "cladogram", cex = 0.6,
     edge.color = "steelblue", edge.width = 2, edge.lty = 2,
     tip.color = "steelblue")









