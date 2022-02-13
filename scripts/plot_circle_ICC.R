setwd("/Users/andrewac/Documents/GitHub/andy1764-Experiments/Collaborations/Aki-FCdists")

library(circlize)
library(ComplexHeatmap)
library(tidyr)

#### Dense ####
dense <- read.csv("results/dense_mean_icc.csv")
ndist <- length(unique(dense$Distance))

heat <- NULL
dists <- c()
rois <- c()
for (d in unique(dense$Distance)) {
  df <- dense[dense$Distance == d,]
  heat_sub <- pivot_wider(df[,-1], 
                      names_from = c("Data_length"), 
                      values_from = Mean_ICC)
  heat <- rbind(heat, heat_sub[,-1])
  dists <- c(dists, rep(d, nrow(heat_sub)))
  rois <- c(rois, heat_sub$Roi)
}
heat <- as.matrix(heat)
rownames(heat) <- rois

# adapted from https://jokergoo.github.io/2020/05/21/make-circular-heatmaps/
png("figures/dense_ICC_circle.png", units="in", width=5, height=5, res=300)

circos.clear()
circos.par(gap.after = c(rep(2, ndist-1), 20))

col_fun <- colorRamp2(c(min(dense$Mean_ICC), max(dense$Mean_ICC)), c("white", "blue"))
text(0, 0, 'rownames.side = "inside"')
circos.heatmap(heat[,10:1], split = as.factor(dists), col = col_fun,
               track.height = 0.4,
               cluster = FALSE,
               rownames.side = "inside",
               show.sector.labels = TRUE)

# draw y-axis labels
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == ndist) { # the last sector
    cn = colnames(heat)
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                1:n - 0.5, cn, 
                cex = 0.5, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)

# draw legend
lgd = Legend(title = "ICC", col_fun = col_fun)
grid.draw(lgd)

dev.off()

#### Sparse ####
sparse <- read.csv("results/sparse_mean_icc.csv")
ndist <- length(unique(sparse$Distance))

heat <- NULL
dists <- c()
rois <- c()
for (d in unique(sparse$Distance)) {
  df <- sparse[sparse$Distance == d,]
  heat_sub <- pivot_wider(df[,-1], 
                          names_from = c("Data_length"), 
                          values_from = Mean_ICC)
  heat <- rbind(heat, heat_sub[,-1])
  dists <- c(dists, rep(d, nrow(heat_sub)))
  rois <- c(rois, heat_sub$Roi)
}
heat <- as.matrix(heat)
rownames(heat) <- rois

# adapted from https://jokergoo.github.io/2020/05/21/make-circular-heatmaps/
png("figures/sparse_ICC_circle.png", units="in", width=5, height=5, res=300)

circos.clear()
circos.par(gap.after = c(rep(2, ndist-1), 20))

col_fun <- colorRamp2(c(min(dense$Mean_ICC), max(dense$Mean_ICC)), c("white", "blue"))
text(0, 0, 'rownames.side = "inside"')
circos.heatmap(heat[,10:1], split = as.factor(dists), col = col_fun,
               track.height = 0.4,
               cluster = FALSE,
               rownames.side = "inside",
               show.sector.labels = TRUE)

# draw y-axis labels
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == ndist) { # the last sector
    cn = colnames(heat)
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                1:n - 0.5, cn, 
                cex = 0.5, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)

# draw legend
lgd = Legend(title = "ICC", col_fun = col_fun)
grid.draw(lgd)

dev.off()
