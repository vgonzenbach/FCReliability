library(circlize)
library(ComplexHeatmap)
library(tidyr)

#### Dense Plots ####
#* ICC ----
dense <- read.csv("results/mean_icc_dense.csv")
ndist <- length(unique(dense$Distance))
dense$Roi <- gsub("ROI", "", dense$Roi)

heat <- NULL
dists <- c()
rois <- c()

for (d in unique(dense$Distance)) {
  df <- dense[dense$Distance == d,]
  heat_sub <- pivot_wider(df[,-1], 
                      names_from = c("Data_length..min."), 
                      values_from = Mean_ICC)
  heat <- rbind(heat, heat_sub[,-1])
  dists <- c(dists, rep(d, nrow(heat_sub)))
  rois <- c(rois, heat_sub$Roi)
}
heat <- as.matrix(heat)
rownames(heat) <- rois
dists = factor(dists, levels = rev(c("correlation", "cosine", "covariance", "euclidean", "sqeuclidean", "tangent", "partial_correlation")))

# adapted from https://jokergoo.github.io/2020/05/21/make-circular-heatmaps/
png("figures/dense_ICC_circle.png", units="in", width=5, height=5, res=300)

circos.clear()
circos.par(gap.after = c(rep(2, ndist-1), 20))

col_fun <- colorRamp2(seq(0,1,.1), RColorBrewer::brewer.pal(11, 'Spectral'))
##text(0, 0, 'rownames.side = "inside"')
circos.heatmap(heat[,10:1], split = dists, col = col_fun,
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
lgd = Legend(title = "ICC", title_position = 'topcenter', col_fun = col_fun)
grid.draw(lgd)

dev.off()

#* Discriminability ----
dense <- read.csv("results/dis_idr_sub.csv")
ndist <- length(unique(dense$Distance))

heat <- NULL
dists <- c()
rois <- c()

for (d in unique(dense$Distance)) {
  df <- dense[dense$Distance == d,]
  heat_sub <- pivot_wider(df[,2:4], 
                          names_from = Session, 
                          values_from = Sub_Discriminability)
  heat <- rbind(heat, heat_sub[,-1])
  dists <- c(dists, rep(d, nrow(heat_sub)))
  rois <- c(rois, heat_sub$ROI)
}
heat <- as.matrix(heat)
rownames(heat) <- rois
dists = factor(dists, levels = rev(c("correlation", "cosine", "covariance", "euclidean", "sqeuclidean", "tangent", "partial_correlation")))

# Plot
png("figures/dense_discrim_circle.png", units="in", width=5, height=5, res=300)

circos.clear()
circos.par(gap.after = c(rep(2, ndist-1), 20))

col_fun <- colorRamp2(seq(0,1,.1), RColorBrewer::brewer.pal(11, 'Spectral'))
#text(0, 0, 'rownames.side = "inside"')
circos.heatmap(heat[,10:1], split = dists, col = col_fun,
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
lgd = Legend(title = "Discriminability", title_position = 'topcenter', col_fun = col_fun)
grid.draw(lgd)

dev.off()

#* Identification ----
dense <- read.csv("results/dis_idr_sub.csv")
ndist <- length(unique(dense$Distance))

heat <- NULL
dists <- c()
rois <- c()
for (d in unique(dense$Distance)) {
  df <- dense[dense$Distance == d,]
  heat_sub <- pivot_wider(df[,c(2:3, 5)], 
                          names_from = Session, 
                          values_from = Sub_Identification)
  heat <- rbind(heat, heat_sub[,-1])
  dists <- c(dists, rep(d, nrow(heat_sub)))
  rois <- c(rois, heat_sub$ROI)
}
heat <- as.matrix(heat)
rownames(heat) <- rois
dists = factor(dists, levels = rev(c("correlation", "cosine", "covariance", "euclidean", "sqeuclidean", "tangent", "partial_correlation")))

# Plot
png("figures/dense_identify_circle.png", units="in", width=5, height=5, res=300)

circos.clear()
circos.par(gap.after = c(rep(2, ndist-1), 20))

col_fun <- colorRamp2(seq(0,1,.1), RColorBrewer::brewer.pal(11, 'Spectral'))
#text(0, 0, 'rownames.side = "inside"')
circos.heatmap(heat[,10:1], split = dists, col = col_fun,
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
lgd = Legend(title = "Identifiability", title_position = 'topcenter', col_fun = col_fun)
grid.draw(lgd)

dev.off()

#### Sparse Plots ####
#* ICC ----
dense <- read.csv("results/mean_icc_sparse.csv")
ndist <- length(unique(dense$Distance))
dense$Roi <- gsub("ROI", "", dense$Roi)

heat <- NULL
dists <- c()
rois <- c()
for (d in unique(dense$Distance)) {
  df <- dense[dense$Distance == d,]
  heat_sub <- pivot_wider(df[,-1], 
                          names_from = c("Data_length..min."), 
                          values_from = Mean_ICC)
  heat <- rbind(heat, heat_sub[,-1])
  dists <- c(dists, rep(d, nrow(heat_sub)))
  rois <- c(rois, heat_sub$Roi)
}
heat <- as.matrix(heat)
rownames(heat) <- rois
dists = factor(dists, levels = rev(c("correlation", "cosine", "covariance", "euclidean", "sqeuclidean", "tangent", "partial_correlation")))

# adapted from https://jokergoo.github.io/2020/05/21/make-circular-heatmaps/
png("figures/sparse_ICC_circle.png", units="in", width=5, height=5, res=300)

circos.clear()
circos.par(gap.after = c(rep(2, ndist-1), 20))

col_fun <- colorRamp2(seq(0,1,.1), RColorBrewer::brewer.pal(11, 'Spectral'))
#text(0, 0, 'rownames.side = "inside"')
circos.heatmap(heat[,10:1], split = dists, col = col_fun,
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
lgd = Legend(title = "ICC", title_position = 'topcenter', col_fun = col_fun)
grid.draw(lgd)

dev.off()

#* Discriminability ----
dense <- read.csv("results/dis_idr_sub_sparse.csv")
ndist <- length(unique(dense$Distance))

heat <- NULL
dists <- c()
rois <- c()
for (d in unique(dense$Distance)) {
  df <- dense[dense$Distance == d,]
  heat_sub <- pivot_wider(df[,2:4], 
                          names_from = Session, 
                          values_from = Sub_Discriminability)
  heat <- rbind(heat, heat_sub[,-1])
  dists <- c(dists, rep(d, nrow(heat_sub)))
  rois <- c(rois, heat_sub$ROI)
}
heat <- as.matrix(heat)
rownames(heat) <- rois
dists = factor(dists, levels = rev(c("correlation", "cosine", "covariance", "euclidean", "sqeuclidean", "tangent", "partial_correlation")))

# Plot
png("figures/sparse_discrim_circle.png", units="in", width=5, height=5, res=300)

circos.clear()
circos.par(gap.after = c(rep(2, ndist-1), 20))

col_fun <- colorRamp2(seq(0,1,.1), RColorBrewer::brewer.pal(11, 'Spectral'))
#text(0, 0, 'rownames.side = "inside"')
circos.heatmap(heat[,10:1], split = dists, col = col_fun,
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
lgd = Legend(title = "Discriminability", title_position = 'topcenter', col_fun = col_fun)
grid.draw(lgd)

dev.off()

#* Identification ----
dense <- read.csv("results/dis_idr_sub_sparse.csv")
ndist <- length(unique(dense$Distance))

heat <- NULL
dists <- c()
rois <- c()
for (d in unique(dense$Distance)) {
  df <- dense[dense$Distance == d,]
  heat_sub <- pivot_wider(df[,c(2:3, 5)], 
                          names_from = Session, 
                          values_from = Sub_Identification)
  heat <- rbind(heat, heat_sub[,-1])
  dists <- c(dists, rep(d, nrow(heat_sub)))
  rois <- c(rois, heat_sub$ROI)
}
heat <- as.matrix(heat)
rownames(heat) <- rois
dists = factor(dists, levels = rev(c("correlation", "cosine", "covariance", "euclidean", "sqeuclidean", "tangent", "partial_correlation")))

# Plot
png("figures/sparse_identify_circle.png", units="in", width=5, height=5, res=300)

circos.clear()
circos.par(gap.after = c(rep(2, ndist-1), 20))

col_fun <- colorRamp2(seq(0,1,.1), RColorBrewer::brewer.pal(11, 'Spectral'))
#text(0, 0, 'rownames.side = "inside"')
circos.heatmap(heat[,10:1], split = dists, col = col_fun,
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
lgd = Legend(title = "Identifiability", title_position = 'topcenter', col_fun = col_fun)
grid.draw(lgd)

dev.off()

