# https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/

dendro_data_k <- function(dendro, k) {
    #adapted from https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
  #dendro is an ape::dendrogram object
  #k is an integer
  dendr <- ggdendro::dendro_data(dendro, type = "rectangle")
  
  seg = dendr$segments
  labclust = dendextend::cutree_1k.dendrogram(dendro, k = k, order_clusters_as_data = F)
  segclust = rep(0L, nrow(seg))
  heights = sort(dendextend::get_nodes_attr(dendro, "height"), decreasing = TRUE)
  height =  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)

  for (i in 1:k) {
      xi      <-  dendr$labels$x[labclust == i]
      idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
      idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
      idx3    <-  seg$yend < height
      idx     <-  idx1 & idx2 & idx3
      segclust[idx] <- i
  }

  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  dendr$segments$clust  <-  segclust
  dendr$segments$line   <-  as.integer(segclust < 1L)
  dendr$labels$clust    <-  labclust

  return(dendr)
}
