cluster_groups <- function(unit.compare, cut.level = .2, minclus = 2, group.label = "Association"){
  compared <- unit.compare 
  dis.matrix <- bec_dist_matrix(compared) 
  ss_clst <- agnes(dis.matrix,
                   diss = TRUE, stand = TRUE,
                   method = "average")
  dendro_hc <- as.hclust(ss_clst)
  groups <- cutree(dendro_hc, h = cut.level) %>% as.data.frame() %>% rename(!!group.label := 1)
  # dendro_hc.dend <- as.dendrogram(dendro_hc)
  # dendro_hc.dend <- cut(dendro_hc.dend, h = cut.level)
  # 
  # dend.co <- stats::cophenetic(dendro_hc)
  # dend.dis <- as.dist(dis.matrix)
  # cluster_grps <- cutreeHybrid(dendro_hc, distM = dis.matrix, cutHeight = cut.level,
  #                              minClusterSize = minclus, deepSplit = 0)
  #groups <- cluster_grps$labels
  su_grps <- cbind(as.data.frame(groups), as.data.frame(row.names(dis.matrix))) %>% rename(SiteUnit = 2 )
  return( su_grps)
}
