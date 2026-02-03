create_dendro_bybgc <- function(bgc.choose, unit.compare,
                                threshold.low = .1, threshold.high = .2){
  
  # Filter to chosen BGC
  compared <- unit.compare %>%
    filter(bgc1 == bgc2, bgc1 == bgc.choose)
  
  # Distance matrix
  dis.matrix <- bec_dist_matrix(compared)
  
  # Clustering
  ss_clst <- agnes(dis.matrix, diss = TRUE, stand = TRUE, method = "average")
  dendro_hc <- as.hclust(ss_clst)
  
  # Cophenetic correlation
  dend.co <- stats::cophenetic(dendro_hc)
  dend.dis <- as.dist(dis.matrix)
  cophenetic <- round(cor(dend.dis, dend.co), 2)
  
  # Annotation
  coph_annotation <- data.frame(
    x = 9, y = .85,
    label = paste0("Cophenetic:", cophenetic)
  )
  
  # Dendrogram data
  hcdata <- dendro_data(dendro_hc, type = "rectangle")
  
  # ---- Similarity table ----
  sim_table <- as.data.frame(as.table(dis.matrix)) %>%
    rename(unit1 = Var1, unit2 = Var2, dissimilarity = Freq) %>%
    filter(unit1 != unit2) %>%
    filter(dissimilarity < threshold.high)
  
  # Groups from cutree
  groups <- cutree(dendro_hc, h = threshold.high) %>%
    tibble(SiteUnit = names(.),
           working_unit = paste0(bgc.choose, "_", as.integer(.)))
  
  # Add group info to similarity table
  sim_table <- sim_table %>%
    left_join(groups, by = c("unit1" = "SiteUnit")) %>%
    rename(group_unit1 = working_unit) %>%
    left_join(groups, by = c("unit2" = "SiteUnit")) %>%
    rename(group_unit2 = working_unit)
  

  
  
  # ---- Plot ----
  yy <- ggplot() +
    geom_segment(data = segment(hcdata),
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = label(hcdata),
              aes(x = x, y = y, label = label, hjust = 0),
              size = 3) +
    geom_hline(yintercept = threshold.low, linetype = "dashed", color = "red") +
    geom_hline(yintercept = threshold.high, linetype = "dashed", color = "darkgreen") +
    geom_text(aes(x = 0, y = threshold.low + 0.02,
                  label = paste0(threshold.low * 100, "%"), hjust = 0),
              angle = 90, color = "red", size = 3) +
    geom_text(aes(x = 0, y = threshold.high + 0.02,
                  label = paste0(threshold.high * 100, "%"), hjust = 0),
              angle = 90, color = "darkgreen", size = 3) +
    geom_text(data = coph_annotation,
              aes(x = x, y = y, label = label),
              color = "black", size = 3, fontface = "bold") +
    coord_flip() +
    scale_y_reverse(limits = c(1, -.3)) +
    labs(x = "", y = "Dissimilarity") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle(paste0("Cluster Dendrogram of ", bgc.choose, " Site Series"))
  
  # Print plot
  print(yy)
  
  # Return both plot and tables
  list(
    #plot = yy,
    similar_units = sim_table,
    groups = groups
  )
}
