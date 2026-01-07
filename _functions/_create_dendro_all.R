create_dendro_all <- function(unit.compare, threshold.low = .1, cut.level = .2, method = "average", parameter = 1){
  compared <- unit.compare 
  dis.matrix <- bec_dist_matrix(compared) 
  ss_clst <- agnes(dis.matrix,
                   diss = TRUE, stand = TRUE,
                   method = method, par.method = parameter)
                   #method = "average")
  dendro_hc <- as.hclust(ss_clst)
  dend.co <- stats::cophenetic(dendro_hc)
  dend.dis <- as.dist(dis.matrix)
  cophenetic <- cor(dend.dis, dend.co) %>% round(2)
  cophenetic ## 0.891 this value shows how well the clusters align with the data
  # ### >0.7 is considered good
  coph_annotation <- data.frame(x = 9, y=.85, label = paste0("Cophenetic:",cophenetic))
  hcdata <- ggdendro::dendro_data(dendro_hc, type = "rectangle")
  yy <- ggplot() +
    geom_segment(data = segment(hcdata), 
                 aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    geom_text(data = label(hcdata), 
              aes(x = x, y = y, label = label, hjust = 0), 
              size = 3
    ) +
    geom_hline(yintercept = threshold.low, linetype = "dashed", color = "red")+
    
    
    geom_hline(yintercept = cut.level, linetype = "dashed", color = "darkgreen")+
    # add value label to hline
    geom_text(aes(x = 0, y = (threshold.low+0.02), label = paste0((threshold.low*100),"%"), hjust = 0), angle = 90,color = "red", size = 3)+
    geom_text(aes(x = 0, y = (cut.level+0.02), label = paste0((cut.level*100),"%"), hjust = 0), angle = 90, color = "darkgreen", size = 3)+
    
    # add label to graph 
    # geom_text(data=coph_annotation, aes( x=x, y=y, label=label), 
    #           color="black", 
    #           size=3 , angle=0, fontface="bold" ) +
    labs(title = paste0("Cluster Dendrogram of Site Units: "),
         subtitle = paste0("Cophenetic Correlation: ", round(cophenetic, 3), 
                           ", Agglomerative Coefficient: ", round(ss_clst$ac, 3))) +
    coord_flip()+
    scale_y_reverse(limits = c(1, -.3))+
    labs(x = "", y = "Difference")+
    theme_minimal()+
    theme(axis.text.y=element_blank(), axis.title.y = element_blank())#+
    #ggtitle(paste0("Cluster Dendrogram of Site Units") )
  print(yy)
  return(dendro_hc)
  
}
