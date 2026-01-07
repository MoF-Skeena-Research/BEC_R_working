
draw_dendro_split <- function(unit.compare, subass.level = .1, assoc.level = .2, alliance.level = .4, cut.level = .5){
  singles.count = 0
  singles.list = data.frame(SiteUnit = character(), stringsAsFactors = FALSE)
  new.unit <- data.frame(SiteUnit = character(), stringsAsFactors = FALSE)
  compared <- unit.compare 
  dis.matrix <- bec_dist_matrix(compared) 
  ss_clst <- agnes(dis.matrix,
                   diss = TRUE, stand = TRUE,
                   method = "average")
  dendro_hc <- as.hclust(ss_clst)
  dendro_hc.dend <- as.dendrogram(dendro_hc)
  if (!is.null(cut)){
    dendro_hc.dend <- cut(dendro_hc.dend, h = cut.level)
  }
  n <- length(dendro_hc.dend$lower)
  for (i in 1:n){
    dendro_hc.dend.i <- dendro_hc.dend$lower[[i]]
    n2 <- length(dendro_hc.dend.i)
    if (n2 == 1) {singles.count = singles.count + 1}
    if (n2 == 1) {singles.list <- rbind(singles.list, x = setNames(as.data.frame(partition_leaves(dendro_hc.dend.i)), names(singles.list)))}
    if (n2 == 1) next
    # plot(dendro_hc.dend.i, main = paste0("Cluster Dendrogram of Site Units: branch ", i, " at hcut = ", h))
    # }
    #   dend.co <- stats::cophenetic(dendro_hc.i)
    # dend.dis <- as.dist(dis.matrix)
    # cophenetic <- cor(dend.dis, dend.co) %>% round(2)
    # cophenetic ## shows how well the clusters align with the data >0.7 is considered good
    # coph_annotation <- data.frame(x = 9, y=.85, label = paste0("Cophenetic:",cophenetic))
    hcdata <- dendro_data(dendro_hc.dend.i, type = "rectangle")
    yy <- ggplot() +
      # Draw cluster segments
      geom_segment(data = segment(hcdata), 
                   aes(x = x, y = y, xend = xend, yend = yend)) +
      
      # Label leaves
      geom_text(data = label(hcdata), 
                aes(x = x, y = y, label = label, hjust = 0), 
                size = 3) +
      
      # Add shaded band between 0.07 and 0.10
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = subass.level - .03, ymax = subass.level,
               fill = "red", alpha = 0.15) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = assoc.level - .02, ymax = assoc.level + .02,
               fill = "darkgreen", alpha = 0.15) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = alliance.level - .02, ymax = alliance.level + .02,
               fill = "purple", alpha = 0.15) +
      # Additional horizontal lines
      # geom_hline(yintercept = assoc.level + .02, linetype = "dashed", color = "green") +
      # geom_hline(yintercept = assoc.level - .02, linetype = "dashed", color = "green") +
      # geom_hline(yintercept = subass.level - .03, linetype = "dashed", color = "red") +
      # geom_hline(yintercept = subass.level, linetype = "dashed", color = "red") +
      # geom_hline(yintercept = alliance.level +2, linetype = "dashed", color = "purple") +
      # geom_hline(yintercept = alliance.level -2, linetype = "dashed", color = "purple") +
      # Add labels
      geom_text(aes(x = 0, y = subass.level - .015, label = "Subassociation", hjust = 0),
                angle = 90, color = "grey30", size = 3) +
      geom_text(aes(x = 0, y = assoc.level, label = "Association", hjust = 0),
                angle = 90, color = "grey30", size = 3) +
      geom_text(aes(x = 0, y = alliance.level, label = "Alliance", hjust = 0),
                angle = 90, color = "grey30", size = 3)+
      # add label to graph 
      # geom_text(data=coph_annotation, aes( x=x, y=y, label=label), 
      #           color="black", 
      #           size=3 , angle=0, fontface="bold" ) +
      # annotate(cophenetic, x = .1, y = .1,
      #      label = "Cophonetic" , color="orange",
      #       size=7 , angle=0, fontface="bold")+
      coord_flip()+
      scale_y_reverse(limits = c(1, -.3))+
      labs(x = "", y = "Difference")+
      theme_minimal()+
      theme(axis.text.y=element_blank(), axis.title.y = element_blank())+
      ggtitle(paste0("Cluster Dendrogram of Site Units: branch ", i, " at hcut = ", cut.level))
    print(yy)
  }
  #print(paste0("The total number of site units is ", ss.count))
  print(paste0("The number of singles at hcut ", cut.level, " is ", singles.count))
 return(singles.list)
}
