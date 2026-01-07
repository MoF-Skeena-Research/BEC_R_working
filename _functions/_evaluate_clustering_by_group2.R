evaluate_clustering_by_group <- function(unit.compare,
                                         clst_method = "average",
                                         grouping_table = NULL,
                                         siteunit_col = "SiteUnit",
                                         analgroup_col = "Order") {
  library(dplyr)
  library(cluster)
  library(data.table)
  
  # Helper: compute distance matrix from pairwise similarity
  bec_dist_matrix <- function(compare) {
    units <- unique(c(compare$Unit1, compare$Unit2))
    mat <- matrix(1, nrow = length(units), ncol = length(units),
                  dimnames = list(units, units))
    for (i in seq_len(nrow(compare))) {
      u1 <- compare$Unit1[i]
      u2 <- compare$Unit2[i]
      sim <- compare$BEC.sim.min[i]
      mat[u1, u2] <- sim
      mat[u2, u1] <- sim
    }
    return(1 - mat)  # convert similarity to dissimilarity
  }
  
  # Helper: evaluate clustering
  evaluate_clustering <- function(compare, clst_method = "average") {
    dis.matrix <- bec_dist_matrix(compare)
    dis.dis <- as.dist(dis.matrix)
    ss_clst <- agnes(dis.matrix, diss = TRUE, stand = TRUE, method = clst_method)
    CCC <- cor(dis.dis, cophenetic(ss_clst))
    AC <- round(ss_clst$ac, 3)
    data.frame(
      clst_method = clst_method,
      AC = AC,
      CCC = round(CCC, 3),
      stringsAsFactors = FALSE
    )
  }
  
  clustering_stats <- list()
  
  if (is.null(grouping_table)) {
    # Evaluate clustering across all units
    all_units <- union(unit.compare$Unit1, unit.compare$Unit2)
    if (length(all_units) < 2) return(NULL)
    
    clust_eval <- tryCatch({
      #evaluate_clustering(compare = unit.compare, clst_method = clst_method) %>%
        clust_eval <- evaluate_clustering(compare = unit.compare, clst_method = clst_method)
        clust_eval$Spp.group <- "All Units"
        clust_eval$n_SiteUnits<- length(all_units)
        clust_eval$n_comparisons <- nrow(unit.compare)
    }, error = function(e) {
      message(sprintf("⚠️ Skipping All_Units due to error: %s", e$message))
      NULL
    })
    
    if (!is.null(clust_eval)) clustering_stats[["All_Units"]] <- clust_eval
    
  } else {
    # Clean SiteUnit names
    grouping_table <- grouping_table %>%
      mutate(!!siteunit_col := gsub("/", "_", gsub(" ", "", .data[[siteunit_col]])))
    
    spp.groups <- unique(grouping_table[[analgroup_col]])
    grp = "MH"
    for (grp in spp.groups) {
      units.choose <- grouping_table %>%
        filter(.data[[analgroup_col]] == grp) %>%
        pull(.data[[siteunit_col]]) %>%
        unique()
      
      #units.in.compare <- unit.choose[unit.choose %in% union(unit.compare$Unit1, unit.compare$Unit2)]
      
      if (length(units.choose) < 2) next
      
      unit_subset <- unit.compare %>%
        filter(Unit1 %in% units.choose) %>% filter(Unit2 %in% units.choose)
      
      if (nrow(unit_subset) == 0) next
      
      clust_eval <- tryCatch({
        clust_eval <- evaluate_clustering(compare = unit_subset, clst_method = clst_method) %>% 
          mutate(Spp.group = grp, 
                 n_SiteUnits = length(units.choose),
                 n_comparisons = nrow(unit_subset))
        # clust_eval$Spp.group <- grp
        # clust_eval$n_SiteUnits <- length(units.choose)
        # clust_eval$n_comparisons <- nrow(unit_subset)
      }, error = function(e) {
        message(sprintf("⚠️ Skipping %s due to error: %s", grp, e$message))
        NULL
      })
      
      if (!is.null(clust_eval)) clustering_stats[[grp]] <- clust_eval
    }
  }
  
  result <- do.call(rbind, clustering_stats)# %>% arrange(desc(CCC))
  return(result)
}
