evaluate_clustering_by_group <- function(unit.compare,
                                         anal.groups = NULL,
                                         clst_method = "average",
                                         siteunit_col = "SiteUnit",
                                         sppgroup_col = "Spp.group") {
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
  
  if (is.null(anal.groups)) {
    # Evaluate clustering across all units
    all_units <- union(unit.compare$Unit1, unit.compare$Unit2)
    if (length(all_units) < 2) return(NULL)
    
    clust_eval <- tryCatch({
      evaluate_clustering(compare = unit.compare, clst_method = clst_method) %>%
        mutate(Spp.group = "All_Units",
               n_SiteUnits_in_compare = length(all_units),
               n_comparisons = nrow(unit.compare))
    }, error = function(e) {
      message(sprintf("⚠️ Skipping All_Units due to error: %s", e$message))
      NULL
    })
    
    if (!is.null(clust_eval)) clustering_stats[["All_Units"]] <- clust_eval
    
  } else {
    # Clean SiteUnit names
    anal.groups <- anal.groups %>%
      mutate(!!siteunit_col := gsub("/", "_", gsub(" ", "", .data[[siteunit_col]])))
    
    spp.groups <- unique(anal.groups[[sppgroup_col]])
    
    for (grp in spp.groups) {
      unit.choose <- anal.groups %>%
        filter(.data[[sppgroup_col]] == grp) %>%
        pull(.data[[siteunit_col]]) %>%
        unique()
      
      units.in.compare <- unit.choose[unit.choose %in% union(unit.compare$Unit1, unit.compare$Unit2)]
      
      if (length(units.in.compare) < 2) next
      
      unit_subset <- unit.compare %>%
        filter(Unit1 %in% units.in.compare & Unit2 %in% units.in.compare)
      
      if (nrow(unit_subset) == 0) next
      
      clust_eval <- tryCatch({
        evaluate_clustering(compare = unit_subset, clst_method = clst_method) %>%
          mutate(Spp.group = grp,
                 n_SiteUnits_in_compare = length(units.in.compare),
                 n_comparisons = nrow(unit_subset))
      }, error = function(e) {
        message(sprintf("⚠️ Skipping %s due to error: %s", grp, e$message))
        NULL
      })
      
      if (!is.null(clust_eval)) clustering_stats[[grp]] <- clust_eval
    }
  }
  
  result <- do.call(rbind, clustering_stats) %>% arrange(desc(CCC))
  return(result)
}
