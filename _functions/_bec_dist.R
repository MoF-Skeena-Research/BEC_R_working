bec_dist <- function(vegsum.pairs) {
  dis.matrix1 <- as.data.frame(vegsum.pairs) |>
    dplyr::mutate(diss = 1 - diff.ratio) |>
    dplyr::mutate(Unit1 = as.character(Unit1)) |>
    dplyr::select(Unit1, Unit2, diss) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(Unit1), !is.na(Unit2)) |>
    dplyr::arrange(Unit1)
  ## go from one way to two way matrix for clustering
  dis.matrix2 <- dis.matrix1 |>
    dplyr::rename(Unit1 = Unit2, Unit2 = Unit1) |>
    dplyr::select(Unit1, Unit2, diss) |>
    dplyr::arrange(Unit1)
  dis.matrix <- rbind(dis.matrix1, dis.matrix2) |>
    dplyr::arrange(Unit1, Unit2) |>
    tidyr::pivot_wider(id_cols = Unit1, names_from = Unit2, values_from = diss,
                names_sort = TRUE) |>
    tibble::column_to_rownames("Unit1") |>
    dplyr::mutate_all(~ replace(., is.na(.), 1)) |>
    as.matrix()
  bec_dist <- as.dist(dis.matrix)
  attr(bec_dist, "method") <- "bec"
  return(bec_dist)
}
