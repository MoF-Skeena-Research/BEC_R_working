## pair-wise comparison of vegetation summary data

# su = su2
# minimportance = 0
# minconstancy = 60
# noiseconstancy = 10
# minplots = 5
# minor = 1
# use.ksi = TRUE; ksi = key.site.indicators; ksi.value = 1.2
# reduce.lifeform = TRUE; reduced.lifeforms = reduced.lifeforms; reduction = .1
# reduced.exceptions = NULL

do_pairwise <- function(veg.dat, su, minimportance = 0, minconstancy = 60,
                        noiseconstancy = 10,
                        minplots = 5, 
                        use.ksi = FALSE, ksi = NULL, ksi.value = 1.5,
                        reduce.lifeform = FALSE, reduced.lifeforms = NULL,
                        reduction = 1, reduced.exceptions = NULL,
                        minor = 1.0) {
  ### ---Create vegetation summary
  tic()
  su.choice <- su %>% select(SiteUnit)
  vegdat <- as.data.table(veg.dat)
  vegdat[su, SiteUnit := i.SiteUnit, on = "PlotNumber"] ## limit data to those listed in SiteUnit
  vegdat <- vegdat[!is.na(SiteUnit) & SiteUnit != "", ]
  vegdat <- unique(vegdat[!is.na(SiteUnit) & SiteUnit != "", ])
  vegdat3 <- vegdat[, if (.N > 1) .SD, by = .(SiteUnit, Species)]
  vegdat3[, nplots := length(unique(PlotNumber)), by = .(SiteUnit)]
  vegsum <- vegdat3[, .(MeanCov = sum(Cover, na.rm = TRUE) / nplots[1],
                        Constancy = (.N / nplots[1]) * 100, nplots = nplots[1]),
                    by = .(SiteUnit, Species)]

  ## --------Create the analysis set
  vegsum <- as.data.frame(vegsum)
  vegsum$MeanCov[vegsum$MeanCov > 100] <- 100
  vegsum$spp_importance <- vegsum$MeanCov^1 / 2
  vegsum <- vegsum %>%
    dplyr::mutate(spp_importance = spp_importance * (Constancy / 100)) %>%
    dplyr::filter(spp_importance > minimportance)
  ### remove species where the maximum constancy for species is < contancy
  speciesmax <- vegsum %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(maxcons = max(Constancy)) %>%
    dplyr::filter(maxcons >= minconstancy)
  ### remove species with no constancy > minconstancy in any unit, 
  #where a unit has less than minplots, or
  #any unit species constancy is < noiseconstancy
  vegsum <- vegsum %>%
    dplyr::filter(Species %in% speciesmax$Species) %>%
    dplyr::filter(nplots >= minplots) %>%
    filter(Constancy > noiseconstancy)

  # vegsum <- create_diagnostic_veg(vegsum)
  ## --------assign potential diagnostics before pairwise-----------------
  vegsum <- vegsum %>%
    rowwise() %>%
    mutate(constant_type = ifelse((Constancy >= minconstancy & MeanCov >= 10), "cd",
      ifelse((Constancy >= minconstancy & MeanCov <= minor), "cm",
        ifelse(Constancy >= minconstancy, "c", NA)
      )
    ))
  ### Calculate differential potential
  vegsum <- vegsum %>%
    mutate(d.potential = ifelse(Constancy > 50 & Constancy < minconstancy,
              (Constancy^(1 / 3)) * (Constancy / 100) * ((Constancy - 50) / 10),
                  ifelse(Constancy >= minconstancy & MeanCov > minor, (Constancy^(1 / 3) * (Constancy / 100)), 
                         ifelse(Constancy >= minconstancy & MeanCov <= minor, (Constancy^(1 / 3) * (Constancy / 100)*0.5),0)
      #ifelse(constant_type %in% c("c", "cd"), (Constancy^(1 / 3) * (Constancy / 100)),
        #ifelse(constant_type %in% c("cm"), (Constancy^(1 / 3) * (Constancy / 100) * 0.5), 0)
      )
    ))

  vegsum <- vegsum %>% mutate(d.potential = ifelse(d.potential < .6, 0,
    ifelse(d.potential > 4, 4, d.potential)
  ))

  ### Calculate dominant differential potential
  vegsum <- vegsum %>% 
    mutate(dd.potential = ifelse(MeanCov >= 10, (MeanCov^(1 / 3)),
    ifelse((MeanCov >= 5 & MeanCov < 10), (MeanCov^(1 / 3) * ((10 - MeanCov) / 10)), 0)
  ))
  vegsum <- vegsum %>% mutate(dd.potential = ifelse(dd.potential > 4, 4,
    ifelse(dd.potential < 0.6, 0, dd.potential)
  ))

  vegsum <- vegsum %>%
    mutate(diagnostic.potential = ifelse((d.potential > 0 & dd.potential > 0),
      ((d.potential + dd.potential) * 1.25),
      (d.potential + dd.potential)
    ))
  setDT(vegsum)
  vegsum <- merge(vegsum, taxon.lifeform, by.x = c("Species"), by.y = c("Code"))
  # Adjust points for reduced.lifeforms (moss layer)
  if (isTRUE(reduce.lifeform)) {
    vegsum[,  diagnostic.potential := ifelse(Lifeform %in% reduced.lifeforms & !Species %in% reduced.exceptions, 
                          diagnostic.potential * reduction, diagnostic.potential)
           ]
  }
  
  ## adjust points for key indicators
  if (isTRUE(use.ksi)) {
    vegsum[, diagnostic.potential := ifelse(Species %in% ksi, diagnostic.potential * ksi.value, diagnostic.potential)
           ]
  }
 vegsum <- vegsum %>% dplyr::select(-ScientificName, -EnglishName, -Lifeform)
   ### unit summary of diagnostic potential
  vegsum <- vegsum %>%
    dplyr::group_by(SiteUnit) %>%
    mutate(unit.diag.sum = sum(diagnostic.potential))
  vegsum <- vegsum %>%
    dplyr::group_by(SiteUnit) %>%
    mutate(n_constants = sum(!is.na(constant_type)))


  ### ________________DO PAIR-WISE Comparison_______

  ## build pairs
  pairs <- unique(vegsum$SiteUnit) %>%
    combn(m = 2) %>%
    t() %>%
    data.frame() %>%
    dplyr::rename(Unit1 = 1, Unit2 = 2) %>%
    arrange(Unit1)
  setDT(pairs)
  setDT(vegsum)
  setDT(taxon.lifeform)

  vegsum.pairs1 <- pairs[vegsum, on = c("Unit1" = "SiteUnit"), allow.cartesian = TRUE]
  vegsum.pairs2 <- pairs[vegsum, on = c("Unit2" = "SiteUnit"), allow.cartesian = TRUE]

  # vegsum.pairs both ways merged
  vegsum.pairs <- merge(vegsum.pairs1, vegsum.pairs2, by = c("Unit1", "Unit2", "Species"), all = TRUE)
  vegsum.pairs <- vegsum.pairs %>% mutate_if(is.numeric, replace_na, replace = 0)
  vegsum.pairs <- vegsum.pairs %>% filter(!is.na(Unit1), !is.na(Unit2))
  # add taxon.lifeform
  vegsum.pairs <- merge(vegsum.pairs, taxon.lifeform, by.x = c("Species"), by.y = c("Code"))

  setDT(vegsum.pairs)[, c("cov.diff", "const.diff") := .(MeanCov.x - MeanCov.y, Constancy.x - Constancy.y)]
  setkey(vegsum.pairs, "Unit1", "Unit2", "Species")
  ## calculate shared diagnostic potential
  vegsum.pairs[, `:=`(
    shared.diag = pmin(diagnostic.potential.x, diagnostic.potential.y, na.rm = TRUE)
  ), by = .(Unit1, Unit2, Species)]
  # Differential
  ## assigns differential type based on cut-off ranges
  const_cut <- c(37, 60, 80, 101)
  const_cut2 <- c(-101, -80, -60, -37)
  const_labels <- c("d3", "d2", "d1")
  const_labels2 <- c("d1", "d2", "d3")
  vegsum.temp <- vegsum.pairs
  vegsum.pairs <- vegsum.temp
  vegsum.pairs[, d.type.x := cut(const.diff, breaks = const_cut, labels = const_labels)]
  vegsum.pairs[, d.type.y := cut(const.diff, breaks = const_cut2, labels = const_labels2)]

  ## if minimum constancy is not met sets to NA
  vegsum.pairs[, d.type.x := ifelse(Constancy.x < minconstancy, NA, paste0("", d.type.x))]
  vegsum.pairs[, d.type.y := ifelse(Constancy.y < minconstancy, NA, paste0("", d.type.y))]

  ## calculate d points
  vegsum.pairs[, d.points.x := ifelse(((Constancy.x >= minconstancy) & (const.diff > 30)), (const.diff ^(1 / 3)*(const.diff / 100)),
              ifelse(((Constancy.x <minconstancy) & (Constancy.x >(minconstancy-10)) & (const.diff>30)), ((const.diff ^(1 / 3)*(const.diff / 100)*((abs(10-(60.01-Constancy.x)))))/10),0))]

  vegsum.pairs[, d.points.y := ifelse(((Constancy.y >= minconstancy) & (const.diff < (-30))), (abs(const.diff)^(1 / 3) * (abs(const.diff) / 100)),
                                    ifelse(((Constancy.y < minconstancy) & (Constancy.y >(minconstancy-10)) & (const.diff < (-30))), ((((abs(const.diff)) ^(1 / 3)) * (abs(const.diff))/ 100) * ((abs(10-(60.01- Constancy.y)))/10)),0))]
  # vegsum.pairs[, const2.x := ifelse(Constancy.x >=60, (Constancy.x^(1 / 3))*(Constancy.x / 100), (Constancy.x^(1 / 3)*(Constancy.x / 200)))]
  # vegsum.pairs[, const2.y := ifelse(Constancy.y >=60, (Constancy.y^(1 / 3))*(Constancy.y / 100), (Constancy.y^(1 / 3)*(Constancy.y / 200)))]
  # 
  # vegsum.pairs[, const2.x := ifelse(const2.x>4, 4, const2.x)]
  # vegsum.pairs[, const2.y := ifelse(const2.y>4, 4, const2.y)]
  # 
  # vegsum.pairs[, const2.diff := const2.x - const2.y]

  ### adjust d points for constancy of minor species
  vegsum.pairs[, c("d.points.x", "d.points.y") := .(
    ifelse((MeanCov.x < minor), (d.points.x / 2), d.points.x),
      #ifelse(const.diff >=40, const2.diff, NA_real_)
   # ),
    ifelse((MeanCov.y < minor), (d.points.y / 2), d.points.y)
    #  ifelse(const.diff <= -40,  (0 - const2.diff), NA_real_))
  )]
  # vegsum.pairs[, "d.points.x" := .(
  #   ifelse(Constancy.x < 60, d.points.x * (Constancy.x / 100), d.points.x)
  # )]
  # vegsum.pairs[, "d.points.y" := .(
  #   ifelse(Constancy.y < 60, d.points.y * (Constancy.y / 100), d.points.y)
  # )]

  ## remove non significant differences and limit to max 4 points
  vegsum.pairs[, "d.points.x" := .(
    ifelse(d.points.x < 0.6, NA_real_,
      ifelse(d.points.x > 4, 4, d.points.x)
    )
  )]
  vegsum.pairs[, "d.points.y" := .(
    ifelse(d.points.y < 0.6, NA_real_,
      ifelse(d.points.y > 4, 4, d.points.y)
    )
  )]

  ## calculate dd points
  # vegsum.pairs[, `:=`(cover2.x = ifelse(MeanCov.x >= 10, MeanCov.x^(1 / 3), (MeanCov.x^(1 / 3) * (1 - abs((MeanCov.x - 10) / 30)))))]
  # vegsum.pairs[, `:=`(cover2.y = ifelse(MeanCov.y >= 10, MeanCov.y^(1 / 3), (MeanCov.y^(1 / 3) * (1 - abs((MeanCov.y - 10) / 30)))))]
  # vegsum.pairs[, cover2.diff := cover2.x - cover2.y]
  vegsum.pairs[, dd.points.x := ifelse(MeanCov.x > MeanCov.y & MeanCov.x >= 10, (MeanCov.x^(1 / 3) - (MeanCov.y ^(1 / 3))),
                                       ifelse(MeanCov.x > MeanCov.y & MeanCov.x < 10 & MeanCov.x > 5, ((MeanCov.x^(1 / 3) - (MeanCov.y ^(1 / 3)))* ((10 - MeanCov.x) / 10)),0))]
  vegsum.pairs[, dd.points.y := ifelse(MeanCov.y > MeanCov.x & MeanCov.y >= 10, (MeanCov.y^(1 / 3) - (MeanCov.x ^(1 / 3))),
                                       ifelse(MeanCov.y > MeanCov.x & MeanCov.y < 10 & MeanCov.y > 5, ((MeanCov.y^(1 / 3) - (MeanCov.x ^(1 / 3)))* ((10 - MeanCov.y) / 10)),0))]
  
  # vegsum.pairs[, `:=`(dd.points.y = ifelse(MeanCov.y >= 10, (MeanCov.y^(1 / 3) - (MeanCov.y ^(1 / 3))),0))]
  # vegsum.pairs[, `:=`(dd.points.x = ifelse((MeanCov.x < 10 & MeanCov.x >5), ((MeanCov.x^(1 / 3) - (MeanCov.y ^(1 / 3)))* ((10 - MeanCov.x) / 10)),0))]
  # vegsum.pairs[, `:=`(dd.points.y = ifelse((MeanCov.y < 10 & MeanCov.y >5), ((MeanCov.y^(1 / 3) - (MeanCov.x ^(1 / 3)))* ((10 - MeanCov.y) / 10)),0))]
 
  setDT(vegsum.pairs)

  # Eliminate dd points that fall below or above thresholds
  vegsum.pairs[, `:=`(
    dd.points.x = fifelse(dd.points.x < 0.6, 0, dd.points.x),
    dd.points.y = fifelse(dd.points.y < 0.6, 0, dd.points.y)
  )][, `:=`(
    dd.points.x = fifelse(dd.points.x > 4, 4, dd.points.x),
    dd.points.y = fifelse(dd.points.y > 4, 4, dd.points.y)
  )]

  ## assign dd type
  vegsum.pairs[, `:=`(
    dd.type.x = ifelse((dd.points.x >= 3) & (constant_type.x == "cd"), "dd1",
      ifelse((dd.points.x>= 2 & dd.points.x < 3) & (constant_type.x == "cd"), "dd2",
        ifelse((dd.points.x >= 1 & dd.points.x < 2) & (constant_type.x == "cd"), "dd3",
          ifelse((dd.points.x >= 0.6 & dd.points.x < 1) & (constant_type.x == "cd"), "dd4", NA)
        )
      )
    ),
    dd.type.y = ifelse((dd.points.y >= 3) & (constant_type.y == "cd"), "dd1",
                       ifelse((dd.points.y>= 2 & dd.points.y < 3) & (constant_type.y == "cd"), "dd2",
                              ifelse((dd.points.y >= 1 & dd.points.y < 2) & (constant_type.y == "cd"), "dd3",
                                     ifelse((dd.points.y >= 0.6 & dd.points.y < 1) & (constant_type.y == "cd"), "dd4", NA)
        )
      )
    )
  )]

  # Adjust points for reduced.lifeforms (moss layer)
  if (isTRUE(reduce.lifeform)) {
    vegsum.pairs[, `:=`(
      d.points.x = ifelse(Lifeform %in% reduced.lifeforms & !Species %in% reduced.exceptions, d.points.x * 0.1, d.points.x),
      d.points.y = ifelse(Lifeform %in% reduced.lifeforms & !Species %in% reduced.exceptions, d.points.y * 0.1, d.points.y),
      dd.points.x = ifelse(Lifeform %in% reduced.lifeforms & !Species %in% reduced.exceptions, dd.points.x * 0.1, dd.points.x),
      dd.points.y = ifelse(Lifeform %in% reduced.lifeforms & !Species %in% reduced.exceptions, dd.points.y * 0.1, dd.points.y)
    )]
  }
  # Sum sum of diagnostic differentials by species
  setkey(vegsum.pairs, "Unit1", "Unit2", "Species")

  ### sum total points
  vegsum.pairs[, `:=`(
    diff.pts.x = sum(d.points.x, dd.points.x, na.rm = TRUE),
    diff.pts.y = sum(d.points.y, dd.points.y, na.rm = TRUE)
  ), by = .(Unit1, Unit2, Species)]

  ## adjust points for being d as well as dd
  vegsum.pairs[, `:=`(
    diff.pts.x = ifelse((!is.na(d.type.x) & !is.na(dd.type.x)), (diff.pts.x * 1.25), diff.pts.x),
    diff.pts.y = ifelse((!is.na(d.type.y) & !is.na(dd.type.y)), (diff.pts.y * 1.25), diff.pts.y)
  ), by = .(Unit1, Unit2, Species)]

  ## adjust points for key indicators
  if (isTRUE(use.ksi)) {
    vegsum.pairs[, `:=`(
      diff.pts.x = ifelse(Species %in% ksi, (diff.pts.x * ksi.value), diff.pts.x),
      diff.pts.y = ifelse(Species %in% ksi, (diff.pts.y * ksi.value), diff.pts.y)
    ), by = .(Unit1, Unit2, Species)]
  }

  # Sums by pair
  vegsum.pairs[, `:=`(
    sum.shared.diag = sum(shared.diag),
    diff.tot.x = sum(diff.pts.x, na.rm = TRUE),
    diff.tot.y = sum(diff.pts.y, na.rm = TRUE)
  ), by = .(Unit1, Unit2)]

  vegsum.pairs[, diff.tot := rowSums(.SD, na.rm = TRUE), .SDcols = c("diff.tot.x", "diff.tot.y"), by = .(Unit1, Unit2)]

  vegsum.pairs[, `:=`(
    diff.ratio = (sum.shared.diag * 2) / (diff.tot + (sum.shared.diag * 2)),
    diff.ratio.x = sum.shared.diag / (max(diff.tot.x) + sum.shared.diag),
    diff.ratio.y = sum.shared.diag / (max(diff.tot.y) + sum.shared.diag)
  ), by = .(Unit1, Unit2)]

  vegsum.pairs[, `:=`(
    BEC.sim.mean = (diff.ratio.x + diff.ratio.y) / 2,
    BEC.sim.max = pmax(diff.ratio.x, diff.ratio.y),
    BEC.sim.min = pmin(diff.ratio.x, diff.ratio.y)
  ), by = .(Unit1, Unit2)]
  # Ungroup the data.table
  vegsum.pairs <- vegsum.pairs[, .SD, .SDcols = names(vegsum.pairs)]
  vegsum.pairs <- vegsum.pairs %>%
    # dplyr::select(-EnglishName.x, -EnglishName.y, -ScientificName.x, -ScientificName.y) %>%
    dplyr::select(
      Unit1, Unit2, BEC.sim.min, BEC.sim.mean, BEC.sim.max, diff.ratio,
      unit.diag.sum.x, unit.diag.sum.y, everything()
    )
 toc()
   return(vegsum.pairs)
}
