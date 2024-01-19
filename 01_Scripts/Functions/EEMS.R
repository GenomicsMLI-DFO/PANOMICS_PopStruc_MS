# Use the "pairwise.complete.obs" method to compute pairwise dissimilarities
# This straightforward implementation
# uses a double loop, so would be slow if the sample size is large.
bed2diffs_v1 <- function(genotypes) {
  
  nIndiv <- nrow(genotypes)
  nSites <- ncol(genotypes)
  diffs <- matrix(0, nIndiv, nIndiv)
  
  for (i in seq(nIndiv - 1)) {
    for (j in seq(i + 1, nIndiv)) {
      x <- genotypes[i, ]
      y <- genotypes[j, ]
      diffs[i, j] <- mean((x - y)^2, na.rm = TRUE)
      diffs[j, i] <- diffs[i, j]
    }
  }
  
  diffs
}



# Compute the diffs matrix using the "mean allele frequency"
# imputation method
bed2diffs_v2 <- function(genotypes) {
  
  nIndiv <- nrow(genotypes)
  nSites <- ncol(genotypes)
  missing <- is.na(genotypes)
  
  ## Impute NAs with the column means (= twice the allele frequencies)
  geno_means <- colMeans(genotypes, na.rm = TRUE)
  # nIndiv rows of genotype means
  geno_means <- matrix(geno_means, nrow = nIndiv, ncol = nSites, byrow = TRUE)
  
  ## Set the means which correspond to observed genotypes to 0
  geno_means[missing == FALSE] <- 0
  ## Set the missing genotypes to 0 (used to be NA)
  genotypes[missing == TRUE] <- 0
  genotypes <- genotypes + geno_means
  
  similarities <- genotypes %*% t(genotypes) / nSites
  self_similarities <- diag(similarities)
  vector1s <- rep(1, nIndiv)
  
  diffs <-
    self_similarities %*% t(vector1s) +
    vector1s %*% t(self_similarities) - 2 * similarities
  diffs
}



extract.EEMS.raster <- function(mcmcpath, is.mrates = TRUE){

    # BELOW AN EXAMPLE TO RECREATE A RASTER
  
  #mcmcpath <-  file.path(here::here(),"./02_Results/01_Overall_PopGen/06_EEMS", paste0("nDemes", 1000, "-chain",1))
  
  plot.params <- list(
    plot.width = 7,
    plot.height = 7,
    out.png = TRUE,
    res = 600,
    xpd = TRUE,
    col.grid = "gray80",
    lwd.grid = 1,
    add.demes = FALSE,
    col.demes = "black",
    pch.demes = 19,
    min.cex.demes = 1,
    max.cex.demes = 3,
    add.outline = FALSE,
    col.outline = "gray90",
    lwd.outline = 2,
    add.grid = T,
    projection.in = NULL,
    projection.out = NULL,
    add.map = FALSE,
    col.map = "gray60",
    lwd.map = 2,
    eems.colors = NULL,
    prob.levels = c(0.9, 0.95),
    add.colbar = TRUE,
    m.colscale = NULL,
    q.colscale = NULL,
    remove.singletons = TRUE,
    add.abline = FALSE,
    add.r.squared = FALSE,
    add.title = TRUE,
    m.plot.xy = NULL,
    q.plot.xy = NULL,
    xy.coords = NULL
  )
  
  #plot.params <-  rEEMSplots:::check.plot.params(plot.params)
  
  dimns <- rEEMSplots:::read.dimns(mcmcpath[1], longlat = TRUE, coord = NULL)
  
  # average.eems.contours
  #zero_mean <- plot.params$m.zero_mean
  #log_scale <- plot.params$m.log_scale
  
  Zmean <- matrix(0, dimns$nxmrks, dimns$nymrks)
  
  #PrGT0 <- matrix(0, dimns$nxmrks, dimns$nymrks)
  #PrLT0 <- matrix(0, dimns$nxmrks, dimns$nymrks)
  #coord_extra <- nrow(dimns$coord)
  #Zmean_extra <- numeric(coord_extra)
  
  niters <- 0
  
  for (path in mcmcpath) {
  
  voronoi <- rEEMSplots:::read.voronoi(path, longlat = TRUE, is.mrates = is.mrates, log_scale = TRUE)
  
  tiles <- voronoi$tiles
  rates <- voronoi$rates
  xseed <- voronoi$xseed
  yseed <- voronoi$yseed
  
  rslt <- rEEMSplots:::transform.rates(dimns, tiles, rates, xseed, yseed, zero_mean = TRUE)
  Zmean <- Zmean + rslt$Zvals
  #Zmean_extra <- Zmean_extra + rslt$Zvals_extra
  
  niters <- niters + rslt$niters
  
  #if (zero_mean) {
  #  PrGT0 <- PrGT0 + rslt$PrGT0
  #  PrLT0 <- PrLT0 + rslt$PrLT0
  #}
  }
  
  #PrGT0 <- PrGT0 + rslt$PrGT0
  #PrLT0 <- PrLT0 + rslt$PrLT0
  
  Zmean <- Zmean / niters
  #PrGT0 <- PrGT0 / niters
  #PrLT0 <- PrLT0 / niters
  #Zmean_extra <- Zmean_extra / niters
  
  # THIS is the RASTER we want to export - and cut at the righ format with QGIS
  ras.EEMS <- rEEMSplots:::points_to_raster(Zmean, dimns, plot.params)
  
  return(ras.EEMS)
  #writeRaster(ras.EEMS,'test.tif',options=c('TFW=YES'))
  
}



extract.EEMS.prob.raster <- function(mcmcpath, is.mrates = TRUE){
  
  # BELOW AN EXAMPLE TO RECREATE A RASTER
  
  #mcmcpath <-  file.path(here::here(),"./02_Results/01_Overall_PopGen/06_EEMS", paste0("nDemes", 1000, "-chain",1))
  
  plot.params <- list(
    plot.width = 7,
    plot.height = 7,
    out.png = TRUE,
    res = 600,
    xpd = TRUE,
    col.grid = "gray80",
    lwd.grid = 1,
    add.demes = FALSE,
    col.demes = "black",
    pch.demes = 19,
    min.cex.demes = 1,
    max.cex.demes = 3,
    add.outline = FALSE,
    col.outline = "gray90",
    lwd.outline = 2,
    add.grid = T,
    projection.in = NULL,
    projection.out = NULL,
    add.map = FALSE,
    col.map = "gray60",
    lwd.map = 2,
    eems.colors = NULL,
    prob.levels = c(0.9, 0.95),
    add.colbar = TRUE,
    m.colscale = NULL,
    q.colscale = NULL,
    remove.singletons = TRUE,
    add.abline = FALSE,
    add.r.squared = FALSE,
    add.title = TRUE,
    m.plot.xy = NULL,
    q.plot.xy = NULL,
    xy.coords = NULL
  )
  
  #plot.params <-  rEEMSplots:::check.plot.params(plot.params)
  
  dimns <- rEEMSplots:::read.dimns(mcmcpath[1], longlat = TRUE, coord = NULL)
  
  # average.eems.contours
  #zero_mean <- plot.params$m.zero_mean
  #log_scale <- plot.params$m.log_scale
  
  #Zmean <- matrix(0, dimns$nxmrks, dimns$nymrks)
  
  PrGT0 <- matrix(0, dimns$nxmrks, dimns$nymrks)
  PrLT0 <- matrix(0, dimns$nxmrks, dimns$nymrks)
  #coord_extra <- nrow(dimns$coord)
  #Zmean_extra <- numeric(coord_extra)
  
  niters <- 0
  for (path in mcmcpath) {
    
  voronoi <- rEEMSplots:::read.voronoi(path, longlat = TRUE, is.mrates = is.mrates, log_scale = TRUE)
  
  tiles <- voronoi$tiles
  rates <- voronoi$rates
  xseed <- voronoi$xseed
  yseed <- voronoi$yseed
  
  rslt <- rEEMSplots:::transform.rates(dimns, tiles, rates, xseed, yseed, zero_mean = TRUE)
  #Zmean <- Zmean + rslt$Zvals
  #Zmean_extra <- Zmean_extra + rslt$Zvals_extra
  
  niters <- niters + rslt$niters
  
  PrGT0 <- PrGT0 + rslt$PrGT0
  PrLT0 <- PrLT0 + rslt$PrLT0
  
  }
  
  #Zmean <- Zmean / niters
  PrGT0 <- PrGT0 / niters
  PrLT0 <- PrLT0 / niters
  #Zmean_extra <- Zmean_extra / niters
  
  # THIS is the RASTER we want to export - and cut at the righ format with QGIS

  ras.EEMS <- rEEMSplots:::points_to_raster(PrGT0 - PrLT0, dimns, plot.params)
  # 2 sides tests
  ras.EEMS <-  (ras.EEMS + 1) / 2
  
  return(ras.EEMS)
  #writeRaster(ras.EEMS,'test.tif',options=c('TFW=YES'))
  
}



