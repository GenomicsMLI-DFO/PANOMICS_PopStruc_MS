# Source Capblancq

rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}



adaptive_index <- function(RDA, K, env_pres, range = NULL, method = "loadings", scale_env, center_env){
  
  # Formatting environmental rasters for projection
  var_env_proj_pres <- as.data.frame(rasterToPoints(env_pres[[row.names(RDA$CCA$biplot)]]))
  
  # Standardization of the environmental variables
  var_env_proj_RDA <- as.data.frame(scale(var_env_proj_pres[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))
  
  # Predicting pixels genetic component based on RDA axes
  Proj_pres <- list()
  if(method == "loadings"){
    for(i in 1:K){
      ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(apply(var_env_proj_RDA[,names(RDA$CCA$biplot[,i])], 1, function(x) sum( x * RDA$CCA$biplot[,i])))), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Prediction with RDA model and linear combinations
  if(method == "predict"){ 
    pred <- predict(RDA, var_env_proj_RDA[,names(RDA$CCA$biplot[,i])], type = "lc")
    for(i in 1:K){
      ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred[,i])), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Mask with the range if supplied
  if(!is.null(range)){
    Proj_pres <- lapply(Proj_pres, function(x) mask(x, range))
  }
  
  # Returning projections for current climates for each RDA axis
  return(Proj_pres = Proj_pres)
}

