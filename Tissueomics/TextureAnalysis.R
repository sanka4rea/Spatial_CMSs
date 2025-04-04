# texture of whole slides------------------------------------------------
texture_features <- function(LocMatrix=Loc_matrix){ # label=Labels,
  require(radiomics)

  texture <- mclapply(1:length(LocMatrix), function(x){
    mat <- LocMatrix[[x]]
    first_or <- calc_features(mat)
    #GLCMs
    glcm1 <- glcm(mat, angle=0, d=1, n_grey = length(unique(as.numeric(mat))))
    glcm_or <- calc_features(glcm1)
    # Grey Level Run Length Matrix (GLRLM)
    glrlm_1 <- glrlm(mat, angle=0, verbose=FALSE,n_grey = length(unique(as.numeric(mat))))
    glrlm_or <- calc_features(glrlm_1)
    # Grey Level Size Zone Matrix (GLSZM) and Multiple-GLSZM
    glszm_1 <- glszm(mat, n_grey = length(unique(as.numeric(mat))))  
    glszm_or <- calc_features(glszm_1)
    
    all <- c(as.numeric(first_or),as.numeric(glcm_or), as.numeric(glrlm_or), as.numeric(glszm_or))
    names(all) <- c(colnames(first_or),colnames(glcm_or),colnames(glrlm_or),colnames(glszm_or))
    return(all)
  },mc.cores=20)
  texture <- do.call(rbind,texture)
  rownames(texture) <- names(LocMatrix)
  
  return(t(texture))
}
