# mask function
GetMask <- function(kernel1 = 5,kernel2 = 3,isFill = TRUE,locmatrix = EQ_loc_select,tumorLabel=8,minishape=64){
  #cl <- makeCluster(10)
  #clusterExport(cl,varlist = c("locmatrix","kernel1","kernel2","isFill")) 
  density_strtum_aver <- mclapply(1:length(locmatrix),function(i){
    tmp_mat <- locmatrix[[i]]
    mask <- tmp_mat==tumorLabel
    mask[which(mask==TRUE)] <- 1
    mask[which(mask==FALSE)] <- 0
    
    kern = EBImage::makeBrush(kernel1, shape='box')  # 3
    mask_1 <- EBImage::dilate(mask,kern = kern)
    kern1 = EBImage::makeBrush(kernel2, shape='box') 
    mask_1 <- EBImage::erode(mask_1,kern = kern1)
    if(isFill){
      mask_1 = EBImage::fillHull(mask_1)
    }
    return(mask_1)
  },mc.cores = 20)
  names(density_strtum_aver) <-names(locmatrix)
  mask_k53_ErQi <- density_strtum_aver
  
  mask_k53del_ErQi <- lapply(1:length(mask_k53_ErQi),function(x){
    tmp <- mask_k53_ErQi[[x]]
    # display(tmp)
    imgS=EBImage::bwlabel(tmp)
    imgSd=as.vector(EBImage::imageData(imgS)+1)
    hc=tabulate(imgSd)
    imgSdN=EBImage::imageData(imgS)+1
    rm(imgS)
    a=array(hc[imgSdN] < minishape,dim(imgSdN))
    imgSdN[a]=-1
    failures=which(imgSdN==-1)
    tmp[failures] <- 0
    # display(tmp)
    return(tmp)
  })
  names(mask_k53del_ErQi) <- names(mask_k53_ErQi)
  
  return(mask_k53del_ErQi)
}