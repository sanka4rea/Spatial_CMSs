##  calculate for each patients: average value of slides------------------------------------------------
whole_slides_infiltration_R <- function(CountMatrix=whole_slides_matrix,SlideNamesIndex1=1,SlideNamesIndex2=9){ # label=Labels,
  totalnum <- rowSums(CountMatrix)
  CountMatrix <- CountMatrix/totalnum
  #namesDen <- substr(rownames(CountMatrix),SlideNamesIndex1,SlideNamesIndex2)
  #namesDen_uni <- unique(namesDen)
  ## whole_slides_w64_t50_ErQi[which(whole_slides_w64_t50_ErQi==-1)] <- 0
  #density_uni <- sapply(1:length(namesDen_uni), function(x){
    #index <- which(namesDen==namesDen_uni[x])
    #if(length(index)>1){
      #tmp <- apply(CountMatrix[index,], 2, mean)
    #}else{
      #tmp <- CountMatrix[index,]
    #}
    #tmp <- as.numeric(tmp)
    #return(tmp)
  #})
  #colnames(density_uni) <- namesDen_uni
  #rownames(density_uni) <- colnames(CountMatrix)
  
  #label=Labels,
  #if(!is.null(label)){
  #  density_uni <- as.data.frame(t(density_uni))
  #  EQ_index_fm <- cbind(density_uni[match(names(label),rownames(density_uni)),],label)
  #  colnames(EQ_index_fm) <- c(colnames(CountMatrix),'Subtype')
  #  EQ_index_fm$Subtype <- factor(EQ_index_fm$Subtype)
  #}
  colnames(CountMatrix) <- paste0("WS_infil_ratio",colnames(CountMatrix))
  return(CountMatrix)
}

## calculate for each patients: average value of slides------------------------------------------------
whole_slides_intratumor_R <- function(CountMatrix=whole_slides_matrix,SlideNamesIndex1=1,SlideNamesIndex2=9,tumorindex=7){ # label=Labels,
  CountMatrix <- CountMatrix[,-tumorindex]/CountMatrix[,tumorindex]
  #namesDen <- substr(rownames(CountMatrix),SlideNamesIndex1,SlideNamesIndex2)
  #namesDen_uni <- unique(namesDen)
  ## whole_slides_w64_t50_ErQi[which(whole_slides_w64_t50_ErQi==-1)] <- 0
  #density_uni <- sapply(1:length(namesDen_uni), function(x){
    #index <- which(namesDen==namesDen_uni[x])
    #if(length(index)>1){
      #tmp <- apply(CountMatrix[index,], 2, mean)
    #}else{
      #tmp <- CountMatrix[index,]
    #}
    #tmp <- as.numeric(tmp)
    #return(tmp)
  #})
  #colnames(density_uni) <- namesDen_uni
  #rownames(density_uni) <- colnames(CountMatrix)
  
  #label=Labels,
  #if(!is.null(label)){
  #  density_uni <- as.data.frame(t(density_uni))
  #  EQ_index_fm <- cbind(density_uni[match(names(label),rownames(density_uni)),],label)
  #  colnames(EQ_index_fm) <- c(colnames(CountMatrix),'Subtype')
  #  EQ_index_fm$Subtype <- factor(EQ_index_fm$Subtype)
  #}
  colnames(CountMatrix) <- paste0("WS_intra_ratio",colnames(CountMatrix))
  return(CountMatrix)
}

# calculate all/specific entropy------------------------------------------------
whole_slides_entropy <- function(CountMatrix=whole_slides_matrix,specificLabels=NULL, # specificLabels=c("lym","str","tum")
                                 SlideNamesIndex1=1,SlideNamesIndex2=9){
  if(is.null(specificLabels)){
    allentropy <- apply(CountMatrix,1,entropy::entropy) 
  }else{
    allentropy <- apply(CountMatrix[,specificLabels],1,entropy::entropy)  
  }

  #namesDen <- substr(rownames(CountMatrix),SlideNamesIndex1,SlideNamesIndex2)
  #namesDen_uni <- unique(namesDen)
  ## whole_slides_w64_t50_ErQi[which(whole_slides_w64_t50_ErQi==-1)] <- 0
  #density_uni <- sapply(1:length(namesDen_uni), function(x){
    #index <- which(namesDen==namesDen_uni[x])
    #tmpallentropy <- mean(allentropy[index])
    #return(tmpallentropy)
  #})
  #names(density_uni) <- namesDen_uni
  return(allentropy)
}

# whole-slides within sum-of-squares for cluster------------------------------------------------
whole_slides_centroid <- function(LocMatrix=Loc_matrix,tissueLabels=c(2:8),tissueNames=NULL){
  
  Level1 <- mclapply(1:length(LocMatrix),function(i){
    mat <- LocMatrix[[i]]
    # reference location matrix
    ref_mat <- matrix(0,nrow = nrow(mat) ,ncol = ncol(mat))
    for (x in 1:nrow(ref_mat)) {
      for(y in 1:ncol(ref_mat)){
        ref_mat[x,y] <- paste0(x,"_",y)
      }
    }
    covert2mat <- function(index_n){
      mat <- ref_mat[index_n]
      tmpmat <- matrix(0,nrow = length(mat),ncol = 2)
      for (i in 1:length(mat)) {
        tmpmat[i,1] <- as.numeric(str_split(mat[i],"_")[[1]][1]) 
        tmpmat[i,2] <- as.numeric(str_split(mat[i],"_")[[1]][2]) 
      }
      return(tmpmat)
    }
    
    # each class
    eachClass <- lapply(tissueLabels, function(j){
      
      if(length(which(mat==j))> 1){ # at least two
        con_dis <- covert2mat(which(mat==j))
        centroid <- c(mean(as.numeric(con_dis[,1])), mean(as.numeric(con_dis[,2])) )
        
        # centroid distance
        cen_con_mat <- rbind(centroid,con_dis)
        cen_con_dis <- dist(cen_con_mat)
        cen_con_dis <- cen_con_dis[1:nrow(con_dis)]
        
        maxCendis <- max(cen_con_dis)
        averCendis <- mean(cen_con_dis)
        sdCendis <- sd(cen_con_dis)
        
        # the within-cluster variance
        wcv <- sum(as.matrix(dist(con_dis)^2)) / (2 * nrow(con_dis))
        aver_wcv <- wcv/nrow(con_dis)
        
        three_con_dis <- c(maxCendis,averCendis,sdCendis,aver_wcv)
      }else{
        three_con_dis <- c(NA,NA,NA,NA)
      }
      return(three_con_dis)
    })
    eachClass <- do.call(rbind,eachClass)
    #rownames(eachClass) <- tissueLabels
    #colnames(eachClass) <- c("WS_maxCendis","WS_averCendis","WS_aver_wcv")
    eachClass <- as.vector(eachClass)
    names(eachClass) <- c(paste0(tissueNames,"_","WS_maxCendis"),paste0(tissueNames,"_","WS_averCendis"),
                          paste0(tissueNames,"_","WS_sdCendis"),paste0(tissueNames,"_","WS_aver_wcv"))
    return(eachClass)
  },mc.cores = 40)
  
  Level1 <- do.call(rbind,Level1)
  rownames(Level1) <-names(LocMatrix)
  return(Level1) 
}

# whole-slides between distance for clusters------------------------------------------------
whole_slides_between <- function(LocMatrix=Loc_matrix,tissueLabels=c(2:8),tissueNames=NULL){
  
  Level1 <- mclapply(1:length(LocMatrix),function(i){
    mat <- LocMatrix[[i]]
    # reference location matrix
    ref_mat <- matrix(0,nrow = nrow(mat) ,ncol = ncol(mat))
    for (x in 1:nrow(ref_mat)) {
      for(y in 1:ncol(ref_mat)){
        ref_mat[x,y] <- paste0(x,"_",y)
      }
    }
    covert2mat <- function(index_n){
      mat <- ref_mat[index_n]
      tmpmat <- matrix(0,nrow = length(mat),ncol = 2)
      for (i in 1:length(mat)) {
        tmpmat[i,1] <- as.numeric(str_split(mat[i],"_")[[1]][1]) 
        tmpmat[i,2] <- as.numeric(str_split(mat[i],"_")[[1]][2]) 
      }
      return(tmpmat)
    }
    
    # each class centroid
    eachCentroid <- sapply(length(tissueLabels):1, function(j){
      if(length(which(mat==tissueLabels[j]))> 1){ # at least two
        con_dis <- covert2mat(which(mat==tissueLabels[j]))
        centroid <- c(mean(as.numeric(con_dis[,1])), mean(as.numeric(con_dis[,2])) )
      }else{
        centroid <- c(NA,NA)
      }
      return(centroid)
    })
    eachCentroid <- t(eachCentroid) # col=x,y
    rownames(eachCentroid) <- tissueNames[length(tissueLabels):1]
    
    dis_betw <- as.numeric(dist(eachCentroid))
    name <- lapply(1:(nrow(eachCentroid)-1), function(z){
      tmp <- paste0("WS_betwDis_",rownames(eachCentroid)[z],"_",rownames(eachCentroid)[(z+1):nrow(eachCentroid)])
      return(tmp)
    })
    names(dis_betw) <- do.call(c,name) 
    return(dis_betw)
    
  },mc.cores = 40)
  
  Level1 <- do.call(rbind,Level1)
  rownames(Level1) <-names(LocMatrix)
  return(Level1) 
}


# gaussian clustering based distance to centroid------------------------------------------------
whole_slides_gaussian4 <- function(LocMatrix=Loc_matrix,tissueLabels=7,tumorLabel=8,delBackLabel=1,tissueName="STR"){
  require(mclust)
  # obtain all scores
  allscores <- mclapply(1:length(LocMatrix),function(i){
    mat <- LocMatrix[[i]]
    ref_mat <- matrix(0,nrow = nrow(mat) ,ncol = ncol(mat))
    for (x in 1:nrow(ref_mat)) {
      for(y in 1:ncol(ref_mat)){
        ref_mat[x,y] <- paste0(x,"_",y)
      }
    }
    
    covert2mat <- function(index_n){
      mat <- ref_mat[index_n]
      tmpmat <- matrix(0,nrow = length(mat),ncol = 2)
      for (i in 1:length(mat)) {
        tmpmat[i,1] <- as.numeric(str_split(mat[i],"_")[[1]][1]) 
        tmpmat[i,2] <- as.numeric(str_split(mat[i],"_")[[1]][2]) 
      }
      return(tmpmat)
    }
    
    if(length(which(mat==tissueLabels))>0){
      tissue.loc <- covert2mat(which(mat==tissueLabels))
      cancer.loc <- covert2mat(which(mat==tumorLabel))
      cancer.centroid <- c(mean(as.numeric(cancer.loc[,1])), mean(as.numeric(cancer.loc[,2])) )
      
      score <- dist(rbind(cancer.centroid,tissue.loc))[1:nrow(tissue.loc)]
    }else{
      score <- NA
    }
    return(score)
    
  },mc.cores = 40)
  
  # mclust
  allscores4mc <- do.call(c,allscores)
  allscores4mc <- allscores4mc[!is.na(allscores4mc)]
  set.seed(10)
  x <- sample(allscores4mc,length(allscores4mc)*0.1)
  res <- Mclust(x, G=3)
  
  #plot(hist(x,breaks=500, plot=FALSE),col="grey",border="grey",freq=FALSE,
       #xlab="Lymphocyte proximity to cancer",main="")
  #lines(density(x),lty=2)
  #plot.normal.components.mclust(x, res)
  
  th1 <- max(res$data[which(res$classification==1)])
  th2 <- max(res$data[which(res$classification==2)])
  
  # classification all scores and calculate the ratio in each slides
  Ratios <- mclapply(1:length(allscores),function(zz){
    scores <- allscores[[zz]]
    mat <- LocMatrix[[zz]]
    
    Len_infil <- sum(scores<=th1)
    Len_adjac <- sum(th1<scores & scores<=th2)
    Len_distal <- sum(scores>th2)
    
    infil_infil_ratio <- Len_infil/length(intersect(which(mat!=0),which(mat!=delBackLabel))) 
    infil_adjac_ratio <- Len_adjac/length(intersect(which(mat!=0),which(mat!=delBackLabel)))
    infil_distal_ratio <- Len_distal/length(intersect(which(mat!=0),which(mat!=delBackLabel)))
    
    intra_infil_ratio <- Len_infil/length(which(mat==tumorLabel)) 
    intra_adjac_ratio <- Len_adjac/length(which(mat==tumorLabel)) 
    intra_distal_ratio <- Len_distal/length(which(mat==tumorLabel)) 
    
    return(c(infil_infil_ratio,infil_adjac_ratio,infil_distal_ratio,
             intra_infil_ratio,intra_adjac_ratio,intra_distal_ratio))
  },mc.cores=20)
  Ratios <- do.call(rbind,Ratios)
  rownames(Ratios) <- names(LocMatrix)
  colnames(Ratios) <- paste0(tissueName,c("infil_infil_ratio","infil_adjac_ratio","infil_distal_ratio",
                        "intra_infil_ratio","intra_adjac_ratio","intra_distal_ratio"))
  
  return(Ratios)
}

#  plot.normal.components.mclust
plot.normal.components.mclust <- function(x, mixture,...){
  # Function to add Gaussian mixture components for mclust object, vertically scaled, to the current plot
  # Presumes the mixture object has the structure used by mixtools
  sapply(1:mixture$G,function(y)
    curve(mixture$parameters$pro[y] * dnorm(x,mean=mixture$parameter$mean[y], 
                                            sd=(ifelse(length(mixture$parameter$variance$sigmasq)==1, 
                                                       mixture$parameter$variance$sigmasq, mixture$parameter$variance$sigmasq[y]))^(1/2)), 
          add=TRUE, ...))
  
}


# gaussian clustering based proximity------------------------------------------------
whole_slides_gaussian3 <- function(LocMatrix=Loc_matrix,tissueLabels=7,tumorLabels=8){
  require(EBImage)
  require(splancs)
  # select h
  MSE <- NULL
  set.seed(10)
  ffs <- sample(1:length(LocMatrix), 20)
  for (ff in ffs){
    res <- LocMatrix[[ff]]
    ref_mat <- matrix(0,nrow = nrow(mat) ,ncol = ncol(mat))
    for (x in 1:nrow(ref_mat)) {
      for(y in 1:ncol(ref_mat)){
        ref_mat[x,y] <- paste0(x,"_",y)
      }
    }
    covert2mat <- function(index_n){
      mat <- ref_mat[index_n]
      tmpmat <- matrix(0,nrow = length(mat),ncol = 2)
      for (i in 1:length(mat)) {
        tmpmat[i,1] <- as.numeric(str_split(mat[i],"_")[[1]][1]) 
        tmpmat[i,2] <- as.numeric(str_split(mat[i],"_")[[1]][2]) 
      }
      return(tmpmat)
    }
    cancer.loc <- covert2mat(which(mat==tumorLabels))
    cv <- mse2d(as.points(cancer.loc), poly=cbind(c(0,0,nrow(mat),nrow(mat)),c(0,ncol(mat),ncol(mat),0)),nsmse=50,range=10)
    MSE <- rbind(MSE, cv$mse)
  }
  h=ceiling(which.min(apply(MSE, 2, mean))*0.2)
  
  # obtain all scores
  allscores <- mclapply(1:length(LocMatrix),function(i){
    mat <- LocMatrix[[i]]
    ref_mat <- matrix(0,nrow = nrow(mat) ,ncol = ncol(mat))
    for (x in 1:nrow(ref_mat)) {
      for(y in 1:ncol(ref_mat)){
        ref_mat[x,y] <- paste0(x,"_",y)
      }
    }
    
    covert2mat <- function(index_n){
      mat <- ref_mat[index_n]
      tmpmat <- matrix(0,nrow = length(mat),ncol = 2)
      for (i in 1:length(mat)) {
        tmpmat[i,1] <- as.numeric(str_split(mat[i],"_")[[1]][1]) 
        tmpmat[i,2] <- as.numeric(str_split(mat[i],"_")[[1]][2]) 
      }
      return(tmpmat)
    }
    
    tissue.loc <- as.data.frame(covert2mat(which(mat==tissueLabels))) 
    cancer.loc <- covert2mat(which(mat==tumorLabels))
    
    # get mask 
    # get location of the first non-zero row and col
    
    
    res <- kernel2d(as.points(cancer.loc),poly = cbind(c(non-zero-row,non-zero-row,max-zero-row,nrow(mask)),c(non-zero-col,max-zero-col,ncol(mat),max-zero-col)),h=h,nx=nrow(mat),ny=ncol(mat))
    # res <- kernel2d(as.points(cancer.loc),poly = cbind(c(30,30,120,120),c(30,130,130,30)),h=h,nx=130-30,ny=120-30)
    z.l <- unlist(sapply(1:length(tissue.loc$V1), function(x) res$z[tissue.loc$V1[x], tissue.loc$V2[x]]))
    
    return(z.l)
    
  },mc.cores = 20)
  
  # mclust
   ......
}
