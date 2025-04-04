# get Morisita-Horn similarity index----------------------------------------------------------------------------
MH_Index <- function(windows=8,tumorRatio=0.25,LocMatrix=EQ_loc_select,maskMatrix=mask_k53del_ErQi,
                              delBackLabel=1,tissueLabels=2:8,tumorLabel=8){
  require(spaa)
  Level1 <- lapply(1:length(LocMatrix),function(i){
    locmat <- LocMatrix[[i]]
    mask <- maskMatrix[[i]]
    
    ref_mat <- matrix(0,nrow = nrow(locmat) ,ncol = ncol(locmat))
    for (x in 1:nrow(ref_mat)) {
      for(y in 1:ncol(ref_mat)){
        ref_mat[x,y] <- paste0(x,"_",y)
      }
    }
    covert2mat <- function(ref_4,mat_4){
      tmpmat <- matrix(0,nrow = length(ref_4),ncol = 3)
      for (i in 1:length(ref_4)) {
        tmpmat[i,1] <- as.numeric(str_split(ref_4[i],"_")[[1]][1]) 
        tmpmat[i,2] <- as.numeric(str_split(ref_4[i],"_")[[1]][2]) 
      }
      tmpmat[,3] <-as.numeric(mat_4) 
      return(tmpmat)
    }
    
    squnR <- seq(1,nrow(mask)-round(windows/2),round(windows/2))
    lenR <- length(squnR)
    row_inter <- sapply(squnR[-lenR],function(row){
      # count <- vector("numeric",length=1)
      squn <- seq(1,ncol(mask)-round(windows/2),round(windows/2))
      len <- length(squn)
      all_1000 <- sapply(squn[-len],function(col){  # step = 4, if not enough for the last window, just discard it.
        mat_mask <- mask[row:(row+(windows-1)),col:(col+(windows-1))] 
        mat_4 <- locmat[row:(row+(windows-1)),col:(col+(windows-1))] 
        ref_4 <- ref_mat[row:(row+(windows-1)),col:(col+(windows-1))] 
        
        if( length(which(mat_4==0 | mat_4==delBackLabel))<=round(windows*windows/2) & length(which(mat_mask==0))< round(windows*windows/2) & 
            length(which(mat_mask==1))>= round(windows*windows*tumorRatio) & length(which(mat_4==tumorLabel))>=round(windows*windows*tumorRatio) ){ 
          # mat_mask ==1 means tumor regions
          all_count <- sapply(tissueLabels, function(zzz){
            len <- length(which(mat_4==zzz))
            return(len)
          })
          return(all_count)
        }else{
          return(as.numeric(rep(-1,length(tissueLabels))))
        }
      }) # all_1000
      
      all_1000 <- all_1000[,-which(all_1000[1,]==-1)]
      if(length(dim(all_1000))>0 ){
        if(dim(all_1000)[2]>0){
          aver <- all_1000
        }else{
          aver <-  all_1000
        }
      }else{
        aver <-  as.numeric(rep(0,length(tissueLabels)))
      }
      return(aver)
    }) # row_inter
    
    row_inter <- row_inter[lapply(row_inter,sum)>0]  # !!
    
    name <- NULL
    for (ii in length(tissueLabels):2) { # from tumor to con
      for (jj in (ii-1):1) {
        name <- c(name,paste0("MH_",tissueLabels[ii],"_",tissueLabels[jj]))
      }
    } 
    
    if(length((row_inter))>1){
      row_inter <- do.call(cbind,row_inter)
      row_inter <- t(row_inter)
      colnames(row_inter) <- tissueLabels
      # MH index 
      MH_index <- NULL
      #name <- NULL
      for (ii in ncol(row_inter):2) { # from tumor to con
        for (jj in (ii-1):1) {
          tmpp <- niche.overlap(row_inter[,c(ii,jj)],method = "morisita")
          MH_index <- c(MH_index,tmpp)
          #name <- c(name,paste0("MH_",colnames(row_inter)[ii],"_",colnames(row_inter)[jj]))
        }
      }
      
    }else{
      MH_index <- rep(NA,length(name))
    }
        
    names(MH_index) <- name
    return(MH_index)
  })
  Level1 <- do.call(rbind,Level1)
  rownames(Level1) <-names(LocMatrix)
  return(Level1) 
}


# get pearson correlation----------------------------------------------------------------------------
Cor_Index <- function(windows=8,tumorRatio=0.25,LocMatrix=EQ_loc_select,maskMatrix=mask_k53del_ErQi,
                    delBackLabel=1,tissueLabels=2:8,tumorLabel=8){
  require(spaa)
  Level1 <- mclapply(1:length(LocMatrix),function(i){
    tmp_mat <- LocMatrix[[i]]
    mask <- maskMatrix[[i]]
    
    squnR <- seq(1,nrow(mask)-round(windows/2),round(windows/2))
    lenR <- length(squnR)
    row_inter <- sapply(squnR[-lenR],function(row){
      # count <- vector("numeric",length=1)
      squn <- seq(1,ncol(mask)-round(windows/2),round(windows/2))
      len <- length(squn)
      all_1000 <- sapply(squn[-len],function(col){  # step = 4, if not enough for the last window, just discard it.
        mat_mask <- mask[row:(row+(windows-1)),col:(col+(windows-1))] 
        mat_4 <- tmp_mat[row:(row+(windows-1)),col:(col+(windows-1))] 
        
        if( length(which(mat_4==0 | mat_4==delBackLabel))<=round(windows*windows/2) & length(which(mat_mask==0))< round(windows*windows/2) & 
            length(which(mat_mask==1))>= round(windows*windows*tumorRatio) & length(which(mat_4==tumorLabel))>=round(windows*windows*tumorRatio) ){ 
          # mat_mask ==1 means tumor regions
          all_count <- sapply(tissueLabels, function(zzz){
            len <- length(which(mat_4==zzz))
            return(len)
          })
          return(all_count)
        }else{
          return(as.numeric(rep(-1,length(tissueLabels))))
        }
      }) # all_1000
      
      all_1000 <- all_1000[,-which(all_1000[1,]==-1)]
      if(length(dim(all_1000))>0 ){
        if(dim(all_1000)[2]>0){
          aver <- all_1000
        }else{
          aver <-  all_1000
        }
      }else{
        aver <-  as.numeric(rep(0,length(tissueLabels)))
      }
      return(aver)
    }) # row_inter
    
    row_inter <- row_inter[lapply(row_inter,sum)>0]  # !!
    name <- NULL
    for (ii in length(tissueLabels):2) { # from tumor to con
      for (jj in (ii-1):1) {
        name <- c(name,paste0("Cor_",tissueLabels[ii],"_",tissueLabels[jj]))
      }
    } 
    
    if(length((row_inter))>1){
      row_inter <- do.call(cbind,row_inter)
      row_inter <- t(row_inter)
      colnames(row_inter) <- tissueLabels # sample*7
      # cor 
      Cor_index <- NULL
      #name <- NULL
      for (ii in ncol(row_inter):2) { # from tumor to con
        for (jj in (ii-1):1) {
          tmpp <- cor(row_inter[,ii],row_inter[,jj],method = "pearson")
          Cor_index <- c(Cor_index,tmpp)
          #name <- c(name,paste0("Cor_",colnames(row_inter)[ii],"_",colnames(row_inter)[jj]))
        }
      }
      #names(Cor_index) <- name
    }else{
      Cor_index <- rep(NA,length(name))
    }
    names(Cor_index) <- name
    return(Cor_index)
  },mc.cores=20)
  Level1 <- do.call(rbind,Level1)
  rownames(Level1) <-names(LocMatrix)
  return(Level1) 
}

# get ROI K function----------------------------------------------------------------------------
ROI_KFunction <- function(windows=8,tumorRatio=0.25,LocMatrix=EQ_loc_select,maskMatrix=mask_k53del_ErQi,
                    specificLabels=2:8,tumorLabel=8){
  # cl <- makeCluster(10)
  # clusterExport(cl,varlist = c("windows","tumorRatio","LocMatrix","maskMatrix"))  
  require(spatstat)
  KLDiver_aver <- mclapply(1:length(LocMatrix),function(i){
    tmp_mat <- LocMatrix[[i]]
    mask <- maskMatrix[[i]]
    
    squnR <- seq(1,nrow(tmp_mat)-round(windows/2),round(windows/2))
    lenR <- length(squnR)
    row_inter <- sapply(squnR[-lenR],function(row){
      # count <- vector("numeric",length=1)
      squn <- seq(1,ncol(tmp_mat)-round(windows/2),round(windows/2))
      len <- length(squn)
      all_1000 <- sapply(squn[-len],function(col){  # step = 4, if not enough for the last window, just discard it.
        mat_4 <- tmp_mat[row:(row+(windows-1)),col:(col+(windows-1))] 
        mat_mask <- mask[row:(row+(windows-1)),col:(col+(windows-1))] 
        if( length(which(mat_4==0 | mat_4==1))<=round(windows*windows/2) & length(which(mat_mask==0))< round(windows*windows/2) & 
            length(which(mat_mask==1))>=round(windows*windows*tumorRatio) & length(which(mat_4==tumorLabel))>=round(windows*windows*tumorRatio) 
        ){  #  & length(which(mat_4==9))>=round(windows*windows*tumorRatio) 
          all_count <- sapply(specificLabels, function(zzz){
            len <- length(which(mat_4==zzz))
            return(len)
          })
          prop <- all_count/sum(all_count)
          return(prop)
        }else{
          return(as.numeric(rep(-111,length(specificLabels))) )
        }
      }) # all_1000
      
      all_1000 <- all_1000[,-which(all_1000[length(specificLabels),]==-111)]
      
      if(length(dim(all_1000))>1){
        rownames(all_1000) <- specificLabels
      }else{
        names(all_1000) <- specificLabels
      }
      return(all_1000)
    }) # row_inter
    row_inter <- row_inter[lapply(row_inter,sum)>0]
    row_inter <- do.call(cbind,row_inter)
    
    aver_kld <- apply(row_inter, 1, mean)
    KLD <- sapply(1:ncol(row_inter), function(xxx){
      tmp <- KL.plugin(row_inter[,xxx],aver_kld,unit="log2")
      return(tmp)
    })
    KLD_aver <- mean(KLD)
    return(KLD_aver)
    
  },mc.cores = 20)
  
  names(KLDiver_aver) <-names(LocMatrix)
  KLDiver_aver <- do.call(c,KLDiver_aver)
  
  return(KLDiver_aver) 
}

