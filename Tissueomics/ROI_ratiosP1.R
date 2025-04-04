# get ROI infil ratio----------------------------------------------------------------------------
ROI_infiltration_R <- function(windows=8,tumorRatio=0.25,LocMatrix=EQ_loc_select,maskMatrix=mask_k53del_ErQi,
                               delBackLabel=1,tissueLabels=2:8,tumorLabel=8,
                               tissueNames = c("con","epi","gla","lym","mus","str","tum")){
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
          all <- length(intersect(which(mat_4!=0),which(mat_4!=delBackLabel)))   # change from | to intersect# length(which(mat_4!=0 | mat_4!=1))
          all_count <- sapply(tissueLabels, function(zzz){
            len <- length(which(mat_4==zzz))
            return(len)
          })
          re <- as.numeric(all_count)/all
          return(re)
        }else{
          return(as.numeric(rep(-1,length(tissueLabels))))
        }
      }) # all_1000
      
      all_1000 <- all_1000[,-which(all_1000[1,]==-1)]
      if(length(dim(all_1000))>0 ){
        if(dim(all_1000)[2]>0){
          aver <- t(apply(all_1000, 1, mean))
        }else{
          aver <-  as.numeric(rep(0,length(tissueLabels)))
        }
      }else{
        aver <-  as.numeric(rep(0,length(tissueLabels)))
      }
      return(aver)
    }) # row_inter
    row_inter <- as.data.frame(row_inter)
    row_inter <- row_inter[,row_inter[length(tissueLabels),]!=0] # total row = 7, 7=tumor
    if(length(dim(row_inter))>0 ){
      aver_row <- t(apply(row_inter, 1, mean))
    }else{
      aver_row <-  as.numeric(rep(0,length(tissueLabels)))
    }
    return(aver_row)
  },mc.cores=20)
  Level1 <- do.call(rbind,Level1)
  rownames(Level1) <-names(LocMatrix)
  colnames(Level1) <- paste0("ROI_infil_",tissueNames)
  
  return(Level1) 
}

# get ROI intra ratio----------------------------------------------------------------------------
ROI_intratumor_R <- function(windows=8,tumorRatio=0.25,LocMatrix=EQ_loc_select,maskMatrix=mask_k53del_ErQi,
                             delBackLabel=1,tissueLabels=2:7,tumorLabel=8,
                             tissueNames = c("con","epi","gla","lym","mus","str")){
  Level1 <- lapply(1:length(LocMatrix),function(i){
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
        if( length(which(mat_4==0 | mat_4==delBackLabel))<= round(windows*windows/2) & length(which(mat_mask==0))< round(windows*windows/2) & 
            length(which(mat_mask==1))>= round(windows*windows*tumorRatio) & length(which(mat_4==tumorLabel))>=round(windows*windows*tumorRatio) ){ 
          
          all_count <- sapply(tissueLabels, function(zzz){
            len <- length(which(mat_4==zzz))
            return(len)
          })
          re <- as.numeric(all_count)/length(which(mat_4==tumorLabel))
          return(re)
        }else{
          return(as.numeric(rep(-1,length(tissueLabels))))
        }
      }) # all_1000
      
      all_1000 <- all_1000[,-which(all_1000[1,]==-1)]
      if(length(dim(all_1000))>0 ){
        if(dim(all_1000)[2]>0){
          aver <- t(apply(all_1000, 1, mean))
        }else{
          aver <-  as.numeric(rep(0,length(tissueLabels)))
        }
      }else{
        aver <-  as.numeric(rep(0,length(tissueLabels)))
      }
      return(aver)
    }) # row_inter
    
    row_inter <- as.data.frame(row_inter)
    row_inter <- row_inter[,row_inter[length(tissueLabels),]!=0]
    # row_inter <- row_inter[,rowSums(row_inter,na.rm = T)!=0] 
    if(length(dim(row_inter))>0 ){
      aver_row <- t(apply(row_inter, 1, mean))
    }else{
      aver_row <-  as.numeric(rep(0,length(tissueLabels)))
    }
    return(aver_row)
  })
  
  Level1 <- do.call(rbind,Level1)
  rownames(Level1) <-names(LocMatrix)
  colnames(Level1) <- paste0("ROI_intra_",tissueNames) # c( 'con', 'epi', 'gla', 'lym', 'mus', 'str', 'tum')
  
  return(Level1) 
}

# get ROI entropy----------------------------------------------------------------------------
ROIEntropy <- function(windows=8,tumorRatio=0.25,LocMatrix=EQ_loc_select,maskMatrix=mask_k53del_ErQi,
                       specificLabels=2:8,tumorLabel=8, delBackLabel=1){
  # cl <- makeCluster(10)
  # clusterExport(cl,varlist = c("windows","tumorRatio","LocMatrix","maskMatrix"))  
  entropy_aver <- lapply(1:length(LocMatrix),function(i){
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
        if( length(which(mat_4==0 | mat_4==delBackLabel))<=round(windows*windows/2) & length(which(mat_mask==0))< round(windows*windows/2) & 
            length(which(mat_mask==1))>=round(windows*windows*tumorRatio) & length(which(mat_4==tumorLabel))>=round(windows*windows*tumorRatio) 
        ){  #  & length(which(mat_4==9))>=round(windows*windows*tumorRatio) 
          all_count <- sapply(specificLabels, function(zzz){
            len <- length(which(mat_4==zzz))
            return(len)
          })
          entro <- entropy::entropy(all_count)
          return(entro)
        }else{
          return(c(-111))
        }
      }) # all_1000
      
      all_1000 <- all_1000[-which(all_1000==-111)]
      if(length(all_1000)>0 ){
        aver <- mean(all_1000)
      }else{
        aver <-  0
      }
      return(aver)
    }) # row_inter
    
    row_inter <- row_inter[row_inter!=0]
    if(length(row_inter)>0 ){
      aver_row <- mean(row_inter)
    }else{
      aver_row <-  0
    }
    return(aver_row)
  })
  names(entropy_aver) <-paste0("Entropy_",names(LocMatrix))
  entropy_aver <- do.call(c,entropy_aver)
  
  return(entropy_aver) 
}


# get ROI KL divergence----------------------------------------------------------------------------
ROI_KLD <- function(windows=8,tumorRatio=0.25,LocMatrix=EQ_loc_select,maskMatrix=mask_k53del_ErQi,
                       specificLabels=2:8,tumorLabel=8, delBackLabel=1){
  # cl <- makeCluster(10)
  # clusterExport(cl,varlist = c("windows","tumorRatio","LocMatrix","maskMatrix"))  
  require(entropy)
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
        if( length(which(mat_4==0 | mat_4==delBackLabel))<=round(windows*windows/2) & length(which(mat_mask==0))< round(windows*windows/2) & 
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
    
    if(length(dim(row_inter))>1){
      aver_kld <- apply(row_inter, 1, mean)
      KLD <- sapply(1:ncol(row_inter), function(xxx){
        tmp <- KL.plugin(row_inter[,xxx],aver_kld,unit="log2")
        return(tmp)
      })
      KLD_aver <- mean(KLD,na.rm=T)
    }else{
      aver_kld <- row_inter
      KLD_aver <- KL.plugin(row_inter,aver_kld,unit="log2")
    }
    return(KLD_aver)
    
  },mc.cores = 30)
  
  names(KLDiver_aver) <-paste0("KLDiver_",names(LocMatrix))  
  KLDiver_aver <- do.call(c,KLDiver_aver)
  
  return(KLDiver_aver) 
}


# get ROI within----------------------------------------------------------------------------
ROIWithin <- function(windows=8,tumorRatio=0.25,LocMatrix=EQ_loc_select,maskMatrix=mask_k53del_ErQi,
                      delBackLabel=1,tissueLabels=2:7,tumorLabel=8,
                      tissueNames = c("con","epi","gla","lym","mus","str") ){
  Level1 <- mclapply(1:length(LocMatrix),function(i){
    ref_mat <- matrix(0,windows,windows)
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
    
    tmp_mat <- LocMatrix[[i]]
    mask <- maskMatrix[[i]]
    
    squnR <- seq(1,nrow(tmp_mat)-round(windows/2),round(windows/2))
    lenR <- length(squnR)
    row_inter <- sapply(squnR[-lenR],function(row){
      # count <- vector("numeric",length=1)
      squn <- seq(1,ncol(tmp_mat)-round(windows/2),round(windows/2))
      len <- length(squn)
      all_1100 <- sapply(squn[-len],function(col){  # step = 4, if not enough for the last window, just discard it.
        mat_4 <- tmp_mat[row:(row+(windows-1)),col:(col+(windows-1))] 
        mat_mask <- mask[row:(row+(windows-1)),col:(col+(windows-1))] 
        if( length(which(mat_4==0 | mat_4==delBackLabel))<=round(windows*windows/2) & length(which(mat_mask==0))< round(windows*windows/2) & 
            length(which(mat_mask==1))>=round(windows*windows*tumorRatio) &  length(which(mat_4==tumorLabel))>=round(windows*windows*tumorRatio) ){ 
          
          all_count <- sapply(tissueLabels, function(zzz){
            if(length(which(mat_4==zzz))> 1){ # at least two
              con_dis <- covert2mat(which(mat_4==zzz))
              average_con_dis <- sum(dist(con_dis))/length(dist(con_dis))
            }else{
              average_con_dis <- -1
            }
            return(average_con_dis)
          })
          all_count <- as.numeric(all_count)
          return(all_count)
          
        }else{
          return(as.numeric(rep(-1,length(tissueLabels)))) # ,-1,-1
        }
      }) # all_1000
      
      if(length(which(all_1100[length(tissueLabels),]==-1))!=ncol(all_1100)){  # whether any windows meeting requires
        all_1000 <- all_1100[,all_1100[length(tissueLabels),]!=-1]
        if(length(dim(all_1000))>0){ # only one or more
          aver <- t(apply(all_1000, 1, function(x){
            if(length(which(x!=-1))!=0){
              tmp <- x[which(x!=-1)]
              tmp <- mean(tmp)
            }else{
              tmp <- -1
            }
            return(tmp)
          }))
        }else{
          aver <-  as.numeric(rep(-1,length(tissueLabels)))
        }         
        
      }else{
        aver <-  as.numeric(rep(-1,length(tissueLabels)))
      }
      return(aver)
    }) # row_inter
    
    if(length(which(row_inter[length(tissueLabels),]==-1))!=ncol(row_inter)){  # whether any windows meeting requires
      row_inter <- row_inter[,row_inter[length(tissueLabels),]!=-1] 
      if(length(dim(row_inter))>0){ # only one or more
        aver_row <- t(apply(row_inter, 1, function(x){
          if(length(which(x!=-1))!=0){
            tmp <- x[which(x!=-1)]
            tmp <- mean(tmp)
          }else{
            tmp <- NA
          }
          return(tmp)
        }))
      }else{
        aver_row <-  row_inter
      }  
    }else{
      aver_row <- as.numeric(rep(NA,length(tissueLabels)))
    }
  },mc.cores=20)
  
  Level1 <- do.call(rbind,Level1)
  rownames(Level1) <-names(LocMatrix)
  colnames(Level1) <- paste0("ROI_Within_",tissueNames) # "MUC", "NORM",
  return(Level1) 
}

# get ROI interaction (distance based: between cluster distance)--------------------------------
ROIInter <- function(windows=8,tumorRatio=0.25,LocMatrix=EQ_loc_select,maskMatrix=mask_k53del_ErQi,
                    delBackLabel=1,tissueLabels=2:8,tumorLabel=8,
                     tissueNames = c("con","epi","gla","lym","mus","str","tum")){
  distance <- mclapply(1:length(LocMatrix),function(i){
    
    ref_mat <- matrix(0,windows,windows)
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
    
    tmp_mat <- LocMatrix[[i]]
    mask <- maskMatrix[[i]]
    squnR <- seq(1,nrow(tmp_mat)-round(windows/2),round(windows/2))
    lenR <- length(squnR)
    row_inter <- sapply(squnR[-lenR],function(row){
      squn <- seq(1,ncol(tmp_mat)-round(windows/2),round(windows/2))
      len <- length(squn)
      all_1100 <- sapply(squn[-len],function(col){
        mat_4 <- tmp_mat[row:(row+(windows-1)),col:(col+(windows-1))] 
        mat_mask <- mask[row:(row+(windows-1)),col:(col+(windows-1))] 
        all_aver <- NULL
        
        if( length(which(mat_4==0 | mat_4==delBackLabel))<=round(windows*windows/2) & length(which(mat_mask==0))< round(windows*windows/2) & 
            length(which(mat_mask==1))>=round(windows*windows*tumorRatio) &  length(which(mat_4==tumorLabel))>=round(windows*windows*tumorRatio) ){
          for (j in length(tissueLabels):2){
            for (zz in (j-1):1){
              if(sum(mat_4==tissueLabels[j])> 1){
                con_dis <- covert2mat(which(mat_4==tissueLabels[j]))
                if(length(which(mat_4==tissueLabels[zz]))>0){
                  con_tmp <- covert2mat(which(mat_4==tissueLabels[zz]))
                  dist_j <- sapply(1:nrow(con_dis), function(mm){
                    tmp <- dist(rbind(con_dis[mm,],con_tmp))
                    tmp <- tmp[1:nrow(con_tmp)]
                    tmp <- mean(tmp)
                    return(tmp)
                  })
                  dist_j <- mean(as.numeric(dist_j))
                  all_aver <- c(all_aver,dist_j)
                }else{all_aver <- c(all_aver,NA)}
              }else{all_aver <- c(all_aver,NA)}
              #name <- c(name,paste0("inter_",tissueLabels[j],"_",tissueLabels[zz]))
            }
          }
          #names(all_aver) <- name
        }else{
          all_aver <- rep(NA,(length(tissueLabels)*(length(tissueLabels)-1)/2))
        }
        
        return(all_aver)
      })
      all_1100 <- t(all_1100)
      #colnames(all_1100) <- name  
      all_1100 <- all_1100[-which(rowSums(all_1100,na.rm =T)==0),]
      return(all_1100)
    })
    row_inter <- row_inter[lapply(row_inter, sum, na.rm=T)>0]
    row_inter <- do.call(rbind,row_inter)
    
    if(length(dim(row_inter))>1){
      row_inter <- apply(row_inter, 2, mean,na.rm=T)
    }else{
      row_inter <- row_inter
    }
    
    return(row_inter)
  },mc.cores = 40)
  
  name <- NULL
  for (j in length(tissueNames):2){
    for (zz in (j-1):1){
      name <- c(name,paste0("inter_",tissueNames[j],"_",tissueNames[zz]))
    }
  }
  names_P <- names(LocMatrix)[lapply(distance,length)!=0]
  distance <- do.call(rbind,distance)
  colnames(distance) <- name
  rownames(distance) <- names_P
  return(distance)
}