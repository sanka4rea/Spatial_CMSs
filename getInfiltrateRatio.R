getInfiltrateRatio <- function(windows=6,overlap=0.5,locationMatrix){
  
  
  Ratio_aver <- sapply(1:length(locationMatrix),function(i){
    tmp_mat <- locationMatrix[[i]]
    
    squnR <- seq(1,nrow(tmp_mat)-windows*overlap,windows*overlap)
    lenR <- length(squnR)
    row_inter <- sapply(squnR[-lenR],function(row){
      squn <- seq(1,ncol(tmp_mat)-windows*overlap,windows*overlap)
      len <- length(squn)
      Ratio <- sapply(squn[-len],function(col){  
        mat_4 <- tmp_mat[row:(row+(windows-1)),col:(col+(windows-1))] 
        if( length(which(mat_4==0 | mat_4==2))<= 18 & length(which(mat_4==9))>=10 ){ 
          
          all <- length(which(mat_4!=0 | mat_4!=2))
          den1tr <- length(which(mat_4==1))
          # den2tr <- length(which(mat_4==2))/all  # no back
          den3tr <- length(which(mat_4==3))
          denLYMtr <- length(which(mat_4==4))
          den5tr <- length(which(mat_4==5))
          den6tr <- length(which(mat_4==6))
          den7tr <- length(which(mat_4==7))
          denSTRtr <- length(which(mat_4==8))
          denTUMtr <- length(which(mat_4==9))
          re <- c(den1tr,den3tr,denLYMtr,den5tr,den6tr,den7tr,denSTRtr,denTUMtr)/all
          return(re)
        }else{
          return(c(-1,-1,-1,-1,-1,-1,-1,-1))
        }
      }) 
      
      Ratio <- Ratio[,-which(Ratio[1,]==-1)]
      if(sum(Ratio)>0 ){
        if(length(dim(Ratio))>0){
          aver <- t(apply(Ratio, 1, mean))
        }else{
          aver <- Ratio
        }
      }else{
        aver <-  c(0,0,0,0,0,0,0,0)
      }
      return(aver)
    })
    
    row_inter <- row_inter[,row_inter[8,]!=0]
    if(sum(row_inter)>0 ){
      if(length(dim(row_inter))>0){
        aver_row <- t(apply(row_inter, 1, mean))
      }else{
        aver_row <-  row_inter
      }
    }else{
      aver_row <-  c(0,0,0,0,0,0,0,0)
    }
    return(aver_row)
  })
  
  colnames(Ratio_aver) <-names(CRC_locationMatrix)
  rownames(Ratio_aver) <- c("ADI","DEB","LYM","MUC","MUS","NORM","STR","TUM")
  
  Ratio_w6_t50_tcga <- Ratio_aver[,Ratio_aver[8,]!=0]
  return(Ratio_w6_t50_tcga)
}