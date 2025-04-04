# check whether one slides contains two biopsy tissues (not support >2)
# cause the loc matrix has delete last right and last bottom, so we calculated margin from right to left, or bottom to top
CheckMultiSlides <- function(LocMatrix=EQ_loc_select,delBackLabel=1,minishape=2000){
  
  new_loc_matrix <- mclapply(1:length(LocMatrix), function(x){# length(LocMatrix) 
    mat <- LocMatrix[[x]]
    
    roworcol <- nrow(mat)>ncol(mat) # true: row false:col
    if(roworcol){ # larger row
      
      # arow0 <- which(rowSums(mat)==0)
      arow0 <- apply(mat, 1, function(y){ # row
        tmp <- length(which(y==0 | y==delBackLabel)) == ncol(mat)
        return(tmp)
      } )
      arow0 <- which(arow0==TRUE)
      arow <- sort(arow0,decreasing = T) # right to left
      diffrow <- diff(arow)
      indexrow <-which(diffrow< (-5))
      indexrow <- c(arow[indexrow],arow[indexrow+1])
      indexrow <- sort(c(nrow(mat),indexrow,1),decreasing = T)  
      
      validtmp <- sapply(1:(length(indexrow)-1), function(y){
        tmp <- mat[indexrow[y+1]:indexrow[y],]
        tmp <- length(which(tmp!=0 & tmp!=delBackLabel))
        return(tmp)
      })
      selec <- which(validtmp>minishape)
      
      if(length(selec)==0){
        return(0)
      }else{
        selecmat <- lapply(selec, function(y){
          tmp_mat <- mat[indexrow[y+1]:indexrow[y],]
          return(tmp_mat)
        })
        name <- as.vector(1:length(selecmat))
        names(selecmat) <- paste0(names(LocMatrix)[x],"part",name)
        return(selecmat)
      }
      
    }else{ # larger col
      
      acol0 <- apply(mat, 2, function(y){ # col
        tmp <- length(which(y==0 | y==delBackLabel)) == nrow(mat)
        return(tmp)
      } )
      acol0 <- which(acol0==TRUE)
      acol <- sort(acol0,decreasing = T) # right to left
      diffcol <- diff(acol)
      indexcol <-which(diffcol< (-5))
      indexcol <- c(acol[indexcol],acol[indexcol+1])
      indexcol <- sort(c(ncol(mat),indexcol,1),decreasing = T)  
      
      validtmp <- sapply(1:(length(indexcol)-1), function(y){
        tmp <- mat[,indexcol[y+1]:indexcol[y]]
        tmp <- length(which(tmp!=0 & tmp!=delBackLabel))
        return(tmp)
      })
      selec <- which(validtmp>minishape)
      
      if(length(selec)==0){
        return(0)
      }else{
        selecmat <- lapply(selec, function(y){
          tmp_mat <- mat[,indexcol[y+1]:indexcol[y]]
          return(tmp_mat)
        })
        name <- as.vector(1:length(selecmat))
        names(selecmat) <- paste0(names(LocMatrix)[x],"part",name)
        return(selecmat)
      }
    }
    
    
  },mc.cores = 20) # 
  
  new_loc_matrix <- do.call(c,new_loc_matrix)
  return(new_loc_matrix)
}