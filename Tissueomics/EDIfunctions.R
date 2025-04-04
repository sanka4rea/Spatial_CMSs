getCellFiles <- function(ff){
  cells <- scan(file=ff, what=c('numeric', 'character', rep('character', 105)), fill=TRUE, quiet=TRUE)
  n <- length(cells)
  cells <- matrix(cells, ncol=105, byrow=TRUE)
  if (n%%105 != 0 )
    cells <- cells[-nrow(cells),]
  colname <- cells[1,-c(1:2)]
  rowname <- cells[-1,2]
  cells <- cells[-1, -c(1:2)]
  cells <- matrix(as.numeric(cells), ncol=103, dimnames=list(rowname,colname)) 
  cells
}

getRegionalComposition <- function(ff, imgDir, outputDir='./', n=400){
# This function computes the number of cells for each cell type in tumour regions based on CRImage result
# ff: sample folder name
# imgDir: directory where CRImage output is stored
# outputDir: directory where output of this function is stored  
# n: region size in number of pixels   
  m <- 2000/n
  ff1 <- dir(paste(imgDir,'/', ff,  '/Cells/', sep=''))
  res <- NULL
  for (cc in ff1){
    tab <- try(getCellFiles(paste(imgDir, '/', ff, '/Cells/', cc, sep='')))
    if(class(tab)!='try-error')
      if (nrow(tab)>5){
        x <- tab[,'g.x']%/%n+1
        y <- tab[,'g.y']%/%n+1
        pos <- factor(m*(x-1)+y, levels=1:m^2)
        cls <- factor(rownames(tab), levels=c('a','c', 'l','o'))
        cellDen <- table(cls, pos)
        colnames(cellDen) <- rep(cc, ncol(cellDen))
        res <- rbind(res, t(cellDen))
      }}
  write.table(res, paste(outputDir,'/', ff, '.txt', sep=''), quote=F, sep='\t')
}

getEDI <- function(ff, regStatDir, n=400, q=0.8, minCell=NULL, maxK=5, sampling=NULL){
# This function computes EDI score based on data extracted by function getRegionalComposition
# ff: sample folder name
# regStatDir: directory where getRegionalComposition output is stored
# n: region size in number of pixels   
# q: quantile of data to compute the score
# minCell: minimal number of cells per region; by default it is n/25
# maxK: maximum number of clusters expected
# sampling: if specified this should be a percentage such as 0.9 for sampling purpose
  require(mclust)
  require(vegan)
  if (is.null(minCell))
    minCell <- n/25
  res <- try(read.table(paste(regStatDir, ff, '.txt', sep=''), as.is=T, sep='\t',row.names=NULL))
  
  if(class(res)!='try-error'){
    if (!is.null(sampling))
      res <- res[sample(1:nrow(res), nrow(res)*sampling), ]
    res <- res[rowSums(res[,3:5]) >= minCell, ]
    dvs <- diversity(res[,3:5], index = "shannon", MARGIN = 1)
    x <- dvs[dvs<quantile(dvs, q, na.rm=T)] # 0-80
    x <- x[x>quantile(dvs, 1-q, na.rm=T)] # 20-80
    x <- x[!is.na(x)]
    if (length(x)>3){
        s <- Mclust(x, G=1:maxK)$G
      }else{
        s <- NA
      }
  }else{
    s <- NA
  }
  s
}

group2 <- function(x, th=.5)
  (x > quantile(x, th, na.rm=T))*1+1



color.code <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000"))

group3 <- function(x, q=c(.25, .75)){
  q <- quantile(x, q, na.rm=TRUE)
  (x > q[1])*1+(x>q[2])*1+1
}

Boxplot <- function(f, main=NULL, subset=NULL, ...){
  ks <- do.call("aov", list(formula=f, subset=subset, ...))
  p <- signif(summary(ks)[[1]][1,5],2)
  do.call("boxplot", list(formula=f, main=paste(main,'p=', p), varwidth=TRUE, subset=subset,...))
}
Boxplotkw <- function(f, xlab=NULL, ylab=NULL, main=NULL, subset=NULL, ...){
  ks <- do.call("kruskal.test", list(formula=f, subset=subset, ...))
  p <- signif(ks$p.value,2)
  do.call("boxplot", list(formula=f, main=paste(main,'p=', p), ylab=ylab, xlab=xlab, varwidth=TRUE, subset=subset,...))
}

replace.vector <-function(x, tochange=unique(x), toreplace=1:length(unique(x))){
	for (y in 1:length(tochange))  x=replace(x, list=(x==tochange[y]), toreplace[y])
	x
}


plotSurv <- function(S, clusters, name='Cluster', type='', file='Survival', fileType='', rho=0, width=NULL, height=NULL, col=NULL, lwd=2, lty=1, legendpos='bottomleft', SurvType='Survival probability', ...){
  library(survival)
  if(is.null(col)){
    col <- color.code(20)
    set.seed(1)
    col <- sample(col)
  }
  dat <- data.frame(x=clusters, S=S)
  fit <- survfit(S ~ x,data=dat)
  test <- survdiff(S ~ x, data=dat, rho=rho)
  p.val <- 1 - pchisq(test$chisq, length(test$n) - 1)
  if (fileType=='png')
    png(paste(file, '.', fileType, sep=''), width=ifelse(is.null(width), 600, width), height=ifelse(is.null(height), 600, height))
  if (fileType=='pdf')
    pdf(paste(file, '.', fileType, sep=''), width=ifelse(is.null(width), 6, width), height=ifelse(is.null(height), 6, height))
  plot(fit, lwd=lwd, lty=lty, main=paste(name, type, ' p=', signif(p.val,2), sep=''), xlab='Months', ylab=SurvType, cex.lab=1.5, col=col, ...)
  if (name!='')
      legend(legendpos, legend = paste(name, sapply(names(test$n), function(x) substr(x, 2, nchar(x))), ': ', test$n, '(', test$obs, ')', sep=''),  lty=lty, lwd=lwd, inset=c(0.01, 0.01), col=col)
  if (name=='')
        legend(legendpos, legend = paste(sapply(names(test$n), function(x) strsplit(x, split='=', fixed=T)[[1]][2]), ': ', test$n, '(', test$obs, ')', sep=''),  lty=lty, lwd=lwd, inset=c(0.01, 0.01), col=col)

  if (fileType%in%c('png', 'pdf')) dev.off()
  p.val
}



plot.normal.components.mclust <- function(x, mixture,...) 
# Function to add Gaussian mixture components for mclust object, vertically scaled, to the current plot
# Presumes the mixture object has the structure used by mixtools
  sapply(1:mixture$G,function(y)
         curve(mixture$parameters$pro[y] * dnorm(x,mean=mixture$parameter$mean[y], sd=(ifelse(length(mixture$parameter$variance$sigmasq)==1, mixture$parameter$variance$sigmasq, mixture$parameter$variance$sigmasq[y]))^(1/2)), add=TRUE, ...))


colorBar <- function(x, f='colorbar.pdf', col=NULL, digit=1, br=NULL, anno=NULL, brief=TRUE, ...){
  pdf(f, height=3, width=2)
  par(mar=c(.3,6,.3,3))
  if (is.null(br))
    br <- quantile(x, prob=seq(0,1,len=21), na.rm=TRUE)
  n <- length(br)
  if (is.null(col))
    col <- color.code(n-1)

  if(brief){
    br <- br[c(seq(1, n, by=3),n)]
    col <- col[seq(1, n, by=3)]
    n <- length(br)
  }
  image(matrix(1:(n-1), nrow=1), col=col, axes=FALSE)
  axis(4, at=(1:n -1.5)/(n-2), paste(signif(br,digit), sep=''), las=2, cex=.7)
  if (!is.null(anno))
    axis(2, anno, at=(1:(n-1) -1)/(n-2), anno, las=2, cex=.7)
  if(f!='')
    dev.off()
}


plotRegionDiv <- function(imgDir, outputDir, ff, compositionDir, plotStats=FALSE, plotGridded=FALSE, q=0.8, h=40){
#imgDir: input directory; usually is a folder generated by the getCellPos function; it should contain a "CellPosAndMask" folder where cell position files are stored; a "OutputImage" folder where scaled down H&E images are stored.
#outputDir: output directory
#compositionDir: tumour region composition folder with files generated by getRegionalComposition function
#ff: file name
#h: scaling the image with scale factor; the larger the factor the smaller the output image 
#q: same as
#plotStats: if true will plot result statistics of clustering
#plotGridded: if true will plot a gridded heatmap and H&E in pdf
require(vegan)
require(mclust)
require(EBImage)
     
    tmp <- try(load(paste0(imgDir, '/CellPosAndMask/', ff, '.rdata')))
    if(class(tmp)!='try-error'){
        n <- nrow(Mask)/h
        m <- ncol(Mask)/h
            
        if(! paste0(ff, '.rdata') %in% dir(outputDir)){
            mat <- matrix(-2, nrow=n, ncol=m)
            CellPos <- CellPos[CellPos[,1]!='a',]
            for (i in 1:n)
                for (j in 1:m){
                    idx <- ( CellPos$x %in% ((i-1)*h+1): (i*h) )  & ( CellPos$y %in% ((j-1)*h+1): (j*h) )
                    if (sum(idx)>16)
                        mat[i,j] <- diversity(table(CellPos[idx,1]), index = "shannon", MARGIN = 1, base = exp(1))
                    if(sum(idx)>0 & sum(idx)<=16)
                        mat[i,j] <- -1
                }
            save(mat, file=paste0(outputDir, '/', ff, '.rdata'))
        }
        
        load(paste0(outputDir, '/', ff, '.rdata'))
        res <- try(read.table(paste0(compositionDir, '/', ff, '.txt'),as.is=T, sep='\t',row.names=NULL))
        res <- data.frame(res, n=rep(1:25, times=nrow(res)%/%25))
        dvs <- diversity(res[,3:5], index = "shannon", MARGIN = 1, base = exp(1))
        dvs[rowSums(res[,3:5]) < 16] <- -1
        res <- cbind(res,dvs)
        
        x <- dvs[dvs!=-1]
        x1 <- x[x<quantile(x, q, na.rm=T) & x>quantile(x, 1-q, na.rm=T)]  
        cl <- Mclust(x1, G=1:5)
        k <- rep(NA, length(dvs))
        k[dvs!=-1][x<quantile(x, q, na.rm=T) & x>quantile(x, 1-q, na.rm=T)] <- cl$classification
        k[dvs!=-1 & is.na(k)] <- 0
        res$k <- k
        write.csv(res, file=paste0(outputDir, '/', ff, '.csv'))   
        
        img <- readImage(paste0(imgDir,'/OutputImage/',ff,'OrginialImage.jpg'))
        cols <- colorRampPalette(c('yellow', 'orange', 'red', 'darkred'))(max(cl$classification))
        col <- c('grey','white', 'white', cols, 'white')
        br <- sort(aggregate(x1, list(cl$classification), min)[,2])
        br <- c(-3,-1.5, -.5, br, quantile(x, q, na.rm=T),2)
        if(plotStats){
            pdf(paste0(outputDir, '/', ff,'stats.pdf'), width=n/8, height=m/8)
            plot(cl, 'density')
            plot(cl, 'classification')
            plot(hist(x,breaks=30, plot=F),col="grey",border="grey",freq=FALSE,
                 xlab="Microenvironmental diversity",main="")
            lines(density(x1),lty=2)
            plot.normal.components.mclust(x1, cl)
            par(mar=c(3,3,3,3))
            par(mfrow=c(1,2))
            plot(rowSums(res[res$dvs!=-1,3:5]), res$dvs[res$dvs!=-1], main=paste('Before filtering cor=', signif(cor(rowSums(res[res$dvs!=-1,3:5]), res$dvs[res$dvs!=-1]),2)), pch=19, cex=.5)
            plot(rowSums(res[res$dvs!=-1 & res$k!=0,3:5]), res$dvs[res$dvs!=-1 & res$k!=0], main=paste('After filtering cor=', signif(cor(rowSums(res[res$dvs!=-1 & res$k!=0,3:5]), res$dvs[res$dvs!=-1 & res$k!=0]),2)), pch=19, cex=.5)
            dev.off()
        }
        if(plotGridded){
            pdf(paste0(outputDir, '/', ff,'Gridded.pdf'), width=n/8, height=m/8)
            par(mar=c(0,0,0,0))
            plot(x=c(0,1), y=c(0,1), pch='', xaxs = "i", yaxs = "i")
            image(mat[,ncol(mat):1], breaks=br,  col=col, add=T)
            par(mar=c(0,0,0,0))
            plot(x=c(1,n), y=c(1,m), pch='', xaxs = "i", yaxs = "i")
            rasterImage(aperm(img, c(2,1,3)), xleft=1,xright=n,ybottom=1,ytop=m)         
            grid(nrow(mat), ncol(mat))
            
            mask <- mat > br[4] & mat < br[9]
            mask <- channel(mask, 'gray')
            mask <- resize(mask, nrow(img), ncol(img))
            img1 <- paintObjects(mask, img,col='black')

            plot(x=c(0,n), y=c(0,m), pch='', xaxs = "i", yaxs = "i")
            rasterImage(aperm(img1, c(2,1,3)), xleft=1,xright=n,ybottom=1,ytop=m)            
            dev.off()
        }
        
        png(paste0(outputDir, '/', ff, '.png'), , width=n*10, height=m*20)
        par(mfrow=c(2,1))
        par(mar=c(0,0,0,0))
        plot(x=c(0,1), y=c(0,1), pch='', xaxs = "i", yaxs = "i")
        image(mat[,ncol(mat):1], breaks=br,  col=col, add=T)
        par(mar=c(0,0,0,0))
        plot(x=c(1,n), y=c(1,m), pch='', xaxs = "i", yaxs = "i")
        rasterImage(aperm(img, c(2,1,3)), xleft=1,xright=n,ybottom=1,ytop=m)
        dev.off()
        
        colorBar(x, f=paste0(outputDir, '/', ff,'colorbar.pdf'), br=br, digit=2, anno=c('No tissue', 'Few cells', 'Lower 20%', paste0('K', 1:max(cl$classification)), 'Higher 20%'), col=col, brief=FALSE)
    }
}

plotRegionDiv3D <- function(imgDir, outputDir, ff, q=0.8){
    require(rgl)
    require(sm)
    load(paste0(outputDir, '/', ff, '.rdata'))
    
    i <- which(mat>q)
    x <- i%%nrow(mat)+1
    y <- ncol(mat)- i%%ncol(mat)

    data <- data.frame(x=x, y=y)
    data$y <- data$y*10    
    res <- sm.density(cbind(data$x, data$y), h=c(10, 10), display='contour', eval.points=cbind(seq(1, max(data$x), len=200), seq(1, max(data$y), len=200)))
    points(data$x, data$y, cex=.5)
    z <- res$estimate


    z <- mat
    z[z== -1] <- 0
    z[z<0.1] <- 0
    z <- 2^(2^z)    # Exaggerate the relief
    z[z==2] <- 0
    x <- 10 * (1:nrow(z)) # 10x spacing (S to N)
    y <- 10 * (1:ncol(z)) # 10x spacing (E to W)
    zlim <- range(z)
    zlen <- zlim[2] - zlim[1] + 1
    colorlut <- colorRampPalette(c('white', 'lightblue', 'blue'))(zlen) #terrain.colors(zlen,alpha=0)
    col <- colorlut[ z-zlim[1]+1 ] # assign colors to heights for each point
    img <- readImage(paste0(imgDir,'/OutputImage/',ff,'OrginialImage.jpg'))
    writeImage(img, file=paste0(outputDir,'/', ff,'OrginialImage.png'))
    
    open3d()
    z[z!=0]=rnorm(length(z[z!=0]), 10, 5)
    rgl.surface(x, y, z*0.000001-50, texture=paste0(outputDir, '/',ff,'OrginialImage.png'), alpha=.8,   lit=T, add=F)
    rgl.surface(x, y, 2*z, texture=paste0(outputDir, '/',ff,'OrginialImage.png'), color=col, alpha=.6,  lit=T, add=T)
    rgl.snapshot(paste0(outputDir, '/', ff, "3D.png"), "png")
}
