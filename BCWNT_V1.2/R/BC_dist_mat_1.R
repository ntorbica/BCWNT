
# Function to evaluate corrected data with bhattacharyya distances.
# Outputs a diagonal matrix containing all of the possible distances,
# Intended use with principal components.

BCW.dist.mat <- function(d, heatmap = TRUE, PCA.plots = TRUE, PCA.c = 5, PDF = FALSE, csv = FALSE, t = ''){
  # d               - Batch-Corrected data filename, containing information columns and peaks (BCW output)
  # heatmap         - logical, defines whether a heatmap is plotted
  # PCA.plots       - logical, gives a distance plot with means of principal components, which is another visualization next
  #                   to the heatmap
  # PCA.c           - Count of PCs to plot (PCA.plots argument)
  # PDF             - logical, defines whether a PDF should be produced
  # csv             - logical, defines whether a csv of the distances should be produced
  # t               - title for the plots

  suppressMessages(require(BatchCorrMetabolomics))
  suppressMessages(require(RUVSeq))
  suppressMessages(require(ChemometricsWithR))
  suppressMessages(require(fpc))
  suppressMessages(require(gplots))

  cat('Preparing data...\n')

  # Check input type
  type <- class(d)

  if(type == 'character'){
    d <- read.csv(d,row.names = 1)
  }

  if(type == 'list'){
    d <- as.data.frame(d)
  }

  info <- d[,c('SeqNr','SCode','Batch')]
  peaks <- d[,-(which(colnames(d)%in%colnames(info)))]

  nb <- nlevels(info$Batch)

  # check which indices apply to QCs and remove them from the set
  noref.idx <- which(info$SCode != "ref")
  Xsample <- peaks[noref.idx,]
  YSample <- info[noref.idx,]

  Xsample <- Xsample[, apply(Xsample, 2, function(x) !all(is.na(x)))]
  Xsample <- Xsample[, apply(Xsample, 2, sd, na.rm = TRUE) > 0]

  # replace NAs with  column means, given there are NAs
  for (i in 1:ncol(Xsample)){
    Xsample[is.na(Xsample[,i]),i] <- mean(Xsample[,i], na.rm = TRUE)
  }

  # get principal components
  Xsample <- scale(Xsample)
  X.PCA <- PCA(Xsample)

  # Get the number of total PCs
  npc <-ncol(X.PCA$scores)

  ## make a separate object from the score-component of X.PCA
  Xscores <- scores.PCA(X.PCA)[, 1:npc, drop = FALSE]

  ## These are Wehrens loops that give the batch mean vectors
  ## and covariance matrices.
  batch.means <-
    lapply(levels(YSample$Batch),
           function(btch)
             colMeans(Xscores[which(YSample$Batch == btch),,drop=FALSE]))

  batch.covs <-
    lapply(levels(YSample$Batch),
           function(btch)
             cov(Xscores[which(YSample$Batch == btch),,drop=FALSE]))

  ## So here we account for the fact that bhattacharyya.dist uses
  ## mahalanobis distance computation, which uses solve(cov,...),
  ## which is a (reasonable) bitch when doing matrix transposition
  ## and throws singular matrices due to small values.
  noCov.idx <- which(sapply(batch.covs, function(x) all(x < 1e-8)))

  b1m <- t(matrix(unlist(batch.means[[1]]),1,npc))
  rownames(b1m) <- sprintf('PC%d', 1:npc)
  colnames(b1m) <- sprintf('Batch %d', 1)

  b2m <- t(matrix(unlist(batch.means[[2]]),1,npc))
  rownames(b2m) <- sprintf('PC%d', 1:npc)
  colnames(b2m) <- sprintf('Batch %d', 2)


  b1c <- matrix(0,npc,npc)
  b1c.v <- unlist(batch.covs[1])
  for(i in 1:length(b1c.v)){
    b1c[i] <- b1c.v[i]
  }
  colnames(b1c) <- sprintf('PC%d',1:npc)
  rownames(b1c) <- sprintf('PC%d',1:npc)

  b2c <- matrix(0,npc,npc)
  b2c.v <- unlist(batch.covs[2])
  for(i in 1:length(b2c.v)){
    b2c[i] <- b2c.v[i]
  }
  colnames(b2c) <- sprintf('PC%d',1:npc)
  rownames(b2c) <- sprintf('PC%d',1:npc)


  ## EDIT: In order to fulfill the chosen approach, I will have to make
  ## groupings of the covariance matrices, as it has 24x24 format, while
  ## we only can use a 2x2 matrix for the comparison of two batches
  ## due to matrix multiplication logic (cols of the first matrix have
  ## to be equal to the rows of the second).ov.pairlist1 <- list()
  cov.pairlist1 <- list()
  cov.pairlist2 <- list()
  x <- 1
  for(i in 1:(npc-1)){
    for(j in (i+1):npc){
      #print(b1c[c(i,j),c(i,j)])
      cov.pairlist1[[x]] <- matrix(b1c[c(i,j),c(i,j)],2,2)
      cov.pairlist2[[x]] <- matrix(b2c[c(i,j),c(i,j)],2,2)
      x <- x+1
    }
  }

  ## The final preparation steps include:

  # wrap up the mean-matrix
  bm <- cbind(b1m,b2m)

  # wrap up the covariance-array
  sigmaarray <- array(c(unlist(cov.pairlist1),unlist(cov.pairlist2)),dim=c(2,2,length(bm)))

  # compute the distances for the given number of PCs
  cat('Computing distances...\n')
  dist.mat <- bhattacharyya.matrix(t(bm),sigmaarray,ipairs = 'all')

  if(csv){
    filename <- readline('Please enter a filename for the distance matrix, ending in .csv:\n\n')
    write.csv(dist.mat,filename)
    cat(filename,'saved in',getwd(),'\n\n')
  }

  if(PDF){
    pdfname <- readline('Please enter a filename for the plot file, ending in .pdf:\n\n')
    pdf(pdfname)
    if(heatmap){
      heatmap.2(t(dist.mat),dendrogram = 'none',Rowv = F,Colv = F,trace='none', main = paste(t,'\nBhattacharyya Distances PCA'))
    }
    if(PCA.plots){
      b <- info$Batch[noref.idx]
      pairs(Xscores[,1:PCA.c], col=b,pch='.', main = paste(t, "\nPCA distances"))
    }
    suppressMessages(dev.off())
    cat(pdfname,'saved in',getwd(),'\n\n')
  }else{
    if(heatmap){
      heatmap.2(t(dist.mat),dendrogram = 'none',Rowv = F,Colv = F,trace='none',main = paste(t,'\nBhattacharyya Distances PCA'))
    }
    if(PCA.plots){
      s <- info$Batch[-noref.idx]
      b <- info$Batch[noref.idx]
      pairs(Xscores[,1:PCA.c], col=b,pch='.', main = paste(t, "\nPCA combinations"))
    }
  }
  cat('Computation complete!\n\n')
  return(dist.mat)
}
