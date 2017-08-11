
BCW.plot.loading <- function(dat, n.PC = 5, PDF = FALSE, t = ''){
  # dat:        - filename of data table
  # n.PC:       - number of PCs to plot
  # PDF:        - logical, should PDF file be generated
  # t           - title for plots

  suppressMessages(require(BatchCorrMetabolomics))
  suppressMessages(require(RUVSeq))
  suppressMessages(require(ChemometricsWithR))
  suppressMessages(require(ggplot2))
  suppressMessages(require(reshape))
  suppressMessages(require(reshape2))

  # Read the data
  cat('Preparing data...\n')
  type <- class(dat)

  if(type == 'character'){
    dat <- as.data.frame(read.csv(dat,row.names = 1, sep = ';'))
  }
  if(type == 'list'){
    dat <- as.data.frame(dat)
  }


  D <- dat
  info <- D[,c('SeqNr','SCode','Batch')]
  peaks <- D[,-(which(colnames(D)%in%colnames(info)))]

  nb <- nlevels(info$Batch)

  noref.idx <- which(info$SCode != "ref")
  Xsample <- peaks[noref.idx,]
  YSample <- info[noref.idx,]

  Xsample <- Xsample[, apply(Xsample, 2, function(x) !all(is.na(x)))]
  Xsample <- Xsample[, apply(Xsample, 2, sd, na.rm = TRUE) > 0]

  for (i in 1:ncol(Xsample)){
    Xsample[is.na(Xsample[,i]),i] <- mean(Xsample[,i], na.rm = TRUE)
  }

  Xsample <- scale(Xsample)
  X.PCA <- PCA(Xsample)

  Xload <- as.data.frame(abs(X.PCA$loadings))
  f.names <- rownames(Xload)

  if(PDF){
    pdfname <- readline('Please enter a name for the plot file, ending in .pdf:\n\n')
    pdf(pdfname)
    for(i in 1:n.PC){
      p <- ggplot(Xload,
                     aes((1:length(f.names)),Xload[,i])) + geom_col() + theme(axis.text.x = element_text(angle = 90,size = 10)) + labs(x = 'Features', y = 'Loading') + ggtitle(label = t, subtitle = paste('Loadings PC',i))
      print(p)
    }
    dev.off()
  }else{
    cat('Plotting...\n')
    for(i in 1:n.PC){
      p <- ggplot(Xload,
                  aes((1:length(f.names)),Xload[,i])) + geom_col() + theme(axis.text.x = element_text(angle = 90,size = 10)) + labs(x = 'Features', y = 'Loading') + ggtitle(label = t, subtitle = paste('Loadings PC',i))
      print(p)

    }
  }
  cat('Computations complete!\n\n')
}
