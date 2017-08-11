#################################################################
# This function serves to convert raw LC-MS output  to the      #
# format defined by R. Wehrens for batch correction. The        #
# initial assumption is that there are two files containing the #
# peak information and clinical information about the measured  #
# samples. These will be wrapped to several tables intended to  #
# 1. conserve the complete information and 2. to generate the   #
# input for Wehrens batch correction.                           #
#################################################################


BCW.r2w <- function(peaks, clin.dat, ts, long = TRUE, i.cols){
  ## INPUTS:
  #   peaks     - The peak measurement filename, given as string
  #   clin.dat  - The clinical information filename, given as string
  #   ts        - table separator as string
  #   long      - defines if the peak values are given in long or wide format
  #   i.cols    - vector with column names within clin.dat to choose

  cat('Preparing data...\n')
  type.p <- class(peaks)
  type.i <- class(clin.dat)

  if(type.p == 'list'){
    p.tab <- as.data.frame(peaks)
  }
  if(type.p == 'character'){
    p.tab <- read.csv(peaks,sep = ts)
  }

  if(type.i == 'list'){
    c.tab <- as.matrix(as.data.frame(clin.dat))
  }
  if(type.i == 'character'){
    c.tab <- as.matrix(read.csv(clin.dat,sep = ts))
  }

  zrow <- rowCheck(p.tab,20)
  if(!is.null(zrow)){
    p.tab <- p.tab[-zrow,]
  }

  # if peaks in long format (feature names in rows) transpose to wide
  if(long){
    p.tab <- t(p.tab)
  }
  # find the 'border' between peaks and information within table
  border <- which(rownames(p.tab) == 'Raw.abundance')

  # alternative search in case naming varies
  if(length(border) == 0){
    pat <- '^[0-9]{8}_'
    match <- grep(pat, p.tab[,1])
    border <- match[1]
  }

  # Save the information part within the peak table in a separate object
  p.i.tab <- p.tab[1:border-1,]
  p.tab <- p.tab[-(1:border-1),]

  cat('Ordering data. May take a while...\n')

  p.n <- p.tab[,1]
  c.n <- c.tab[,1]

  m <- max(length(p.n), length(c.n))
  length(p.n) <- m
  length(c.n) <- m

  t <- cbind(c.n, p.n)
  t.id <- which(t[,1] %in% t[,2])

  if(length(t.id) == 0){
    c.n <- c.tab[,2]
  }

  eq.id1 <- which(c.n %in% p.n)
  c.tab.eq <- c.tab[as.numeric(eq.id1),]

  eq.id2 <- which(p.n %in% c.n)
  p.tab.eq <- p.tab[as.numeric(eq.id2),]
  check.tab <- as.data.frame(p.tab.eq)   # Backup table for mismatch checking

  # Get information table with specified info
  ip.tab <- c.tab[eq.id1,]
  ip.tab <- ip.tab[,i.cols]

  # Get order information from sample names and add them to table
  if(length(grep('[0-9]{8}', ip.tab[,1])) > 0){
    SeqNr <- getOrd(ip.tab[,1])
  } else {
    SeqNr <- as.numeric(ip.tab[,1])
  }
  ip.tab <- as.data.frame(ip.tab)
  ip.tab$SeqNr <- SeqNr

  ## Different approach on sorting the peaks. eqN seems to mess it up and twist
  ## intensities such that it looks like wrong labelling.
  p.seq <- getOrd(p.tab.eq[,1])
  p.tab.eq <- as.data.frame(p.tab.eq)
  p.tab.eq$SeqNr <- p.seq

  p.tab.eq <- align2(ip.tab, p.tab.eq, check.tab, ini = 1)

  cat('Finalizing...\n\n')
  # Slightly modify the feature names for easier handling with regular
  # expressions (add a 'X' at the beginning of each feature)
  cnames <- list()
  p.i.dim <- dim(p.i.tab)
  if(!is.null(p.i.dim)){
    for(i in 2:length(p.i.tab[1,1:length(p.i.tab[1,])])){
      n <- p.i.tab[1,i]
      cnames[i] <- gsub(' ','',paste('X',n))
    }
  } else {
    for(i in 2:length(p.i.tab)){
      n <- p.i.tab[i]
      cnames[i] <- gsub(' ','',paste('X',n))
    }
  }


  colnames(p.tab.eq) <- cnames
  p.tab.eq <- p.tab.eq[,-1]

  # Finally, generate output
  complete <- cbind(ip.tab, p.tab.eq)
  n <- grep('name', colnames(complete), ignore.case = T)
  rownames(complete) <- complete[,n]
  complete <- complete[,-n]


  filename <- sprintf('R2W_%s', peaks)
  cat('Writing file with name:\n', filename, '\n')
  write.csv(complete, filename)
  cat('Done!\n\n')
  return(complete)
}

