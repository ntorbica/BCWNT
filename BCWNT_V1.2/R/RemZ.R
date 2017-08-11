BCW.remZ <- function(tab,Info,maxZ){
  ################################################################################
  # Zero-column filter. This function takes a dataset (only the actual values,   #
  # remove descriptive information for now) and a ratio, consisting of the zero  #
  # count per column divided by the total row count, which is the decisive       #
  # boundary for the removal of a column. The information table as defined for   #
  # Wehrens approach is needed to account for QCs.                               #
  ################################################################################

  if(missing(maxZ)){
    readline('Please enter the maximum zero content within a single column to be allowed (0-1): ')
  }
  cat("Removing columns with more than", maxZ*100, "% zero content.\n")
  # Backup the data table to re-insert the QCs later on
  BU <- tab

  if(class(tab)=="list"){
    tab <- as.data.frame(tab)
  }
  if(class(tab)=="data.frame"){
    tab <- as.matrix(tab)
  }

  if(class(Info)=="list"){
    Info <- as.data.frame(Info)
  }
  if(class(Info)=="data.frame"){
    Info <- as.matrix(Info)
  }

  fullCol <- nrow(Info)

  # Check if first column is sample names
  tab.t <- matrix(as.numeric(tab), dim(tab)[1], dim(tab)[2])
  if(is.na(tab.t[,1])){
    tab <- tab[,-1]
  }
  # Count the zeroes in each column and save them in a list per column index
  Zcount <- list()
  for(i in 1:ncol(tab)){
    x <- 0
    for(j in 1:fullCol){
      if(as.numeric(tab[,i][j]) == 0){
        x <- x+1
      }
    }
    Zcount[[i]] <- x
  }
  l.z <- length(Zcount)

  # Create a binary vector with the decision concerning zero content
  # (1 = too many zeroes, 0 = keep the column)
  id.dec <- c()
  x <- 1
  C <- 0
  while(x < l.z+1){
    ratio <- Zcount[[x]]/fullCol
    d <- ifelse(ratio >= maxZ, 1, 0)
    if(d != 0){
      C <- C+1
    }
    id.dec[x] <- d
    x <- x+1
  }
  id.rem <- which(id.dec == 1)
  if(C != 0){
    cat("Removed", C, "columns exceeding", maxZ*100, "% zero content.\n")
    #cat("Following columns (index only) are to be removed:\n", sprintf('%s',id.rem))

  }else{
    cat("No columns found that exceed", maxZ*100, "% zero content.\n")
  }

  # Finally remove the decided columns and return a list containing
  # the filtered table, as well as an index-list of all removed columns
  tab <- list(tab[,-id.rem], id.rem)
}
