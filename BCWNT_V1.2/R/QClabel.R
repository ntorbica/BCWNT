QClabel <- function(tab, QC.n, QC.l, G, B, csv = F, ts, rem = F){
  # tab           Data table
  # QC.n          QC label
  # QC.l          new QC batch label
  # B             column name of batch label column
  # csv           output a table file
  if(!is.data.frame(tab)){
    if(is.character(tab)){
      n <- gsub('.csv','',tab)
      tab <- read.csv(tab, row.names = 1, sep = ts)
    }
    if(is.list(tab)){
      tab <- as.data.frame(tab)
    }
  }

  if(is.data.frame(tab)){
    qc.id <- which(tab[,G] == QC.n)

    if(rem){
      tab <- tab[-qc.id,]
    } else {
    tab[,B] <- as.character(tab[,B])
    tab[qc.id,B] <- as.character(QC.l)
    tab[,B] <-  as.factor(tab[,B])
    }

    if(csv){
      write.csv(x = tab, file = ifelse(is.character(n), gsub(' ','', paste(n,'_noQC.csv')), 'Info_noQCs.csv'))
    }
  }

  return(tab)
}
