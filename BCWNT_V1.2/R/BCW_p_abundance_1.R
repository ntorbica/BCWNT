## Function to get p-value abundances from given data

BCW.p.abundance <- function(d, PLOT = TRUE, PDF = FALSE, csv = FALSE, t = ''){
  # d:          - filename of data to evaluate
  # PLOT:       - logical, defines if histograms of the computed p-values should be made
  # PDF:        - logical, defines if the plots should be saved as .pdf
  # t:          - title for plots

  ## I'll just be applying what I got from Chris, to keep everything traceable to the initial script

  cat('Preparing data...\n\n')
  type <- class(d)

  if(type == 'character'){
    d <- read.csv(d,row.names = 1)
  }
  if(type == 'list'){
    d <- as.data.frame(d)
  }

  # Prepare data

  info <- d[,c('SeqNr','SCode','Batch')]
  peaks <- d[,-(which(colnames(d)%in%colnames(info)))]

  n.c <- ncol(d)
  n.p <- ncol(peaks)

  # intialize objects for peaks
  sitep<-rep(0,n.p)
  timep<-rep(0,n.p)
  patientp<-rep(0,n.p)

  if(class(d$SeqNr) != 'integer'){
    d$SeqNr <- as.integer(d$SeqNr)
  }
  # get p-values
  cat('Computing p-values for patients, injection sequence and batch. This may take a while...\n')
  for (i in 4:n.c){
    crudelm <- aov(d[,i]~factor(SCode)+Batch+SeqNr,data=d)
    viasumunl<-unlist(summary(crudelm))
    patientp[i-3]<-viasumunl[["Pr(>F)1"]]
    sitep[i-3]<-viasumunl[["Pr(>F)2"]]
    timep[i-3]<-viasumunl[["Pr(>F)3"]]
  }

  p.abundances <- matrix()
  p.abundances <- cbind(patientp,sitep,timep)

  allpadjcb <- as.data.frame(apply(p.abundances, 2, p.adjust , method = "fdr", n = nrow(p.abundances)))


  p.ALL <- cbind(p.abundances,'sep'=rep(0,nrow(p.abundances)),allpadjcb)

  if(PLOT){
    if(PDF){
      pdfname <- readline('Please enter a name for the plotfile, ending in .pdf:\n\n')

      pdf(pdfname)
      par(mfrow=c(1,3))
      hist(patientp, main = paste(t, "Patient"))
      hist(sitep, main = paste(t, "Site"))
      hist(timep, main = paste(t, "Order"))

      par(mfrow=c(1,3))
      hist(allpadjcb$patientp, main = paste(t, "Patient adj"))
      hist(allpadjcb$sitep, main = paste(t, "Site adj"))
      hist(allpadjcb$timep, main = paste(t, "Order adj"))
      suppressMessages(dev.off())
    }else{
      par(mfrow=c(1,3))
      hist(patientp, main=paste(t,"Patient"))
      hist(sitep, main=paste(t, "Site"))
      hist(timep, main=paste(t, "Order"))

      par(mfrow=c(1,3))
      hist(allpadjcb$patientp, main=paste(t, "Patient adj"))
      hist(allpadjcb$sitep, main=paste(t, "Site adj"))
      hist(allpadjcb$timep, main=paste(t, "Order adj"))
    }
  }

  if(csv){
    csvname <- readline('Please enter a name for the table file, ending in .csv:\n\n')
    write.csv(p.ALL,csvname)
  }
  cat('Computations complete!\n')
  return(p.ALL)
}

