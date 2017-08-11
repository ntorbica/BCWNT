###################################################################################
# Master function to coordinate batch correction by R.Wehrens.                    #
# The main components are:                                                        #
#                                                                                 #
#       - BCW.r2w:              Function to reformat RAW output tables            #
#                               in such a way that they are in wide format        #
#                               and only contain the needed information           #
#                                                                                 #
#       - BCW.format:           Function to draw only the necessary info          #
#                               and the peak values into separate tables          #
#                               for Wehrens doBC function                         #
#                                                                                 #
#       - BCW.correct:          The batch correction step                         #
#                                                                                 #
#       - BCW.dist.mat:         Function to compute pairwise principal            #
#                               component bhattacharyya distances                 #
#                                                                                 #
#       - BCW.plot.loading      Function to generate loading plots of the         #
#                               PC-transformed data                               #
#                                                                                 #
#       - BCW.p.abundance       Function to compute the pvalues for each          #
#                               feature with respect to patient id, injection     #
#                               squence and batch labels                          #
#                                                                                 #
# This function is basically a wrapper and is intended for simplified use         #
# of the the listed functions.                                                    #
#                                                                                 #
###################################################################################

BCW <- function(RAW = TRUE,
                ts.r,
                i.cols,
                ts.f,
                QC = TRUE,
                QC.rem = FALSE,
                QC.n, G, B, S,
                condition = c('', '0', '1', '2', 'c'),
                method = c('lm', 'rlm'),
                k,
                maxZ = 0.1,
                PDF = TRUE,
                tc = TRUE,
                loading = TRUE){

  suppressMessages(require(BatchCorrMetabolomics))
  suppressMessages(require(RUVSeq))
  suppressMessages(require(ChemometricsWithR))
  suppressMessages(require(fpc))
  suppressMessages(require(gplots))
  suppressMessages(require(ggplot2))
  suppressMessages(require(reshape))
  suppressMessages(require(reshape2))


  if(RAW){
    cat("Welcome to the batch correction pipeline using the correction procedure by R.Wehrens!\n")
    cat("\nPlease make sure that the raw data files (RAW output) are in your working directory and that they are labelled with \n'_raw.csv' and '_ClinData.csv' at the end of the filename.\n")
    readline("Press enter to continue")

    if(missing(ts.r)){
      ts.r <- readline("Please input a table separator for the raw table:\n\n")
    }

    if(missing(i.cols)){
      i.cols <- c()
      i.cols <- c(i.cols, readline("Please check the input table and enter the name of column containing the sample description / names (can be the sample descrpition): "))
      i.cols <- c(i.cols, readline("Please check the input table and enter the name of column containing the sample grouping: "))
      i.cols <- c(i.cols, readline("Please check the input table and enter the name of column containing the batch labelling: "))
    }

    if(missing(QC)){
      QC <- readline("Please check the input table and enter the label for quality controls: ")
    }
    if(missing(S)){
      S <- readline("Please check the input table and enter the column name for the injection sequence.\nIf none is available, enter 'SeqNr': ")
    }
    if(missing(G)){
      G <- i.cols[2]
    }
    if(missing(B)){
      B <- i.cols[3]
    }
    if(missing(condition)){
      condition <- readline("Please enter a condition for correction ('', '0', '1', '2' or 'c', without the quote signs): ")
    }

    if(condition != 'c'){
      if(missing(method)){
        method <- readline("Please enter a method for correction ('lm' or 'rlm', without the quote signs)\n\n")
      }
    }
    # Only check for the raw files
    raw.lab <- '_raw.csv$'
    clin.lab <- '_ClinData.csv$'


    raw.data <-  list('raw' = grep(raw.lab, dir(),value = T), 'clin.data' = grep(clin.lab,dir(),value = T))

    if(length(raw.data)==0){
      stop('No compatible files found. Please make sure they have the required labels at the end (_raw.csv & _ClinData.csv)')
    }


    cat('STEP 1: Converting raw data (BCW.R2W)\n')
    r2w <- BCW.r2w(raw.data$raw, raw.data$clin.data, i.cols = i.cols ,ts = ts.r)


    cat('STEP 2: Formatting for correction (BCW.format)\n')
    dat.n <- grep('^R2W_', dir(), value = T)
    if(QC.rem){
      QClabel(tab = dat.n, QC.n = QC.n,QC.l = 'QC', B = B, G = G, csv = T, ts = ',', rem = T)
      dat.n <- grep('_noQC.csv$', dir(), value = T)
      dat <- BCW.format(dat.n, QC = QC.n ,G = G, B = B, S = S, mkdir = TRUE)
    }else{
      dat <- BCW.format(dat.n, QC = QC.n ,G = G, B = B, S = S, mkdir = TRUE, ts = ',')
    }

    cat('STEP 3: Removing columns with high zero content (BCW.remZ))\n')
    peaks.n0 <- BCW.remZ(dat[[1]],dat[[2]], maxZ)

    p.dim <- dim(peaks.n0[[1]])
    if(p.dim[[2]] > 1){
      p.dim <- dim(peaks.n0[[1]])
      peaks <- matrix(unlist(peaks.n0[[1]]),p.dim[[1]],p.dim[[2]])
      colnames(peaks) <- colnames(peaks.n0[[1]])
      rownames(peaks) <- peaks[,1]
      peaks <- peaks[,-1]
      peaks <- matrix(as.numeric(peaks), p.dim[[1]], p.dim[[2]]-1)

    } else {
      p.dim <- dim(dat[[1]][[1]])
      peaks <- matrix(unlist(dat[[1]][[1]]),p.dim[[1]], p.dim[[2]])
      colnames(peaks) <- colnames(dat[[1]][[1]])
      rownames(peaks) <- peaks[,1]
      peaks <- peaks[,-1]
      peaks <- matrix(as.numeric(peaks), p.dim[[1]], p.dim[[2]]-1)
    }
    info <- dat[[2]]
    rownames(info) <- rownames(peaks)

    cat('Writing file containing correction input...\n\n')
    write.csv(peaks, './Wehrens Input/Peaks_no0.csv')

   if(PDF){
      pdf(sprintf('BatchCorrection_Wehrens_%s.pdf', gsub(' ', '', paste(method,condition))), paper = 'a4')
      if(tc){
        textplot(paste('Content:       page 4..........Repeatabilities\n',
                       ifelse(QC, '              page 5..........PCA plots\n','              page 4..........PCA plots\n'),
                       ifelse(QC, '              page 6..........p-values uncorr.\n', '              page 5..........p-values uncorr.\n'),
                       ifelse(QC, '              page 7..........FDR adj. p-values uncorr.\n', '              page 6..........FDR adj. p-values uncorr.\n'),
                       ifelse(QC, '              page 8..........p-values ANCOVA\n', '              page 7..........p-values ANCOVA\n'),
                       ifelse(QC, '              page 9..........FDR adj. p-values ANCOVA\n','              page 8..........FDR adj. p-values ANCOVA\n'),
                       ifelse(QC, '              page 10.........p-values RUV\n', gsub(' ','',paste(''))),
                       ifelse(QC, '              page 11.........FDR adj. p-values RUV\n',gsub(' ','',paste(''))),
                       ifelse(QC, '              page 12.........p-values ANCOVA+RUV\n', gsub(' ','',paste(''))),
                       ifelse(QC, '              page 13.........FDR adj. p-values ANCOVA+RUV\n', gsub(' ','',paste(''))),
                       ifelse(QC, '              page 14.........Bhattacharyya Distance Matrix uncorr.\n','          page 9..........Bhattacharyya Distance Matrix uncorr.\n'),
                       ifelse(QC, '              page 15.........PCA Combination plots uncorr.\n','              page 10.........PCA Combination plots uncorr.\n'),
                       ifelse(QC, '              page 16.........Bhattacharyya Distance Matrix ANCOVA\n','              page 11.........Bhattacharyya Distance Matrix ANCOVA\n'),
                       ifelse(QC, '              page 17.........PCA Combination plots ANCOVA\n','              page 12.........PCA Combination plots ANCOVA\n'),
                       ifelse(QC, '              page 18.........Bhattacharyya Distance Matrix RUV\n', gsub(' ','',paste(''))),
                       ifelse(QC, '              page 19.........PCA Combination plots RUV\n', gsub(' ','',paste(''))),
                       ifelse(QC, '              page 20.........Bhattacharyya Distance Matrix ANCOVA+RUV\n',gsub(' ','',paste(''))),
                       ifelse(QC, '              page 21.........PCA Combination plots ANCOVA+RUV\n', gsub(' ','',paste(''))),
                       ifelse(loading, ifelse(QC, '              pages 22-26.....Loading plots ANCOVA\n', '          pages 13-17.....Loading plots ANCOVA\n'),gsub(' ','',paste(''))),
                       ifelse(loading, ifelse(QC, '              pages 27-31.....Loading plots RUV\n', '\n'), gsub(' ','',paste(''))),
                       ifelse(loading, ifelse(QC, '              pages 32-36.....Loading plots ANCOVA+RUV\n', gsub(' ','',paste(''))),gsub(' ','',paste('')))),
                 halign = 'left', valign = 'top')
      }
    }


    cat('STEP 4: Correcting for batch effect (BCW.correct)\n')

    bc <- BCW.correct(peaks = peaks,
                      info = info,
                      condition = condition,
                      method = method,
                      k = k,
                      QC = QC,
                      PDF = FALSE,
                      csv = TRUE,
                      PLOT = TRUE,
                      nA = 'ANCOVA.csv',
                      nR = 'RUV.csv',
                      nAR = 'ANCOVA_RUV.csv')

    p.bcA.dim <- dim(bc$BC_Data_ANCOVA)

    if(QC){
      p.bcR.dim <- dim(bc$BC_Data_RUV)
      p.bcAR.dim <- dim(bc$BC_Data_AR)
    }

    peaks.bcA <- matrix(as.numeric(unlist(bc$BC_Data_ANCOVA)),p.bcA.dim[[1]],p.bcA.dim[[2]])
    if(QC){
      peaks.bcR <- matrix(as.numeric(unlist(bc$BC_Data_RUV)),p.bcR.dim[[1]],p.bcR.dim[[2]])
      peaks.bcAR <- matrix(as.numeric(unlist(bc$BC_Data_AR)),p.bcAR.dim[[1]],p.bcAR.dim[[2]])
    }


    cat('STEP 5: Computing p-value abundances before and after correction (BCW.p.abuncance)\n')
    pre <- cbind(info,peaks)
    pval.pre <- BCW.p.abundance(pre, PLOT = T)

    post.A <- cbind(info,peaks.bcA)
    if(QC){
      post.R <- cbind(info,peaks.bcR)
      post.AR <- cbind(info,peaks.bcAR)
    }

    pval.post.A <- BCW.p.abundance(post.A, PLOT = T, t = 'ANCOVA')
    if(QC){
      pval.post.R <- BCW.p.abundance(post.R, PLOT = T, t = 'RUV')
      pval.post.AR <- BCW.p.abundance(post.AR, PLOT = T, t = 'ANCOVA + RUV')
    }

    cat('STEP 6: Computing bhattacharyya distances before and after correction (BCW.plot.loading)\n')
    # Some graphics tweaking
    old.par <- par(mar = c(0, 0, 0, 0))

    dist.pre <- BCW.dist.mat(pre, PCA.plots = T, PCA.c = 5, t = 'raw')
    par(old.par)
    dist.post.A <- BCW.dist.mat(post.A, PCA.plots = T, PCA.c = 5, t = 'ANCOVA')
    par(old.par)
    if(QC){
      dist.post.R <- BCW.dist.mat(post.R, PCA.plots = T, PCA.c = 5, t = 'RUV')
      par(old.par)
      dist.post.AR <- BCW.dist.mat(post.AR, PCA.plots = T, PCA.c = 5, t = 'ANCOVA + RUV')
      par(old.par)
    }

    if(loading){
      cat('STEP 7: Produce loading plots for identification of features of high contribution to variability\n')
      BCW.plot.loading(post.A, t = 'ANCOVA')
      if(QC){
        BCW.plot.loading(post.R, t = 'RUV')
        BCW.plot.loading(post.AR, t = 'ANCOVA + RUV')
      }
    }
    cat('End of the (pipe) line!\nFind your corrected data in the directory with the input files.\n\n')
    dev.off()



  } else {
    cat("Welcome to the batch correction pipeline using the correction procedure by R.Wehrens!\n\nPlease make sure that the raw data files (already preformatted from RAW format) are in your working directory.\nNote that you should rename the file such that it contains 'R2W_' at the beginning of the filename to be recognizable.\n")
    readline("Press enter to continue")

    if(missing(ts.f)){
      ts.f <- readline("Please enter a table separator for the preformatted table: ")
    }

    if(missing(QC)){
      QC <- readline("Please check the input table and enter the label for quality controls: ")
    }
    if(missing(S)){
      S <- readline("Please check the input table and enter the column name for the injection sequence.\nIf none is available, enter 'SeqNr'. ")
    }
    if(missing(G)){
      G <- readline("Please check the input table and enter the column name for the sample grouping, should contain the quality control labels: ")
    }
    if(missing(B)){
      B <- readline("Please check the input table and enter the column name for the batch labels: ")
    }
    if(missing(condition)){
      condition <- readline("Please enter a condition for correction ('', '0', '1', '2' or 'c', without the quote signs): ")
    }
    if(condition != 'c'){
      if(missing(method)){
        method <- readline("Please enter a method for correction ('lm' or 'rlm', without the quote signs)\n\n")
      }
    }

    cat('STEP 1: Formatting for correction (BCW.format)\n')
    dat.n <- grep('^R2W_',dir(),value = T)
    if(QC.rem){
      QClabel(tab = dat.n, QC.n = QC.n ,QC.l = 'QC', B = B, G = G, csv = T, ts = ts.f, rem = T)
      dat.n <- grep('_noQC.csv$', dir(), value = T)
      dat <- BCW.format(dat.n, QC = QC.n ,G = G, B = B, S = S, mkdir = TRUE)
    } else {
    dat <- BCW.format(dat.n,QC = QC.n, G = G, B = B, S = S, mkdir = TRUE, ts = ts.f)
    }

    cat('STEP 2: Removing columns with high zero content (BCW.remZ)\n')
    peaks.n0 <- BCW.remZ(dat[[1]],dat[[2]], maxZ = maxZ)

    p.dim <- dim(peaks.n0[[1]])
    if(p.dim[[2]] > 1){
      p.dim <- dim(peaks.n0[[1]])
      peaks <- matrix(unlist(peaks.n0[[1]]),p.dim[[1]],p.dim[[2]])
      colnames(peaks) <- colnames(peaks.n0[[1]])
      rownames(peaks) <- peaks[,1]
      peaks <- peaks[,-1]
      peaks <- matrix(as.numeric(peaks), p.dim[[1]], p.dim[[2]]-1)

    } else {
      p.dim <- dim(dat[[1]][[1]])
      peaks <- matrix(unlist(dat[[1]][[1]]),p.dim[[1]], p.dim[[2]])
      colnames(peaks) <- colnames(dat[[1]][[1]])
      rownames(peaks) <- peaks[,1]
      peaks <- peaks[,-1]
      peaks <- matrix(as.numeric(peaks), p.dim[[1]], p.dim[[2]]-1)
    }
    info <- dat[[2]]
    rownames(info) <- rownames(peaks)


    cat('Writing file containing correction input...\n\n')
    write.csv(peaks, './Wehrens Input/Peaks_no0.csv')

    if(PDF){
      pdf('BatchCorrection_Wehrens.pdf', paper = 'a4')
      if(tc){
        textplot(paste('Content:       page 4..........Repeatabilities\n',
                       ifelse(QC, '              page 5..........PCA plots\n','              page 4..........PCA plots\n'),
                       ifelse(QC, '              page 6..........p-values uncorr.\n', '              page 5..........p-values uncorr.\n'),
                       ifelse(QC, '              page 7..........FDR adj. p-values uncorr.\n', '              page 6..........FDR adj. p-values uncorr.\n'),
                       ifelse(QC, '              page 8..........p-values ANCOVA\n', '              page 7..........p-values ANCOVA\n'),
                       ifelse(QC, '              page 9..........FDR adj. p-values ANCOVA\n','              page 8..........FDR adj. p-values ANCOVA\n'),
                       ifelse(QC, '              page 10.........p-values RUV\n', gsub(' ','',paste(''))),
                       ifelse(QC, '              page 11.........FDR adj. p-values RUV\n',gsub(' ','',paste(''))),
                       ifelse(QC, '              page 12.........p-values ANCOVA+RUV\n', gsub(' ','',paste(''))),
                       ifelse(QC, '              page 13.........FDR adj. p-values ANCOVA+RUV\n', gsub(' ','',paste(''))),
                       ifelse(QC, '              page 14.........Bhattacharyya Distance Matrix uncorr.\n','          page 9..........Bhattacharyya Distance Matrix uncorr.\n'),
                       ifelse(QC, '              page 15.........PCA Combination plots uncorr.\n','              page 10.........PCA Combination plots uncorr.\n'),
                       ifelse(QC, '              page 16.........Bhattacharyya Distance Matrix ANCOVA\n','              page 11.........Bhattacharyya Distance Matrix ANCOVA\n'),
                       ifelse(QC, '              page 17.........PCA Combination plots ANCOVA\n','              page 12.........PCA Combination plots ANCOVA\n'),
                       ifelse(QC, '              page 18.........Bhattacharyya Distance Matrix RUV\n', gsub(' ','',paste(''))),
                       ifelse(QC, '              page 19.........PCA Combination plots RUV\n', gsub(' ','',paste(''))),
                       ifelse(QC, '              page 20.........Bhattacharyya Distance Matrix ANCOVA+RUV\n',gsub(' ','',paste(''))),
                       ifelse(QC, '              page 21.........PCA Combination plots ANCOVA+RUV\n', gsub(' ','',paste(''))),
                       ifelse(loading, ifelse(QC, '              pages 22-26.....Loading plots ANCOVA\n', '          pages 13-17.....Loading plots ANCOVA\n'),gsub(' ','',paste(''))),
                       ifelse(loading, ifelse(QC, '              pages 27-31.....Loading plots RUV\n', '\n'), gsub(' ','',paste(''))),
                       ifelse(loading, ifelse(QC, '              pages 32-36.....Loading plots ANCOVA+RUV\n', gsub(' ','',paste(''))),gsub(' ','',paste('')))),
                 halign = 'left', valign = 'top')
      }

    }
    cat('STEP 3: Correcting for batch effect (BCW.correct)\n')

    bc <- BCW.correct(peaks = peaks, info = info, condition = condition, method = method, k = k, QC = QC, PDF = FALSE, csv = TRUE, PLOT = TRUE, nA = 'ANCOVA.csv', nR = 'RUV.csv', nAR = 'ANCOVA_RUV.csv')

    p.bcA.dim <- dim(bc$BC_Data_ANCOVA)
    if(QC){
      p.bcR.dim <- dim(bc$BC_Data_RUV)
      p.bcAR.dim <- dim(bc$BC_Data_AR)
    }

    peaks.bcA <- matrix(as.numeric(unlist(bc$BC_Data_ANCOVA)),p.bcA.dim[[1]],p.bcA.dim[[2]])
    if(QC){
      peaks.bcR <- matrix(as.numeric(unlist(bc$BC_Data_RUV)),p.bcR.dim[[1]],p.bcR.dim[[2]])
      peaks.bcAR <- matrix(as.numeric(unlist(bc$BC_Data_AR)),p.bcAR.dim[[1]],p.bcAR.dim[[2]])
    }

    cat('STEP 4: Computing p-value abundances before and after correction (BCW.p.abundance)\n')
    pre <- cbind(info,peaks)
    pval.pre <- BCW.p.abundance(pre, PLOT = T)

    post.A <- cbind(info,peaks.bcA)
    if(QC){
      post.R <- cbind(info,peaks.bcR)
      post.AR <- cbind(info,peaks.bcAR)
    }

    pval.post.A <- BCW.p.abundance(post.A, PLOT = T, t = 'ANCOVA')
    if(QC){
      pval.post.R <- BCW.p.abundance(post.R, PLOT = T, t = 'RUV')
      pval.post.AR <- BCW.p.abundance(post.AR, PLOT = T, t = 'ANCOVA + RUV')
    }

    cat('STEP 5: Computing bhattacharyya distances before and after correction (BCW.dist.mat)\n')
    # Some graphics tweaking
    old.par <- par(mar = c(0, 0, 0, 0))

    dist.pre <- BCW.dist.mat(pre, PCA.plots = T, PCA.c = 5)
    par(old.par)

    dist.post.A <- BCW.dist.mat(post.A, PCA.plots = T, PCA.c = 5, t = 'ANCOVA')
    par(old.par)
    if(QC){
      dist.post.R <- BCW.dist.mat(post.R, PCA.plots = T, PCA.c = 5, t = 'RUV')
      par(old.par)
      dist.post.AR <- BCW.dist.mat(post.AR, PCA.plots = T, PCA.c = 5, t = 'ANCOVA + RUV')
      par(old.par)
    }

    if(loading){
      cat('STEP 6: Produce loading plots for identification of features of high contribution to variability (BCW.plot.loading)\n')
      BCW.plot.loading(post.A, t = 'ANCOVA')
      if(QC){
        BCW.plot.loading(post.R, t = 'RUV')
        BCW.plot.loading(post.AR, t = 'ANCOVA + RUV')
      }
    }
    cat('\nEnd of the (pipe) line!\nFind your corrected data in the directory with the input files.\n\n')
    dev.off()
  }
}
