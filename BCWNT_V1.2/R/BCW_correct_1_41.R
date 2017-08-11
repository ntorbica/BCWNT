###############################################################################
# This the 'functionalized' script for batch correction using R. Wehrens      #
# approach. The prerequisites are prepared by BCW.r2w and BCW.format          #
###############################################################################

BCW.correct <- function(peaks, info,
                minQC = 4,
                condition = c('','0','1','2','c'),
                method = c('lm','rlm'),
                k,
                QC = TRUE,
                PCA = TRUE,
                Duplo = TRUE,
                PLOT = TRUE,
                PDF = FALSE,
                csv = FALSE,
                nA = 'ANCOVA.csv',
                nR = 'RUV.csv',
                nAR = 'ANCOVA_RUV.csv',
                ...)
                {
  ## INPUTS:
  # peaks                   - WTabs output filename, table containing only peaks
  # info                    - WTabs output filename, table containing information about peaks (injection order, sample group, batch)
  # minBatchOccurrence.Line - minimum quality control sample count per batch
  # condition               - string defining the approach for imputing NA values in the data, as well as the method used
  #                           to correct the data (either linear regression for "" to "2" or censored regression (tobit) for "c")
  # method                  - The model to apply with the imputed value
  # k                       - Number of PCs to chose for RUVs
  # QC                      - The presence of quality controls. Actually a little obsolete, as the preparation steps give
  #                           insight about this fact. I'll leave it in for the beginning.
  # PCA                     - logical, if TRUE, output a PCA comparison
  # Duplo                   - logical, if TRUE, output a duplicate analysis
  # PLOT                    - Plot the PCA and Duplo results for comparison
  # PDF                     - Depends on PLOT, generates PDF with plots when TRUE
  # csv                     - Output a .csv file with the corrected data.
  # nA, nR, nAR             - File names for outputting the corrected peaks


  #Test lines
  # peaks <- read.csv('DATA/2017_28_03_SET3/Wehrens/Wehrens Input/Peaks.csv')
  # info <- read.csv('DATA/2017_28_03_SET3/Wehrens/Wehrens Input/Info.csv')
  # minBatchOccurrence.Ave = 2
  # minBatchOccurrence.Line = 4


  ## Error catching
  if(missing(peaks)){
    stop("Missing Data! No peak table given.")
  }
  if(missing(info)){
    stop("Missing Data! No information on peaks given.")
  }

  # c <- c('','0','1','2','c')
  # condition <- readline("Please choose a condition for imputed Values ('','0','1','2' or 'c').\n")
  # if(condition %in% c == FALSE){
  #   stop('Invalid input, rerun and input one of the available conditions')
  # }

  if(condition == ''){
    cat("You chose to ignore NA measurements.\n")
    impV <- NA
  }
  if(condition == '0'){
    cat("You chose to impute NA measurements with 0.\n")
  }
  if(condition == '1'){
    cat("You chose to impute NA measurements with half the limit of detection.\n")
  }
  if(condition == '2'){
    cat("You chose to impute NA measurements with the limit of detection (lowest measurement).\n")
  }
  # if(condition != 'c'){
  #   method <- readline("Please enter the method of correction (lm or rlm).\n")
  # }
  if(condition == 'c'){
    cat("You chose to impute Values for a censored regression. Note that this will require 'tobit' as method.\n")
    method <- 'tobit'
  }
  cat('Preparing data...\n\n')
  suppressMessages(require(BatchCorrMetabolomics))
  suppressMessages(require(RUVSeq))


  ## This step serves as class conversion, as the function works best with a matrix for peaks and data.frame
  ## for the information
  type.p <- class(peaks)

  type.i <- class(info)

  if(type.p != 'matrix'){
    if(type.p == 'character'){
      peaks <- as.matrix(read.csv(peaks,row.names = 1))
      type.p <- class(peaks)
    }
    if(type.p == 'list'){
      peaks <- as.matrix(as.data.frame(peaks))
      type.p <- class(peaks)
    }
    if(type.p == 'data.frame'){
      peaks <- as.matrix(peaks)
      type.p <- class(peaks)
    }
  }

  if(type.i != 'data.frame'){
    if(type.i == 'character'){
      info <- as.data.frame(read.csv(info,row.names = 1))
      type.i <- class(info)
    }
    if(type.i == 'list'){
    info <- as.data.frame(info)
    type.i <- class(info)
    }
  }


  if(type.p == 'matrix' & type.i == 'data.frame'){
    set.1 <- peaks
    set.1.Y <- info


    set.1.lod <- as.numeric(min(set.1[!is.na(set.1)]))
    if(condition == '0'){
      impV <- 0
    }
    if(condition == '1'){
      impV <- set.1.lod/2
    }
    if(condition == '2'){
      impV <- set.1.lod
    }
    if(condition == 'c'){
      if(set.1.lod == 0){
        impV <- set.1.lod + 0.01
      }else{
        impV <- set.1.lod - 0.01
      }
    }

    LC.ngenotypes <- nlevels(set.1.Y$SCode)-1             # Get the data specific number of samples
    if(QC){
      LC.nref <- sum(set.1.Y$SCode == "ref")              # Get the number of reference samples/QC measures
    }
    LC.nNA <- apply(set.1, 2, function(x) sum(is.na(x)))  # Get the number of NAs to account for per column


    # get the QC / sample indices
    refSamples <- list("Q" = which(set.1.Y$SCode == "ref"),
                       "S" = which(set.1.Y$SCode != "ref"))


    # with the complete set of possible procedures, we compare the user input and chose the appropriate
    # one in order to avoid computing all of the corrections.
    if(QC){
      actual <- list(str = "Q", exp = condition, method = method, impV = impV)
    }else{
      actual <- list(str = "S", exp = condition, method = method, impV = impV)
    }
    cat('Corrections applied as follows:\n
        correction model:\t',actual$method,'\n',
        '\timputed values:\t\t', actual$impV, '\n',
        '\tstrategy:\t\t', actual$str, '\n\n')
    if(!QC){
      cat('Note that RUV correction is not possible without QCs\n')
    }

    # Now we are technically set for batch correction.
    # ANCOVA type correction:
    CP <- list()

    CP$ANCOVA <- apply(set.1, 2, doBC,
                           ref.idx = as.numeric(refSamples[[actual$str]]),
                           batch.idx = set.1.Y$Batch,
                           minBsamp = minQC,
                           seq.idx = as.numeric(set.1.Y$SeqNr),
                           method = actual$method,
                           imputeVal = actual$impV)

    # RUV normalization
    idx <- which(set.1.Y$SCode == "ref")
    replicates.ind <- matrix(-1, nrow(set.1) - length(idx) + 1, length(idx))
    replicates.ind[1,] <- idx
    replicates.ind[-1,1] <- (1:nrow(set.1))[-idx]
    nColumns <- ncol(set.1)

    if(QC){
      CP$RUVs <- t(RUVs(t(set.1),
                        1:nColumns,
                        k = k,
                        replicates.ind,
                        round = FALSE,
                        isLog = TRUE)$normalizedCounts)

      CP$AR <- t(RUVs(t(CP$ANCOVA),
                      1:nColumns,
                      k = 3,
                      replicates.ind,
                      round = FALSE,
                      isLog = TRUE)$normalizedCounts)
    }

    if(csv){
      if(!dir.exists('CORRECTED_PEAKS/')){
        dir.create('CORRECTED_PEAKS')
      } else {
        Sys.chmod('CORRECTED_PEAKS/')
        cat('Removing old directory...\n')
        unlink('CORRECTED_PEAKS/', recursive = TRUE)
      }

      if(missing(nA)){
        nA <- readline("Enter a filename for ANCOVA corrected data output, ending in .csv:\n")
      }
      if(QC){
        if(missing(nR)){
          nR <- readline("Enter a filename for RUV normlized data output, ending in .csv:\n")
        }
        if(missing(nAR)){
          nAR <- readline("Enter a filename for ANCOVA corrected + RUV normlized data output, ending in .csv:\n")
        }
      }

      cat('Creating new directory...\n')
      Sys.chmod('CORRECTED_PEAKS/')
      write.csv(cbind(set.1.Y,CP$ANCOVA), gsub(' ', '',paste('./CORRECTED_PEAKS/',nA)))
      if(QC){
        write.csv(cbind(set.1.Y,CP$RUVs), gsub(' ', '',paste('./CORRECTED_PEAKS/',nR)))
        write.csv(cbind(set.1.Y,CP$AR), gsub(' ', '',paste('./CORRECTED_PEAKS/',nAR)))
      }
    }
    ## Output the duplicate and PCA evaluation
    RES <- list()
    RES[['BC_Data_ANCOVA']] <- CP$ANCOVA
    if(QC){
      RES[['BC_Data_RUV']] <- CP$RUVs
      RES[['BC_Data_AR']] <- CP$AR
    }

    if(PLOT){
      cat('Plotting...\n')
      par(oma=c(1,1,1,1))
      # par(mfrow=c(1,1))
      textplot(paste('Batch Correctio Wehrens - Pipeline\n\n',
               'Nr. of Samples: ', length(refSamples$S), '\n',
               'Nr. of QCs      ', length(refSamples$Q), '\n\n',
               'Correction:\n\n ',
               'Strategy:       ', ifelse(QC, 'Q', 'S'), '\n',
               ' Method:         ', ifelse(condition == 'c', 'tobit', method), '\n',
               ' Condition:      ', condition, '\n',
               ' Imputed value:  ', actual$impV, '\n\n',
               ifelse(QC, paste('RUV:\n\n nr. of PCs:     ', k, '\n\n'),'')),
               halign = 'left', valign = 'top')

      if(PDF){
        name <- readline('Please enter a filename for the plot file, ending with .pdf\n')

        pdf(name)
        par(mfrow=c(2,2))

        if(Duplo){
          RES['Duplo_Raw'] <- evaluateCorrection(set.1, set.1.Y, what = "duplo", plot = PLOT)
          title(main = paste("Repeatabilities:", round(as.numeric(RES[['Duplo_Raw']]), 3)))

          RES['Duplo_C_ANCOVA'] <- evaluateCorrection(CP$ANCOVA, set.1.Y, what = "duplo", plot = PLOT)
          title(main = paste("Q: Repeatabilities\nANCOVA-corr.:", round(as.numeric(RES[['Duplo_C_ANCOVA']]), 3)))

          if(QC){
            RES['Duplo_C_RUV'] <- evaluateCorrection(CP$RUVs, set.1.Y, what = "duplo", plot = PLOT)
            title(main = paste("Q: Repeatabilities\nRUV-corr.:", round(as.numeric(RES[['Duplo_C_RUV']]), 3)))

            RES['Duplo_C_AR'] <- evaluateCorrection(CP$AR, set.1.Y, what = "duplo", plot = PLOT)
            title(main = paste("Q: Repeatabilities\nRUV-corr.:", round(as.numeric(RES[['Duplo_C_AR']]), 3)))
          }
        }
        if(PCA){

          RES['PCA_Raw'] <- evaluateCorrection(set.1, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright')
          title(main = paste("Interbatch distance:", round(as.numeric(RES[['PCA_Raw']]), 3)))

          RES['PCA_C_ANCOVA'] <- evaluateCorrection(CP$ANCOVA, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright')
          title(main = paste("Q: Interbatch distance\nANCOVA-corr.:", round(as.numeric(RES[['PCA_C_ANCOVA']]), 3)))

          if(QC){
            RES['PCA_C_RUV'] <- evaluateCorrection(CP$RUVs, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright')
            title(main = paste("Q: Interbatch distance\nRUV-corr.:", round(as.numeric(RES[['PCA_C_RUV']]), 3)))

            RES['PCA_C_AR'] <- evaluateCorrection(CP$AR, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright')
            title(main = paste("Q: Interbatch distance\nRUV-corr.:", round(as.numeric(RES[['PCA_C_AR']]), 3)))
          }
        }
        dev.off()
        cat(name, "saved in", getwd(),"\n")

      }else{
        if(Duplo){
          par(mfrow=c(2,2))
          # par(mfrow=c(1,1))

          RES['Duplo_Raw'] <- evaluateCorrection(set.1, set.1.Y, what = "duplo", plot = PLOT)
          title(main = paste("Repeatabilities:", round(as.numeric(RES[['Duplo_Raw']]), 3)))

          RES['Duplo_C_ANCOVA'] <- evaluateCorrection(CP$ANCOVA, set.1.Y, what = "duplo", plot = PLOT)
          title(main = paste("Q: Repeatabilities\nANCOVA-corr.:", round(as.numeric(RES[['Duplo_C_ANCOVA']]), 3)))

          if(QC){
            RES['Duplo_C_RUV'] <- evaluateCorrection(CP$RUVs, set.1.Y, what = "duplo", plot = PLOT)
            title(main = paste("Q: Repeatabilities\nRUV-corr.:", round(as.numeric(RES[['Duplo_C_RUV']]), 3)))

            RES['Duplo_C_AR'] <- evaluateCorrection(CP$AR, set.1.Y, what = "duplo", plot = PLOT)
            title(main = paste("Q: Repeatabilities\nANCOVA + RUV-corr.:", round(as.numeric(RES[['Duplo_C_AR']]), 3)))
          }
        }

        if(PCA){
          RES['PCA_Raw'] <- evaluateCorrection(set.1, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright', legend.col = 2)
          title(main = paste("Interbatch distance:", round(as.numeric(RES[['PCA_Raw']]), 3)))

          RES['PCA_C_ANCOVA'] <- evaluateCorrection(CP$ANCOVA, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright', legend.col = 2)
          title(main = paste("Q: Interbatch distance\nANCOVA-corr.:", round(as.numeric(RES[['PCA_C_ANCOVA']]), 3)))

          if(QC){
            RES['PCA_C_RUV'] <- evaluateCorrection(CP$RUVs, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright', legend.col = 2)
            title(main = paste("Q: Interbatch distance\nRUV-corr.:", round(as.numeric(RES[['PCA_C_RUV']]), 3)))

            RES['PCA_C_AR'] <- evaluateCorrection(CP$AR, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright', legend.col = 2)
            title(main = paste("Q: Interbatch distance\nANCOVA + RUV-corr.:", round(as.numeric(RES[['PCA_C_AR']]), 3)))
          }

        }
      }
    }else{
      if(Duplo){
        RES[['Duplo_Raw']] <- evaluateCorrection(set.1, set.1.Y, what = "duplo", plot = PLOT)
        RES[['Duplo_C_ANCOVA']] <- evaluateCorrection(CP$ANCOVA, set.1.Y, what = "duplo", plot = PLOT)
        if(QC){
          RES[['Duplo_C_RUV']] <- evaluateCorrection(CP$RUVs, set.1.Y, what = "duplo", plot = PLOT)
          RES['Duplo_C_AR'] <- evaluateCorrection(CP$AR, set.1.Y, what = "duplo", plot = PLOT)
        }
      }
      if(PCA){
        RES[['PCA_Raw']] <- evaluateCorrection(set.1, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright', legend.col = 1)
        RES[['PCA_C_ANCOVA']] <- evaluateCorrection(CP$ANCOVA, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright', legend.col = 2)
        if(QC){
          RES[['PCA_C_RUV']] <- evaluateCorrection(CP$RUVs, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright', legend.col = 2)
          RES['PCA_C_AR'] <- evaluateCorrection(CP$AR, set.1.Y, what = "PCA", plot = PLOT, legend.loc = 'bottomright', legend.col = 2)
        }
      }
    }
    cat('Computations complete!\n\n')

  }else{
    stop('Wrong format')
  }
  return(RES)
}
