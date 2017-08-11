
## Read Function to get the peaks and infos from a LC-MS data table
BCW.format <- function(file, n.mass = FALSE, QC, G, B, S, mkdir = FALSE, ts = ','){
  # INPUTS:   file    - name of the Datafile, has to be .csv format
  #           n.mass  - logical, states whether there are neutral mass annotations
  #                     for peaks which originate as SIMCA output
  #           QC      - Naming for quality control samples
  #           G       - Naming of column with group labels; includes the QC-labels
  #           B       - Naming of column with batch labels
  #           S       - Naming of column containg injection sequence
  #           mkdir   - logical, states if a directory with output files should be made
  #           ts      - table separator for read-in

  # The mkdir argument will remove existing directories with same name,
  # thus throw a warning to the user
  if(mkdir){
    if(dir.exists('Wehrens Input/')== TRUE){
      X <- readline("\n\n\t\t\t WARNING:   \n\nthe mkdir argument removes existing directory with name 'Wehrens Input' with its content.
                 \n\t\t\tContinue? (y/n)\n\n")
    }
    else{
      X <- 'y'
    }

    if(X != 'y' & X != 'n'){
      stop('Please rerun and give the answer with either y or n')
    }

    if(X == 'n'){
      options(show.error.messages = FALSE)
      stop("\n\nRerun when ready, good Madam / Sir.")
    }
    else{
      cat('Building Tables...\n')
    }
  }
  # Initialize output objects
  OUT <- list()
  PEAKS <- list()
  INFO <- list()

  # Require a QC naming:
  if(missing(QC)){
    stop("\n\n\tPlease enter a name for quality control/reference samples in your dataset (argument QC).\n")
  }

  # Require a group naming
  if(missing(G)){
    stop("\n\n\tPlease enter the column name for the grouping labels in your dataset (argument G).\n")
  }

  # Require a batch naming
  if(missing(B)){
    stop("\n\n\tPlease enter the column name for the batch labels (argument B).\n")
  }

  # Require a sequence naming
  if(missing(S)){
    stop("\n\n\tPlease enter the column name for the injection sequence (argument S).\n")
  }

  # Check input type
  type <- class(file)

  # Read input
  if(type == 'character'){
    tab <- as.matrix(read.csv(file, sep=ts, row.names = 1))
  }
  if(type == 'list'){
    tab <- as.matrix(as.data.frame(file))
  }
  if(type == 'data.frame'){
    tab <- as.matrix(file)
  }



  # match all the peaks and put them into a matrix
  s <- '^X'
  peaks <- gregexpr(s, colnames(tab))

  # if the matching turns out negative, I add the 'X' to the colname, rather than just alternatively
  # searching for the names. Just for integrity.
  if(length(peaks)-length(which(peaks==-1)) == 0){
    alt.s <- '[0-9]+.[0-9]{2}_'
    peaks <- grep(alt.s, colnames(tab))
    cnames <- colnames(tab)
    for(i in peaks){
      cnames[i] <- gsub(' ', '', paste('X', cnames[i]))
    }
    peak.tab <- tab[,peaks]
  }else{
    peak.tab <- tab[,which(peaks!=-1)]
  }

  PEAKS[[1]] <- as.matrix(peak.tab)
  # names(PEAKS) <- 'all'

  if(n.mass){
    ##################################################################
    # Additionally, output a separate table for neutral mass values  #
    # and m/z values for the features:                               #
    ##################################################################

    # match samples with neutral mass values and save indices to
    # a list. The remaining indices are m/z values
    n.s <- 'n$'
    peaks.n <- gregexpr(n.s, colnames(peak.tab))
    peaks.mz <- as.numeric(which(peaks.n == -1))
    peaks.n <- as.numeric(which(peaks.n != -1))

    peak.tab.mz <- peak.tab[,peaks.mz]
    peak.tab.n <- peak.tab[,peaks.n]

    PEAKS[[2]] <- as.matrix(peak.tab[,peaks.mz])
    PEAKS[[3]] <- as.matrix(peak.tab[,peaks.n])
    names(PEAKS) <- c('all', 'mz', 'n')
  }

  OUT[[1]] <- PEAKS

  # As the peaks are already defined, the information table is easy
  # to generate, requiring injection sequence, grouping information
  # including QC labels, and batch information.
  if(length(peaks)-length(which(peaks==-1)) == 0){
    info <- tab[,peaks]
  }else{
    info <- tab[,which(peaks==-1)]
  }

  # # Compare order columns and take the one with the lowest start value.
  # ord.comp <- matrix(0,nrow(tab),length(info.o))
  # for(n in 1:length(info.o)){
  #   ord.comp[,n] <- info[,info.order[n]]
  # }
  # ord <- which(ord.comp[1,] == min(ord.comp[1,]))

  # Get order info according to user input.
  info.s <- grep(S, colnames(info), ignore.case = TRUE)
  # Get the grouping information according to user input
  # info.g <- unique(grep(paste(G, collapse = '|'), colnames(info), ignore.case = TRUE))
  info.g <- grep(G, colnames(info), ignore.case = TRUE)

  # Get the batch informatioin according to user input
  info.b <- grep(B, colnames(info),  ignore.case = TRUE)

  # Construct Information table
  rownames(info) <- rownames(tab)
  qc.id <- which(info[,info.g] == QC)
  info[qc.id,info.g] <- as.character('ref')
  INFO[[1]] <- info[,c(info.s, info.g, info.b)]
  colnames(INFO[[1]]) <- c('SeqNr', 'SCode', 'Batch')

  OUT[[2]] <- as.data.frame(INFO)


  # Create directory with .csv output if specified
  if(mkdir){
    Sys.chmod('.')

    if(dir.exists('Wehrens Input/')){
      Sys.chmod('./Wehrens Input/')
      cat('Removing old directory...\n')
      unlink('./Wehrens Input', recursive = TRUE)
    }

    dir.create('./Wehrens Input')

    if(n.mass){
      m <- c('all','mz', 'n')
      for(i in m){
        n <- sprintf('./Wehrens Input/Peaks_%s.csv', i)
        write.csv(OUT[[1]][i], n)
        write.csv(OUT[[2]], './Wehrens Input/Info.csv')
      }
    }else{
      Sys.chmod('.')
      write.csv(OUT[1], './Wehrens Input/Peaks.csv')
      write.csv(OUT[2], './Wehrens Input/Info.csv')
    }
    cat("Tables are prepared! Find them in the 'Wehrens Input' directory.\n\n")
  }

  invisible(OUT)

}

#################################################################
# The question remains how to treat the values termed 'neutral' #
# For now, I'll just stick to what we have done before, thus no #
# accounting for this circumstance at all. Also, the function   #
# currently outputs list objects containing the data if not     #
# written as .csv. Will have to rework the output.              #
#################################################################

