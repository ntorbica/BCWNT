\name{BCW.correct}
\alias{BCW.correct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Batch effect correction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Set of functions to correct for batch effects. Performs corrections based on linear models and removal of unwanted variation. Correction algorithm written by R. Wehrens. Removal of unwanted variation algorithm from bioconductor package 'RUV'.
}
\usage{
BCW.correct(peaks, info, minBatchOccurrence.Ave = 2, minBatchOccurrence.Line = 4, condition = c("", "0", "1", "2", "c"), method = c("lm", "rlm"), k, QC = TRUE, PCA = TRUE, Duplo = TRUE, PLOT = TRUE, PDF = FALSE, csv = FALSE, nA = "ANCOVA.csv", nR = "RUV.csv", nAR = "ANCOVA_RUV.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{
%%     ~~Describe \code{peaks} here~~
Peak file outputted by BCW.format. Can be a list, matrix or data.frame object as well
}
  \item{info}{
%%     ~~Describe \code{info} here~~
Information file outputted by BCW.format. Can be a list, matrix or data.frame object as well
}
  \item{minQC}{
%%     ~~Describe \code{minBatchOccurrence.Line} here~~
Minimum quality control samples per batch
}
  \item{condition}{
%%     ~~Describe \code{condition} here~~
The condition defining the imputation for NAs, can be set to '' (ignore NAs), '0' (replace NAs with 0), '1' (replace NAs with the limit of detection), '2' (replace NAs with half of the detection limit) and 'c' (set censoring threshold to the limit of detection - 0.01, + 0.01 if limit of detection is 0). Note that condition 'c' will force set the method to 'tobit'.
}
  \item{method}{
%%     ~~Describe \code{method} here~~
The (univariate) correction method applied. Can be set to either 'lm' (linear regression), 'rlm' (robust linear regression). The 'tobit' method will be forced if the condition argument is set to 'c'
}
  \item{k}{
%%     ~~Describe \code{k} here~~
Number of principal components to model subspace for removal of unwanted variation, based on the presence of quality control samples.
}
  \item{QC}{
%%     ~~Describe \code{QC} here~~
Logical; specifies whether there are quality control samples in the dataset. If set to FALSE, no removal of unwanted variation can be performed.
}
  \item{PCA}{
%%     ~~Describe \code{PCA} here~~
Logical; specifies whether principal component analysis followed by computation of bhattacharyya distances of present batches (within the first two principal components).
}
  \item{Duplo}{
%%     ~~Describe \code{Duplo} here~~
Logical; specifies whether duplicate analysis in order to assess the repeatability / replicability.
}
  \item{PLOT}{
%%     ~~Describe \code{PLOT} here~~
Logical; specifies whether results should be plotted.
}
  \item{PDF}{
%%     ~~Describe \code{PDF} here~~
Logical; specifies whether a .pdf file of the correction plots is outputted.
}
  \item{csv}{
%%     ~~Describe \code{csv} here~~
Logical; specifies whether a .csv file of the corrected peaks is outputted.
}
  \item{nA}{
%%     ~~Describe \code{nA} here~~
Filename for outputting .csv file containing linear model corrected peaks
}
  \item{nR}{
%%     ~~Describe \code{nR} here~~
Filename for outputting .csv file containing removal of variation corrected peaks. Only necessary if QC = TRUE.
}
  \item{nAR}{
%%     ~~Describe \code{nAR} here~~
Filename for outputting .csv file containing linear model + removal of variation corrected peaks. Only necessary if QC = TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The execution of the function requires two input files, optimally in a separate directory (as created by BCW.format). Note that removal of unwanted variation correction will only take place if there are quality control samples available.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Generates .csv files with corrected data and a list object containing following components:

\item{BC_Data_ANCOVA}{
Linear regression corrected data matrix.
}
\item{BC_Data_RUV}{
Removal of unwanted variation corrected data matrix. Only generated if QC = TRUE.
}
\item{BC_Data_AR}{
Linear regression + removal of unwanted variation corrected data matrix. Only generated if QC = TRUE.
}

}
\references{
%% ~put references to the literature/web site here ~
"@ArticleWehrens2016, author =	 Ron Wehrens and Jos.~A.~Hageman and Fred~van~Eeuwijk and Rik~Kooke and P\'adraic~J.~Flood and Erik Wijnker and Joost~J.B.~Keurentjes and Arjen~Lommen and Henri\"ette~D.L.M.~van~Eekelen and Robert~D.~Hall and Roland~Mumm and Ric~C.H.~de~Vos, title =	 Improved batch correction in untargeted MS-based metabolomics, journal =	 Metabolomics, year =	 2016, volume =	 12, DOI =	 10.1007/s11306-016-1015-8, pages =	 1–12 "
}
\author{
%%  ~~who you are~~
Nenad Torbica, Ron Wehrens, (Bioconductor)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (peaks, info, minBatchOccurrence.Ave = 2, minBatchOccurrence.Line = 4,
    condition = c("", "0", "1", "2", "c"), method = c("lm", "rlm"),
    k, QC = TRUE, PCA = TRUE, Duplo = TRUE, PLOT = TRUE, PDF = FALSE,
    csv = FALSE, nA = "ANCOVA.csv", nR = "RUV.csv", nAR = "ANCOVA_RUV.csv")
{
    if (missing(peaks)) {
        stop("Missing Data! No peak table given.")
    }
    if (missing(info)) {
        stop("Missing Data! No information on peaks given.")
    }
    if (condition == "") {
        cat("You chose to ignore NA measurements.\n")
        impV <- NA
    }
    if (condition == "0") {
        cat("You chose to impute NA measurements with 0.\n")
    }
    if (condition == "1") {
        cat("You chose to impute NA measurements with half the limit of detection.\n")
    }
    if (condition == "2") {
        cat("You chose to impute NA measurements with the limit of detection (lowest measurement).\n")
    }
    if (condition == "c") {
        cat("You chose to impute Values for a censored regression. Note that this will require 'tobit' as method.\n")
        method <- "tobit"
    }
    cat("Preparing data...\n\n")
    suppressMessages(require(BatchCorrMetabolomics))
    suppressMessages(require(RUVSeq))
    type.p <- class(peaks)
    type.i <- class(info)
    if (type.p != "matrix") {
        if (type.p == "character") {
            peaks <- as.matrix(read.csv(peaks, row.names = 1,
                sep = ","))
            type.p <- class(peaks)
        }
        if (type.p == "list") {
            peaks <- as.matrix(as.data.frame(peaks))
            type.p <- class(peaks)
        }
        if (type.p == "data.frame") {
            peaks <- as.matrix(peaks)
            type.p <- class(peaks)
        }
    }
    if (type.i != "data.frame") {
        if (type.i == "character") {
            info <- as.data.frame(read.csv(info, row.names = 1))
            type.i <- class(info)
        }
        if (type.i == "list") {
            info <- as.data.frame(info)
            type.i <- class(info)
        }
    }
    if (type.p == "matrix" & type.i == "data.frame") {
        set.1 <- peaks
        set.1.Y <- info
        set.1.lod <- as.numeric(min(set.1[!is.na(set.1)]))
        if (condition == "0") {
            impV <- 0
        }
        if (condition == "1") {
            impV <- set.1.lod/2
        }
        if (condition == "2") {
            impV <- set.1.lod
        }
        if (condition == "c") {
            if (set.1.lod == 0) {
                impV <- set.1.lod + 0.01
            }
            else {
                impV <- set.1.lod - 0.01
            }
        }
        LC.ngenotypes <- nlevels(set.1.Y$SCode) - 1
        if (QC) {
            LC.nref <- sum(set.1.Y$SCode == "ref")
        }
        LC.nNA <- apply(set.1, 2, function(x) sum(is.na(x)))
        refSamples <- list(Q = which(set.1.Y$SCode == "ref"),
            S = which(set.1.Y$SCode != "ref"))
        if (QC) {
            actual <- list(str = "Q", exp = condition, method = method,
                impV = impV)
        }
        else {
            actual <- list(str = "S", exp = condition, method = method,
                impV = impV)
        }
        cat("Corrections applied as follows:\n\n        correction model:\t",
            actual$method, "\n", "\timputed values:\t\t", actual$impV,
            "\n", "\tstrategy:\t\t", actual$str, "\n\n")
        if (!QC) {
            cat("Note that RUV correction is not possible without QCs\n")
        }
        CP <- list()
        CP$ANCOVA <- apply(set.1, 2, doBC, ref.idx = as.numeric(refSamples[[actual$str]]),
            batch.idx = set.1.Y$Batch, minBsamp = minBatchOccurrence.Line,
            seq.idx = as.numeric(set.1.Y$SeqNr), method = actual$method,
            imputeVal = actual$impV)
        idx <- which(set.1.Y$SCode == "ref")
        replicates.ind <- matrix(-1, nrow(set.1) - length(idx) +
            1, length(idx))
        replicates.ind[1, ] <- idx
        replicates.ind[-1, 1] <- (1:nrow(set.1))[-idx]
        nColumns <- ncol(set.1)
        if (QC) {
            CP$RUVs <- t(RUVs(t(set.1), 1:nColumns, k = k, replicates.ind,
                round = FALSE, isLog = TRUE)$normalizedCounts)
            CP$AR <- t(RUVs(t(CP$ANCOVA), 1:nColumns, k = 3,
                replicates.ind, round = FALSE, isLog = TRUE)$normalizedCounts)
        }
        if (csv) {
            if (!dir.exists("CORRECTED_PEAKS/")) {
                dir.create("CORRECTED_PEAKS")
            }
            else {
                Sys.chmod("CORRECTED_PEAKS/")
                cat("Removing old directory...\n")
                unlink("CORRECTED_PEAKS/", recursive = TRUE)
            }
            if (missing(nA)) {
                nA <- readline("Enter a filename for ANCOVA corrected data output, ending in .csv:\n")
            }
            if (QC) {
                if (missing(nR)) {
                  nR <- readline("Enter a filename for RUV normlized data output, ending in .csv:\n")
                }
                if (missing(nAR)) {
                  nAR <- readline("Enter a filename for ANCOVA corrected + RUV normlized data output, ending in .csv:\n")
                }
            }
            cat("Creating new directory...\n")
            Sys.chmod("CORRECTED_PEAKS/")
            write.csv(cbind(set.1.Y, CP$ANCOVA), gsub(" ", "",
                paste("./CORRECTED_PEAKS/", nA)))
            if (QC) {
                write.csv(cbind(set.1.Y, CP$RUVs), gsub(" ",
                  "", paste("./CORRECTED_PEAKS/", nR)))
                write.csv(cbind(set.1.Y, CP$AR), gsub(" ", "",
                  paste("./CORRECTED_PEAKS/", nAR)))
            }
        }
        RES <- list()
        RES[["BC_Data_ANCOVA"]] <- CP$ANCOVA
        if (QC) {
            RES[["BC_Data_RUV"]] <- CP$RUVs
            RES[["BC_Data_AR"]] <- CP$AR
        }
        if (PLOT) {
            cat("Plotting...\n")
            plot.new()
            textplot(paste("Batch Correctio Wehrens - Pipeline\n\n",
                "Nr. of Samples: ", length(refSamples$S), "\n",
                "Nr. of QCs      ", length(refSamples$Q), "\n\n",
                "Correction:\n\n ", "Strategy:       ", ifelse(QC,
                  "Q", "S"), "\n", " Method:         ", ifelse(condition ==
                  "c", "tobit", method), "\n", " Condition:      ",
                condition, "\n", " Imputed value:  ", actual$impV,
                "\n\n", ifelse(QC, paste("RUV:\n\n nr. of PCs:     ",
                  k, "\n\n"), "")), halign = "left", valign = "top")
            if (PDF) {
                name <- readline("Please enter a filename for the plot file, ending with .pdf\n")
                pdf(name)
                par(mfrow = c(2, 2))
                if (Duplo) {
                  RES["Duplo_Raw"] <- evaluateCorrection(set.1,
                    set.1.Y, what = "duplo", plot = PLOT)
                  title(main = paste("Repeatabilities:", round(as.numeric(RES[["Duplo_Raw"]]),
                    3)))
                  RES["Duplo_C_ANCOVA"] <- evaluateCorrection(CP$ANCOVA,
                    set.1.Y, what = "duplo", plot = PLOT)
                  title(main = paste("Q: Repeatabilities\nANCOVA-corr.:",
                    round(as.numeric(RES[["Duplo_C_ANCOVA"]]),
                      3)))
                  if (QC) {
                    RES["Duplo_C_RUV"] <- evaluateCorrection(CP$RUVs,
                      set.1.Y, what = "duplo", plot = PLOT)
                    title(main = paste("Q: Repeatabilities\nRUV-corr.:",
                      round(as.numeric(RES[["Duplo_C_RUV"]]),
                        3)))
                    RES["Duplo_C_AR"] <- evaluateCorrection(CP$AR,
                      set.1.Y, what = "duplo", plot = PLOT)
                    title(main = paste("Q: Repeatabilities\nRUV-corr.:",
                      round(as.numeric(RES[["Duplo_C_AR"]]),
                        3)))
                  }
                }
                if (PCA) {
                  RES["PCA_Raw"] <- evaluateCorrection(set.1,
                    set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                    legend.col = 2)
                  title(main = paste("Interbatch distance:",
                    round(as.numeric(RES[["PCA_Raw"]]), 3)))
                  RES["PCA_C_ANCOVA"] <- evaluateCorrection(CP$ANCOVA,
                    set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                    legend.col = 2)
                  title(main = paste("Q: Interbatch distance\nANCOVA-corr.:",
                    round(as.numeric(RES[["PCA_C_ANCOVA"]]),
                      3)))
                  if (QC) {
                    RES["PCA_C_RUV"] <- evaluateCorrection(CP$RUVs,
                      set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                      legend.col = 2)
                    title(main = paste("Q: Interbatch distance\nRUV-corr.:",
                      round(as.numeric(RES[["PCA_C_RUV"]]), 3)))
                    RES["PCA_C_AR"] <- evaluateCorrection(CP$AR,
                      set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                      legend.col = 2)
                    title(main = paste("Q: Interbatch distance\nRUV-corr.:",
                      round(as.numeric(RES[["PCA_C_AR"]]), 3)))
                  }
                }
                dev.off()
                cat(name, "saved in", getwd(), "\n")
            }
            else {
                if (Duplo) {
                  par(mfrow = c(2, 2))
                  RES["Duplo_Raw"] <- evaluateCorrection(set.1,
                    set.1.Y, what = "duplo", plot = PLOT)
                  title(main = paste("Repeatabilities:", round(as.numeric(RES[["Duplo_Raw"]]),
                    3)))
                  RES["Duplo_C_ANCOVA"] <- evaluateCorrection(CP$ANCOVA,
                    set.1.Y, what = "duplo", plot = PLOT)
                  title(main = paste("Q: Repeatabilities\nANCOVA-corr.:",
                    round(as.numeric(RES[["Duplo_C_ANCOVA"]]),
                      3)))
                  if (QC) {
                    RES["Duplo_C_RUV"] <- evaluateCorrection(CP$RUVs,
                      set.1.Y, what = "duplo", plot = PLOT)
                    title(main = paste("Q: Repeatabilities\nRUV-corr.:",
                      round(as.numeric(RES[["Duplo_C_RUV"]]),
                        3)))
                    RES["Duplo_C_AR"] <- evaluateCorrection(CP$AR,
                      set.1.Y, what = "duplo", plot = PLOT)
                    title(main = paste("Q: Repeatabilities\nANCOVA + RUV-corr.:",
                      round(as.numeric(RES[["Duplo_C_AR"]]),
                        3)))
                  }
                }
                if (PCA) {
                  RES["PCA_Raw"] <- evaluateCorrection(set.1,
                    set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                    legend.col = 2)
                  title(main = paste("Interbatch distance:",
                    round(as.numeric(RES[["PCA_Raw"]]), 3)))
                  RES["PCA_C_ANCOVA"] <- evaluateCorrection(CP$ANCOVA,
                    set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                    legend.col = 2)
                  title(main = paste("Q: Interbatch distance\nANCOVA-corr.:",
                    round(as.numeric(RES[["PCA_C_ANCOVA"]]),
                      3)))
                  if (QC) {
                    RES["PCA_C_RUV"] <- evaluateCorrection(CP$RUVs,
                      set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                      legend.col = 2)
                    title(main = paste("Q: Interbatch distance\nRUV-corr.:",
                      round(as.numeric(RES[["PCA_C_RUV"]]), 3)))
                    RES["PCA_C_AR"] <- evaluateCorrection(CP$AR,
                      set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                      legend.col = 2)
                    title(main = paste("Q: Interbatch distance\nANCOVA + RUV-corr.:",
                      round(as.numeric(RES[["PCA_C_AR"]]), 3)))
                  }
                }
            }
        }
        else {
            if (Duplo) {
                RES[["Duplo_Raw"]] <- evaluateCorrection(set.1,
                  set.1.Y, what = "duplo", plot = PLOT)
                RES[["Duplo_C_ANCOVA"]] <- evaluateCorrection(CP$ANCOVA,
                  set.1.Y, what = "duplo", plot = PLOT)
                if (QC) {
                  RES[["Duplo_C_RUV"]] <- evaluateCorrection(CP$RUVs,
                    set.1.Y, what = "duplo", plot = PLOT)
                  RES["Duplo_C_AR"] <- evaluateCorrection(CP$AR,
                    set.1.Y, what = "duplo", plot = PLOT)
                }
            }
            if (PCA) {
                RES[["PCA_Raw"]] <- evaluateCorrection(set.1,
                  set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                  legend.col = 1)
                RES[["PCA_C_ANCOVA"]] <- evaluateCorrection(CP$ANCOVA,
                  set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                  legend.col = 2)
                if (QC) {
                  RES[["PCA_C_RUV"]] <- evaluateCorrection(CP$RUVs,
                    set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                    legend.col = 2)
                  RES["PCA_C_AR"] <- evaluateCorrection(CP$AR,
                    set.1.Y, what = "PCA", plot = PLOT, legend.loc = "bottomright",
                    legend.col = 2)
                }
            }
        }
        cat("Computations complete!\n\n")
    }
    else {
        stop("Wrong format")
    }
    return(RES)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
