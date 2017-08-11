\name{BCW.format}
\alias{BCW.format}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
File formatting for batch correction pipline BCW()
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generate input file with preformatted LC-MS .csv file containing information on injection sequence, sample grouping and batch information, as well as the peaks. Such a file is outputted by BCW.s2w from raw format output files.
}
\usage{
BCW.format(file, n.mass = FALSE, QC, G, B, S, mkdir = FALSE, ts = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
File name of the preformatted, compact format file, its filename beginning with 'S2W_'.
}
  \item{n.mass}{
%%     ~~Describe \code{n.mass} here~~
Logical; specifies whether neutral masses are present in the peak data (represented by 'n' instead of 'mz' in the feature-names).
}
  \item{QC}{
%%     ~~Describe \code{QC} here~~
Quality control label as found in the sample grouping within the input file.
}
  \item{G}{
%%     ~~Describe \code{G} here~~
Column name within input file containing the Grouping information.
}
  \item{B}{
%%     ~~Describe \code{B} here~~
Column name within input file containing the Batch labels.
}
  \item{S}{
%%     ~~Describe \code{S} here~~
Column name within input file containing the injection Sequence.
}
  \item{mkdir}{
%%     ~~Describe \code{mkdir} here~~
Logical; specifies whether a directory with output files should be generated.
}
  \item{ts}{
%%     ~~Describe \code{ts} here~~
Table separator for the input file.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function is initially executed as formatting step in BCW() if SIMCA = FALSE. It is the second step within the pipeline if SIMCA = TRUE.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list object:

\item{peaks}{
The peak table without other information.
}
\item{info}{
The info table without any peaks.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Nenad Torbica
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
function (file, n.mass = FALSE, QC, G, B, S, mkdir = FALSE, ts = ",")
{
    if (mkdir) {
        if (dir.exists("Wehrens Input/") == TRUE) {
            X <- readline("\n\n\t\t\t WARNING:   \n\nthe mkdir argument removes existing directory with name 'Wehrens Input' with its content.\n                 \n\t\t\tContinue? (y/n)\n\n")
        }
        else {
            X <- "y"
        }
        if (X != "y" & X != "n") {
            stop("Please rerun and give the answer with either y or n")
        }
        if (X == "n") {
            options(show.error.messages = FALSE)
            stop("\n\nRerun when ready, good Madam / Sir.")
        }
        else {
            cat("Building Tables...\n")
        }
    }
    OUT <- list()
    PEAKS <- list()
    INFO <- list()
    if (missing(QC)) {
        stop("\n\n\tPlease enter a name for quality control/reference samples in your dataset (argument QC).\n")
    }
    if (missing(G)) {
        stop("\n\n\tPlease enter the column name for the grouping labels in your dataset (argument G).\n")
    }
    if (missing(B)) {
        stop("\n\n\tPlease enter the column name for the batch labels (argument B).\n")
    }
    if (missing(S)) {
        stop("\n\n\tPlease enter the column name for the injection sequence (argument S).\n")
    }
    type <- class(file)
    if (type == "character") {
        tab <- as.matrix(read.csv(file, sep = ts, row.names = 1))
    }
    if (type == "list") {
        tab <- as.matrix(as.data.frame(file))
    }
    if (type == "data.frame") {
        tab <- as.matrix(file)
    }
    s <- "^X"
    peaks <- gregexpr(s, colnames(tab))
    if (length(peaks) - length(which(peaks == -1)) == 0) {
        alt.s <- "[0-9]+.[0-9]{2}_"
        peaks <- grep(alt.s, colnames(tab))
        cnames <- colnames(tab)
        for (i in peaks) {
            cnames[i] <- gsub(" ", "", paste("X", cnames[i]))
        }
        peak.tab <- tab[, peaks]
    }
    else {
        peak.tab <- tab[, which(peaks != -1)]
    }
    PEAKS[[1]] <- as.matrix(peak.tab)
    if (n.mass) {
        n.s <- "n$"
        peaks.n <- gregexpr(n.s, colnames(peak.tab))
        peaks.mz <- as.numeric(which(peaks.n == -1))
        peaks.n <- as.numeric(which(peaks.n != -1))
        peak.tab.mz <- peak.tab[, peaks.mz]
        peak.tab.n <- peak.tab[, peaks.n]
        PEAKS[[2]] <- as.matrix(peak.tab[, peaks.mz])
        PEAKS[[3]] <- as.matrix(peak.tab[, peaks.n])
        names(PEAKS) <- c("all", "mz", "n")
    }
    OUT[[1]] <- PEAKS
    if (length(peaks) - length(which(peaks == -1)) == 0) {
        info <- tab[, peaks]
    }
    else {
        info <- tab[, which(peaks == -1)]
    }
    info.s <- grep(S, colnames(info), ignore.case = TRUE)
    info.g <- grep(G, colnames(info), ignore.case = TRUE)
    info.b <- grep(B, colnames(info), ignore.case = TRUE)
    rownames(info) <- rownames(tab)
    qc.id <- which(info[, info.g] == QC)
    info[qc.id, info.g] <- as.character("ref")
    INFO[[1]] <- info[, c(info.s, info.g, info.b)]
    colnames(INFO[[1]]) <- c("SeqNr", "SCode", "Batch")
    OUT[[2]] <- as.data.frame(INFO)
    if (mkdir) {
        Sys.chmod(".")
        if (dir.exists("Wehrens Input/")) {
            Sys.chmod("./Wehrens Input/")
            cat("Removing old directory...\n")
            unlink("./Wehrens Input", recursive = TRUE)
        }
        dir.create("./Wehrens Input")
        if (n.mass) {
            m <- c("all", "mz", "n")
            for (i in m) {
                n <- sprintf("./Wehrens Input/Peaks_\%s.csv",
                  i)
                write.csv(OUT[[1]][i], n)
                write.csv(OUT[[2]], "./Wehrens Input/Info.csv")
            }
        }
        else {
            Sys.chmod(".")
            write.csv(OUT[1], "./Wehrens Input/Peaks.csv")
            write.csv(OUT[2], "./Wehrens Input/Info.csv")
        }
        cat("Tables are prepared! Find them in the 'Wehrens Input' directory.\n\n")
    }
    invisible(OUT)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
