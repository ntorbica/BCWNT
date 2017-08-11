\name{BCW.plot.loading}
\alias{BCW.plot.loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot loading vectors of principal components.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Produces plots of the loadings of n principal components.
}
\usage{
BCW.plot.loading(dat, n.PC = 5, PDF = FALSE, t = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
Data to evaluate. Can be a filename, list, or data.frame object.
}
  \item{n.PC}{
%%     ~~Describe \code{n.PC} here~~
Number of principal components to plot.
}
  \item{PDF}{
%%     ~~Describe \code{PDF} here~~
Logical; specifies whether a .pdf file should be generated.
}
  \item{t}{
%%     ~~Describe \code{t} here~~
Title for plots.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a set of plots or a .pdf. Does not output any objects.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
function (dat, n.PC = 5, PDF = FALSE, t = "")
{
    suppressMessages(require(BatchCorrMetabolomics))
    suppressMessages(require(RUVSeq))
    suppressMessages(require(ChemometricsWithR))
    suppressMessages(require(ggplot2))
    suppressMessages(require(reshape))
    suppressMessages(require(reshape2))
    cat("Preparing data...\n")
    type <- class(dat)
    if (type == "character") {
        dat <- as.data.frame(read.csv(dat, row.names = 1, sep = ";"))
    }
    if (type == "list") {
        dat <- as.data.frame(dat)
    }
    D <- dat
    info <- D[, c("SeqNr", "SCode", "Batch")]
    peaks <- D[, -(which(colnames(D) \%in\% colnames(info)))]
    nb <- nlevels(info$Batch)
    noref.idx <- which(info$SCode != "ref")
    Xsample <- peaks[noref.idx, ]
    YSample <- info[noref.idx, ]
    Xsample <- Xsample[, apply(Xsample, 2, function(x) !all(is.na(x)))]
    Xsample <- Xsample[, apply(Xsample, 2, sd, na.rm = TRUE) >
        0]
    for (i in 1:ncol(Xsample)) {
        Xsample[is.na(Xsample[, i]), i] <- mean(Xsample[, i],
            na.rm = TRUE)
    }
    Xsample <- scale(Xsample)
    X.PCA <- PCA(Xsample)
    Xload <- as.data.frame(abs(X.PCA$loadings))
    f.names <- rownames(Xload)
    if (PDF) {
        pdfname <- readline("Please enter a name for the plot file, ending in .pdf:\n\n")
        pdf(pdfname)
        for (i in 1:n.PC) {
            p <- ggplot(Xload, aes((1:length(f.names)), Xload[,
                i])) + geom_col() + theme(axis.text.x = element_text(angle = 90,
                size = 10)) + labs(x = "Features", y = "Loading") +
                ggtitle(label = t, subtitle = paste("Loadings PC",
                  i))
            print(p)
        }
        dev.off()
    }
    else {
        cat("Plotting...\n")
        for (i in 1:n.PC) {
            p <- ggplot(Xload, aes((1:length(f.names)), Xload[,
                i])) + geom_col() + theme(axis.text.x = element_text(angle = 90,
                size = 10)) + labs(x = "Features", y = "Loading") +
                ggtitle(label = t, subtitle = paste("Loadings PC",
                  i))
            print(p)
        }
    }
    cat("Computations complete!\n\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
