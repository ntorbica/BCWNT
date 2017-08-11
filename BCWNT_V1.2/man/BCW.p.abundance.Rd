\name{BCW.p.abundance}
\alias{BCW.p.abundance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute p-values by anova and generate histogram of unadjusted and adjusted p-values
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to evaluate multivariate data on p-value abundances, and plot the resulting abundances. Computes p-values for 'Order' (injection sequence), 'Site' (Batch) and 'Patient' (sample grouping) as response variable.
}
\usage{
BCW.p.abundance(d, PLOT = TRUE, PDF = FALSE, csv = FALSE, t = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
%%     ~~Describe \code{d} here~~
Data file as found in compact format input for correction.
}
  \item{PLOT}{
%%     ~~Describe \code{PLOT} here~~
Logical; specifies whether histograms of the abundances should be outputted.
}
  \item{PDF}{
%%     ~~Describe \code{PDF} here~~
Logical; specifies whether a .pdf file with the result histograms should be outputted.
}
  \item{csv}{
%%     ~~Describe \code{csv} here~~
Logical; specifies whether a .csv file with the p-values (unadjusted separaded from adjusted values by a column containing only 0).
}
  \item{t}{
%%     ~~Describe \code{t} here~~
Title for the plots
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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Chris Nakas, Nenad Torbica
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
function (d, PLOT = TRUE, PDF = FALSE, csv = FALSE, t = "")
{
    cat("Preparing data...\n")
    type <- class(d)
    if (type == "character") {
        d <- read.csv(d, row.names = 1)
    }
    if (type == "list") {
        d <- as.data.frame(d)
    }
    info <- d[, c("SeqNr", "SCode", "Batch")]
    peaks <- d[, -(which(colnames(d) \%in\% colnames(info)))]
    n.c <- ncol(d)
    n.p <- ncol(peaks)
    sitep <- rep(0, n.p)
    timep <- rep(0, n.p)
    patientp <- rep(0, n.p)
    if (class(d$SeqNr) != "integer") {
        d$SeqNr <- as.integer(d$SeqNr)
    }
    cat("Computing p-values for patients, injection sequence and batch. This may take a while...\n")
    for (i in 4:n.c) {
        crudelm <- aov(d[, i] ~ factor(SCode) + Batch + SeqNr,
            data = d)
        viasumunl <- unlist(summary(crudelm))
        patientp[i - 3] <- viasumunl[["Pr(>F)1"]]
        sitep[i - 3] <- viasumunl[["Pr(>F)2"]]
        timep[i - 3] <- viasumunl[["Pr(>F)3"]]
    }
    p.abundances <- matrix()
    p.abundances <- cbind(patientp, sitep, timep)
    allpadjcb <- as.data.frame(apply(p.abundances, 2, p.adjust,
        method = "fdr", n = nrow(p.abundances)))
    p.ALL <- cbind(p.abundances, sep = rep(0, nrow(p.abundances)),
        allpadjcb)
    if (PLOT) {
        if (PDF) {
            pdfname <- readline("Please enter a name for the plotfile, ending in .pdf:\n\n")
            pdf(pdfname)
            par(mfrow = c(1, 3))
            hist(patientp, main = paste(t, "Patient"), xlab = "Patient p-values")
            hist(sitep, main = paste(t, "Site"), xlab = "Site p-values")
            hist(timep, main = paste(t, "Order"), xlab = "Order p-values")
            par(mfrow = c(1, 3))
            hist(allpadjcb$patientp, main = paste(t, "Patient adj"),
                xlab = "Patient p-values")
            hist(allpadjcb$sitep, main = paste(t, "Site adj"),
                xlab = "Site p-values")
            hist(allpadjcb$timep, main = paste(t, "Order adj"),
                xlab = "Order p-values")
            suppressMessages(dev.off())
        }
        else {
            par(mfrow = c(1, 3))
            hist(patientp, main = paste(t, "Patient"), xlab = "Patient p-values")
            hist(sitep, main = paste(t, "Site"), xlab = "Site p-values")
            hist(timep, main = paste(t, "Order"), xlab = "Order p-values")
            par(mfrow = c(1, 3))
            hist(allpadjcb$patientp, main = paste(t, "Patient adj"),
                xlab = "Patient p-values")
            hist(allpadjcb$sitep, main = paste(t, "Site adj"),
                xlab = "Site p-values")
            hist(allpadjcb$timep, main = paste(t, "Order adj"),
                xlab = "Order p-values")
        }
    }
    if (csv) {
        csvname <- readline("Please enter a name for the table file, ending in .csv:\n\n")
        write.csv(p.ALL, csvname)
    }
    cat("Computations complete!\n\n")
    return(p.ALL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
