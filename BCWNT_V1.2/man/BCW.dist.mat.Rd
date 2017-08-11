\name{BCW.dist.mat}
\alias{BCW.dist.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Bhattacharyya distances for all possible principal component combinations.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes bhattacharyya distances for each principal component combination, generating a diagonal matrix with as many principal components as there are samples.
}
\usage{
BCW.dist.mat(d, heatmap = TRUE, PCA.plots = TRUE, PCA.c = 5, PDF = FALSE, csv = FALSE, t = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
%%     ~~Describe \code{d} here~~
Data as found in compact format file (three information columns followed by all the peaks).
}
  \item{heatmap}{
%%     ~~Describe \code{heatmap} here~~
Logical; specifies whether a heatmap of the distance matrix should be plotted.
}
  \item{PCA.plots}{
%%     ~~Describe \code{PCA.plots} here~~
Logical; specifies whether plots of the first PC.c (count) principal combinations should be plotted. Basically shows whether the batch effect is within the first few component.
}
  \item{PCA.c}{
%%     ~~Describe \code{PCA.c} here~~
Principal component count for plotting PC-combinations.
}
  \item{PDF}{
%%     ~~Describe \code{PDF} here~~
Logical; specifies whether a .pdf file should be outputted containing all the plots.
}
  \item{csv}{
%%     ~~Describe \code{csv} here~~
Logical; specifies whether the computed distance matrix should be outputted as .csv file
}
  \item{t}{
%%     ~~Describe \code{t} here~~
Title for the plots.
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
Returns a matrix object and a file containing the diagonal matrix with batch distances for each principal component combination.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
C. Hennig, Nenad Torbica
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
function (d, heatmap = TRUE, PCA.plots = TRUE, PCA.c = 5, PDF = FALSE,
    csv = FALSE, t = "")
{
    suppressMessages(require(BatchCorrMetabolomics))
    suppressMessages(require(RUVSeq))
    suppressMessages(require(ChemometricsWithR))
    suppressMessages(require(fpc))
    suppressMessages(require(gplots))
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
    npc <- ncol(X.PCA$scores)
    Xscores <- scores.PCA(X.PCA)[, 1:npc, drop = FALSE]
    batch.means <- lapply(levels(YSample$Batch), function(btch) colMeans(Xscores[which(YSample$Batch ==
        btch), , drop = FALSE]))
    batch.covs <- lapply(levels(YSample$Batch), function(btch) cov(Xscores[which(YSample$Batch ==
        btch), , drop = FALSE]))
    noCov.idx <- which(sapply(batch.covs, function(x) all(x <
        1e-08)))
    b1m <- t(matrix(unlist(batch.means[[1]]), 1, npc))
    rownames(b1m) <- sprintf("PC\%d", 1:npc)
    colnames(b1m) <- sprintf("Batch \%d", 1)
    b2m <- t(matrix(unlist(batch.means[[2]]), 1, npc))
    rownames(b2m) <- sprintf("PC\%d", 1:npc)
    colnames(b2m) <- sprintf("Batch \%d", 2)
    b1c <- matrix(0, npc, npc)
    b1c.v <- unlist(batch.covs[1])
    for (i in 1:length(b1c.v)) {
        b1c[i] <- b1c.v[i]
    }
    colnames(b1c) <- sprintf("PC\%d", 1:npc)
    rownames(b1c) <- sprintf("PC\%d", 1:npc)
    b2c <- matrix(0, npc, npc)
    b2c.v <- unlist(batch.covs[2])
    for (i in 1:length(b2c.v)) {
        b2c[i] <- b2c.v[i]
    }
    colnames(b2c) <- sprintf("PC\%d", 1:npc)
    rownames(b2c) <- sprintf("PC\%d", 1:npc)
    cov.pairlist1 <- list()
    cov.pairlist2 <- list()
    x <- 1
    for (i in 1:(npc - 1)) {
        for (j in (i + 1):npc) {
            cov.pairlist1[[x]] <- matrix(b1c[c(i, j), c(i, j)],
                2, 2)
            cov.pairlist2[[x]] <- matrix(b2c[c(i, j), c(i, j)],
                2, 2)
            x <- x + 1
        }
    }
    bm <- cbind(b1m, b2m)
    sigmaarray <- array(c(unlist(cov.pairlist1), unlist(cov.pairlist2)),
        dim = c(2, 2, length(bm)))
    cat("Computing distances...\n")
    dist.mat <- bhattacharyya.matrix(t(bm), sigmaarray, ipairs = "all")
    if (csv) {
        filename <- readline("Please enter a filename for the distance matrix, ending in .csv:\n\n")
        write.csv(dist.mat, filename)
        cat(filename, "saved in", getwd(), "\n\n")
    }
    if (PDF) {
        pdfname <- readline("Please enter a filename for the plot file, ending in .pdf:\n\n")
        pdf(pdfname)
        if (heatmap) {
            heatmap.2(t(dist.mat), dendrogram = "none", Rowv = F,
                Colv = F, trace = "none", main = paste(t, "\nBhattacharyya Distances PCA"))
        }
        if (PCA.plots) {
            b <- info$Batch[noref.idx]
            pairs(Xscores[, 1:PCA.c], col = b, pch = ".", main = paste(t,
                "\nPCA distances"))
        }
        suppressMessages(dev.off())
        cat(pdfname, "saved in", getwd(), "\n\n")
    }
    else {
        if (heatmap) {
            heatmap.2(t(dist.mat), dendrogram = "none", Rowv = F,
                Colv = F, trace = "none", main = paste(t, "\nBhattacharyya Distances PCA"))
        }
        if (PCA.plots) {
            s <- info$Batch[-noref.idx]
            b <- info$Batch[noref.idx]
            pairs(Xscores[, 1:PCA.c], col = b, pch = ".", main = paste(t,
                "\nPCA combinations"))
        }
    }
    cat("Computation complete!\n\n")
    return(dist.mat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
