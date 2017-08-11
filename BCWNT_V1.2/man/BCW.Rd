\name{BCW}
\alias{BCW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Batch correction pipeline.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Performs formatting steps, correction, and evaluation. Inputs are searched in working directory, necessary arguments will be read if none are given. Outputs formatted files, corrected data and evaluations into the working directory.
}
\usage{
BCW(RAW = TRUE, ts.sw, i.cols, ts.f, QC = TRUE, QC.rem = FALSE, QC.n, G, B, S, condition = c("", "0", "1", "2", "c"), method = c("lm", "rlm"), k, PDF = TRUE, tc = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RAW}{
%%     ~~Describe \code{RAW} here~~
Logical; has to be specified. If TRUE, the function will search for two files ending with "_raw.csv" (raw peak file), "_ClinData.csv" (clinical data file), otherwise it looks for one compact format file beginning with "S2W_". This defines whether the conversion from raw to compact will happen.
}
  \item{ts.sw}{
%%     ~~Describe \code{ts.sw} here~~
Table separator for raw files, for .csv files usually ',' or ';'. Make sure they are the same, the function will not be capable of proccessing the information properly. Only needs to be specify RAW = TRUE
}
  \item{i.cols}{
%%     ~~Describe \code{i.cols} here~~
Vector of three column names within the clinical data file, containing the sample naming, the sample grouping and the batch information. Only has to be specified if RAW = T
}
  \item{ts.f}{
%%     ~~Describe \code{ts.f} here~~
Table separator for compact format file, only needs to be specified if RAW = FALSE
}
  \item{QC}{
%%     ~~Describe \code{QC} here~~
Logical; specifies whether the dataset contains quality controls.
}
  \item{QC.rem}{
%%     ~~Describe \code{QC.rem} here~~
Logical; specifies whether present quality controls should be removed
}
  \item{QC.n}{
%%     ~~Describe \code{QC.n} here~~
Specifies the label for quality controls, such as "QC-ALL" or "QCpool-All". Should be present within the sample grouping column specified in i.cols or G of the clinical data file.
}
  \item{G}{
%%     ~~Describe \code{G} here~~
Specifies the column name within the compact file containing the sample grouping.
}
  \item{B}{
%%     ~~Describe \code{B} here~~
Specifies the column name within the compact file containing the batch information.
}
  \item{S}{
%%     ~~Describe \code{S} here~~
Specifies the column name within the compact file containing the injection sequence.
}
  \item{condition}{
%%     ~~Describe \code{condition} here~~
Decides which imputation for NA measurements is chosen. Can be set to '', '0','1','2' or 'c'. The conditions result in the following imputed values:

'': NAs are ignored

0: NAs are set to 0

1: NAs are set to the limit of detection

2: NAs are set to half of the limit of detection

c: NAs are set to the limit of detection - 0.01. Note that this will force 'tobit' as method.
}
  \item{method}{
%%     ~~Describe \code{method} here~~
Defines the method for correction. Can be set to 'lm' (linear model) and 'rlm' (robust linear model), though the method will be set to 'tobit' (censored regression) if condition is set to 'c'. Note that the corrections defined by the method  argument are univariate.
}
  \item{k}{
%%     ~~Describe \code{k} here~~
Defines the number of prinipal components used by RUVs (Bioconductor package RUV).
}
  \item{PDF}{
%%     ~~Describe \code{PDF} here~~
Logical; specifies if a pdf containing evaluation plots will be outputted.
}
  \item{tc}{
%%     ~~Describe \code{tc} here~~
Logical; specifies if the pdf should contain a table of content
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
In order to run the function, the working directory should optimally contain only the input files needed. The pipeline will come to a halt if multiple files are found as input.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Several objects are produced during the course of the function, the main output is intended to be as new files. The value of each of the objects produced by the pipeline functions can be found in the corresponding documentation.
}
\references{
%% ~put references to the literature/web site here ~
"@ArticleWehrens2016, author =	 Ron Wehrens and Jos.~A.~Hageman and Fred~van~Eeuwijk and Rik~Kooke and P\'adraic~J.~Flood and Erik Wijnker and Joost~J.B.~Keurentjes and Arjen~Lommen and Henri\"ette~D.L.M.~van~Eekelen and Robert~D.~Hall and Roland~Mumm and Ric~C.H.~de~Vos, title =	 Improved batch correction in untargeted MS-based metabolomics, journal =	 Metabolomics, year =	 2016, volume =	 12, DOI =	 10.1007/s11306-016-1015-8, pages =	 1–12"
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
###############################################################################
#                                                                             #
# Test script: example for BCWNT correction pipeline, function BCW            #
#                                                                             #
###############################################################################



data("Test_raw")                                # load example data (reduced LC-MS dataset)
data("Test_ClinData")
R <- Test_raw
CD <- Test_ClinData

######################################################################################################################

dir.create('TEST_raw_QC')                       # Create and set demo-directory (correction with quality controls)
setwd('TEST_raw_QC')
write.csv(R, "Test_raw.csv")
write.csv(CD, "Test_ClinData.csv")


# Correct raw data with quality control samples by linear regression
BCW(RAW = T,
    ts.sw = ',',
    i.cols = c('File.Name', 'Sample.group', 'Site'),
    S = 'SeqNr',
    QC = TRUE,
    QC.rem = FALSE,
    QC.n = 'QC-ALL',
    condition = 0,
    method = 'lm',
    k = 5,
    PDF = TRUE,
    tc = TRUE)

#####################################################################################################################

dir.create('TEST_raw_noQC')                     # Create and set demo-directory (corr. without quality controls)
setwd('TEST_raw_noQC')
write.csv(R, "Test_raw.csv")
write.csv(CD, "Test_ClinData.csv")

# Correct raw data without quality control samples by tobit regression
BCW(RAW = T,
    ts.sw = ',',
    i.cols = c('File.Name', 'Sample.group', 'Site'),
    S = 'SeqNr',
    QC = FALSE,
    QC.rem = TRUE,
    QC.n = 'QC-ALL',
    condition = 'c',
    k = 5,
    PDF = TRUE,
    tc = TRUE)

######################################################################################################################

dir.create('TEST_comp_QC')                      # Create and set demo-directory (corr. with quality controls)
setwd('TEST_comp_QC')
C <- data("Test_compact")
write.csv('S2W_Test.csv')

# Correct compact (preformatted data) with quality control by robust linear regression
BCW(RAW = FALSE,
    ts.f = ',',
    QC = TRUE,
    QC.n = 'QC-ALL',
    QC.rem = FALSE,
    B = 'Site',
    G = 'Sample.group',
    S = 'SeqNr',
    condition = 0,
    method = 'rlm',
    k = 5,
    PDF = TRUE,
    tc = TRUE)

######################################################################################################################

dir.create('TEST_comp_noQC')                    # Create and set demo-directory (corr. with quality controls)
setwd('TEST_comp_noQC')
C <- data("Test_compact")
write.csv('S2W_Test.csv')

# Correct compact (preformatted data) with quality control by robust linear regression
BCW(RAW = FALSE,
    ts.f = ',',
    QC = TRUE,
    QC.n = 'QC-ALL',
    QC.rem = TRUE,
    B = 'Site',
    G = 'Sample.group',
    S = 'SeqNr',
    condition = 0,
    method = 'rlm',
    k = 5,
    PDF = TRUE,
    tc = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
