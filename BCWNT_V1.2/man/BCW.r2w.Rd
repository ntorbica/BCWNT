\name{BCW.r2w}
\alias{BCW.r2w}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Conversion from raw format to compact format conversion
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes Raw output, which is comprised of two .csv-files, one with peak-data and one with grouping information. The two tables are wrapped together to a handier format for further manipulation.
}
\usage{
BCW.r2w(peaks, clin.dat, ts, long = TRUE, i.cols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{
%%     ~~Describe \code{peaks} here~~
The file containing peaks, ending with '_raw.csv' in the file name.
}
  \item{clin.dat}{
%%     ~~Describe \code{clin.dat} here~~
The file containing the grouping information, ending with '_ClinData.csv' in the file name.
}
  \item{ts}{
%%     ~~Describe \code{ts} here~~
Table separator for the input files, make sure it is the same for both.
}
  \item{long}{
%%     ~~Describe \code{long} here~~
Logical; specifies whether the peak files orientation (long = samples in columns, wide = samples in rows).
}
  \item{i.cols}{
%%     ~~Describe \code{i.cols} here~~
Vector with column names within the information file. Should contain the names of sample description, sample grouping and batch information. Note that the sample grouping column to chose should contain the labels for quality control samples.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function takes the Raw output files and wraps them up to one compact file with three necessary information colums from the information file and the feature columns from the peak file. The process also orders the data by the injection sequence usually given in the sample description. The function may fail at this point if sample naming does not contain any sequence information.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a matrix object containing the wrapped table. Also outputs a file with 'R2W_' at the beginning of the original filename.
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
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
