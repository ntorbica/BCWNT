\name{BCWNT-package}
\alias{BCWNT-package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Batch Correction Wehrens by Nenad Torbica
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Package containing an algorithm to perform batch correction on LC-MS Data with a 8-step procedure.
Corrections implemented by R.Wehrens, pipeline built by Nenad Torbica.
Example data is coming soon.
}

%- maybe also 'usage' for other objects documented here.

\details{
%%  ~~ If necessary, more details than the description above ~~
Content (ordered as in pipeline):

1. BCW()              - Batch correction pipeline

2. BCW.r2w()          - Raw format (implemented for SIMCA output) to compact format

3. BCW.format()       - Compact format to input format for correction function

4. BCW.remZ()         - Removal tool for features richt in zero content

5. BCW.correct()      - Correction containing functions by R. Wehrens

6. BCW.p.abundance()  - P-value computation for uncorrected and corrected data

7. BCW.dist.mat()     - Bhattacharyya distance computation for all possible principal component combinations

8. BCW.plot.loading() - Loading plots of principal components


}

\references{
%% ~put references to the literature/web site here ~
"@ArticleWehrens2016, author =	 Ron Wehrens and Jos.~A.~Hageman and Fred~van~Eeuwijk and Rik~Kooke and P\'adraic~J.~Flood and Erik Wijnker and Joost~J.B.~Keurentjes and Arjen~Lommen and Henri\"ette~D.L.M.~van~Eekelen and Robert~D.~Hall and Roland~Mumm and Ric~C.H.~de~Vos, title =	 Improved batch correction in untargeted MS-based metabolomics, journal =	 Metabolomics, year =	 2016, volume =	 12, DOI =	 10.1007/s11306-016-1015-8, pages =	 1–12

"
}
\author{
%%  ~~who you are~~
Nenad Torbica
}
