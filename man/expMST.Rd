\name{expMST}
\alias{expMST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Minimum spanning tree of exponential dissimilarities
}
\description{
Construct a tree that links all features together such that the likelihood of obtaining any of the linking correlations from the null distribution is minimized. Produce clusters from cutting this tree at links above a certain p-value threshold.
}
\usage{
expMST(gm, cutoff = 0.05, min.size = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gm}{
    the output of the \code{gelMatrix} function.
}
  \item{cutoff}{
    the threshold p-value beyond which links will be severed.
}
  \item{min.size}{
    minimum cluster size to be output.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a numeric vector of cluster labels. \code{0} indicates
  "unclustered"
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Anooj Arkatkar
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(foldChange)
data(traitData)
gm <- gelMatrix(foldChange, traitData, squared=TRUE, output="orig")

# Total number of links in the MST for Bonferroni correction
N <- dim(gm)[1]-1

clust <- expMST(gm, cutoff=0.05/N)
table(clust)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
