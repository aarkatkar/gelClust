\name{gammaTest}
\alias{gammaTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test gamma clustering parameters}
\description{
Test a range of values for the parameter, \code{n}, which represents an upper bound on the number of clusters to find. Optimal clustering is achieved at the parameter that minimizes the modified partition entropy.
}
\usage{
gammaTest(start, end, gm, maxiter = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
    the lower bound of the tested range.
}
  \item{end}{
    the upper bound of the tested range.
}
  \item{gm}{
    the output of a \code{gelMatrix} call.
}
  \item{maxiter}{
    the maximum iterations of the clustering algorithm for each
    test.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A plot where a blue star indicates the minimum. A red "x" denotes a test where convergence was not reached before the maximum iterations.
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
\code{\link{gammaCluster}}
}
\examples{
  data(foldChange)
  gm <- gelMatrix(foldChange)
  gammaTest(2, 15, gm)
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
