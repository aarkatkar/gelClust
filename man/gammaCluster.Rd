\name{gammaCluster}
\alias{gammaCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perform gamma clustering}
\description{
  Cluster features (such as genes) from the output of
  \code{gelMatrix} under the assumption that inter-cluster
  associations are obtained from the null, gamma distribution.
}
\usage{
gammaCluster(gm, n, maxiter = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gm}{a matrix of p-values from the \code{gelMatrix}
  function.}
  \item{n}{the maximum number of clusters to produce.}
  \item{maxiter}{the maximum number of iterations for the
  algorithm.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{labels}{a numeric vector of the cluster assigned to each
    feature.}
  \item{membership}{the n_clusters x n_features membership
    probability matrix.}
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
  \code{\link{gammaTest}}.
}
\examples{
data(foldChange)
gm <- gelMatrix(foldChange)
clust <- gammaCluster(gm, 7)
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
