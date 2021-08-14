\name{eigenGenes}
\alias{eigenGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Eigengenes of clusters
}
\description{
Obtain eigengenes representative of each cluster found by \code{gammaCluster}. This is useful for creating heatmaps and simplifying the network.
}
\usage{
eigenGenes(x, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    The original feature dataset used to construct the
    \code{gelMatrix}
}
  \item{m}{
    A matrix of membership probabilities obtained from the
    \code{membership} attribute of \code{gammaCluster}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A matrix of dimensions (n_eigengenes, n_samples) that contains
  representations of each cluster.
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
clust <- gammaCluster(gm, 7)
egenes <- eigenGenes(foldChange, clust$membership)
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
