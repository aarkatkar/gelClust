\name{gammaPlot}
\alias{gammaPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
UMAP visualization of clusters
}
\description{
Plot features using Uniform Manifold Approximation and Projection and colored according to the \code{gammaCluster} function. This is useful for verifying that the gamma clustering model accurately characterizes the underlying structure of the correlations.
}
\usage{
gammaPlot(gm, clusters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gm}{
    the output of the \code{gelMatrix} function
}
  \item{clusters}{
    a numeric vector of cluster labels, usually obtained from
    \code{labels} attribute of the \code{gammaCluster} function
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a plot of features placed according to UMAP and colored according to gamma clustering.
}
\references{
McInnes, L., Healy, J., & Melville, J. (2018). Umap: Uniform manifold approximation and projection for dimension reduction. arXiv preprint arXiv:1802.03426.
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
gammaPlot(gm, clust$labels)
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
