\name{minPath}
\alias{minPath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Minimal path
}
\description{
  Obtain the path of links through the feature space that links
  two traits together such that the probability of obtaining any
  link from the null distribution is minimized. For the nodeNames
  vector, note that trait names come before feature names.
}
\usage{
minPath(s, t, gm, clusters, r = 1, plotpath = TRUE, nodeNames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
    a string (if nodeNames is supplied) or an index (otherwise)
    corresponding to the starting node of the path
}
  \item{t}{
    a string (if nodeNames is supplied) or an index (otherwise)
    corresponding to the ending node of the path
}
  \item{gm}{
    the output of \code{gelMatrix}.
}
  \item{clusters}{
    a numeric vector of cluster labels, usually obtained from the
    \code{labels} attribute of \code{gammaCluster}
}
  \item{r}{
    a regularization term. More info in the details section.
}
  \item{plotpath}{
    whether to plot the minimal path.
}
  \item{nodeNames}{
    a character vector of names for the nodes in the network.
    Trait names come before feature names
}
}
\details{
For regularization, consider the following three options, where N is the total number of nodes:

  1. r = 1 (no regularization)

  2. r = 1 - 2 / (N*(N-1)) (the smallest meaningful p-value is the probability of drawing the single greatest association out of the entire network)

  3. r = 1 - 1/(N-1) (the smallest meaningful p-value is the probability of drawing a given feature's closest neighbor out of the entire network)
}
\value{
  \item{pvals}{the individual p-values of each link from the first node to the last node. The sum of these p-values is the total distance}
  \item{clusters}{cluster labels for each feature in the path}
  \item{names}{the names of the nodes in the path}
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
\code{\link{distBoot}}, \code{\link{pBoot}}
}
\examples{
# Load data
data(foldChange)
data(traitData)

# Obtain clusters
c.gm <- gelMatrix(foldChange)
clust <- gammaCluster(c.gm, 7)

# Set up the path data
path.gm <- gelMatrix(foldChange, phen=traitData, squared=TRUE, output="orig")
nodeNames <- c(rownames(traitData), rownames(foldChange))

# Find the minimal path
mp <- minPath("TCDD", "Fat", path.gm, clust$labels, nodeNames=nodeNames)
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
