\name{pBoot}
\alias{pBoot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap path p-value
}
\description{
Obtain a p-value for a \code{minPath} by comparison to a bootstrap distribution obtained by \code{distBoot}.
}
\usage{
pBoot(d, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
    Output of the \code{distBoot} function
}
  \item{z}{
    The minimal path length to be tested
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A bootstrap p-value
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
\code{\link{minPath}}, \code{\link{distBoot}}
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

# Obtain the p-value
bdist <- distBoot(foldChange)
p <- pBoot(bdist, sum(mp$pvals))
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
