\name{distBoot}
\alias{distBoot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap distance distribution
}
\description{
Create a distribution of minimal distances between features by shuffling each row of the original data independently. Used for computed p-values of minimal paths with \code{pBoot}.
}
\usage{
distBoot(x, partial = FALSE, squared = TRUE, output = "orig")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Original dataset of features to be shuffled
}
  \item{partial}{
    Whether to compute partial correlations
}
  \item{squared}{
    Whether to square correlations
}
  \item{output}{
    Output of the internal \code{gelMatrix} call. Either "orig" or     "unif" for forcing a uniform distribution of p-values.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A n_features x n_features matrix of the minimal distance between   each pair of features in the bootstrap distribution.
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
\code{\link{minPath}}, \code{\link{pBoot}}
}
\examples{
# Load data
data(foldChange)

# Plot the bootstrap distribution of minimal paths
bdist <- distBoot(foldChange)
hist(bdist)
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
