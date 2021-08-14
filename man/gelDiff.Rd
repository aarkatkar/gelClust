\name{gelDiff}
\alias{gelDiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Laplace p-value matrix
}
\description{
  Obtain a p-value matrix of the laplace distributed difference
  between two different outputs of the \code{gelMatrix} function.
  If \code{twoTailed == FALSE}, the p-value describes correlations
  that are more significant in the first \code{gelMatrix} than the
  second.
}
\usage{
gelDiff(gm1, gm2, twoTailed = FALSE, unifOut = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gm1}{
    The first \code{gelMatrix}
}
  \item{gm2}{
    The second \code{gelMatrix}
}
  \item{twoTailed}{
    Whether to make the resulting p-values two-tailed
}
  \item{unifOut}{
    Whether to force the output to the uniform distribution
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A p-value matrix of the same dimensions as the each input
  \code{gelMatrix}. This matrix may be used in place of
  \code{gelMatrix} for any function in this package.
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
\code{\link{gelMatrix}}
}
\examples{
data(foldChange)

# Compare the networks between two sets of samples
gm1 <- gelMatrix(foldChange[,1:7])
gm2 <- gelMatrix(foldChange[,8:14])

gd <- gelDiff(gm1, gm2, twoTailed=TRUE)
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
