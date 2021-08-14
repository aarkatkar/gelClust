\name{gelMatrix}
\alias{gelMatrix}
\title{Obtain a correlation/p-value matrix}
\description{
  Generate a matrix of p-values ("orig"), p-values forced to a
  uniform distribution ("unif"), or correlations ("cor").
}
\usage{
gelMatrix(x, phen = NULL, z = NULL, squared = FALSE, partial = FALSE, output = "unif")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A matrix of dimensions (n_features, n_samples)
}
  \item{phen}{
    A phenotypic trait matrix of dimensions (n_features, n_samples)
}
  \item{z}{
    A covariate matrix of dimensions (n_cov, n_smaples)
}
  \item{squared}{
    Square correlations to ignore their sign
}
  \item{partial}{
    Generate a partial correlation matrix
}
  \item{output}{
    The type of matrix to be calculated ("unif", "orig", "cor")
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A matrix of dimensions (n_phen + n_feat, n_phen + n_feat)
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
gm <- gelMatrix(foldChange)
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
