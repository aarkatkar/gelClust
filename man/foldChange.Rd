\name{foldChange}
\alias{foldChange}
\docType{data}
\title{
 TCDD-induced gene expression
}
\description{
  The log2(fold-change) of 932 DEGs identified across 14 mouse
  strains exposed to TCDD.
}
\usage{data("foldChange")}
\format{
  A data frame with 932 genes within the following 14 strains.
  \describe{
    \item{\code{A.J}}{a numeric vector}
    \item{\code{C57Bl.6J}}{a numeric vector}
    \item{\code{CBA.J}}{a numeric vector}
    \item{\code{DBA.J}}{a numeric vector}
    \item{\code{BXD100}}{a numeric vector}
    \item{\code{FVB.nJ}}{a numeric vector}
    \item{\code{C3Heb.FeJ}}{a numeric vector}
    \item{\code{Balb.cj}}{a numeric vector}
    \item{\code{NOD.ShiLtJ}}{a numeric vector}
    \item{\code{NZO.HilTJ}}{a numeric vector}
    \item{\code{BXD91}}{a numeric vector}
    \item{\code{X129S1.SvlmJ}}{a numeric vector}
    \item{\code{BXD40}}{a numeric vector}
    \item{\code{CC019}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(foldChange)
hist(as.matrix(foldChange))
}
\keyword{datasets}
