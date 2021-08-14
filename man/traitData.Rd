\name{traitData}
\alias{traitData}
\docType{data}
\title{
 Phenotypic traits
}
\description{
  Phenotypic traits of 14 mouse strains exposed to TCDD. The traits
  are hepatic TCDD accumulation, change in fat percentage, and
  change in body weight gain.
}
\usage{data("traitData")}
\format{
  A data frame with 3 traits within the following 14 strains.
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
data(traitData)
barplot(as.matrix(traitData["TCDD",]), las=2)
}
\keyword{datasets}
