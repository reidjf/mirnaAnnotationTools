\name{@ANNOBJPREFIX@MIRNA}
\alias{@ANNOBJPREFIX@MIRNA}
\title{Annotation of miRBase identifiers to associated data}
\description{
  This data set gives mappings between miRBase identifiers and their
  respective associated data.
}
\details{
  This is an R object containing key and value pairs. Keys are miRBase
  identifiers and values are their associated data. Values are simple
  named lists of the following members:
  \describe{
    \item{MiRBase.ID}{The miRBase identifier, it starts with a three
      letter species code.}
    \item{MiRBase.Accession}{Accession number of the miRNA.}
    \item{Seed.m8}{miRNA seed, seven bases from the 5' end og the miRNA.}
    \item{Species}{Scientific name of the species.}
    \item{Mature.sequence}{Mature sequence of the miRNA.}
    \item{Family.conservation}{Level of conservation of the miRNA
      family.}
  }

  The package contains the data for \emph{all} miRBase identifiers, the
  ones for other species as well.
}
\references{
  \url{http://www.targetscan.org} for TargetScan and
  \url{http://microrna.sanger.ac.uk/sequences/} for miRBase.
}
\examples{
## Query some random miRBase identifiers
fams <- sample(ls(@ANNOBJPREFIX@MIRNA), 3)
mget(fams, @ANNOBJPREFIX@MIRNA)
}
\keyword{datasets}
