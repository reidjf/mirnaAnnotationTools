\name{@ANNOBJPREFIX@FAMILY2MIRBASE}
\alias{@ANNOBJPREFIX@FAMILY2MIRBASE}
\title{An annotation data object that maps miRNA families to miRBase
  identifiers}
\description{
  @ANNOBJPREFIX@FAMILY2MIRBASE maps miRNA family names to miRBase
  identifiers. The mappings are taken from the TargetScan database.
}
\details{
  This is an R object containing key and value pairs. Keys are miRNA
  family names and values are miRBase identifiers. Values are vectors of
  various lengths, as a single miRNA family usually corresponds to many
  miRBase miRNAs.
}
\references{
  \url{http://www.targetscan.org} for TargetScan and
  \url{http://www.mirbase.org/} for miRBase.
}
\seealso{
 The miRBase database is also provided as a custom annotation package in
 Bioconductor: mirbase.db.
}
\examples{
## Map some random miRNA families
fams <- sample(ls(@ANNOBJPREFIX@FAMILY2MIRBASE), 3)
mget(fams, @ANNOBJPREFIX@FAMILY2MIRBASE)
}
\keyword{datasets}
