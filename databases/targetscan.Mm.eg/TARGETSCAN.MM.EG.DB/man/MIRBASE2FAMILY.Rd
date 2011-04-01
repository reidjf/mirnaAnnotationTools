\name{@ANNOBJPREFIX@MIRBASE2FAMILY}
\alias{@ANNOBJPREFIX@MIRBASE2FAMILY}
\title{An annotation data object that maps miRBase identifiers to miRNA
  family names}
\description{
  @ANNOBJPREFIX@MIRBASE2FAMILY maps miRBase identifiers to miRNA family
  names. The mappings are taken from the TargetScan database. 
}
\details{
  This is an R object containing key and value pairs. Keys are miRBase
  identifieers and values are miRNA family names. Values are
  vectors of length 1, as a miRBase miRNA belongs to a single miRNA
  family.
}
\references{
  \url{http://www.targetscan.org} for TargetScan and
  \url{http://www.mirbase.org} for miRBase.
}
\seealso{
 The miRBase database is also provided as a custom annotation package in
 Bioconductor: mirbase.db.
}
\examples{
## Map some random miRBase miRNAs
fams <- sample(ls(@ANNOBJPREFIX@MIRBASE2FAMILY), 3)
mget(fams, @ANNOBJPREFIX@MIRBASE2FAMILY)
}
\keyword{datasets}
