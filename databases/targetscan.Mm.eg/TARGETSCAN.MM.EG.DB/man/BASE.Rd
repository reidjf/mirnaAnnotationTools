\name{@ANNOBJPREFIX@.db}
\alias{@ANNOBJPREFIX@.db}
\alias{@ANNOBJPREFIX@}
\title{Bioconductor annotation data package}
\description{

  Welcome to the @ANNOBJPREFIX@.db annotation Package.  The purpose of
  this package is to provide detailed information about the latest
  version of the TargetScan miRNA target prediction database. This
  package is updated  biannually. 
  
  You can learn what objects this package supports with the following command:

  \code{ls("package:@ANNOBJPREFIX@.db")} 
  
  Each of these objects has their own manual page detailing where
  relevant data was obtained along with some examples of how to use it.

  There are currently two annotation packages for TargetScan, one for
  human and another one for mouse. These correspond to the
  TargetScanHuman and the TargetScanMouse databases, see the TargetScan
  website for details.

  Please note, that the package currently contains _only_ the prediction
  for conserved miRNA targets, separately for each miRNA family.
  These are the predicted targets that are returned are the 'conserved
  across mammals' ones which are listed in the table
  'Predicted Conserved Targets Info' from the download page.

}
\references{
  \url{http://www.targetscan.org}
}
\examples{
  ls("package:@ANNOBJPREFIX@.db")
}
\keyword{datasets}

