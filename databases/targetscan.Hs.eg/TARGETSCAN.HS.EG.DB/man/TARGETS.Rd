\name{@ANNOBJPREFIX@TARGETS}
\alias{@ANNOBJPREFIX@TARGETS}
\title{Annotation data for TargetScane miRNA target predictions}
\description{
  @ANNOBJPREFIX@TARGETS maps Entrez gene idenifiers to miRNA
  families. The mappings are taken from the TargetScan database.
}
\details{
  TargetScan predicts biological targets of miRNAs by searching for a
  presence of conserved 8mer and 7mer sites that match the seed region
  of each miRNA.

  There are currently two annotation packages for TargetScan, one for
  human and another one for mouse. These correspond to the
  TargetScanHuman and the TargetScanMouse databases, see the TargetScan
  website for details.
  
  @ANNOBJPREFIX@TARGETS is an R object containing key and value
  pairs. Keys are Entrez gene identifiers and values are miRNA
  families. Values are vectors of various lengths, as some genes as
  predicted to be potentially regulated by many miRNA families. Note
  that a given miRNA family can come up many times in the value vector,
  if the miRNA seed matches the UTR of gene at multiple positions.
}
\references{
  \url{http://www.targetscan.org} for TargetScan and
  \url{http://www.mirbase.org/} for miRBase.
}
\seealso{
  \link{@ANNOBJPREFIX@TARGETSFULL} for more information about the miRNA
  targets.
  The miRBase database is also provided as a custom annotation package
  in Bioconductor: mirbase.db.
}
\examples{
## Get the miRNA families that regulate these genes
genes <- sample(ls(@ANNOBJPREFIX@TARGETS), 3)
mget( genes, @ANNOBJPREFIX@TARGETS )

## Get all targets of a given miRNA family
mget("miR-10", revmap(@ANNOBJPREFIX@TARGETS))
}
\keyword{datasets}
