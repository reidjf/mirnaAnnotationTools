\name{@ANNOBJPREFIX@TARGETSFULL}
\alias{@ANNOBJPREFIX@TARGETSFULL}
\title{Annotation data for TargetScane miRNA target predictions}
\description{
  @ANNOBJPREFIX@TARGETSFULL maps Entrez gene idenifiers to miRNA
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
  
  @ANNOBJPREFIX@TARGETSFULL is an R object containing key and value
  pairs. Keys are Entrez gene identifiers and values are lists, each
  list member corresponds to a predicted targeting miRNA family. Each
  list member is itself a list with the following members:
  \describe{
    \item{miR.Family}{Name of the miRNA family.}
    \item{UTR.start}{Start position on UTR.}
    \item{UTR.end}{End position on UTR.}
    \item{MSA.start}{Start of multiple sequence alignment, with gaps.}
    \item{MSA.end}{End of multiple sequence alignment.}
    \item{Seed.match}{Type of the matching. Possible values are
      \sQuote{8mer}, which means that all eight bases match, \sQuote{m8},
      which means an exact match to positions 2-8 of the mature miRNA, and
      \sQuote{1A}, which means that positions 2-7 of the mature miRnA match,
      followed by an 'A'.}
    \item{PCT}{The probability of conserved targeting, see Friedman et
      al. in the references section below for details.}
  }
}
\references{
  Friedman RC, Farh KK, Burge CB, Bartel DP.: Most mammalian mRNAs
  are conserved targets of microRNAs. Genome Res. 2009 Jan;19(1):92-105
  
  \url{http://www.targetscan.org} for TargetScan and
  \url{http://www.mirbase.org/} for miRBase.
}
\seealso{
 The miRBase database is also provided as a custom annotation package in
 Bioconductor: mirbase.db.
}
\examples{
## Get the miRNA families that regulate these genes
genes <- sample(ls(@ANNOBJPREFIX@TARGETSFULL), 3)
mget(genes, @ANNOBJPREFIX@TARGETSFULL)
}
\keyword{datasets}
