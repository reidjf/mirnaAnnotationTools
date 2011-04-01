\name{@ANNOBJPREFIX@MFE}
\alias{@ANNOBJPREFIX@MFE}

\title{MicroRNA IDs to Minimum Fold Energy}

\description{
  @ANNOBJPREFIX@MFE is an R object that provides mappings between
  microRNA identifiers and the Minimum Fold Energy of the folded
  precursor (stem-loop sequence).
}

\details{
  Each microRNA identifier maps to a unique numeric value
  representing the Minimum Fold Energy of the folded
  precursor (stem-loop sequence) computed by the RNAfold program from
  the ViennaRNA suite \url{http://www.tbi.univie.ac.at/~ivo/RNA/}.
  
  A graphical representation of the folded sequence can be found by
  using \code{\link{@ANNOBJPREFIX@MFE}}.
  
  % this data is not in the db, it was processed ad-hoc
  % using ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz

  Source: @MFESOURCE@

}

\seealso{\code{\link{@ANNOBJPREFIX@HAIRPIN}}}

\references{
  Hofacker IL, Stadler PF.
  Memory efficient folding algorithms for circular RNA secondary
  structures. Bioinformatics. 2006 May 15; 22(10):1172-6.
  \url{http://www.ncbi.nlm.nih.gov/pubmed/16452114}
}

\examples{
    x <- @ANNOBJPREFIX@MFE
    # Get the microRNA identifiers that are mapped to a MFE
    mapped_keys <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_keys])
    if(length(xx) > 0) {
        # Get the MFE for the first five entries
        xx[1:5]
    }
}

\keyword{data}

