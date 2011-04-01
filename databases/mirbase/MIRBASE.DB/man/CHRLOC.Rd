\name{@ANNOBJPREFIX@CHRLOC}
\alias{@ANNOBJPREFIX@CHRLOC}
\alias{@ANNOBJPREFIX@CHRLOCEND}

\title{MicroRNA IDs to Chromosomal Location}

\description{
  @ANNOBJPREFIX@CHRLOC is an R object that maps microRNA identifiers to
  the starting position of the microRNA.
  The position of a microRNA is measured as the number of base pairs.

  The CHRLOCEND mapping is the same as the CHRLOC mapping except that it
  specifies the ending base of a microRNA instead of the start.
}
\details{
  Each microRNA identifier maps to a named vector of chromosomal locations,
  where the name indicates the chromosome.

  Chromosomal locations on both the sense and antisense strands are
  measured as the number of base pairs from the p (5' end of the sense
  strand) to q (3' end of the sense strand) arms.
  Chromosomal locations on the antisense strand have a leading "-" sign
  (e. g. -1234567).

  Since some microRNAs have multiple start sites, this field can map to
  multiple locations.

  Source: @CHRLOCSOURCE@
}

\examples{
    x <- @ANNOBJPREFIX@CHRLOC
    # Get the microRNA identifiers that are mapped to chromosome locations
    mapped_keys <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_keys])
    if(length(xx) > 0) {
        # Get the CHRLOC for the first five entries
        xx[1:5]
    }
}

\keyword{data}

