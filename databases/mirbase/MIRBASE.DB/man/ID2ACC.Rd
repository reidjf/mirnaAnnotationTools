\name{@ANNOBJPREFIX@ID2ACC}
\alias{@ANNOBJPREFIX@ID2ACC}
\alias{@ANNOBJPREFIX@ACC2ID}

\title{MicroRNA IDs to Accessions}

\description{
  @ANNOBJPREFIX@ID2ACC is an R object that provides mappings between
  microRNA identifiers and their Accession numbers.
}

\details{
  Each microRNA identifier maps to a unique Accession number.

  Source: @ID2ACCSOURCE@
}

\examples{
    x <- @ANNOBJPREFIX@ID2ACC
    # Get the microRNA identifiers that are mapped to an Accession
    mapped_keys <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_keys])
    if(length(xx) > 0) {
      # Get the ID2ACC for the first five entries
      xx[1:5]
     }
}

\keyword{data}
