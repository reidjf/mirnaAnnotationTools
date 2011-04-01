\name{@ANNOBJPREFIX@CHR}
\alias{@ANNOBJPREFIX@CHR}

\title{MicroRNA IDs to Chromosomes}

\description{
  @ANNOBJPREFIX@CHR is an R object that provides mappings between
  microRNA identifiers and the chromosome that contains the microRNA of
  interest.
}

\details{
  Each microRNA identifier maps to a vector of character strings
  representing possibly multiple chromosomes.

  Source: @CHRSOURCE@
}

\examples{
    x <- @ANNOBJPREFIX@CHR
    # Get the microRNA identifiers that are mapped to a chromosome
    mapped_keys <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_keys])
    if(length(xx) > 0) {
        # Get the CHR for the first five entries
        xx[1:5]
    }
}

\keyword{data}
