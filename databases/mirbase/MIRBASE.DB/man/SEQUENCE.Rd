\name{@ANNOBJPREFIX@SEQUENCE}
\alias{@ANNOBJPREFIX@SEQUENCE}

\title{MicroRNA IDs to Sequence}

\description{
  @ANNOBJPREFIX@SEQUENCE is an R object that provides mappings between
  microRNA identifiers and their precursor sequence (stem-loop).
}

\details{
  Each microRNA identifier maps to a unique character string
  representing the precursor (stem-loop) sequence of the microRNA.
  A graphical representation of the folded sequence can be found by
  using \code{\link{@ANNOBJPREFIX@HAIRPIN}}.
  
  Source: @SEQUENCESOURCE@
}

\seealso{\code{\link{@ANNOBJPREFIX@HAIRPIN}}}

\examples{
    x <- @ANNOBJPREFIX@SEQUENCE
    # Get the microRNA identifiers that are mapped to a SEQUENCE
    mapped_keys <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_keys])
    if(length(xx) > 0) {
        # Get the SEQUENCE for the first five entries
        xx[1:5]
    }
}

\keyword{data}

