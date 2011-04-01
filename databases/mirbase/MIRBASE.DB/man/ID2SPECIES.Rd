\name{@ANNOBJPREFIX@ID2SPECIES}
\alias{@ANNOBJPREFIX@ID2SPECIES}
\alias{@ANNOBJPREFIX@SPECIES2ID}

\title{MicroRNA IDs to Species}

\description{
  @ANNOBJPREFIX@ID2SPECIES is an R object that provides mappings between
  microRNA identifiers and the species they belong to.
}

\details{
  Each microRNA identifier maps to a unique character string which is
  an abbreviated name of the species.
  Further information concerning the species can be found by using
  the name with \code{\link{@ANNOBJPREFIX@SPECIES}}.

  Source: @ID2SPECIESSOURCE@
}

\examples{
    x <- @ANNOBJPREFIX@ID2SPECIES
    # Get the microRNA identifiers that are mapped to a species
    mapped_keys <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_keys])
    if(length(xx) > 0) {
      # Get the ID2SPECIES for the first five entries
      xx[1:5]
     }
}

\keyword{data}
