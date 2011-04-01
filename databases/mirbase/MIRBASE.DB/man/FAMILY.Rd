\name{@ANNOBJPREFIX@FAMILY}
\alias{@ANNOBJPREFIX@FAMILY}

\title{MicroRNA IDs to Family}

\description{
  @ANNOBJPREFIX@FAMILY is an R object that provides mappings between
  microRNA identifiers and their family.
}

\details{
  Each microRNA identifier maps to a character string which is the
  miRNA gene family ID.
  The name of the returned value is the family name.

  Source: @FAMILYSOURCE@
}

\examples{
    x <- @ANNOBJPREFIX@FAMILY
    # first 3 microRNA with a family.
    mget(mappedkeys(x)[1:3], x)
}

\keyword{data}
