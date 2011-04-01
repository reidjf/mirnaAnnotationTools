\name{@ANNOBJPREFIX@DESCRIPTION}
\alias{@ANNOBJPREFIX@DESCRIPTION}

\title{MicroRNA IDs to Descriptions}

\description{
  @ANNOBJPREFIX@DESCRIPTION is an R object that provides mappings between
  microRNA identifiers and their full names.
}

\details{
  Each microRNA identifier maps to a character string.

  Source: @DESCRIPTIONSOURCE@
}

\examples{
    x <- @ANNOBJPREFIX@DESCRIPTION
    # first 3 entries with a description.
    mget(mappedkeys(x)[1:3], x)
}

\keyword{data}
