\name{@ANNOBJPREFIX@LINKS}
\alias{@ANNOBJPREFIX@LINKS}

\title{MicroRNA IDs to External Database Links}

\description{
  @ANNOBJPREFIX@LINKS is an R object that maps microRNA identifiers to a
  set external database identifiers.
}

\details{
  Each microRNA identifier maps to a \code{mirnaLinks} object that has 3
  slots:
   linksDbLink: Accession number
   linksDbId: Database
   linksDbSecondary: Second accession number (or name)

  Source: @LINKSSOURCE@
}

\seealso{\code{\link{mirnaLinks}}}

\examples{
    x <- @ANNOBJPREFIX@LINKS
    mapped_keys <- mappedkeys(x)
    # Get the LINKS for the first element of xx
    get(mapped_keys[1], x)
}

\keyword{data}
