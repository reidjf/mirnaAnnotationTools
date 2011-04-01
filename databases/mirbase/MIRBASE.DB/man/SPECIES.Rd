\name{@ANNOBJPREFIX@SPECIES}
\alias{@ANNOBJPREFIX@SPECIES}

\title{Organism Acronym to Species}

\description{
  @ANNOBJPREFIX@SPECIES is an R object that provides mappings between
  organism acronyms and the species specifications.
}

\details{
  Each organism acronym to a unique character string which is
  the full name of the species.
  Further information concerning the species can be found by using
  the \code{\link{toTable}} function which reports the 'division',
  'taxonomy', 'genome_assembly' and the 'ensembl_db' used.

  The organism acronym of a microRNA is returned by using
  \code{\link{@ANNOBJPREFIX@ID2SPECIES}}.
  
  Source: @SPECIESSOURCE@
}

\seealso{\code{\link{@ANNOBJPREFIX@ID2SPECIES}}}

\examples{
    x <- @ANNOBJPREFIX@ID2SPECIES
    y <- @ANNOBJPREFIX@SPECIES
    # get full species information for the first microRNA
    mir1 <- mappedkeys(x)[1]
    mir1Species <- get(mir1, x)
    get(mir1Species, y)
    toTable(y[mir1Species])
}

\keyword{data}
