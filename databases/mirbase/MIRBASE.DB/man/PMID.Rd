\name{@ANNOBJPREFIX@PMID}
\alias{@ANNOBJPREFIX@PMID}

\title{MicroRNA IDs to References}

\description{
  @ANNOBJPREFIX@PMID is an R object that maps microRNA identifiers to a
  set of references (PubMed IDentifiers).
}

\details{
  Each microRNA identifier maps to a \code{\link{mirnaPmid}} object that
  has 5 slots:
   pmidAuthor: List of authors
   pmidTitle: Title
   pmidJournal: Citation
   pmidMedline: Pubmed identifier
   pmidOrderAdded: Order

  The Order is the same as the one used in
  \code{\link{@ANNOBJPREFIX@COMMENT}} and in
  \code{\link{@ANNOBJPREFIX@MATURE}}.

  Source: @PMIDSOURCE@
}

\seealso{\code{\link{mirnaPmid}}}

\examples{
    x <- @ANNOBJPREFIX@PMID
    mapped_keys <- mappedkeys(x)
    # Get the PMIDs for the first element of xx
    get(mapped_keys[1], x)
}

\keyword{data}
