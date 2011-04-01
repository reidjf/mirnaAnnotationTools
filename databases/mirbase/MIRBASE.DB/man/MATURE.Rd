\name{@ANNOBJPREFIX@MATURE}
\alias{@ANNOBJPREFIX@MATURE}

\title{MicroRNA IDs to Mature MicroRNAs}

\description{
  @ANNOBJPREFIX@MATURE is an R object that maps microRNA identifiers to
  information related to their corresponding mature microRNA(s).
}

\details{
  Each microRNA identifier maps to a \code{mirnaMATURE} object that has 7
  slots:
   matureAccession: Accession number
   matureName: ID (or name)
   matureFrom: Start position in precursor sequence
   matureTo: End position in precursor sequence
   matureEvidence: Experimental evidence
   matureExperiment: Experiment description and citation
   matureSimilarity: Accession number of similar precursor microRNA

  The citations in the Experiment description can be retrieved by using
  \code{\link{@ANNOBJPREFIX@PMID}}.

  Source: @MATURESOURCE@
}

\examples{
    x <- @ANNOBJPREFIX@MATURE
    mapped_keys <- mappedkeys(x)
    # Get the MATURE for the first element of xx
    get(mapped_keys[1], x)
}

\keyword{data}
