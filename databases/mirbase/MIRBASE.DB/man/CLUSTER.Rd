\name{@ANNOBJPREFIX@CLUSTER}
\alias{@ANNOBJPREFIX@CLUSTER}

\title{MicroRNA IDs to Clusters}

\description{
  @ANNOBJPREFIX@CLUSTER is an R object that provides mappings between
  microRNA identifiers and other microRNA identifiers within a 10kb
  window ('genomic cluster').
}

\details{
  Each microRNA identifier maps to a vector of a microRNA identifiers
  belonging to a cluster (if any beyound itself).
  This information was computed ad-hoc as it is not included in the
  original data tables provided by mirbase.org.

  Source: @CLUSTERSOURCE@
}

\examples{
    x <- @ANNOBJPREFIX@CLUSTER
    # Get the first cluster containg more than one element
    clL <- eapply(x, length)
    id <- names(clL[clL > 1])[1]
    get(id, x)
}

\keyword{data}
