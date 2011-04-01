\name{@ANNOBJPREFIX@COMMENT}
\alias{@ANNOBJPREFIX@COMMENT}

\title{MicroRNA IDs to Comments}

\description{
  @ANNOBJPREFIX@COMMENTS is an R object that provides mappings between
  microRNA identifiers and comments related to their description.
}

\details{
  Each microRNA identifier maps to a character string.
  The references cited in the comment can be found by using
  \code{\link{@ANNOBJPREFIX@PMID}}.

  Source: @COMMENTSOURCE@
}

\seealso{\code{\link{@ANNOBJPREFIX@PMID}}}

\examples{
  # first microRNA with a comment that includes a citation.
  mir1 <- intersect(mappedkeys(@ANNOBJPREFIX@COMMENT),
                    mappedkeys(@ANNOBJPREFIX@PMID))[1]
  get(mir1, @ANNOBJPREFIX@COMMENT)
  get(mir1, @ANNOBJPREFIX@PMID)
}

\keyword{data}
