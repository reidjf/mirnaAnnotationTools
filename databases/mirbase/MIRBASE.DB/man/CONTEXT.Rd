\name{@ANNOBJPREFIX@CONTEXT}
\alias{@ANNOBJPREFIX@CONTEXT}

\title{MicroRNA IDs to Genomic Context}

\description{
  @ANNOBJPREFIX@CONTEXT is an R object that maps microRNA identifiers to
  information related to overlapping transcripts.
}

\details{
  Each microRNA identifier maps to a \code{mirnaContext} object that has 6
  slots:
  contextTranscriptID: Transcript identifer
  contextOverlapSense: Strand of transcript
  contextOverlapType: Type of overlap (exon, intron, 3' or 5' UTR)
  contextNumber: Type of overlap number (eg exon 3)
  contextTranscriptSource: Transcript database
  contextTranscriptName: Transcript Name

  Source: @CONTEXTSOURCE@
}

\seealso{\code{\link{mirnaContext}}}

\examples{
    x <- @ANNOBJPREFIX@CONTEXT
    # Get the CONTEXT for the first element mapped 
    get(mappedkeys(x)[1], x)
}

\keyword{data}
