\name{mirnaContext-class}
\docType{class}

\alias{mirnaContext}
\alias{mirnaContext-class}

\alias{mirbaseContextBimap}
\alias{mirbaseContextBimap-class}

\alias{contextTranscriptID}
\alias{contextTranscriptID,mirnaContext-method}
\alias{contextOverlapSense}
\alias{contextOverlapSense,mirnaContext-method}
\alias{contextOverlapType}
\alias{contextOverlapType,mirnaContext-method}
\alias{contextNumber}
\alias{contextNumber,mirnaContext-method}
\alias{contextTranscriptSource}
\alias{contextTranscriptSource,mirnaContext-method}
\alias{contextTranscriptName}
\alias{contextTranscriptName,mirnaContext-method}

\alias{show,mirnaContext-method}

\title{Class "mirnaContext"}

\description{A class to represent the genomic context of a microRNA}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("mirnaContext", contextTranscriptID, contextOverlapSense,
    contextOverlapType, contextNumber, contextTranscriptSource,
    contextTranscriptName)}.
}

\section{Slots}{
  \describe{
    \item{\code{contextTranscriptID}:}{
      Object of class \code{"character"}
      A character string for Transcript identifer overlapping a
      microRNA.
    }
    \item{\code{contextOverlapSense}:}{
      Object of class \code{"character"}
      A character string for the Strand of transcript.
    }
    \item{\code{contextOverlapType}:}{
      Object of class \code{"character"}
      A character string for the of Type of overlap (exon, intron, 3' or
      5' UTR).
    }
    \item{\code{contextNumber}:}{
      Object of class \code{"numeric"}
      A numeric vector for the Type of overlap number (eg exon 3).
    }
    \item{\code{contextTranscriptSource}:}{
      Object of class \code{"character"}
      A character string for the Transcript database.
    }
    \item{\code{contextTranscriptName}:}{
      Object of class \code{"character"}
      A character string for the Transcript Name.
    }
  }
}

\section{Methods}{
  \describe{
    \item{contextTranscriptID}{
      \code{signature(object = "mirnaContext")}:
      The get method for slot contextTranscriptID.
    }
    \item{contextOverlapSense}{
      \code{signature(object = "mirnaContext")}:
      The get method for slot contextOverlapSense.
    }
    \item{contextOverlapType}{
      \code{signature(object = "mirnaContext")}:
      The get method for slot contextOverlapType.
    }
    \item{contextNumber}{
      \code{signature(object = "mirnaContext")}:
      The get method for slot contextNumber.
    }
    \item{contextTranscriptSource}{
      \code{signature(object = "mirnaContext")}:
      The get method for slot contextTranscriptSource.
    }
    \item{contextTranscriptName}{
      \code{signature(object = "mirnaContext")}:
      The get method for slot contextTranscriptName.
    }
    \item{show}{
      \code{signature(x = "mirnaContext")}:
      The method for pretty print.
    }
  }
}

\seealso{\code{\link{@ANNOBJPREFIX@CONTEXT}}}

\keyword{methods}
\keyword{classes}

