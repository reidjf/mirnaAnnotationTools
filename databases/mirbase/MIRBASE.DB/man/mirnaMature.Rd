\name{mirnaMature-class}
\docType{class}

\alias{mirnaMature}
\alias{mirnaMature-class}

\alias{mirbaseMatureBimap}
\alias{mirbaseMatureBimap-class}

\alias{matureAccession}
\alias{matureAccession,mirnaMature-method}
\alias{matureName}
\alias{matureName,mirnaMature-method}
\alias{matureFrom}
\alias{matureFrom,mirnaMature-method}
\alias{matureTo}
\alias{matureTo,mirnaMature-method}
\alias{matureEvidence}
\alias{matureEvidence,mirnaMature-method}
\alias{matureExperiment}
\alias{matureExperiment,mirnaMature-method}
\alias{matureSimilarity}
\alias{matureSimilarity,mirnaMature-method}

\alias{show,mirnaMature-method}

\title{Class "mirnaMature"}

\description{A class to represent mature microRNAs}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("mirnaMature", matureAccession, matureName, matureFrom,
    matureTo, matureEvidence, matureExperiment, matureSimilarity)}.
}

\section{Slots}{
  \describe{
    \item{\code{matureAccession}:}{
      Object of class \code{"character"}
      A character string for the Accession number of a mature microRNA.
    }
    \item{\code{matureName}:}{
      Object of class \code{"character"}
      A character string for the ID (or name) of a mature microRNA.}
    \item{\code{matureFrom}:}{
      Object of class \code{"character"}
      A character string for Start position in precursor sequence of a
      mature microRNA.
    }
    \item{\code{matureTo}:}{
      Object of class \code{"character"}
      A character string for the End position in precursor sequence of a
      mature microRNA.
    }
    \item{\code{matureEvidence}:}{
      Object of class \code{"character"}
      A character string for the Experimental evidence of a mature
      microRNA.
    }
    \item{\code{matureExperiment}:}{
      Object of class \code{"character"}
      A character string for the Experiment description and citation of
      a  mature microRNA.
    }
    \item{\code{matureSimilarity}:}{
      Object of class \code{"character"}
      A character string for the Accession number of similar precursor
      microRNA of a mature microRNA.
    }
  }
}

\section{Methods}{
  \describe{
    \item{matureAccession}{
      \code{signature(object = "mirnaMature")}:
      The get method for slot matureAccession.
    }
    \item{matureName}{
      \code{signature(object = "mirnaMature")}:
      The get method for slot matureName.
    }
    \item{matureFrom}{
      \code{signature(object = "mirnaMature")}:
      The get method for slot matureFrom.
    }
    \item{matureTo}{
      \code{signature(object = "mirnaMature")}:
      The get method for slot matureTo.
    }
    \item{matureEvidence}{
      \code{signature(object = "mirnaMature")}:
      The get method for slot matureEvidence.
    }
    \item{matureExperiment}{
      \code{signature(object = "mirnaMature")}:
      The get method for slot matureExperiment.
    }
    \item{matureSimilarity}{
      \code{signature(object = "mirnaMature")}:
      The get method for slot matureSimilarity.
    }
    \item{show}{
      \code{signature(x = "mirnaMature")}:
      The method for pretty print.
    }
  }
}

\seealso{\code{\link{@ANNOBJPREFIX@MATURE}}}

\keyword{methods}
\keyword{classes}

