\name{mirnaPmid-class}
\docType{class}

\alias{mirnaPmid}
\alias{mirnaPmid-class}

\alias{mirbasePmidBimap}
\alias{mirbasePmidBimap-class}

\alias{pmidAuthor}
\alias{pmidAuthor,mirnaPmid-method}
\alias{pmidTitle}
\alias{pmidTitle,mirnaPmid-method}
\alias{pmidJournal}
\alias{pmidJournal,mirnaPmid-method}
\alias{pmidMedline}
\alias{pmidMedline,mirnaPmid-method}
\alias{pmidOrderAdded}
\alias{pmidOrderAdded,mirnaPmid-method}

\alias{show,mirnaPmid-method}


\title{Class "mirnaPmid"}

\description{A class to represent microRNA references}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("mirnaPmid", pmidAuthor, pmidTitle, pmidJournal,
  pmidMedline, pmidOrderAdded)}.
}

\section{Slots}{
  \describe{
    \item{\code{pmidAuthor}:}{
      Object of class \code{"character"}
      A character string for the authors of a reference.
    }
    \item{\code{pmidTitle}:}{
      Object of class \code{"character"}
      A character string for the title of a reference.
    }
    \item{\code{pmidJournal}:}{
      Object of class \code{"character"}
      A character string for the citation of a reference.
    }
    \item{\code{pmidMedline}:}{
      Object of class \code{"numeric"}
      A numeric vector for the Pubmed identifier of a reference.
    }
    \item{\code{pmidOrderAdded}:}{
      Object of class \code{"numeric"}
      A numeric vector for the order of the reference.
    }
  }
}

\section{Methods}{
  \describe{
    \item{pmidAuthor}{
      \code{signature(object = "mirnaPmid")}:
      The get method for slot pmidAuthor.
    }
    \item{pmidTitle}{
      \code{signature(object = "mirnaPmid")}:
      The get method for slot pmidTitle.
    }
    \item{pmidJournal}{
      \code{signature(object = "mirnaPmid")}:
      The get method for slot pmidJournal.
    }
    \item{pmidMedline}{
      \code{signature(object = "mirnaPmid")}:
      The get method for slot pmidMedline.
    }
    \item{pmidOrderAdded}{
      \code{signature(object = "mirnaPmid")}:
      The get method for slot pmidOrderAdded.
    }
    \item{show}{
      \code{signature(x = "mirnaPmid")}:
      The method for pretty print.
    }
  }
}

\seealso{\code{\link{@ANNOBJPREFIX@PMID}}}

\keyword{methods}
\keyword{classes}

