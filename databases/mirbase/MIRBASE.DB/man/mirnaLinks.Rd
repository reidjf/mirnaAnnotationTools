\name{mirnaLinks-class}
\docType{class}

\alias{mirnaLinks}
\alias{mirnaLinks-class}

\alias{mirbaseLinksBimap}
\alias{mirbaseLinksBimap-class}

\alias{linksDbLink}
\alias{linksDbLink,mirnaLinks-method}
\alias{linksDbId}
\alias{linksDbId,mirnaLinks-method}
\alias{linksDbSecondary}
\alias{linksDbSecondary,mirnaLinks-method}

\alias{show,mirnaLinks-method}

\title{Class "mirnaLinks"}

\description{A class to represent microRNA external database links}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("mirnaLinks", linksDbLink, linksDbId, linksDbSecondary)}.
}

\section{Slots}{
  \describe{
    \item{\code{linksDbLink}:}{
      Object of class \code{"character"}
      A character string for the Accession number of a microRNA.
    }
    \item{\code{linksDbId}:}{
      Object of class \code{"character"}
      A character string for the Database of a microRNA.
    }
    \item{\code{linksDbSecondary}:}{
      Object of class \code{"character"}
      A character string for the Second accession number (or name) of a
      mature microRNA.
    }
  }
}

\section{Methods}{
  \describe{
    \item{linksDbLink}{
      \code{signature(object = "mirnaLinks")}:
      The get method for slot linksDbLink.
    }
    \item{linksDbId}{
      \code{signature(object = "mirnaLinks")}:
      The get method for slot linksDbId.
    }
    \item{linksDbSecondary}{
      \code{signature(object = "mirnaLinks")}:
      The get method for slot linksDbSecondary.
    }
    \item{show}{
      \code{signature(x = "mirnaLinks")}:
      The method for pretty print.}
  }
}

\seealso{\code{\link{@ANNOBJPREFIX@LINKS}}}

\keyword{methods}
\keyword{classes}

