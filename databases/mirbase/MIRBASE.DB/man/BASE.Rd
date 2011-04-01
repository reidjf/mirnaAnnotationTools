\name{@ANNOBJPREFIX@.db}
\alias{@ANNOBJPREFIX@.db}
\alias{@ANNOBJPREFIX@}

\title{Bioconductor Custom Annotation Data Package}

\description{
  Welcome to the @ANNOBJPREFIX@.db custom annotation package.
  This package contains multiple organisms.
  The purpose is to provide detailed information about the mirRBase
  microRNA database (\url{http://www.mirbase.org/}).
  In particular, it covers the searchable database of published microRNA
  sequences and annotation and not the miRBase Registry and Targets
  database.
  As requested by the providers of this resource, the references below
  should be cited when making use of the data.
  Ambros et al. 2003 provides guidelines on microRNA annotation.

  The current version and release date are:
  @CHRSOURCE@.
  This information is also available by using the 
  \code{@ANNOBJPREFIX@()} function.

  This package is updated biannually.

  You can learn what objects this package supports with the following
  command:
  \code{ls("package:@ANNOBJPREFIX@.db")}

  Each of these objects has their own manual page detailing where relevant
  data was obtained along with examples of how to use it.
  Many of these objects also have a reverse map available.
  When this is true, expect to usually find relevant information on the
  same manual page as the forward map.
}

\references{
  Griffiths-Jones S, Saini HK, van Dongen S, Enright AJ.
  miRBase: tools for microRNA genomics. NAR 2008 36(Database
  Issue):D154-D158
  \url{http://nar.oxfordjournals.org/cgi/content/full/36/suppl_1/D154}

  Griffiths-Jones S, Grocock RJ, van Dongen S, Bateman A, Enright AJ.
  miRBase: microRNA sequences, targets and gene nomenclature. NAR 2006
  34(Database Issue):D140-D144
  \url{http://nar.oxfordjournals.org/cgi/content/full/34/suppl_1/D140}

  Griffiths-Jones S.
  The microRNA Registry. NAR 2004 32(Database Issue):D109-D111
  \url{http://nar.oupjournals.org/cgi/content/full/32/suppl_1/D109}

  Ambros V, Bartel B, Bartel DP, Burge CB, Carrington JC, Chen X,
  Dreyfuss G, Eddy SR, Griffiths-Jones S, Marshall M, Matzke M, Ruvkun
  G, Tuschl T.
  A uniform system for microRNA annotation. RNA 2003 9(3):277-279
  \url{http://www.rnajournal.org/cgi/content/full/9/3/277}
}

\examples{
  ls("package:@ANNOBJPREFIX@.db")
}

\keyword{data}
