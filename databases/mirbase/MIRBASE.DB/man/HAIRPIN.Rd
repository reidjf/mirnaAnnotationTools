\name{@ANNOBJPREFIX@HAIRPIN}
\alias{@ANNOBJPREFIX@HAIRPIN}

\title{MicroRNA IDs to Hairpin}

\description{
  @ANNOBJPREFIX@HAIRPIN is an R object that provides mappings between
  microRNA identifiers and an ASCII representation of the the folded
  precursor (stem-loop sequence).
}

\details{
  Each microRNA identifier maps to a character string.
  The mature microRNA(s) sequence are highlighted in capital letters.

  The representation was created the RNAfold program from the ViennaRNA
  suite \url{http://www.tbi.univie.ac.at/~ivo/RNA/}.
  Extra information include the minimum free energy ('MFE') which
  can be found by using \code{\link{@ANNOBJPREFIX@MFE}} and the
  position on the sequence of mature mirna(s) which can be found
  by using \code{\link{@ANNOBJPREFIX@MATURE}}.

  % this data is not in the db, it was processed ad-hoc
  % using ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz

  Source: @HAIRPINSOURCE@

}

\references{
  Hofacker IL, Stadler PF.
  Memory efficient folding algorithms for circular RNA secondary
  structures. Bioinformatics. 2006 May 15; 22(10):1172-6.
  \url{http://www.ncbi.nlm.nih.gov/pubmed/16452114}
}

\examples{
    x <- @ANNOBJPREFIX@HAIRPIN
    # hairpin representation sequences of all microRNAs
    mirnaHairpin <- mget(mappedkeys(x), x)
    # print first one
    cat(mirnaHairpin[[1]], "\n")
}

\keyword{data}

