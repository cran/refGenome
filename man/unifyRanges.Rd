\name{unifyRanges}
\Rdversion{1.1}
%
\alias{unifyRanges}
\alias{unifyRanges-methods}
\alias{unifyRanges,refExons-method}
%
\alias{plotUbs}
\alias{plotUbs-methods}
\alias{plotUbs,unifiedExons-method}
%
\title{Unification of Exonic ranges.}
\description{Unified ranges are intended to be used as preparation of 
exonic coordinates for counting of alignments (reads) inside these regions.
In order to prevent ambiguities, overlaps between exons are removed.
}
\usage{
unifyRanges(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{refExons} object. Contains exonic coordinates data.}
}
\details{Removal of overlaps works as follows: Ranges are sorted
according to their seqid and begin value.
Ranges which are contained in a precedent range are removed.
For ranges with partial overlap with a precedent range, the start position
is shifted to one position after the end of the precedent range.
The size of the shift is contained in the 'ubs' value.
An overview over the size and the number (begin-) shifts can be obtained
using the \code{plotUbs} function.
}
\value{
The \code{unifyRanges} function returns a refExons object. The data table
contains the following columns:
\tabular{ll}{
  id \tab Numeric index retained from the incoming refExons object.\cr
  seqid \tab Chromosome identifier\cr
  begin \tab begin of range\cr
  end \tab end of range\cr
  strand \tab strand for exon (+ or -)\cr
  gene_id \tab Fene identifier\cr
  transcript_id \tab Transcript identifier\cr
  gene_name \tab Gene name\cr
  exon_number \tab Exon numbe\cr
  ubs \tab Unification begin shift. Number by which begin has been increased
    during the unification process\cr
}
}
\author{
Wolfgang Kaisers
}
\examples{
# + + + + + + + + + + + + + + + + + + #
# A) Ensembl
# + + + + + + + + + + + + + + + + + + #
ensfile<-system.file("extdata", 
                        "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ensfile)
enex<-refExons(ens)
saveGenome(enex,"enex.RData", useBasedir=FALSE)
er<-loadGenome("enex.RData")
# Unify exons
uex <- unifyRanges(enex)
plotUbs(uex, xlim=c(0,400), ylim=c(0,5))

# + + + + + + + + + + + + + + + + + + #
# B) UCSC
# + + + + + + + + + + + + + + + + + + #
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(ucfile)
ucex<-refExons(uc)
# Unify exons
uex <- unifyRanges(enex)
plotUbs(uex, xlim=c(0,400), ylim=c(0,5))
}
\keyword{refJunctions}
