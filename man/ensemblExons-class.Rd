\name{ensemblExons-class}
\Rdversion{1.1}
\docType{class}
\alias{ensemblExons-class}
\alias{ensemblExons}
\alias{ensemblExons-methods}
\alias{ensemblExons,ensemblGenome-method}
\alias{getSpliceTable}
\alias{getSpliceTable-methods}
\alias{getSpliceTable,ensemblExons-method}
\alias{load.ensembl.exons}


\title{Class \code{"ensemblExons"}}
\description{
ensemblExons represents ensembl genomic annotation data where exon-features are extracted from primary assembly seqid's. Relative locations of CDS, start_codon and stop_codon features are added.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{ensemblExons(eg)}. 'eg' represents an object of class ensemblGenome.
}
\section{Slots}{
  \describe{
    \item{\code{basedir}:}{Object of class \code{"character"} Directory where SQLite database is written.  }
    \item{\code{ev}:}{Object of class \code{"environment"} Environment that contains data structures. Optionally, there are gtf and attr data.frames. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "ensemblExons")}: Creates a sensible printout. }
    \item{getSpliceTable}{\code{signature(object = "ensemblExons")}: Returns tabled splice sites. }
	 }
}
\references{
http://www.ensembl.org/info/data/ftp/index.html
http://mblab.wustl.edu/GTF22.html#fields
}
\author{
Wolfgang Kaisers
}
\examples{
# +++++++++++++++++++++++++++++++++++++ #
# Loading of a saved genome
ensfile<-system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens<-load.ensembl(ensfile)
enex<-ensemblExons(ens)
# Save and load ensemblExonData
saveGenome(ens,"ensemblExons62.RData",useBasedir=FALSE)
load.ensembl.exons("ensemblExons62.RData")
junc<-getSpliceTable(enex)
# +++++++++++++++++++++++++++++++++++++ #
}
\keyword{classes}
