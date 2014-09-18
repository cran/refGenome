\name{refExons-class}
\Rdversion{1.1}
\docType{class}
\alias{refExons-class}
\alias{ensemblExons-class}
\alias{ucscExons-class}
\alias{refExons}
\alias{refExons-methods}
\alias{refExons,refGenome-method}
\alias{ensemblExons}
\alias{ucscExons}

\title{Class \code{"refExons"}}
\description{
refExons represents genomic annotation data on exon-features for Ensembl
and UCSC genomes.
Relative locations of CDS, start_codon and stop_codon features are added.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{refExons(rg)}.
'rg' represents an object of class \code{ensemblGenome} or \code{ucscGenome}.
}
\section{Slots}{
    \describe{
        \item{\code{basedir}:}{Object of class \code{"character"} Directory
            where SQLite database is written.}
        %
        \item{\code{ev}:}{Object of class \code{"environment"}
            Environment that contains data structures. 
            Optionally, there are gtf and attr data.frames. }
        %
    }
}
\section{Methods}{
    \describe{
        \item{show}{\code{signature(object = "ensemblExons")}: 
            Creates a sensible printout. }
        %
        \item{getSpliceTable}{\code{signature(object = "ensemblExons",
            coding="logical")}: 
            Returns tabled splice sites. 
            When \code{coding}=TRUE only entries with
            \code{gene_biotype=="protein_coding"} are included.}
        %
    }
}
\references{
\describe{
  \item{Ensembl File index}{http://www.ensembl.org/info/data/ftp/index.html}
  \item{GTF Field definitions}{http://mblab.wustl.edu/GTF22.html#fields}
  \item{UCSC home page}{http://genome.ucsc.edu/}
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

# + + + + + + + + + + + + + + + + + + #
# B) UCSC
# + + + + + + + + + + + + + + + + + + #
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(ucfile)
ucex<-refExons(uc)
saveGenome(ucex, "ucex.RData", useBasedir=FALSE)
ur<-loadGenome("ucex.RData")
}
\keyword{classes}
