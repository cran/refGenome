\name{refJunctions-class}
\Rdversion{1.1}
\docType{class}
\alias{refJunctions-class}
\alias{ensemblJunctions-class}
\alias{ucscJunctions-class}
\alias{refJunctions}
\alias{ensemblJunctions}
\alias{ucscJunctions}
\alias{getSpliceTable}
\alias{getSpliceTable-methods}
\alias{getSpliceTable,refGenome-method}
\title{Class \code{"refJunctions"}}
\description{
refJunctions represents ensembl genomic annotation data for splice-junctions. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form getSpliceTable(rg) where 'rg' is
an object of class \code{refGenome} (ensemblGenome or ucscGenome).
}
\section{Slots}{
\describe{
    \item{\code{basedir}:}{Object of class \code{"character"}
        Directory where SQLite database is written.}
    %
    \item{\code{ev}:}{Object of class \code{"environment"} 
        Environment that contains data structures.
        The splice-junction data is stored in 'gtf' named data.frame inside ev.
        Optionally the environment also contains the result tables from 
        \code{unifyJuncs} or \code{getGenePositions} functions. 
        They are automatically included in save and load procedures. }
    %
}}
\section{Methods}{
    \describe{
        \item{show}{Creates a sensible printout.}
        %
        \item{unifyJuncs}{\code{signature(object = "refJunctions")}:
            Calculates unique splice-sites and associates each site with gene-id.
            Adds uid to 'gtf' table and creates a new 'unique junction site'
            (ujs) table. }
        %
        \item{getGenePositions}{\code{signature(object="refJunctions", 
            force="logical")}:
            Extracts table with position data for whole genes (smallest exon
            start position and largest exon end position.
            A copy of the table will be placed inside the internal environment.
            Upon subsequent call only a copy of the contained table is returned
            unless \code{force=TRUE} is given.
            Upon \code{force=TRUE} new gene positions are calculated regardless
            of existing tables.)}
        %
}}
\details{The getSpliceTable function eventually contains a "transcript_biotype"
column (present in Ensembl genome). This coulumn is used by the "unifyJuncs"-
method.
}
\references{
http://www.ensembl.org/info/data/ftp/index.html
http://mblab.wustl.edu/GTF22.html#fields
}
\author{
Wolfgang Kaisers
}
\examples{
# + + + + + + + + + + + + + + + + + + + #
# A) Ensembl
# + + + + + + + + + + + + + + + + + + + #
ef<-system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ef)
enj<-getSpliceTable(ens)
ufe<-unifyJuncs(enj)
saveGenome(enj,"enj.RData",useBasedir=FALSE)
enjr<-loadGenome("enj.RData")

# + + + + + + + + + + + + + + + + + + + #
# B) UCSC
# + + + + + + + + + + + + + + + + + + + #
uf<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(uf)
ucj<-getSpliceTable(uc)
ufu<-unifyJuncs(ucj)
saveGenome(ucj,"ucj.RData",useBasedir=FALSE)
ucjr<-loadGenome("ucj.RData")
}
\keyword{refJunctions}
